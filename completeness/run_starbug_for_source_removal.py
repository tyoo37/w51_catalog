
import glob
import time
import numpy
import crowdsource
import regions
import numpy as np
from functools import cache
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
from astropy.modeling.fitting import LevMarLSQFitter
from astropy import wcs
from astropy import table
from astropy import stats
from astropy import units as u
from astropy.nddata import NDData
from astropy.io import fits
from scipy import ndimage
from astropy.stats import sigma_clipped_stats, sigma_clip
import requests
import requests.exceptions
import urllib3
import urllib3.exceptions
from photutils.detection import DAOStarFinder, IRAFStarFinder
from photutils.psf import IntegratedGaussianPRF, extract_stars, EPSFStars, EPSFModel
try:
    # version >=1.7.0, doesn't work: the PSF is broken (https://github.com/astropy/photutils/issues/1580?)
    from photutils.psf import PSFPhotometry, IterativePSFPhotometry, SourceGrouper
except:
    # version 1.6.0, which works
    from photutils.psf import BasicPSFPhotometry as PSFPhotometry, IterativelySubtractedPSFPhotometry as IterativePSFPhotometry, DAOGroup as SourceGrouper
try:
    from photutils.background import MMMBackground, MADStdBackgroundRMS, MedianBackground, Background2D, LocalBackground
except:
    from photutils.background import MMMBackground, MADStdBackgroundRMS, MedianBackground, Background2D
    from photutils.background import MMMBackground as LocalBackground

from photutils.psf import EPSFBuilder
from photutils.psf import extract_stars

import warnings
from astropy.utils.exceptions import AstropyWarning, AstropyDeprecationWarning
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=AstropyDeprecationWarning)

from crowdsource import crowdsource_base
from crowdsource.crowdsource_base import fit_im, psfmod

from astroquery.svo_fps import SvoFps
from astropy.table import Table, vstack

import pylab as pl
pl.rcParams['figure.facecolor'] = 'w'
pl.rcParams['image.origin'] = 'lower'

import os
print("Importing webbpsf", flush=True)
import stpsf as webbpsf
print(f"Webbpsf version: {webbpsf.__version__}")
from webbpsf.utils import to_griddedpsfmodel
import datetime
import subprocess

print("Done with imports", flush=True)
# step i )  load the catalog and do the quality assessment
# for the quality assessment, some tests are required to see the effect of masking
# load the catalog obtained from each exposure file

nvisits = {'2221': {'brick': 1, 'cloudc': 2},
               '1182': {'brick': 2},
               '6151': {'w51': 1 , 'w51_miri': 2 }
               }
def update_param(filepath, param_name, new_value):
    lines = []
    with open(filepath, "r") as f:
        for line in f:
            # Skip blank lines and comment-only lines
            stripped = line.strip()
            if not stripped or stripped.startswith("//") or stripped.startswith("#"):
                lines.append(line)
                continue
                
            if line.strip().startswith(param_name + " "):  # match beginning of line
                parts = line.split("=")
                if len(parts) >= 2:
                    # keep comments if exist
                    comment = "//" + parts[1].split("//")[1] if "//" in parts[1] else ""
                    line = f"{param_name} = {new_value} {comment}\n"
            lines.append(line)

    with open(filepath, "w") as f:
        f.writelines(lines)
    print(f"Updated {param_name} to {new_value}")
def get_filenames(basepath, filtername, proposal_id, field, each_suffix, module, pupil='clear', visitid='001'):

    # jw01182004002_02101_00012_nrcalong_destreak_o004_crf.fits
    # jw02221001001_07101_00012_nrcalong_destreak_o001_crf.fits
    # jw02221001001_05101_00022_nrcb3_destreak_o001_crf.fits
    glstr = f'{basepath}/{filtername}/pipeline/jw0{proposal_id}*{module}*_{each_suffix}.fits'
    
  
    fglob = glob.glob(glstr)
    for st in fglob:
        
        if 'align' in st or 'uncal' in st:
            print(f"Removing {st} from glob string because it is an alignment file")
            fglob.remove(st)
    if len(fglob) == 0:
        raise ValueError(f"No matches found to {glstr}")
    else:
        return fglob

def create_apfile(image, catalogfile, outputfilename, radius=1.5, sky_in=3, sky_out=4.5, error=None, apcorr=1.0, sig_sky=3):
    from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
    mask=np.isnan(image)
    if error is None:
        error=np.sqrt(image)
    catalog = Table.read(catalogfile)
    pos=[(line["x_fit"],line["y_fit"]) for line in catalog if np.isfinite(line["x_fit"]) and np.isfinite(line["y_fit"])]
    print(f'number of finite elements: {len(pos)}, total elements: {len(catalog)}')
    apertures=CircularAperture(pos,radius)
    smooth_apertures=CircularAperture(pos, min(1.5*radius,sky_in))
    annulus_aperture=CircularAnnulus(pos, r_in=sky_in, r_out=sky_out)

    phot=aperture_photometry(image, (apertures,smooth_apertures), error=error, mask=mask)
    masks=annulus_aperture.to_mask(method="center")
    dat=list(map(lambda a:a.multiply(image),masks))

    try: dat=np.array(dat).astype(float)
    except:
        ## Cases where the array is inhomegenoeus
        ## If annulus reaches the edge of the image, it will create a mask the wrong shape
        ## If for whatever reason the point lies outside the image, it will have None
        ## in the list, this needs to be caught too
        print("Ran into issues with the sky annuli, trying to fix them..\n")
        size=np.max( [np.shape(d) for d in dat if d is not None ])
        for i,d in enumerate(dat):
            if d is None: dat[i]=np.zeros((size,size))
            elif (shape:=np.shape(d))!=(size,size):
                dat[i]=np.zeros((size,size))
                dat[i][:shape[0],:shape[1]]+=d
        dat=np.array(dat)

    mask=(dat>0) & np.isfinite(dat)
    dat[~mask]=np.nan
    dat=sigma_clip(dat.reshape(dat.shape[0],-1), sigma=sig_sky,axis=1)
    sky=np.ma.median(dat,axis=1).filled(fill_value=0)

    tab_ap = Table()
    tab_ap['xcentroid'] = catalog['x_fit'][np.isfinite(catalog['x_fit'])]
    tab_ap['ycentroid'] = catalog['y_fit'][np.isfinite(catalog['y_fit'])]
    tab_ap['sky'] = sky
    tab_ap['flux'] = apcorr*(phot["aperture_sum_0"] - (sky*apertures.area))
    tab_ap["flux"][tab_ap["flux"]==0]=np.nan

    print('tab_ap', tab_ap)
    tab_ap.write(outputfilename, format='fits', overwrite=True)
    print(f"Created aperture file {outputfilename}")
    if not os.path.exists(outputfilename):
        raise ValueError(f"AP file {outputfilename} not found after creation")
    return tab_ap

def main():
    
    proposal_id = '6151'
    filternames = ['F140M', 'F162M', 'F187N', 'F182M', 'F210M', 'F335M', 'F360M', 'F405N',
     'F410M', 'F480M', 'F560W', 'F770W', 'F1000W', 'F1280W', 'F2100W']


    basepath = f'/orange/adamginsburg/jwst/w51'
    #paramdir = '/blue/adamginsburg/t.yoo/from_red/w51/w51_catalog/completeness/starbug_params'
   
    cmd = ["starbug2", "--local-params"]

    result = subprocess.run(cmd, capture_output=True, text=True)  

    if True:
        for filtername in filternames:
            #detector = module # no sub-detectors for long-NIRCAM
            if filtername in ['F140M', 'F162M', 'F187N', 'F182M', 'F210M']:
                modules = ['nrca1', 'nrca2', 'nrca3', 'nrca4', 'nrcb1', 'nrcb2', 'nrcb3', 'nrcb4']
                target = 'w51'
            elif filtername in ['F335M', 'F360M', 'F410M', 'F405N', 'F480M']:
                modules = ['nrcalong', 'nrcblong']
                target = 'w51'
            elif filtername in ['F560W', 'F770W', 'F1000W', 'F1280W', 'F2100W']:
                modules = ['mirimage']
                target = 'w51_miri'
            
            field_to_reg_mapping = {'2221': {'001': 'brick', '002': 'cloudc'},
                                        '1182': {'004': 'brick'},
                                        '6151': {'001': 'w51', '002': 'w51_miri'}}[proposal_id]
            reg_to_field_mapping = {v:k for k,v in field_to_reg_mapping.items()}
            field = reg_to_field_mapping[target]
            for module in modules:
                if True:
                    for visitid in range(1, nvisits[proposal_id][target] + 1):
                        visitid = f'{visitid:03d}'
                        
                        filenames = get_filenames(basepath, filtername, proposal_id,
                                                field, visitid=visitid,
                                                each_suffix='cal',
                                                module=module, pupil='clear')
                        print('filenames:', filenames)
                        for i, filename in enumerate(filenames):
                            if True:
                                exposurenumber = int(filename.split("_")[2])
                                exposure_id = filename.split("_")[2]
                                visit_id = filename.split("_")[0][-3:]
                                vgroup_id = filename.split("_")[1]
                                exposure_ = f'_exp{exposurenumber:05d}'
                                visitid_ = f'_visit{visitid}' if visitid is not None else ''
                                vgroupid_ = f'_vgroup{vgroup_id}' if vgroup_id is not None else ''
                                catalogfile = f"{basepath}/{filtername}/{filtername.lower()}_{module}{visitid_}{vgroupid_}{exposure_}_daophot_basic.fits"
                                print('catalogfile:', catalogfile)
                                #catalog = fits.open(catalogfile)
                                # change the format of the catalog accordingly to be used in starbug2
                                

                                # step ii) using refined catalogs, run starbugs2 to remove background for each exposure file
                                print('create starbug2 param file')
                                img = fits.open(filename)['SCI'].data
                                err = fits.open(filename)['ERR'].data

                                fwhm_tbl = Table.read(f'{basepath}/reduction/fwhm_table.ecsv')
                                row = fwhm_tbl[fwhm_tbl['Filter'] == filtername]
                                fwhm = fwhm_arcsec = float(row['PSF FWHM (arcsec)'][0])
                                fwhm_pix = float(row['PSF FWHM (pixel)'][0])
                                ap_file = catalogfile.replace('.fits','-ap.fits')

                                create_apfile(img, catalogfile, ap_file, error=err, radius= 2.0*fwhm_pix, sky_in=3.0*fwhm_pix, sky_out=4.0*fwhm_pix, apcorr=1.0)
                                if not os.path.exists(ap_file):
                                    raise ValueError(f"AP file {ap_file} not found")
                                update_param('starbug.param', 'AP_FILE', f"{ap_file}")
                                update_param('starbug.param', 'SHARP_LO', '0.3')
                                update_param('starbug.param', 'SHARP_HI', '1.4')
                                update_param('starbug.param', 'SIGSRC', '4.0')
                                update_param('starbug.param', 'PROF_SCALE', '10')
                                update_param('starbug.param', 'BGD_CHECKFILE', f"{filename.replace('.fits','_starbug2-bgdcheck.fits')}")
                                output_file = f"{filename.replace('.fits','_starbug2.fits')}"
                                print('output_file', output_file)
                                update_param('starbug.param', 'OUTPUT', output_file)
                                update_param('starbug.param', 'FILTER', f"{filtername.upper()}")
                                update_param('starbug.param', 'HDUNAME', 'SCI')
                                cmd = ['starbug2', '--update-param']
                                result = subprocess.run(cmd, capture_output=True, text=True)
                    
                        
                                
                                print('make catalog for starbug2')
                                
                                cmd = ["starbug2", "-vd", f"{ap_file}", '-B', f"{filename}"]
                                result = subprocess.run(cmd, capture_output=True, text=True)  
                                print(result.stdout)
                                print(result.stderr)
                                if result.returncode != 0:
                                    print(f"Error running starbug2: {result.stderr}")
                                else:
                                    print(f"starbug2 output: {result.stdout}")

                                outfile = f"{filename.replace('.fits','_starbug2-bgd.fits')}"
                                if not os.path.exists(outfile):
                                    raise ValueError(f"Output file {outfile} not found")

    # step iii) using the same coordinates of the orignal catalog, do PSF photometry on the background subtracted images
if __name__ == "__main__":
    main()




