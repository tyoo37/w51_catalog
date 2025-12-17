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
from astropy.coordinates import SkyCoord, FK5
from regions import PixCoord
print("Done with imports", flush=True)
# step i )  load the catalog and do the quality assessment
# for the quality assessment, some tests are required to see the effect of masking
# load the catalog obtained from each exposure file



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


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--filternames", dest="filternames",
                    default='F140M',
                    help="filter name list", metavar="filternames")
    parser.add_option("--proposal_id", dest="proposal_id",
                    default='6151',
                    help="proposal_id", metavar="proposal_id")
    parser.add_option("--target", dest="target",
                    default='w51',
                    help="target", metavar="target")
    (options, args) = parser.parse_args()

    filternames = options.filternames.split(",")

    proposal_id = options.proposal_id
    target = options.target

    nvisits = {'2221': {'brick': 1, 'cloudc': 2},
                '1182': {'brick': 2},
                '6151': {'w51': 1, 'w51_miri': 2}
                }
    field_to_reg_mapping = {'2221': {'001': 'brick', '002': 'cloudc'},
                            '1182': {'004': 'brick'},
                            '6151': {'001': 'w51', '002':'w51_miri'}}[proposal_id]
    reg_to_field_mapping = {v:k for k,v in field_to_reg_mapping.items()}
    field = reg_to_field_mapping[target]

    basepath = f'/orange/adamginsburg/jwst/w51/'


    nircam_short_filters = ['F140M', 'F162M', 'F182M', 'F187N', 'F200W', 'F210M', 'F277W', ]
    nircam_long_filters = ['F300M', 'F335M', 'F356W', 'F360M','F405N', 'F410M', 'F430M', 'F444W', 'F460M', 'F480M']
    miri_filters = ['F560W', 'F770W', 'F1000W', 'F1130W', 'F1280W', 'F1500W', 'F1800W', 'F2100W', 'F2550W']



    for filtername in filternames:
        if filtername.upper() in nircam_short_filters:
            modules = ['nrca1', 'nrcb1', 'nrca2', 'nrcb2', 'nrca3', 'nrcb3', 'nrca4', 'nrcb4']
            instrument = 'NIRCam'
        elif filtername.upper() in nircam_long_filters:
            modules = ['nrcalong', 'nrcblong']
            instrument = 'NIRCam'
        elif filtername.upper() in miri_filters:
            modules = ['mirimage']
            instrument = 'MIRI'
        else:
            raise ValueError(f"Filter {filtername} not recognized as NIRCam or MIRI")
        
        for module in modules:
            
            for visitid in range(1, nvisits[proposal_id][target] + 1):
                visitid = f'{visitid:03d}'
                
                filenames = get_filenames(basepath, filtername, proposal_id,
                                            field, visitid=visitid,
                                            each_suffix='cal',
                                            module=module, pupil='clear')
                for i, filename in enumerate(filenames):
                    if True:
                        exposurenumber = int(filename.split("_")[2])
                        exposure_id = filename.split("_")[2]
                        visit_id = filename.split("_")[0][-3:]
                        vgroup_id = filename.split("_")[1]
                        exposure_ = f'_exp{exposurenumber:05d}'
                        visitid_ = f'_visit{int(visitid):03d}' if visitid is not None else ''
                        if instrument == 'NIRCam':
                            vgroupid_ = f'_vgroup{int(vgroup_id):05d}' if vgroup_id is not None else ''
                        elif instrument == 'MIRI':
                            vgroupid_ = f'_vgroup{vgroup_id}' if vgroup_id is not None else ''
                            module2 = 'mirimage'
                        detector = filename.split("_")[-2]
                        print(detector, flush=True)
                        #f360m_nrcb_visit001_vgroup3105_exp00008_daophot_basic.fits
                        #/orange/adamginsburg/jwst/w51/F360M/f360m_nrcalong_visit001_vgroup3105_exp00005_daophot_basic.fits
                    
                        wav = int(filtername[1:4])
                        if instrument == 'NIRCam' and wav < 250:
                            catalogfile = f"{basepath}/{filtername}/{filtername.lower()}_{detector}{visitid_}{vgroupid_}{exposure_}_daophot_basic.fits"
                        elif instrument == 'NIRCam' and wav >= 250:
                            catalogfile = f"{basepath}/{filtername}/{filtername.lower()}_{module}{visitid_}{vgroupid_}{exposure_}_daophot_basic.fits"
                        else:
                            catalogfile = f"{basepath}/{filtername}/{filtername.lower()}_{module2}{visitid_}{vgroupid_}{exposure_}_daophot_basic.fits"
                        print('catalogfile', catalogfile)
                        #f560w_mirimage_visit002_vgroup02101_exp00001_daophot_iterative.fits
                        cat = Table.read(catalogfile)
                        roundness1 = cat['roundness1']
                        roundness2 = cat['roundness2']
                        sharpness = cat['sharpness']
                        flux = cat['flux_fit']
                        fluxerr = cat['flux_err']
                        snr = flux/fluxerr
                        qfit = cat['qfit']
                        cfit = cat['cfit']
                        #skycoord = cat['skycoord']
                        #n_match_good = cat['nmatch_good']
                        #std_ra = cat['std_ra']
                        #std_dec = cat['std_dec']

                        if instrument == 'NIRCam':              
                            good_sources = ((roundness1 < 0.8) & (roundness1 > -0.9) & (roundness2 < 0.6) & (roundness2 > -0.6) & (sharpness < 1.2) & (sharpness>0.25) 
                            & (snr > 3) & (qfit < 0.33) & (cfit < 0.2) & (cfit > -0.2) & ~((snr < 20) & (flux > 50)))
                            
                        elif instrument == 'MIRI':
                            good_sources = ((roundness1 < 0.7) & (roundness1 > -0.7) & (roundness2 < 0.5) & (roundness2 > -0.5) & (sharpness < 1.1) & (sharpness>0.4) 
                            & (snr > 3) & (qfit < 0.6) & (cfit < 0.06) & (cfit > -0.06))
                        
                        refined_cat = cat[good_sources]
                        print(f"Filter {filtername} Module {module} Visit {visitid} Exposure {exposure_id}: {len(cat)} sources, {len(refined_cat)} after quality cuts", flush=True)
                        outcatfile = catalogfile.replace('_daophot_basic.fits', '_daophot_refined.fits')
                        refined_cat.write(outcatfile, overwrite=True)

if __name__ == "__main__":
    main()