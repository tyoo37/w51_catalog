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

print("Done with imports", flush=True)
# step i )  load the catalog and do the quality assessment
# for the quality assessment, some tests are required to see the effect of masking
# load the catalog obtained from each exposure file

nvisits = {'2221': {'brick': 1, 'cloudc': 2},
               '1182': {'brick': 2},
               '6151': {'w51': 1  }
               }
def update_param(filepath, param_name, new_value):
    lines = []
    with open(filepath, "r") as f:
        for line in f:
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
    parser.add_option("--modules", dest="modules",
                    default='miri',
                    help="modules", metavar="modules")
                
    (options, args) = parser.parse_args()

    filternames = options.filternames.split(",")

    proposal_id = options.proposal_id
    target = options.target
#modules = ['miri']
#filternames = ['F560W', 'F770W', 'F1000W', 'F1280W',  'F2100W']
#proposal_id = '6151'
#target = 'w51'
#modules = ['miri']

    field_to_reg_mapping = {'2221': {'001': 'brick', '002': 'cloudc'},
                                '1182': {'004': 'brick'},
                                '6151': {'001': 'w51'}}[proposal_id]
    reg_to_field_mapping = {v:k for k,v in field_to_reg_mapping.items()}
    field = reg_to_field_mapping[target]

    basepath = f'/orange/adamginsburg/jwst/w51/'

    cmd = ["starbug2", "--local-param" ]

    result = subprocess.run(cmd, capture_output=True, text=True)  


    if True:
        for module in modules:
            detector = module # no sub-detectors for long-NIRCAM
            for filtername in filternames:
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
                                #catalogfile = f"{basepath}/{filtername}/{filtername.lower()}_{module}{detector}{visitid_}{vgroupid_}{exposure_}_daophot_basic.fits"

                                #catalog = fits.open(catalogfile)
                                # change the format of the catalog accordingly to be used in starbug2

                                # step ii) using refined catalogs, run starbugs2 to remove background for each exposure file
                                print('create starbug2 param file')
                                update_param('starbug.param', 'SHARP_LO', '0.3')
                                update_param('starbug.param', 'SHARP_HI', '1.4')
                                update_param('starbug.param', 'SIGSRC', '4.0')
                                update_param('starbug.param', 'PROF_SCALE', '5')
                                update_param('starbug.param', 'BGD_CHECKFILE', f"{filename.replace('.fits','_starbug2-bgdcheck.fits')}")
                                output_file = f"{filename.replace('.fits','_starbug2.fits')}"
                                print(output_file)
                                update_param('starbug.param', 'OUTPUT', output_file)
                                update_param('starbug.param', 'FILTER', f"{filtername.upper()}")
                                ap_file = output_file.split('.')[0]+'-ap.fits'
                                update_param('starbug.param', 'AP_FILE', f"{ap_file}")
                        

                                print('make catalog for starbug2')
                                cmd = ["starbug2", "-vD", f"{filename}" ]

                                result = subprocess.run(cmd, capture_output=True, text=True)  
                                print(result.stdout)
                                print(result.stderr)
                                if result.returncode != 0:
                                    print(f"Error running starbug2: {result.stderr}")
                                else:
                                    print(f"starbug2 output: {result.stdout}")

                                print('make bkg image')
                                cmd = ["starbug2", "-vd", f"{ap_file}", '-B', f"{filename}"]
                                result = subprocess.run(cmd, capture_output=True, text=True)  
                                print(result.stdout)
                                print(result.stderr)
                                if result.returncode != 0:
                                    print(f"Error running starbug2: {result.stderr}")
                                else:
                                    print(f"starbug2 output: {result.stdout}")

    # step iii) using the same coordinates of the orignal catalog, do PSF photometry on the background subtracted images


