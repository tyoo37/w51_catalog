JWST W51 NIRCam/MIRI cataloging process

1. py/create_catalog.py : It creates catalogs for each individual frame. Starfinder can be chosen between Crowdsource and DAOFinder. For W51 data, DAOfinder was chosen. At the positions found by those starfinders, Photutil fits PSF function based on the model grid provided by STPSF.
   -> outputs: ..._daophot_basic.fits

~~2. refine_indv_frame_cat.py: It refines catalogs from each individual frame using some parameter cuts that were emprically obtained.
   -> outputs: ..._daophot_refined.fits~~ --> updated Mar 14, 2026

2. py/saturated_star_finding.py: Find saturated stars based on DQ array flags.
   -> outputs: ..._satstar_catalog_newnewnewnew.fits

3. py/add_saturated_star.py: Add saturated star catalog into the original catalog.
   -> outputs: ..._daophot_combined_with_satstars.fits
   
** added new refinement by nmatch ** (Mar 14, 2026)
4. notebook/filter_nmatch.ipynb: creates two catalogs with grade A where nmatch == max(nmatch) and grade B where nmatch >=max(nmatch)-1. Here, nmatch is the number of frames matched.
   

5. py/merging_catalog.py: merge the catalog from individual frame. ~~also merge across filters.~~
   -> outputs: for the catalog where individual frame is merged: catalogs/..._indivexp_merged_...
             ~~: for the catalog where all the filters are merged: catalogs/..._photometry_tables_merged_...~~
6. notebook/fix_miri_offset.ipynb: manually adjust the offset between MIRI and NIRCam catalog 
   
7. notebook/merge_nircam_miri.ipynb: merge nrca and nrcb for nircam and merge across filters. In the first stage of merging, F140M, F480M, F560W, F770W, F1000W, F1280W, F2100W catalogs are merged. In this process, any sources that are not spatially matched (0.1 arcsec) are considered as new sources. In the second round, other filter catalogs are just attached to the first merged catalog, now allowing to create new sources. 
   
