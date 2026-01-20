JWST W51 NIRCam/MIRI cataloging process

1. create_catalog.py : It creates catalogs for each individual frame. Starfinder can be chosen between Crowdsource and DAOFinder. For W51 data, DAOfinder was chosen. At the positions found by those starfinders, Photutil fits PSF function based on the model grid provided by STPSF.
   -> outputs: ..._daophot_basic.fits

2. refine_indv_frame_cat.py: It refines catalogs from each individual frame using some parameter cuts that were emprically obtained.
   -> outputs: ..._daophot_refined.fits

3. saturated_star_finding.py: Find saturated stars based on DQ array flags.
   -> outputs: ..._satstar_catalog_newnewnewnew.fits

4. add_saturated_star.py: Add saturated star catalog into the original catalog.
   -> outputs: ..._daophot_combined_with_satstars.fits

5. merging_catalog_nircam.py / merging_catalog_miri.py : merge the catalog from individual frame. also merge across filters.
   -> outputs: for the catalog where individual frame is merged: catalogs/..._indivexp_merged_...
             : for the catalog where all the filters are merged: catalogs/..._photometry_tables_merged_...
6. notebook/fix_miri_offset.ipynb
7. notebook/merge_nircam_miri.ipynb
   
