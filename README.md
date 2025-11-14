JWST W51 NIRCam/MIRI cataloging process

1. create_catalog.py : It creates catalogs for each individual frame. Starfinder can be chosen between Crowdsource and DAOFinder. For W51 data, DAOfinder was chosen. At the positions found by those starfinders, Photutil fits PSF function based on the model grid provided by STPSF.

2. refine_indv_frame_cat.py: It refines catalogs from each individual frame using some parameter cuts that were emprically obtained.

3. saturated_star_finding.py: Find saturated stars based on DQ array flags.

4. add_saturated_star.py: Add saturated star catalog into the original catalog.

5. merging_catalog_nircam.py / merging_catalog_miri.py : merge the catalog from individual frame. also merge across filters.
