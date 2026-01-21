#!/bin/bash

for module in nrca nrcb; do
    #sbatch --array=0-20 --job-name=webb-cat-${module}-singlefields-dao --output=webb-cat-${module}-singlefields-dao_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=1 --nodes=1 --mem=64gb --time=96:00:00 --wrap "python /home/t.yoo/w51/w51_nircam/catalog/merge_catalog.py --merge-singlefields --modules=${module} --indiv-merge-methods=dao --skip-crowdsource"
    sbatch --array=0-31 --job-name=merge-${module}-singlefields-dao --output=merge-${module}-singlefields-dao_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=1 --nodes=1 --mem=64gb --time=96:00:00 --wrap "python /blue/adamginsburg/t.yoo/from_red/w51/w51_catalog/py/merging_catalog_nircam.py --merge-singlefields --modules=${module} --indiv-merge-methods=dao_after_merger --skip-crowdsource"
   # sbatch --array=0-20 --job-name=webb-cat-${module}-singlefields-crowdsource --output=webb-cat-${module}-singlefields-crowdsource_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=1 --nodes=1 --mem=64gb --time=96:00:00 --wrap "python /home/t.yoo/w51/w51_nircam/catalog/merge_catalog.py --merge-singlefields --modules=${module} --indiv-merge-methods=crowdsource --skip-dao"
   # sbatch --array=0-20 --job-name=webb-cat-${module}-singlefields-iterative --output=webb-cat-${module}-singlefields-iterative_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=1 --nodes=1 --mem=64gb --time=96:00:00 --wrap "python /home/t.yoo/w51/w51_nircam/catalog/merge_catalog.py --merge-singlefields --modules=${module} --indiv-merge-methods=iterative --skip-crowdsource"
done

