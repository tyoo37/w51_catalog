#!/bin/bash
mem=64gb
taskname=catalog


export STPSF_PATH=/blue/adamginsburg/t.yoo/from_red/stpsf-data

daoloop=("--daophot --skip-crowdsource" " ")
mem=64gb
#filter=F140M
#module=nrca1
#dao="--daophot --skip-crowdsource"
#sbatch --array=0-1 --job-name=webb-cat-${filter}-${module}-eachexp --output=webb-cat-${filter}-${module}-eachexp_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /home/t.yoo/w51/w51_nircam/catalog/crowdsource_catalogs_long.py --filternames=${filter} --modules=${module} --each-exposure ${dao}"
#or filter in F140M F162M F182M; do
#    for modnum in 1 2 3 4; do
#        for module in nrca${modnum} nrcb${modnum}; do
#            for dao in "${daoloop[@]}"; do
#                sbatch --array=0-23 --job-name=webb-cat-${filter}-${module}-eachexp --output=webb-cat-${filter}-${module}-eachexp_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /blue/adamginsburg/t.yoo/from_red/w51/w51_catalog/py/create_catalog.py --filternames=${filter} --modules=${module} --each-exposure ${dao}"
#            done
#        done
#    done
#done



#for filter in F187N F210M ; do
#    for modnum in 1 2 3 4; do
#        for module in nrca${modnum} nrcb${modnum}; do
#            for dao in "${daoloop[@]}"; do
#                sbatch --array=0-23 --job-name=webb-cat-${filter}-${module}-eachexp --output=webb-cat-${filter}-${module}-eachexp_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /blue/adamginsburg/t.yoo/from_red/w51/w51_catalog/py/create_catalog.py --filternames=${filter} --modules=${module} --each-exposure ${dao}"
#            done
#        done
#    done
#done


for filter in F405N; do
    for module in nrcalong nrcblong; do
        for dao in "${daoloop[@]}"; do
            sbatch --array=0-23 --job-name=webb-cat-${filter}-${module}-eachexp --output=webb-cat-${filter}-${module}-eachexp_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /blue/adamginsburg/t.yoo/from_red/w51/w51_catalog/py/create_catalog.py --filternames=${filter} --modules=${module} --each-exposure ${dao}"
        done
    done
done








#for filter in F140M F150W F162M F182M F187N F210M F335M F360M F405N F410M F480M; do   
#for filter in F140M F150W F210M F405N ; do
#    for module in nrca nrcb merged; do
#    for module in merged; do
#        if [ "${module}" == "merged" ]; then
#            sbatch --job-name=${taskname}-${filter}-${module}-eachexp --output=${taskname}-${filter}-${module}-eachexp_%j-%A.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /home/t.yoo/w51/w51_nircam/catalog/crowdsource_catalogs_long.py --filternames=${filter} --proposal_id=6151 --module=${module} --daophot " 
#        else
#            sbatch --job-name=${taskname}-${filter}-${module}-eachexp --output=${taskname}-${filter}-${module}-eachexp_%j-%A.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /home/t.yoo/w51/w51_nircam/catalog/crowdsource_catalogs_long.py --filternames=${filter} --proposal_id=6151 --module=${module} --each-exposure --daophot "
#        fi
           
#    done  
#done