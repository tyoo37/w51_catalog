#!/bin/bash
mem=64gb
taskname=syn_stars

export STPSF_PATH=/blue/adamginsburg/t.yoo/from_red/stpsf-data

echo $STPSF_PATH
mem=64gb
#filter=F140M
#module=nrca1
#dao="--daophot --skip-crowdsource"
#sbatch --array=0-1 --job-name=webb-cat-${filter}-${module}-eachexp --output=webb-cat-${filter}-${module}-eachexp_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /home/t.yoo/w51/w51_nircam/catalog/crowdsource_catalogs_long.py --filternames=${filter} --modules=${module} --each-exposure ${dao}"


#for filter in F560W; do
for temp in 3000 4000 5000 6000 ; do
     
    sbatch --array=0-30 --job-name=${taskname}_${temp} --output=${taskname}_${temp}_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /blue/adamginsburg/t.yoo/from_red/w51/w51_catalog/syn_star/inject_synthetic_star.py --temp=${temp}"     
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