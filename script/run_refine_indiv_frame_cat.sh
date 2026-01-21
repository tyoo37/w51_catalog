#!/bin/bash
mem=64gb
taskname=refine


mem=64gb
#filter=F140M
#module=nrca1
#dao="--daophot --skip-crowdsource"
#sbatch --array=0-1 --job-name=webb-cat-${filter}-${module}-eachexp --output=webb-cat-${filter}-${module}-eachexp_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /home/t.yoo/w51/w51_nircam/catalog/crowdsource_catalogs_long.py --filternames=${filter} --modules=${module} --each-exposure ${dao}"


#for filter in F560W; do
for filter in F140M F162M F182M F187N F210M F335M F360M F410M F405N F480M; do
    sbatch --job-name=refine-${filter} --output=refine-${filter}_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /blue/adamginsburg/t.yoo/from_red/w51/w51_catalog/py/refine_indv_frame_cat.py --filternames=${filter} --target=w51"
done

#for filter in F335M F360M F410M F405N F480M ; do
    
#    for module in nrcalong nrcblong; do
        
#        sbatch --array=0-32 --job-name=refine-${filter}-${module}-eachexp --output=webb-cat-${filter}-${module}-eachexp_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /blue/adamginsburg/t.yoo/from_red/w51/w51_catalog/py/refine_indv_frame_cat.py --filternames=${filter}  --target=w51"
        
        
#    done
#done

for filter in F560W F770W F1000W F1280W F2100W ; do
    sbatch --job-name=refine-${filter} --output=refine-${filter}_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /blue/adamginsburg/t.yoo/from_red/w51/w51_catalog/py/refine_indv_frame_cat.py --filternames=${filter} --target=w51_miri"
done
