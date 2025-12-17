#!/bin/bash
mem=64gb
taskname=refine


mem=64gb

for filter in F140M F162M F187N F182M F210M ; do
    for module in nrca1 nrca2 nrca3 nrca4 nrcb1 nrcb2 nrcb3 nrcb4; do
        sbatch --job-name=starbug2-${filter} --output=starbug2-${filter}_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /blue/adamginsburg/t.yoo/from_red/w51/w51_catalog/completeness/run_starbug_for_source_removal.py --filternames=${filter} --target=w51 --modules=${module}"
    done
    
done

for filter in F335M F360M F405N F410M F480M ; do
    for module in nrcalong nrcblong; do
        sbatch --job-name=starbug2-${filter} --output=starbug2-${filter}_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /blue/adamginsburg/t.yoo/from_red/w51/w51_catalog/completeness/run_starbug_for_source_removal.py --filternames=${filter} --target=w51 --modules=${module}"
    done
    
done
for filter in F560W F770W F1000W F1280W F2100W ; do
    for module in mirimage; do
        sbatch --job-name=starbug2-${filter} --output=starbug2-${filter}_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /blue/adamginsburg/t.yoo/from_red/w51/w51_catalog/completeness/run_starbug_for_source_removal.py --filternames=${filter} --target=w51_miri --modules=${module}"
    done
    
done
