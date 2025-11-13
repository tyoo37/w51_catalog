#!/bin/bash
mem=64gb
taskname=sat_star_finding


#daoloop=("--daophot --skip-crowdsource" " ")
#mem=32gb
#for filter in F140M; do
for filter in F140M F162M F182M F187N F210M F335M F360M F405N F410M F480M F560W F770W F1000W F1280W F2100W; do
#for filter in F1280W F2100W; do
#for filter in F560W F770W F1000W F1280W F2100W; do
    sbatch --job-name=sat_star_find-${filter} --output=sat_star_find-${filter}_%j-%A_%a.log  --account=astronomy-dept --qos=astronomy-dept-b --ntasks=2 --nodes=1 --mem=${mem} --time=96:00:00 --wrap "python /red/adamginsburg/t.yoo/w51/w51_catalog/py/add_saturated_stars.py --filter=${filter}"
done


