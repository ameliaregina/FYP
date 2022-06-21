#!/bin/bash

# module load anaconda/3
# conda activate py39

cd ~/fyp

mkdir -p temp
python genes_sep.py 'GI_CN.txt' 'compiled_somatic_mutation_170522.txt' 'TP53' 'METABRIC' 'TCGA'

mkdir -p temp2
cp temp/* temp2

mkdir -p res
for FILE in temp/*; 
do Rscript dotplot2.r $FILE 'res2'; 
done

python boxplot.py 'GI_CN.txt' 'compiled_somatic_mutation_170522.txt' 'TP53' 'GID' 'TCGA'

rm -r temp