#!/bin/bash

echo "Find genes using prodigal:"

cd ~/Projects/AD
    
mkdir Annotate

cd Annotate/

python $DESMAN/scripts/LengthFilter.py ../Assembly/final_contigs_c10K.fa -m 1000 >     final_contigs_gt1000_c10K.fa

prodigal -i final_contigs_gt1000_c10K.fa -a final_contigs_gt1000_c10K.faa -d     final_contigs_gt1000_c10K.fna  -f gff -p meta -o final_contigs_gt1000_c10K.gff 

echo "Assign COGs"


export COGSDB_DIR=~/Databases/rpsblast_cog_db

$CONCOCT/scripts/RPSBLAST.sh -f final_contigs_gt1000_c10K.faa -p -c 12 -r 1

cd ../Concoct
python $CONCOCT/scripts/COG_table.py -b ../Annotate/final_contigs_gt1000_c10K.out  -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c clustering_gt1000.csv  --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > clustering_gt1000_scg.tsv


sed '1d' clustering_gt1000.csv > clustering_gt1000_R.csv
cut -f1,3- < clustering_gt1000_scg.tsv | tr "\t" "," > clustering_gt1000_scg.csv
$CONCOCT/scripts/Sort.pl < clustering_gt1000_scg.csv > clustering_gt1000_scg_sort.csv

concoct_refine clustering_gt1000_R.csv original_data_gt1000.csv clustering_gt1000_scg_sort.csv > concoct_ref.out

python $CONCOCT/scripts/COG_table.py -b ../Annotate/final_contigs_gt1000_c10K.out  -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c clustering_refine.csv  --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > clustering_refine_scg.tsv

