#!/bin/bash


#Then we calculate coverage of each cluster/MAG in each sample.

cd ~/Projects/AD/Concoct
sed '1d' clustering_refine.csv > clustering_refineR.csv
python $DESMAN/scripts/ClusterMeanCov.py Coverage.csv clustering_refineR.csv ../Assembly/final_contigs_c10K.fa > clustering_refine_cov.csv
sed 's/Map\///g' clustering_refine_cov.csv > clustering_refine_covR.csv

cd ~/Projects/AD/Annotate
python $DESMAN/scripts/ExtractCogs.py -b final_contigs_gt1000_c10K.out --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv -g final_contigs_gt1000_c10K.gff > final_contigs_gt1000_c10K.cogs

python $DESMAN/scripts/ExtractGenes.py -g final_contigs_gt1000_c10K.gff > final_contigs_gt1000_c10K.genes
cd ..

mkdir Split
cd Split
$DESMAN/scripts/SplitClusters.pl ../Annotate/final_contigs_gt1000_c10K.fa ../Concoct/clustering_refine.csv
SplitCOGs.pl ../Annotate/final_contigs_gt1000_c10K.cogs ../Concoct/clustering_refine.csv
SplitGenes.pl ../Annotate/final_contigs_gt1000_c10K.genes ../Concoct/clustering_refine.csv
SplitFaa.pl ../Annotate/final_contigs_gt1000_c10K.faa ../Concoct/clustering_refine.csv
cd ..


#Kegg ortholog assignment on genes:

cd ~/Projects/AD/Split
python ~/bin/CompleteClusters.py ../Concoct/clustering_refine_scg.tsv > Cluster75.txt


while read line
do 
    file=${line}/${line}.faa
    stub=${file%.faa}
    base=${stub##*/}
    echo $base

    diamond blastp -d $KEGG_DB/genes/fasta/genes.dmnd -q $file -p 12 -o ${stub}.m8
done < Cluster75.txt



#The above maps onto Kegg genes these are then mapped to kegg orthologs by the Perl script:

COUNT=0
for file in Cluster*/*m8
do
	dir=${file%.m8}
	echo $file
	echo $dir
     Assign_KO.pl < $file > ${dir}.hits&
    let COUNT=COUNT+1

    if [ $COUNT -eq 8 ]; then
        wait;
        COUNT=0
    fi
done

#Discussion point, trivial parallelisation using bash.

~/repos/MAGAnalysis/scripts/CollateHits.pl > CollateHits75.csv

python ~/repos/MAGAnalysis/scripts/KO2MODULEclusters2.py -i CollateHits75.csv -o Collate_modules.csv 

awk -F"," 'NR>1{print $1}' Collate_modules.csv | xargs -I {} curl -s http://rest.kegg.jp/find/module/{} > ModNames.txt

export NR_DMD=$HOME/Databases/NR/nr.dmnd



cd ~/Projects/AD/
mkdir AssignTaxa
cd AssignTaxa
cp ../Annotate/final_contigs_gt1000_c10K.faa .



diamond blastp -p 8 -d $NR_DMD -q final_contigs_gt1000_c10K.faa -o final_contigs_gt1000_c10K_nr.m8 > d.out


python $DESMAN/scripts/Lengths.py -i final_contigs_gt1000_c10K.faa > final_contigs_gt1000_c10K.len
python $DESMAN/scripts/ClassifyContigNR.py final_contigs_gt1000_c10K_nr.m8 final_contigs_gt1000_c10K.len -o final_contigs_gt1000_c10K_nr -l /home/ubuntu/Databases/NR/all_taxa_lineage_notnone.tsv -g /home/ubuntu/Databases/NR/gi_taxid_prot.dmp


$DESMAN/scripts/Filter.pl 8 < final_contigs_gt1000_c10K_nr_contigs.csv | grep -v "_6" | grep -v "None" > final_contigs_gt1000_c10K_nr_species.csv

$CONCOCT/scripts/Validate.pl --cfile=../Concoct/clustering_refine.csv --sfile=final_contigs_gt1000_c10K_nr_species.csv --ffile=../Assembly/final_contigs_c10K.fa --ofile=Taxa_Conf.csv

$CONCOCT/scripts/ConfPlot.R -c Taxa_Conf.csv -o Taxa_Conf.pdf

cd ~/Projects/AD/Split
cp ~/repos/MAGAnalysis/cogs.txt .
mkdir SCGs

while read line
do
    cog=$line
    echo $cog
     ~/repos/MAGAnalysis/phyloscripts/SelectCogsSCG.pl ../Concoct/clustering_refine_scg.tsv ../Annotate/final_contigs_gt1000_c10K.fna $cog > SCGs/$cog.ffn
done < cogs.txt



mkdir AlignAll

while read line
do
    cog=$line
    echo $cog
    cat ~/Databases/NCBI/Cogs/All_$cog.ffn SCGs/${cog}.ffn > AlignAll/${cog}_all.ffn
    mafft --thread 12 AlignAll/${cog}_all.ffn > AlignAll/${cog}_all.gffn
done < cogs.txt

for file in  AlignAll/*gffn
do
    echo $stub
    stub=${file%.gffn}
    trimal -in $file -out ${stub}_al.gfa -gt 0.9 -cons 60
done


#The next script requires the IDs of any cluster or taxa that may appear in fasta files, therefore:


cat AlignAll/*gffn | grep ">" | sed 's/_COG.*//' | sort | uniq | sed 's/>//g' > Names.txt


~/repos/MAGAnalysis/phyloscripts/CombineGenes.pl Names.txt AlignAll/COG0*_al.gfa > AlignAll.gfa

~/repos/MAGAnalysis/phyloscripts/MapTI.pl /home/ubuntu/repos/MAGAnalysis/data/TaxaSpeciesR.txt < AlignAll.gfa > AlignAllR.gfa

fasttreeMP -nt -gtr < AlignAllR.gfa 2> SelectR.out > AlignAllR.tree
