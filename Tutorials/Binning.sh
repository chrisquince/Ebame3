#!/bin/bash


cd ~/Projects/AD

ls ReadsSub/*R1.fastq | tr "\n" "," | sed 's/,$//' > R1.csv
ls ReadsSub/*R2.fastq | tr "\n" "," | sed 's/,$//' > R2.csv



megahit -1 $(<R1.csv) -2 $(<R2.csv) -t 12 -o Assembly > megahit.out

contig-stats.pl < Assembly/final.contigs.fa


cd Assembly
bwa index final_contigs_c10K.fa
cd ..


echo "Start mapping"

mkdir Map

for file in ./ReadsSub/*R1.fastq
do 
   
   stub=${file%_R1.fastq}
   name=${stub##*/}
   
   echo $stub

   file2=${stub}_R2.fastq

   bwa mem -t 12 Assembly/final_contigs_c10K.fa $file $file2 > Map/${name}.sam
done

echo "Calculate coverages"

python $DESMAN/scripts/Lengths.py -i Assembly/final_contigs_c10K.fa > Assembly/Lengths.txt

for file in Map/*.sam
do
    stub=${file%.sam}
    stub2=${stub#Map\/}
    echo $stub  
    samtools view -h -b -S $file > ${stub}.bam
    samtools view -b -F 4 ${stub}.bam > ${stub}.mapped.bam
    samtools sort -m 1000000000 ${stub}.mapped.bam -o ${stub}.mapped.sorted.bam
    bedtools genomecov -ibam ${stub}.mapped.sorted.bam -g Assembly/Lengths.txt > ${stub}_cov.txt
done

echo "Delete sam files"

rm Map/*sam

echo "Collate coverages together"

```
for i in Map/*_cov.txt 
do 
   echo $i
   stub=${i%_cov.txt}
   stub=${stub#Map\/}
   echo $stub
   awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' $i > Map/${stub}_cov.csv
done

$DESMAN/scripts/Collate.pl Map > Coverage.csv

echo "Now run concoct"

mkdir Concoct

mv Coverage.csv Concoct

cd Concoct

tr "," "\t" < Coverage.csv > Coverage.tsv

concoct --coverage_file Coverage.tsv --composition_file ../Assembly/final_contigs_c10K.fa -t 8 
