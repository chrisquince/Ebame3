# DESMAN TARA Tutorial

We are going to locate ecotypes in one of the TARA MAGs, a SAR11 relative TARA_PSW_MAG_00074.
We begin by downloading the Ebame3 repo onto our VMs which has some useful scripts.
Log into the computers:

```
ssh ubuntu@myVM
```

```
cd ~/repos
git clone https://github.com/chrisquince/Ebame3.git
```

Create a directory to work in:

```
cd ~/Projects
mkdir DESMANTutorial
cd DESMANTutorial
```

Download the pre-computed frequency file on the single copy core genes. 

```
wget https://desmantutorial.s3.climb.ac.uk/TARA_PSW_MAG_00074_scg.freq
```

How many core COGs did this MAG contain?

## Variant selection

Then we run the variant filter to identify variant positions. We also filter for non-core cogs '-c' and set a sample coverage threshold of '-m 1.0':

```
python $DESMAN/desman/Variant_Filter.py TARA_PSW_MAG_00074_scg.freq -p -o TARA_PSW_MAG_00074_scg -f 25.0 -c -sf 0.80 -t 2.5 -m 1.0
```

The '-p' option uses the slower but more precision minor variant frequency optimisation.

How many COGs survived the filtering? How many variant positions are there?

## Haplotype inference

Then we run the haplotype inference on the variant positions that are selected.

```
stub=TARA_PSW_MAG_00074_scg
varFile=${stub}sel_var.csv

eFile=${stub}tran_df.csv

for g in 1 2 3 4  
do
for r in 0 1 2 3 4
    do
        echo $g:
        desman $varFile -e $eFile -o ${stub}_${g}_${r} -g $g -s $r -m 1.0 -i 100 > ${stub}_${g}_${r}.out&
    done
done
```

Create a deviance plot for this MAG:

```
cat */fit.txt | cut -d"," -f2- > Dev.csv
sed -i '1iH,G,LP,Dev' Dev.csv 
```

Which we can then visualise:

```
$DESMAN/scripts/PlotDev.R -l Dev.csv -o Dev.pdf
```

![Posterior deviance](../Figures/Dev.pdf)

Now we need the TARA sample meta data:

```
wget https://desmantutorial.s3.climb.ac.uk/TARA-samples.csv
```

An