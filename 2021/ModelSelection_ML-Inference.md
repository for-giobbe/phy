### INTRO:

Currently, two divergent systematic methods are commonly applied for inferring species trees: the supermatrix approach (concatenation) and the coalescent approach (gene trees are calculated and then reconciled in a species tree). You will find a couple of interesting papers on the topic at the end of this tutorial, but here we will focus on the supermatrix approach.

Before starting you should have all aligned & filtered MSA, or if you missed the previous lesson you can grab them from [here](https://github.com/for-giobbe/phy/blob/master/2021/Analyses/Precomputed_1-to-1_Alignments.tar.gz)

## CONCATENATION OF MSA

After having aligned and filtere our genes, we will concatenate them using AMAS, a nice Python toolkit. However, ss a preliminar step we need to change a bit our filtered alignments. First of all we want to convert the gblock format into standard fasta. Then we need to change all fasta headers leaving only the species name. This step is mandatory beacuse concatenation is **always** made looking at the headers.

```
for i in Analyses/1-to-1_Alignments/gblock/*gb; do varOG=$( echo "$i" | cut -d"." -f1); cat "$i" | tr -d " " | sed 's/_1//g' > Analyses/1-to-1_Alignments/gblock/"$varOG".t_cofee-gb.fa; rm Analyses/1-to-1_Alignments/gblock/"$i"; done

for i in Analyses/1-to-1_Alignments/gblock/*.fa; do sed -i.old 's/|.*$//g' "$i"; rm Analyses/1-to-1_Alignments/gblock/*.old; done ## "-i" and "rm" is necessary only for macOS users
```

Now let's have a look at the AMAS help:

```
AMAS.py --help
usage: AMAS <command> [<args>]

The AMAS commands are:
  concat      Concatenate input alignments
  convert     Convert to other file format
  replicate   Create replicate data sets for phylogenetic jackknife
  split       Split alignment according to a partitions file
  summary     Write alignment summary
  remove      Remove taxa from alignment
  translate   Translate DNA alignment into protein alignment
  trim        Remove columns from alignment

Use AMAS <command> -h for help with arguments of the command of interest

positional arguments:
  command     Subcommand to run

optional arguments:
  -h, --help  show this help message and exit
```

