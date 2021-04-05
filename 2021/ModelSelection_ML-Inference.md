### INTRO:

Currently, two divergent systematic methods are commonly applied for inferring species trees: the supermatrix approach (concatenation) and the coalescent approach (gene trees are calculated and then reconciled in a species tree). You will find a couple of interesting papers on the topic at the end of this tutorial, but here we will focus on the supermatrix approach.

Before starting you should have all aligned & filtered MSA, or if you missed the previous lesson you can grab them from [here](https://github.com/for-giobbe/phy/blob/master/2021/Analyses/Precomputed_1-to-1_Alignments.tar.gz)

## CONCATENATION OF MSAs

After having aligned and filtere our genes, we will concatenate them using AMAS, a nice Python toolkit. However, as preliminar steps we need to: 
  * Convert the gblock format into standard fasta. 
  * Change all fasta headers leaving only the species name. This step is mandatory beacuse concatenation is **always** made looking at the headers.
  * Create a new working directory

```
makdir mkdir Analyses/IQ-TREE
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

and specifically, to the ```concat``` function:

```
AMAS.py concat --help
usage: AMAS.py [-h] [-p CONCAT_PART] [-t CONCAT_OUT] [-u {fasta,phylip,nexus,phylip-int,nexus-int}] [-y {nexus,raxml,unspecified}] [-e] [-c CORES] -i IN_FILES [IN_FILES ...] -f
               {fasta,phylip,nexus,phylip-int,nexus-int} -d {aa,dna}

Concatenate input alignments

optional arguments:
  -h, --help            show this help message and exit
  -p CONCAT_PART, --concat-part CONCAT_PART
                        File name for the concatenated alignment partitions. Default: 'partitions.txt'
  -t CONCAT_OUT, --concat-out CONCAT_OUT
                        File name for the concatenated alignment. Default: 'concatenated.out'
  -u {fasta,phylip,nexus,phylip-int,nexus-int}, --out-format {fasta,phylip,nexus,phylip-int,nexus-int}
                        File format for the output alignment. Default: fasta
  -y {nexus,raxml,unspecified}, --part-format {nexus,raxml,unspecified}
                        Format of the partitions file. Default: 'unspecified'
  -e, --check-align     Check if input sequences are aligned. Default: no check
  -c CORES, --cores CORES
                        Number of cores used. Default: 1

required arguments:
  -i IN_FILES [IN_FILES ...], --in-files IN_FILES [IN_FILES ...]
                        Alignment files to be taken as input. You can specify multiple files using wildcards (e.g. --in-files *fasta)
  -f {fasta,phylip,nexus,phylip-int,nexus-int}, --in-format {fasta,phylip,nexus,phylip-int,nexus-int}
                        The format of input alignment
  -d {aa,dna}, --data-type {aa,dna}
                        Type of data
```

Now we can concatenate our alignments:

```
AMAS.py concat -i Analyses/1-to-1_Alignments/gblock/*.fa -f fasta -d aa --part-format nexus
mv concatenated.out Analyses/IQ-TREE/My_Concat.fa
mv partitions.txt Analyses/IQ-TREE/
```

As a results you should have two files, a merged MSA and a partition file, where are stored the coordinate of the gene boundaries. Let's have a look at it:

```
head -n 20 partitions.txt 
#NEXUS

Begin sets;
	charset p1_OG0000030 = 1-197;
	charset p2_OG0000031 = 198-348;
	charset p3_OG0000032 = 349-1070;
	charset p4_OG0000033 = 1071-1512;
	charset p5_OG0000034 = 1513-1925;
	charset p6_OG0000035 = 1926-2300;
	charset p7_OG0000036 = 2301-2463;
	charset p8_OG0000037 = 2464-2867;
	charset p9_OG0000038 = 2868-3281;
	charset p10_OG0000039 = 3282-3703;
	charset p11_OG0000040 = 3704-4054;
	charset p12_OG0000041 = 4055-5124;
	charset p13_OG0000042 = 5125-5285;
	charset p14_OG0000043 = 5286-5482;
	charset p15_OG0000044 = 5483-6214;
	charset p16_OG0000045 = 6215-7227;
	charset p17_OG0000046 = 7228-7383;
```
>NOTE: We are working with proteins so our partitions file just store the starting and the ending position of each gene. However, if we were using nucleotide sequences of PCGs, we might also want the coordinates of all the first, second and third codon positions. Indeed, they evolve differently due to the gen code degeneracy and codon usage bias.

Now we are ready to carry out our model selection and tree inference in a Maximum Likelihood framework...

## MODEL SELECTION AND MAXIMUM LIKELIHOOD INFERENCE
