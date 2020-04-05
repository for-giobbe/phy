# Concatenation & Model Selection

Currently, two divergent systematic methods are commonly applied for inferring species trees: the supermatrix approach (concatenation) and the coalescent approach (gene trees are calculated and then reconciled in a species tree). 
You will find an interesting paper on the topic at the end of this tutorial, but here we will focus on the supermatrix approach.


---


## concatenation: 

After having aligned our genes, we will concatenate them using Phyutility: the software has a wide array of functions and can be considere a swiss-knife for phylogeneticists.
As you always should, when using a new software, take a look at its [manual](https://github.com/blackrim/phyutility/blob/master/manual.pdf). 

After adjusting the path to the java executable of phyutility, try this string:

```
java -jar /Applications/bio/phyutility/phyutility.jar -concat -in *.fasta -out concatenation.nxs
```

Then take a look at the output file ```concatenation.nxs```. It's a nexus file which 

* stores sequence information 
* from which we will extract the block which codes informations for gene boundaries. 

We can reformat the nexus line which stores the gene boundaries information,

```
[12S_total.fasta_gene1 1-816 CO1_total.fasta_gene2 817-2353 ND2_total.fasta_gene3 2354-3380 ]
```

into something similar to:

```
DNA, 12S = 1-816
DNA, CO1 = 817-2353
DNA, ND2 = 2354-3380
```

and save it to a file named ```gene.prt``` using an editor as ```nano```. I strongly discourage manual editing of files, but for us it's a nice way to understand the structure of different formats. Btw if you are curious you can also get familiar with the nexus [format](http://informatics.nescent.org/wiki/NEXUS_Specification).
This is a step in phylogenetic pipelines which is often automated, as manual editing of partition files for hundred/thousands of genes is not possible. Later you can explore this in-house [script](https://github.com/for-giobbe/phy/blob/master/scripts/concatenate_partitions.sh) which I made for the purpose,
just remember to adjust the path to phyutility in it. Moreover future versions of IQ-Tree can accept a folder of alignments as input, completely removing the need for the user to concatenate and edit partitions.  

Nonetheless this step is necessary for you to understand its underlying logic! For example, let's edit the partition file to take into account the different codon position (they evolve under different constrains due to the gen code degeneracy).

Let's use ```nano``` to transform our ```gene.prt``` file from:

```
DNA, 12S = 1-816
DNA, CO1 = 817-2353
DNA, ND2 = 2354-3380
```

into:

```
DNA, 12S = 1-816
DNA, CO1st = 817-2353\3
DNA, CO1nd = 818-2353\3
DNA, CO1rd = 819-2353\3
DNA, ND2st = 2354-3380\3
DNA, ND2nd = 2355-3380\3
DNA, ND2rd = 2356-3380\3
```

and save it as ```codon.prt```. As you can notice the ```/3``` notation informs the program to consider each three positions, and the start of the partition needs to be adjusted as well.


At the end of this part we should have 

* ```concatenation.nxs``` - a necus file which contains the concatenation of our loci and
* ```gene_and_codon.prt```: an ["a priori" partitioning scheme](https://github.com/for-giobbe/phy/blob/master/examples/codon.prt), based on a priori biological information.


---


## model of evolution & partitioning scheme selection: 

To carry out the model selection we will use ModelFinder [Kalyaanamoorthy et al., 2017](https://www.nature.com/articles/nmeth.4285).
As for most of the steps which are carried out to build a tree, there is a staggering choice of different tools and is sometimes difficult to decide which to use.
The more widespread tool for model selection is [PartitionFinder2](http://www.robertlanfear.com/partitionfinder/) but we are using ModelFinder as implemented in IQ-Tree for two key reason:

* usability
* speed

Also I think it's a good idea to keep constantly using new and shiny tools. Let's try the string:

```iqtree -s CO1_total.fasta -m MF```

We can read the best model from the standard output or open the relative file by ```zcat CO1_total.fasta.model```.

The MF word stands for ModelFinder, which tells IQ-TREE to perform ModelFinder:
this tool computes the log-likelihoods of an initial parsimony tree for many different models and the Akaike information criterion (AIC), 
corrected Akaike information criterion (AICc), and the Bayesian information criterion (BIC). 
Then ModelFinder chooses the model that minimizes the BIC score (you can also change to AIC or AICc by adding the option -AIC or -AICc, respectively).

The -m flag can also specify a model name to use during the analysis, which can be a priori specified by the user (here's a [list](http://www.iqtree.org/doc/Substitution-Models) of models implemented in ModelFinder if you feel you will nail it better than ModelFinder).

```iqtree -s example.phy -m HKY+I+G```

As you see several additional parameters are possible in order to:

* describe base frequencies, such as ```+F	```, ```+FQ```, ```+FO```
* model rate heterogeneity across sites: ```+I```, ```+G``` and ```+R```

What we've seen until now is the process through which we select the "best" model of evolution for our sequence data, according to a metric of choice.
In a concatenation framework we should carry out the process on the whole concatenation instead of single alignements, but without loosing the information of the single genes boundaries. Let's try:

```iqtree -s concatenation.nxs -sp gene_and_codon.prt -m MFP```

and the we can take a look at the file ```codon.prt.best_scheme.nex```

The previous analysis will result in separate models for each partion. Nonetheless, there are several reasons for which we wanto to merge partitions which can be described by similar models of evolution, which include:

*computational speed
*better estimation of model parameters. 

To carry out simultaneously model of evolution & partitioning scheme selection let's use:

```iqtree -s concatenation.nxs -sp gene_and_codon.prt -m TESTMERGEONLY  -redo```

We will overwrite the previous analysis, as the merging of partition will almost certainly result in a better model for our dataset.

Moreover, IQ-TREE provides edge-linked or edge-unlinked branch lengths between partitions:

* -q   partition_file: all partitions share the same set of branch lengths.
* -spp partition_file: like above but allowing each partition to have its own evolution rate.
* -sp  partition_file: each partition has its own set of branch lengths to account for, e.g. heterotachy (Lopez et al., 2002).

-spp is recommended for typical analysis while -q is unrealistic and -sp is very parameter-rich.
In real scenarios one should perform all three analyses and compare e.g. the BIC scores to determine the best-fit partition model, but for now we will skip this part
and chose the more reasonable assumption of a separate evolutionary rate of each partitions.










---


## further reading: 

[Here](http://www.iqtree.org/doc/Tutorial) you'll great tutorials from the authors themselves on ModelFinder and IQ-Tree.

[Very interesting paper on how concatenation/coalescence impacts mammalian phylogeny](https://onlinelibrary.wiley.com/doi/full/10.1111/cla.12170?casa_token=X0ctrSm4S1AAAAAA%3AgiB9v0MtJDO6vMWOigdvW9JrgYuJTebMen6zYxg9S0nP8MWIi2zA2fwWfi-lJlMCD9Ir1MDCzkBeyVwg).