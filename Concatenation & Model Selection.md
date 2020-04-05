# Concatenation & Model Selection


---


## concatenation: 

Currently, two divergent systematic methods are commonly applied for inferring species trees: the supermatrix approach (concatenation) and the coalescent approach (gene trees are calculated and then reconciled in a species tree). 
You will find an interesting paper on the topic at the end of this tutorial, but here we will focus on the supermatrix approach.
After having aligned our genes, we will concatenate them using phyutility: the software has a wide array of functions and can be considere a swiss-knife for phylogeneticists.
As you always should, when using a new software, take a look at its [manual](https://github.com/blackrim/phyutility/blob/master/manual.pdf). 

After adjusting the path to the java executable of phyutility, try this string:

```java -jar /Applications/bio/phyutility/phyutility.jar -concat -in *.fasta -out concatenation.fasta```

Then take a look at the three outputfiles:

* a fasta file which stores all sequence information. This file is not designed to store the information of gene boundaries.
* a nexus file from which we will extract the nexus block which codes the information for the gene boundaries. 

We can reformat the nexus line which stores the gene boundaries information 

```[12S_total.fasta_gene1 1-816 CO1_total.fasta_gene2 817-2353 ND2_total.fasta_gene3 2354-3380 ]```

into something similar to 

```
DNA, 12S = 1-816
DNA, CO1 = 817-2353
DNA, ND2 = 2354-3380```

using an editor as ```nano```.


You can also get familiar with the nexus [format](http://informatics.nescent.org/wiki/NEXUS_Specification).
This is a step in phylogenetic pipelines which is often automated, as manual editing of partition files for hundred/thousands of genes is not possible: you can explore this in-house script which I made for the purpose. 
Moreover future versions of IQ-Tree can accept a folder of alignments as input, completely removing the need for the user to concatenate and edit partitions. 


Nonetheless this step is necessary for you to understand its underlying logic! For example, let's edit the partition file to take into account the different codon position (they evolve under different constrains due to the gen code degeneracy).
Let's use ```nano``` to transform our partition file from

```
DNA, 12S = 1-816
DNA, CO1 = 817-2353
DNA, ND2 = 2354-3380```

into 

```
DNA, 12S = 1-816
DNA, CO1st = 817-2353\3
DNA, CO1nd = 818-2353\3
DNA, CO1rd = 819-2353\3
DNA, ND2st = 2354-3380\3
DNA, ND2nd = 2355-3380\3
DNA, ND2rd = 2356-3380\3```


---


## model of evolution & partitioning scheme selection : 

To carry out the model selection we will use ModelFinder [Kalyaanamoorthy et al., 2017](https://www.nature.com/articles/nmeth.4285).
As for most of the steps which are carried out to build a tree, there is a staggering choice of different tools and is sometimes difficult to decide which to use.
The more widespread tool for model selection is PartitionFinder2 but we are using Model finder for two key reason:

* usability
* speed

Also I think it's a good idea to keep constantly using new and shiny tools. Let's try:


```iqtree -s CO1_total.fasta -m MFP```

and then we can open the relevant file by:

```gunzip CO1_total.fasta.model.gz; cat CO1_total.fasta.model```

The special MFP key word stands for ModelFinder Plus, which tells IQ-TREE to perform ModelFinder and the remaining analysis using the selected model.
ModelFinder computes the log-likelihoods of an initial parsimony tree for many different models and the Akaike information criterion (AIC), 
corrected Akaike information criterion (AICc), and the Bayesian information criterion (BIC). 
Then ModelFinder chooses the model that minimizes the BIC score (you can also change to AIC or AICc by adding the option -AIC or -AICc, respectively).

Otherwise the -m flag to specify the model name to use during the analysis, which can be a priori specified by the user (herse's a [list](http://www.iqtree.org/doc/Substitution-Models) of models implemented in ModelFinder).

```iqtree -s example.phy -m HKY+I+G```

What we've seen until now is the process we call model selection, for which, given a sequence the best-fit model get chosen (according to a metric of choice).
But in a concatenation framework we should carry out thi process on the whole concatenation instead of single alignements, but without loosing the information of the single genes. Let's try:

```iqtree -s concatenation.fasta -sp gene.prt -m MFP```

This string will result in separate models for each gene or partion. But there are several reasons for which we wanto to merge partitions which can be described by similar models of evolution,
due to several reasons which include computational burdens and better parameters estimation. To carry out simultaneously model of evolution & partitioning scheme selection let's use:

```iqtree -s concatenation.fasta -sp gene.prt -m TESTMERGEONLY```


---


## further reading: 

[Here](http://www.iqtree.org/doc/Tutorial) you'll great tutorials from the authors themselves on ModelFinder and IQ-Tree.

[Very interesting paper on how concatenation/coalescence impact mammalian phylogeny](https://onlinelibrary.wiley.com/doi/full/10.1111/cla.12170?casa_token=X0ctrSm4S1AAAAAA%3AgiB9v0MtJDO6vMWOigdvW9JrgYuJTebMen6zYxg9S0nP8MWIi2zA2fwWfi-lJlMCD9Ir1MDCzkBeyVwg)