# Building a Tree in a ML Framework


---


## intro:

Maximum likelihood is a general statistical method for estimating parameters of a
probability model. A familiar one might be the normal distribution of a population with two parameters: 
the mean and variance. In molecular phylogenetics there are a wide plethora of parameters, which may include:

* rates of transitions between bases
* base composition
* descriptors of rate heterogeneity across sites
* branchlengths
* and most importantly, the tree itself!

Likelihood is defined to be a quantity proportional to the probability of observing the data given
the model: 
```
P(D|M)
```
Thus, if we have a model we can calculate the probability the observations would have actually been observed as a function of the model. 
We then examine this likelihood function to see where it is at its greatest, and the value of the parameter of
interests (usually the tree and/or branch lengths) at that point is the maximum likelihood estimate of the parameter.

In this lesson we are going to compute a phylogenetic tree in a ML framework and explore a bit the relative support metrics,
which can inform us of the confidence relative to a split. At this point you should have a concatenation and a partition file; 
if not you can use mine [here](https://github.com/for-giobbe/phy/tree/master/examples).


## inferring gene trees - unpartitioned analyses


During last lesson we found the best-fit model of evolution without doing tree reconstruction, running the line:

```
iqtree -s ND2_p_aligned.n.gb.fasta -s  -m MFP 
```

But we can perform both the search for the best-fit model and the phylogenetic inference by just using the ```-s``` flag, so that
after ModelFinder, IQ-TREE will immediately start the tree reconstruction under the best-fit partition model. Let's use the line:

```
iqtree -s example.phy -spp example.nex -m TESTNEWMERGE
```


Let's take a look at the results:


Remember that most phylogenetic programs produce unrooted trees as they are not aware about any biological background.




When large amount of loci are available for phylogenetic inference, IQ-TREE provides the -S flag to compute individual loci trees 
given a partition file or a directory:

```
iqtree -s ALN_FILE -S PARTITION_FILE --prefix loci -T AUTO
```

or

```
iqtree -S ALN_DIR --prefix loci -T AUTO
```

In the second case, IQ-TREE automatically detects that ALN_DIR is a directory and will load all alignment files within the directory. 
So -S takes the same argument as -p except that it performs model selection (ModelFinder) and tree inference separately for each 
partition/alignment. The output files are similar to those from a partitioned analysis,
except that loci.treefile now contains a set of trees.



















---


## ingerring species tree - partitioned analyses


---


## inferring nodal support using different metrics

* UFB - ultra fast boootstrap

The standard nonparametric bootstrap is invoked by the -b option:

```
iqtree -s example.phy -m TIM2+I+G -b 100
```

-b specifies the number of bootstrap replicates where 100 is the minimum recommended number. 

[Hoang et al., 2018](https://academic.oup.com/mbe/article/35/2/518/4565479)
-B specifies the number of bootstrap replicates where 1000 is the minimum number recommended. 


 provide a new option -bnni to reduce the risk of overestimating branch supports with UFBoot due to severe model violations. 
 With this option UFBoot will further optimize each bootstrap tree using a hill-climbing nearest neighbor interchange (NNI) search based directly on the corresponding bootstrap alignment.
 Thus, if severe model violations are present in the data set at hand, users are advised to append -bnni to the regular UFBoot command:

```
iqtree -s example.phy -m TIM2+I+G -B 1000 -bnni
```


* SH-like approximate likelihood ratio test 


---