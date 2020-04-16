# Building a Tree in a ML Framework


---


## intro:

Maximum likelihood is a general statistical method for estimating parameters of a
probability model. A familiar one might be the normal distribution of a population with two parameters: 
the mean and variance. In phylogenetics there are many parameters, including:

* rates of transitions between states
* base composition
* descriptors of among sites variation
* ...
* most importantly, the tree itself!

Likelihood is defined to be a quantity proportional to the probability of observing the data given
the model: 
```

P(D|M)
```

Thus, if we have a model we can calculate the probability the observations would have actually been observed as a function of the model. 
We then examine this likelihood function to see where it is at its greatest, and the value of the parameter of
interests (usually the tree and/or branch lengths) at that point is the maximum likelihood estimate of the parameter.

In this lesson we are going to compute a phylogenetic tree in a ML framework and explore a bit the relative support metrics,
which can inform us of the confidence relative to a split.


## inferring gene trees - unpartitioned analyses


During last lesson we found the best-fit partition model without doing tree reconstruction, running the line:

```
iqtree -s concatenation.nxs -spp gene_and_codon.prt -m MF+MERGE  -redo
```

We can also carry out the phylogenetic inference using the ```-m TESTNEWMERGE``` flag.
After ModelFinder found the best partition, IQ-TREE will immediately start the tree reconstruction under the best-fit partition model.

```
iqtree -s example.phy -spp example.nex -m TESTNEWMERGE
```

Most phylogenetic programs produce unrooted trees as they are not aware about any biological background.

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