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
iqtree -s ND2_p_aligned.n.gb.fasta -m MFP 
```

But we can perform both the search for the best-fit model and the phylogenetic inference by just using the ```-s``` flag, so that
after ModelFinder, IQ-TREE will immediately start the tree reconstruction under the best-fit partition model. Let's use the line:

```
iqtree -s ND2_p_aligned.n.gb.fasta -m TESTNEWMERGE
```

As always if one takes the effort to read the standard output of a software, he/she can lean a lot. For example we can notice that
IQ-Tree is using the model selection we previously calculated, as shown by the line:

```
NOTE: Restoring information from model checkpoint file ND2_p_aligned.n.gb.fasta.model.gz
```

and we can see the inference is carried out in three main steps:

```
--------------------------------------------------------------------
|             INITIALIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Generating 98 parsimony trees... 0.173 second
Computing log-likelihood of 98 initial trees ... 0.657 seconds
Current best score: -9537.664

Do NNI search on 20 best initial trees
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 1: -9536.544
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 2: -9536.138
Iteration 10 / LogL: -9536.325 / Time: 0h:0m:2s
Iteration 20 / LogL: -9538.238 / Time: 0h:0m:3s
Finish initializing candidate tree set (7)
Current best tree score: -9536.138 / CPU time: 1.841
Number of iterations: 20
```


```
--------------------------------------------------------------------
|               OPTIMIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Iteration 30 / LogL: -9540.555 / Time: 0h:0m:3s (0h:0m:9s left)
Iteration 40 / LogL: -9540.146 / Time: 0h:0m:4s (0h:0m:6s left)
Iteration 50 / LogL: -9537.749 / Time: 0h:0m:4s (0h:0m:5s left)
Iteration 60 / LogL: -9536.378 / Time: 0h:0m:5s (0h:0m:3s left)
Iteration 70 / LogL: -9539.614 / Time: 0h:0m:6s (0h:0m:2s left)
Iteration 80 / LogL: -9536.957 / Time: 0h:0m:6s (0h:0m:1s left)
Iteration 90 / LogL: -9540.113 / Time: 0h:0m:7s (0h:0m:0s left)
Iteration 100 / LogL: -9536.266 / Time: 0h:0m:7s (0h:0m:0s left)
TREE SEARCH COMPLETED AFTER 103 ITERATIONS / Time: 0h:0m:7s

```


```
--------------------------------------------------------------------
|                    FINALIZING TREE SEARCH                        |
--------------------------------------------------------------------
Performs final model parameters optimization
Estimate model parameters (epsilon = 0.010)
1. Initial log-likelihood: -9536.138
Optimal log-likelihood: -9536.135
Rate parameters:  A-C: 1.00000  A-G: 2.43981  A-T: 1.00000  C-G: 1.00000  C-T: 8.41355  G-T: 1.00000
Base frequencies:  A: 0.448  C: 0.113  G: 0.079  T: 0.361
Proportion of invariable sites: 0.187
Gamma shape alpha: 0.722
Parameters optimization took 1 rounds (0.012 sec)
BEST SCORE FOUND : -9536.135
Total tree length: 7.923

```

Let's take a look at the outputs:

* IQ-TREE report:                .iqtree
* Maximum-likelihood tree:       .treefile
* Likelihood distances:          .mldist
* Screen log file:                .log

Let's take a deeper look at the IQ-TREE report:



We can find some information regarding the model of substitution:

```
SUBSTITUTION PROCESS
--------------------

Model of substitution: TN+F+I+G4

Rate parameter R:

  A-C: 1.0000
  A-G: 2.4398
  A-T: 1.0000
  C-G: 1.0000
  C-T: 8.4136
  G-T: 1.0000

State frequencies: (empirical counts from alignment)

  pi(A) = 0.4475
  pi(C) = 0.1132
  pi(G) = 0.07858
  pi(T) = 0.3607

Rate matrix Q:

  A   -0.4903   0.08342    0.1412    0.2657
  C    0.3297    -2.623   0.05788     2.235
  G    0.8043   0.08342    -1.153    0.2657
  T    0.3297    0.7019   0.05788    -1.089

Model of rate heterogeneity: Invar+Gamma with 4 categories
Proportion of invariable sites: 0.1867
Gamma shape alpha: 0.7223

 Category  Relative_rate  Proportion
  0         0              0.1867
  1         0.09706        0.2033
  2         0.4596         0.2033
  3         1.147          0.2033
  4         3.214          0.2033
Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category.
```

along with the best tree:

```
MAXIMUM LIKELIHOOD TREE
-----------------------

Log-likelihood of the tree: -9536.1354 (s.e. 204.8100)
Unconstrained log-likelihood (without tree): -5497.8470
Number of free parameters (#branches + #model parameters): 40
Akaike information criterion (AIC) score: 19152.2707
Corrected Akaike information criterion (AICc) score: 19155.6108
Bayesian information criterion (BIC) score: 19349.4905

Total tree length (sum of branch lengths): 7.9233
Sum of internal branch lengths: 1.1308 (14.2712% of tree length)

NOTE: Tree is UNROOTED although outgroup taxon 'Mantis_religiosa' is drawn at root
```

Remember that most phylogenetic programs produce unrooted trees, as they are not aware about any biological background.


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


## further reading: 

[Here](http://www.iqtree.org/doc/Tutorial) you'll find great tutorials from the authors themselves on ModelFinder and IQ-Tree.

resources on concordance factors: [paper](https://www.biorxiv.org/content/10.1101/487801v2) & [tutorial](http://www.robertlanfear.com/blog/files/concordance_factors.html)

