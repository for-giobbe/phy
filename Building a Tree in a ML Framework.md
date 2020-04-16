# Building a Tree in a ML Framework


---


## intro:

Maximum likelihood is a general statistical method for estimating parameters of a
probability model: for example a normal distribution can be described by two parameters: 
the mean and variance. Instead, in molecular phylogenetics there are a wide plethora of parameters, which include:

* rates of transitions between bases
* base composition
* descriptors of rate heterogeneity across sites
* branchlengths
* the tree itself!

The likelihood is defined as a quantity proportional to the probability of observing the data given the model: 
```
P(D|M)
```
If we have a model, we can calculate the probability the observations would have actually been observed as a function of the model. 
We then examine this likelihood function to see where it is at its greatest, and the value of the parameter of
interests (usually the tree and branch lengths) at that point is the maximum likelihood estimate of the parameter.

In this lesson we are going to compute a phylogenetic tree in a ML framework and explore a bit the relative support metrics,
which can inform us of the confidence relative to a split. At this point you should have a concatenation and a partition file; 
if not you can use mine [here](https://github.com/for-giobbe/phy/tree/master/examples).


## inferring gene trees - unpartitioned analyses:


During last lesson we found the best-fit model of evolution without doing tree reconstruction, running the line:

```
iqtree -s ND2_p_aligned.n.gb.fasta -m MF 
```

But we can perform both the search for the best-fit model and the phylogenetic inference by just using the ```-m MFP``` flag, so that
after ModelFinder, IQ-TREE will immediately start the tree reconstruction under the best-fit partition model. Let's use the line:

```
iqtree -s ND2_p_aligned.n.gb.fasta -m MFP
```

As always if one takes the effort to read the standard output of a software, he/she can lean a lot. For example we can notice that
IQ-TREE is using the model selection we previously calculated, as shown by the line:

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

* IQ-TREE report (.iqtree)
* Maximum-likelihood tree (.treefile)
* Likelihood distances (.mldist)
* Screen log file (.log)

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
We can root them using our _a priori_ biological knowledge or use approches as the mid point rooting.


---


## inferring species tree - partitioned analyses:

We can carry out the phylogenetic inference on our concatenation just by specifying the relative partition file with the 
```-spp``` flag and adding ```+MERGE``` in the model specification, so that the model tries to find the best-fit partitioning scheme
by possibly merging partitions.

```
iqtree -s concatenation.nxs -spp gene_and_codon.prt -m MFP+MERGE -redo 
```

The outputs generated will be the same as the ones produced by the unpartitioned analysis. 
Among the large number of parameters which can affect the tree search process in IQ-TREE, two of the more decisive are:

```-nstop```  which specify the number of unsuccessful iterations to stop. DEFAULT: 100
```-pers```   which specify perturbation strength (between 0 and 1) for randomized NNI. DEFAULT: 0.5

---


## inferring nodal support using different metrics:


Several metrics of clade support are possible and should be combined to gain more confidence.
Here are the more frequently used in IQ-TREE:

* Parametric & Nonparametric bootstrap
* SH-like approximate likelihood ratio test 
* ...

I really like this explanation of parametric and non-parametric bootstrap:

> Non-parametric bootstrapping was developed by Efron in 1979 as a general statistical method for estimating the parameters 
> of an unknown probability distribution by resampling from an existing sample that was drawn from this distribution. 
> The method was transferred to phylogenetic reconstruction by Felsenstein in [this](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.1985.tb00420.x) paper from 1985.
> Within molecular phylogenetics it works as follows: 
> from an alignment of length n, columns are randomly drawn with replacement n times. 
> The drawn columns are arranged in a new dataset, a bootstrapped alignment of length n. 
> From this bootstrapped alignment, a phylogenetic tree is constructed by following the same method of phylogenetic analysis 
> as was used for the analysis of the original alignment. 
> This process of constructing bootstrap alignments and bootstrap trees is repeated a large number of times,
> and the resulting trees are stored. The percentage with which a certain bipartition of the taxon set is present in the bootstrap trees
> (the bootstrap value) can be taken as a measure of how homogeneously this bipartition of sequences 
> (i.e., the respective branch in the underlying topology) is supported by the data. 
> Bootstrap values are often summarized by constructing the majority-rule consensus from the bootstrap trees
> or by annotating them on the "best" tree.

> In a probabilistic context, there is an alternative way of generating replicate alignments from given data by computer simulation. 
> This approach proposed by Efron in 1985 is model-based, and hence is commonly referred to as parametric bootstrapping. In its first step, 
> a model of DNA substitution and a phylogenetic tree T are estimated from the original alignment X.
> Using this model, replicate alignments Xi are generated, i.e., sequences are simulated along T. 
> Subsequently, phylogenetic trees Ti are computed for each of the alignments Xi, and branch support values are derived 
> (as in non-parametric bootstrapping) by computing the percentage with which a certain branch occurs in the set of generated trees Ti.
> Support values derived by parametric bootstrapping depend to a large extent on the model estimated from the original alignment. 
> For this reason the method can be used for testing the model inferred from the original alignment as a null hypothesis (Goldman 1993).

Along with its non-replacement version (not too popular nowadays):

> The jackknife is a resampling method closely related to non-parametric bootstrapping. 
> It works by randomly deleting a certain percentage of columns from the original alignment. 
> Usually 50 per cent of the columns are deleted (delete-half jackknife proposed by Felsenstein in 1985), 
> which is equivalent to drawing n/2 columns from the original alignment of length n without replacement. 
> As in non-parametric bootstrapping, the resampling is iterated, trees are computed from the jackknife alignments, 
> and branch-support values are derived as the percentage with which a certain branch is present in the jackknife topologies.

Adapted (slightly) from:
_M. David Posada. Bioinformatics for DNA Sequence Analysis, Springer Protocols, pp.113-137 Methods in Molecular Biology._

The functioning of aLRT is quite interesting as well:

> The approximate likelihood ratio test (aLRT) for branches is closely related to a conventional LRT. 
> The standard LRT uses the test statistics 2(L1 −L0), where L1 (alternative hypothesis) is the log-likelihood of the current tree 
> and L0 (null hypothesis) is the log-likelihood of the same tree, but with the branch of interest being collapsed.
> The aLRT approximates this test statistics in a slightly conservative but practical way,
> where L2 corresponds to the second best NNI configuration around the branch of interest.
> Such test is fast because the log-likelihood value L2 is computed by optimising only over
> the branch of interest and the four adjacent branches, while other parameters are fixed at
> their optimal values corresponding to the best ML tree.
> Thus, the aLRT does not account for other possible topologies that would be highly likely but quite different from the current topology. 
> This implies that the aLRT performs well when the data contains a clear phylogenetic
> signal, but not as well in the opposite case, where it tends to give a local view on
> the branch of interest only.
> Note that parametric branch supports are based on the assumption that the evolutionary model used to infer the trees is the correct one. 
> The rational behind the aLRT clearly differs from bootstrap. 
> Basically, while aLRT values are derived from testing hypotheses, the bootstrap proportion is a repeatability measure. 

Adapted from:
_Guindon et al., 2009. Estimating maximum likelihood phylogenies with PhyML._

---

Let's get some hands-on:

* Nonparametric bootstrap

	The standard nonparametric bootstrap is invoked by the ```-b``` option, which also specifies the number 
	of bootstrap replicates (100 is the minimum recommended number).

	```
	iqtree -s ND2_p_aligned.n.gb.fasta -b 100
	```


* Parametric bootstrap:
	
	IQ-TREE implements UFB2 - Ultra Fast Bootstrap 2 described in [Hoang et al., 2018](https://academic.oup.com/mbe/article/35/2/518/4565479)
	The ```-B``` flag specifies the number of replicates where 1000 is the minimum number recommended.  IQ-TREE also has the option to further optimize each bootstrap tree using a hill-climbing nearest neighbor interchange (NNI) search,
	based directly on the corresponding bootstrap alignment. It's specified through the ```-bnni ``` option to reduce the risk of overestimating branch supports with UFBoot due to severe model violations. 

	```
	iqtree -s ND2_p_aligned.n.gb.fasta -bb 1000 -bnni
	```


* SH-like approximate likelihood ratio test:

	IQ-TREE implements a non-parametric approximate likelihood ratio test based on a Shimodaira-Hasegawa-like procedure via the
	flag ```-alrt```.

	```
	iqtree -s ND2_p_aligned.n.gb.fasta -alrt 1000
	```

We can combine the three metrics in the same analysis and have them annotated on the "best" Maximum Likelilhood phylogeny. 
It's then easier to observe wether the different support metrics are giving contrasting results through our phylogenies.

```
iqtree -s ND2_p_aligned.n.gb.fasta -bb 1000 -bnni -b 100 -alrt 1000
```

Let's take a look at the ```.iqtree``` among the other output files, bearing in mind that the values related to different metrics 
should be treated differently: with the non-parametric bootstrap and SH-aLRT you should start to believe in a clade if 
it has >= 80% support, while with UFBoot it should be >= 95%, 

To conclude: there are several metrics of support in phylogenetics which can provide different perspective on the confidence 
of a clade/bipartition. Moreover they can sometimes be informative of biological processes such as ILS (Incomplete Lineage Sorting) 
or adaptive radiations. Aside the traditional ones (which we just went trough) some new ones get proposed and/or implemented 
from time to time. This is the case of gCF and sCF (genes and sites Concordance Factors) for which I left some additional information 
in the further reading paragraph at the end of the lesson. Moreover consider that different frameworks can have different support metrics,
as the Posterior Probabilities (PP) in the Bayesian Inference.


---


## automation

When large amount of loci are available for phylogenetic inference, IQ-TREE provides the -S flag to compute individual loci trees 
given a partition file or a directory:

```
iqtree -s ALN_FILE -S PARTITION_FILE --prefix loci -T AUTO
```

or

```
iqtree -S ALN_DIR --prefix loci -T AUTO

```

IQ-TREE automatically detects that ALN_DIR is a directory and will load all alignments within the directory. 
The -S takes the same argument as -s except that it performs model selection and tree inference separately for each 
partition or alignment. The output files are similar to those from a partitioned analysis,
except that loci.treefile now contains a set of trees.




---




## further reading: 

[Here](http://www.iqtree.org/doc/Tutorial) you'll find great tutorials from the authors themselves on ModelFinder and IQ-TREE.

resources on concordance factors: [paper](https://www.biorxiv.org/content/10.1101/487801v2) & [tutorial](http://www.robertlanfear.com/blog/files/concordance_factors.html)

[Here](http://www.iqtree.org/doc/iqtree-doc.pdf) you'll find the manual for IQ-TREE, it's quite user-friendly and exhaustive.   

[Here](http://www.iqtree.org/doc/Command-Reference) a coprensive commands list for IQ-TREE.