# Building a Tree in a BI Framework



 
---




## intro: 

In this lesson we will use MrBayes [(Ronquist et al., 2012)](https://academic.oup.com/sysbio/article/61/3/539/1674894) to understand
the underlying concepts of Bayesian Inference in phylogenetics. In the meantime we will explore a bit the command-line to convert different
MSA formats and explore CIPRESS, the more comprehensive and outstanding computational infrastructure for phylogenetic




---



## MSA format conversion: 

At this stage we have already used ```.fas``` and ```.nxs```formats for Multiple Sequence Alignment. 
In this lesson we will also need a ```.phy``` - which is known as phylip-formatted file - for a tool called PartitionFinder2.
This format name comes from a very *ancient* piece of software - called of course [PHYLIP](https://en.wikipedia.org/wiki/PHYLIP) - which is turning 40 this year :-O
It is a very bare bone format which has in the firt line the number of species and the number of charachters separated by a space and then 
one entry for each line, with the sequence i.d. as the first element followed by a space and by the sequence itself.
The more format we encounter, the more we are required to convert them from one to another.
This can be done through a ton of programs, either [online](http://sequenceconversion.bugaco.com/converter/biology/sequences/nexus_to_phylip.php), 
with a graphical interface - as Aliview - or _via_ the command-line. We can learn a lot on the shell and
on phylogenetic formats by trying to handle conversions using bash! Here I put two one-liners which you can take a look at, written
not to be efficient or concise, but to use a good number of common bash commands:

from ```.nxs``` to ```.fasta```:

```
awk '/MATRIX/,EOF { print $0 }' concatenation.nxs | tail +2 | tr  \\t ' ' | sed -e 's/^ />/' 
| tr ' ' '\n' | sed '$d' | sed '$d' | sed '$d' > concatenation.fasta
```

from ```.nxs``` to ```.phy```:

```
nxs=concatenation.nxs; awk '/MATRIX/,EOF { print $0 }' $nxs | tail +2 | tr  \\t ' ' | sed '$d' | sed '$d' | sed '$d' > tmp.phy; 
n_sp=$(wc -l tmp.phy); n_car=$(grep -A 1 "MATRIX" $nxs | tail -1 | awk -F "\t" '{print $3}' | wc -c); cat tmp.ph
```

You can rewrite them in a more efficient way and/or add the missing conversion from ```.fas``` to ```.nxs``` and from ```.phy``` to ```.nxs```. 
it's a funny exercise I guess, and I can add your code here so it will be available for us all.




---




## another model selection & CIPRESS: 

Before the phylogenetic Bayesian Inference (BI), we have to carry out model selection again, mainly for two reasons:

* MrBayes supports slightly different models than ModelFinder and

* even if possible t's really painfull to translate / approximate models from ModelFinder to MrBayes.

Thus we will take advantage of this situation to leverage PartitionFinder2 [(Lanfear et al., 2016)]( https://doi.org/10.1093/molbev/msw260) 
which can carry out the model selection specifically for MrBayes. This situation teaches us a lesson on usability and how sometimes phylogenetic 
pipeline consist in a straight process (ModelFinder -> IQ-TREE), but sometimes they result in a quite fragmented workflow with multiple tools
and steps in between them.

The input of PartitionFinder2 are a ```.phy``` and a ```.cfg```.
```.cfg``` stands for configuration and is a quite widespread approach to customise analyses while using several pieces of software;
all the options are found inside the configuration file, including all the inputs required for the analysis. For convenience,
the input alignment is usually placed in the same folder where the configuration file is located, but a path can be specified as well.


You should e already familiar with many of the options, as they are quite similar to the ones of ModelFinder.
Here is how I _configured_ my ```.cfg```:


```
## ALIGNMENT FILE ##
alignment = infile.phy;

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##
models = mrbayes;

# MODEL SELECCTION: AIC | AICc | BIC #
model_selection = bic;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]

12S = 1-745;
ND2st = 746-1768/3;
ND2nd = 747-1768/3;
ND2rd = 748-1768/3;

## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##
[schemes]

search = greedy;
```

As anticipated before, we are running this analysis on CIPRESS For this reason, we need to specify ```infile.phy``` as the alignment, 
independently to whichever name we chose for our alignment: this is due to how CIPRESS handles the processes and files. 
We can then upload the ```.cfg``` and ```.phy``` and launch the analysis. If you need them you can find mine respectively 
[here](https://github.com/for-giobbe/phy/blob/master/examples/concatenation.phy) and 
[here](https://github.com/for-giobbe/phy/blob/master/examples/gene_and_codon_PF2.cfg).



As you can see CIPRESS is providing us with an almost perfect 1:1 trasposition of the program running, 
including the standard output, error and intermediate files! When the analysis is finished download the ```analysis.zip``` folder and 
open the ```best_scheme.txt```. As you can see this file contains similar information to ModelFinder:

```
Best partitioning scheme

Scheme Name       : start_scheme
Scheme lnL        : -14366.09814453125
Scheme BIC        : 29248.1509818
Number of params  : 69
Number of sites   : 1768
Number of subsets : 4

Subset | Best Model | # sites    | subset id                        | Partition names
1      | GTR+G      | 745        | 4c7442519ef19b615b794f84797e2a90 | 12S
2      | GTR+I+G    | 341        | 4405b051e7d135949d8131c87a543e99 | ND2st
3      | HKY+G      | 341        | 94d689124e5f756113ab5b5b286bf9fd | ND2nd
4      | GTR+G      | 341        | 576f3ac022d60c482c4ed21eb1fdaeae | ND2rd
```

In my case the model selection is quite similar to the one carried out by ModelFinder, which was:
```
    GTR+F+G4: 12S,
    TN+F+I+G4: ND2st,
    HKY+F+G4: ND2nd,
    TIM2+F+R2: ND2rd;
```



The ```best_scheme.txt``` also includes the MrBayes block:

```
begin mrbayes;

	charset Subset1 = 1-745;
	charset Subset2 = 746-1768\3;
	charset Subset3 = 747-1768\3;
	charset Subset4 = 748-1768\3;

	partition PartitionFinder = 4:Subset1, Subset2, Subset3, Subset4;
	set partition=PartitionFinder;

	lset applyto=(1) nst=6 rates=gamma;
	lset applyto=(2) nst=6 rates=invgamma;
	lset applyto=(3) nst=2 rates=gamma;
	lset applyto=(4) nst=6 rates=gamma;

	prset applyto=(all) ratepr=variable;
	unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all) tratio=(all);

end;
```



Warning: MrBayes only allows a relatively small collection of models. If any model in your analysis is not one that is included in MrBayes 
(e.g. by setting nst = 1, 2, or 6 for DNA sequences; or is not in the available list of protein models for MrBayes)then this MrBayes block 
will just set that model to nst = 6 for DNA, or 'wag' for Protein. Similarly, the only additional parameters that this 
MrBayes block will include are +I and +G. Other  parameters, such as +F and +X, are ignored. If you want to use this MrBayes block for 
your analysis, please make sure to check it carefully before you use it we've done our best to make it accurate, 
but there may be errors that remain!




---




## Bayesian Inference: 

Now we can proceed to append the model selection to the concatenation ```.nxs``` file - which is the format required by MrBayes.
You should have the one generated by phyutility or you can use [mine](https://github.com/for-giobbe/phy/blob/master/examples/concatenation.nxs).
Finally, we need to put the parameters for the search and run it. Go back to the NEXUS file and add at the end of the MrBayes block:

```
mcmc ngen=500000 printfreq=5000 samplefreq=5000 nruns=2 nchains=8 temp=0.02;
```

This specifies that:

* the analysis must run for 500k generations
* every 5k parameters are printed to standard output
* every 5k parameters are sampled
* two independent runs are carried out
* each run is composed of 8 different chains
* the temperature of the hot chain is 0.02

In the end the nexus file for MrBayes should look like [this](https://github.com/for-giobbe/phy/blob/master/examples/concatenation_mrbayes.nxs)


Then to run MrBayes you should type ```mb``` and then ```exe .nxs``` which in my case is ```exe concatenation_mrbayes.nxs```.


While running, MrBayes will print the (estimated) remaining time at the end of each line, showing the chain parameters, which
are also printed to the .mcmc file in the same folder of the NEXUS file. At the end of the run, you can check for convergence 
by plotting the last column of that file: the burnin must be selected when the standard deviation of average split frequencies starts
to be (i) relatively very low and (ii) stable.






---




## Post-Inference diagnotics: 

Tracer




---

Uncertainty can generally be observed either from scarce nodal support, polytomies  and/or  
sensitivity to the use of independent sampling of species and analytical frameworks (Yuan et al., 2016). 
Moreover, standard measures of clade support, such as posterior probabilities (Lewis et al., 2005) 
and bootstrap proportions (Simmons and Norton, 2014), can support several conflicting hypotheses with high apparent confidence.




---




## further reading: 

[Here](http://mrbayes.sourceforge.net/wiki/index.php/Manual) you will find the manual of MrBayes and some tutorials as well.