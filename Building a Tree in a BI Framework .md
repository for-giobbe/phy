# Building a Tree in a BI Framework



 
---




## intro: 

In this lesson we will use MrBayes [(Ronquist et al., 2012)](https://academic.oup.com/sysbio/article/61/3/539/1674894) to understand
the underlying concepts of Bayesian Inference in phylogenetics. Recall that in the Bayesian framework we estimate parameters from their posterior distribution, 
instead of finding the best point estimates as is done in a Maximum Likelihood framework.
In the meantime we will explore a bit the command-line to convert different
MSA formats and explore CIPRESS, the more comprehensive and outstanding computational infrastructure for phylogenetic





---




## MSA format conversion: 

At this stage we have already used ```.fas``` and ```.nxs```formats for Multiple Sequence Alignment. 
In this lesson we will also need a ```.phy``` - a phylip-formatted file - for a tool called PartitionFinder2.
This format name comes from a very *ancient* piece of software - called of course [PHYLIP](https://en.wikipedia.org/wiki/PHYLIP) - which is turning 40 this year :-O
It is a very bare bone format which has in the first line the number of species and the number of characters separated by a space and then 
one entry (whether species or specimen) for each line, with the sequence i.d. as the first element followed by a space and by the sequence itself.

The more format we encounter, the more we are required to convert them from one to another.
This can be done through a ton of programs, either [online](http://sequenceconversion.bugaco.com/converter/biology/sequences/nexus_to_phylip.php), 
or _via_ several software with a graphical interface - as Aliview. But we can learn a lot on both the shell and
the different MSA formats by trying to handle conversions using the command-line! Here I put two one-liners which you can take a look at, written
not to be efficient or concise, but to use a good number of common bash commands:

from ```.nxs``` to ```.fasta```:

```
awk '/MATRIX/,EOF { print $0 }' concatenation.nxs | tail +2 | tr  \\t ' ' | sed -e 's/^ />/' | tr ' ' '\n' | sed '$d' | sed '$d' | sed '$d' > concatenation.fasta
```

from ```.nxs``` to ```.phy```:

```
nxs=concatenation.nxs;
awk '/MATRIX/,EOF { print $0 }' $nxs | tail +2 | tr  \\t ' ' | sed 's/^ //g' | sed '$d' | sed '$d' | sed '$d' > tmp.phy; 
n_sp=$(wc -l tmp.phy | awk '{print $1}'); 
n_car=$(grep -A 1 "MATRIX" $nxs | tail -1 | awk -F "\t" '{print $3}' | wc -c); 
echo $n_sp $n_car > first_line.tmp; cat first_line.tmp tmp.phy > concatenation.phy; rm *tmp*
```

For now we can just use the latter to convert our concatenation file from ```.nxs``` to ```.fasta```, 
but you can rewrite them in a more efficient way and/or add the missing conversion from ```.fas``` to ```.nxs``` and from ```.phy``` to ```.nxs```. 
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

While the first five parameters are quite straight forward, temp can be a little bit more obscure in the beginning 
(being somehow similar to the perturbance in IQTREE).
The higher the temperature, the more likely the heated chains are to move between isolated peaks in the posterior distribution. 
However, excessive heating may lead to very low acceptance rates for swaps between different chains.

In the end the nexus file for MrBayes should look like [this](https://github.com/for-giobbe/phy/blob/master/examples/concatenation_mrbayes.nxs)


Then to run MrBayes you should type ```mb``` and then ```exe .nxs``` which in my case is ```exe concatenation_mrbayes.nxs```.
The origram will instantly start writing to the standard output:


Here there is a nice _resumee_ of model parameters:

```
      Active parameters: 

                             Partition(s)
         Parameters          1  2  3  4
         ------------------------------
         Tratio              .  .  1  .
         Revmat              2  3  .  4
         Statefreq           5  6  7  8
         Shape               9 10 11 12
         Pinvar              . 13  .  .
         Ratemultiplier     14 14 14 14
         Topology           15 15 15 15
         Brlens             16 16 16 16
         ------------------------------
```

Here are the initial lnL of the eight chains relative to the two runs:

```
      Initial log likelihoods and log prior probs for run 1:
         Chain 1 -- -20685.580631 -- 67.855733
         Chain 2 -- -20960.717367 -- 67.855733
         Chain 3 -- -21280.231843 -- 67.855733
         Chain 4 -- -21439.401204 -- 67.855733
         Chain 5 -- -20776.225714 -- 67.855733
         Chain 6 -- -20936.729510 -- 67.855733
         Chain 7 -- -21255.572979 -- 67.855733
         Chain 8 -- -21218.825526 -- 67.855733

      Initial log likelihoods and log prior probs for run 2:
         Chain 1 -- -20688.742427 -- 67.855733
         Chain 2 -- -20962.169831 -- 67.855733
         Chain 3 -- -21057.290134 -- 67.855733
         Chain 4 -- -21436.439803 -- 67.855733
         Chain 5 -- -21181.932110 -- 67.855733
         Chain 6 -- -21330.400414 -- 67.855733
         Chain 7 -- -21155.010820 -- 67.855733
         Chain 8 -- -21001.961738 -- 67.855733
```

While running, MrBayes will print the the chain parameters, which are also printed to the ```.p``` files. 
Also the average standard deviation of split frequencies are printed: this values is a measure of similarity the tree toplogies sampled 
by the two independent runs and is highly informative of analysis convergence.
It's usually considered that average standard deviation below 0.01 are quite good, values between 0.01 and 0.05 may be adequate 
depending on the purpose of your analysis, while higher value shouldn't be accepted.


At the end of each line there's also the estimated remaining time.


```
190000 -- (-14392.295) (-14392.776) [-14389.627] (-14388.005) (-14383.686) (-14391.499) (-14386.012) (-14393.668) * [...8 more local chains...] (...0 remote chains...) -- 0:29:40

Average standard deviation of split frequencies: 0.048766
```


We will then stop the run, as half a million generations are way out of reach with our computational resources.


---




## Post-Inference diagnotics: 

sump 

sumt

Tracer

autocorrelation


---

Uncertainty can generally be observed either from scarce nodal support, polytomies  and/or  
sensitivity to the use of independent sampling of species and analytical frameworks (Yuan et al., 2016). 
Moreover, standard measures of clade support, such as posterior probabilities (Lewis et al., 2005) 
and bootstrap proportions (Simmons and Norton, 2014), can support several conflicting hypotheses with high apparent confidence.




---




## further reading: 

[Here](http://mrbayes.sourceforge.net/wiki/index.php/Manual) you will find the manual of MrBayes and some tutorials as well.