# Building a Tree in a BI Framework

<br/>
<br/>

## intro: 

In this lesson we will use MrBayes [(Ronquist et al., 2012)](https://academic.oup.com/sysbio/article/61/3/539/1674894) to understand
the underlying concepts of Bayesian Inference in phylogenetics. Remember that in a Bayesian framework we estimate parameters from their posterior distribution, 
instead of finding the best single estimates as in a Maximum Likelihood framework.
In the meantime we will explore a bit the command-line to convert different
MSA formats and explore [CIPRESS](http://www.phylo.org/), the more comprehensive and outstanding computational infrastructure for phylogenetic


---

<br/>
<br/>

## MSA format conversion: 

At this stage we have already used ```.fas``` and ```.nxs``` formats for Multiple Sequence Alignment. 
In this lesson we will also need a ```.phy``` - a phylip-formatted file - which is required by [PartitionFinder2](http://www.robertlanfear.com/partitionfinder/).
This format name comes from a very ancient piece of software - called of course [PHYLIP](https://en.wikipedia.org/wiki/PHYLIP) - which is turning 40 this year :-O
It is a very bare bone format which has in the first line the number of species and the number of characters (nucleotides, proteins, ...) separated by a space and subsequently 
one entry (whether species or specimen) for each line, with the sequence i.d. as the first element followed by a space and then by the sequence itself.

The more format we encounter, the more we are required to convert them from one to another.
This can be done through a ton of programs, either [online](http://sequenceconversion.bugaco.com/converter/biology/sequences/nexus_to_phylip.php), 
or _via_ several software with a graphical interface, as Aliview. But we can learn a lot on both the shell and
the different MSA formats by trying to handle conversions using the command-line! Here I put two one-liners which you can take a look at - written
not to be efficient or concise - but to use a good number of the more usefull bash commands:

from ```.nxs``` to ```.fasta```:

```
awk '/MATRIX/,EOF { print $0 }' concatenation.nxs | tail +2 | tr  \\t ' ' | sed -e 's/^ />/' | tr ' ' '\n' | sed '$d' | sed '$d' | sed '$d' > concatenation.fasta
```

from ```.nxs``` to ```.phy```:

```
nxs=concatenation.nxs;
awk '/MATRIX/,EOF { print $0 }' $nxs | tail +2 | tr  \\t ' ' | sed 's/^ //g' | sed '$d' | sed '$d' | sed '$d' > tmp.phy; 
n_tax=$(head -4 $nxs | tail -1 | awk -F "NTAX=" '{print $2}' | awk '{print $1}');
n_car=$(head -4 $nxs | tail -1 | awk -F "NCHAR=" '{print $2}' | tr -d ";");
echo $n_tax $n_car > first_line.tmp; cat first_line.tmp tmp.phy > concatenation.phy; 
rm *tmp*
```

For now we can just use the latter to convert our concatenation file from ```.nxs``` to ```.phy```, 
but you can rewrite them in a more efficient way and/or add the missing conversion from ```.fas``` to ```.nxs``` and from ```.phy``` to ```.nxs```. 
it's a funny exercise I guess, and I can add your code here so it will be available for us all.


---

<br/>
<br/>

## another model selection & CIPRESS: 

Before the phylogenetic Bayesian Inference (BI), we have to carry out model selection again, mainly for two reasons:

* MrBayes supports slightly different models than ModelFinder (_e.g._ no base frequencies additional parameters)

* even if possible it's really painfull to translate / approximate models from ModelFinder to MrBayes.

We will take advantage of this situation to leverage PartitionFinder2 [(Lanfear et al., 2016)]( https://doi.org/10.1093/molbev/msw260) 
which can carry out the model selection specifically for MrBayes. This situation teaches us a lesson on usability and how sometimes phylogenetic 
pipeline consist in a straight process (ModelFinder -> IQ-TREE), but sometimes they result in a quite fragmented workflow with multiple tools
and steps in between them.

The input of PartitionFinder2 are a ```.phy``` and a ```.cfg```.
The latter stands for configuration and is a quite widespread approach to customise analyses while using several pieces of software:
all the options are found inside the configuration file, including all the inputs required for the analysis. For convenience,
the input alignment is usually placed in the same folder where the configuration file is located, but a path can be specified as well.


You should e already familiar with many of the options, as they are quite similar to the ones of ModelFinder.
Here is how I configured my ```.cfg```:


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

As anticipated before, we will be running this analysis on CIPRESS. For this reason, we need to specify ```infile.phy``` as the alignment, 
independently to whichever name we chose for it: this is due to how CIPRESS handles the processes and files. 
We can then upload the ```.cfg``` and ```.phy``` and launch the analysis. If you need them you can find mine respectively 
[here](https://github.com/for-giobbe/phy/blob/master/examples/concatenation.phy) and 
[here](https://github.com/for-giobbe/phy/blob/master/examples/gene_and_codon_PF2.cfg).

As you can see CIPRESS is providing us with an almost perfect 1:1 trasposition of the program running, 
including the standard output, error and intermediate files! When the analysis is finished download the ```analysis.zip``` folder and 
open the ```best_scheme.txt```. This file contains similar information to ModelFinder:

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

The ```best_scheme.txt``` also includes the MrBayes nexus block:

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

In this way of describing substitution models ```nst=6 rates=invgamma``` means ```GTR+I+G```. 
Here ```nst``` stands for number of substitution type and ```nst=6``` exactly describing a General Time Reversible model.


---

<br/>
<br/>

## Bayesian Inference: 

You already should have the one generated by phyutility or you can use [mine](https://github.com/for-giobbe/phy/blob/master/examples/concatenation.nxs).

Now we need to include the model and the parameter for the inference inside the ```.nxs``` file - which is the format required by MrBayes.
We can specify the parameters using the subsequent line, which will be included in the MrBayes block:

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

While the first five parameters are quite straightforward, ```temp``` can be a little bit more obscure in the beginning 
(being somehow similar to the perturbance in IQTREE).
The higher the temperature, the more likely the heated chains are to move between isolated peaks in the posterior distribution. 
However, excessive heating may lead to very low acceptance rates for swaps between different chains.

In the end the  the MrBayes block looking like:

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
	mcmc ngen=500000 printfreq=5000 samplefreq=5000 nruns=2 nchains=8 temp=0.02;
end;
```

The nexus file for MrBayes should look like [this](https://github.com/for-giobbe/phy/blob/master/examples/concatenation_mrbayes.nxs).


To run MrBayes you should type ```mb``` and then ```exe .nxs``` which in my case is ```exe concatenation_mrbayes.nxs```.
The program will instantly start writing to the standard output, _e.g._ this nice _resumee_ of model parameters:

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

While running, MrBayes will print the chain lnL, which are also written to the ```.p``` files. 
Along with the generation, the average standard deviation of split frequencies are printed: this values are a measure of similarity the tree toplogies sampled 
by the two independent runs and is highly informative of analysis convergence.
It's usually considered that average standard deviation below 0.01 are quite good, values between 0.01 and 0.05 may be adequate 
depending on the purpose of your analysis, while higher value imply that the analysis is far from stationarity.

```
190000 -- (-14392.295) (-14392.776) [-14389.627] (-14388.005) (-14383.686) (-14391.499) (-14386.012) (-14393.668) * [...8 more local chains...] (...0 remote chains...) -- 0:29:40

Average standard deviation of split frequencies: 0.048766
```

We will then stop the run before it finishes, as half a million generations are way out of reach with our computational resources.
The program has generated several files, including:

* ```.nxs.ckp``` a checkpoint file which let us resume analyses
* two ```.t``` files where are stored all topologies sampled during each run
* two ```.p``` files where are stored all the parameters sampled during each run

Here is how a line from a ```.p``` file looks like:
```
Gen	LnL	LnPr	TL{all}	kappa{3}	r(A<->C){1}	r(A<->G){1}	r(A<->T){1}	r(C<->G){1}	r(C<->T){1}	r(G<->T){1}	r(A<->C){2}	r(A<->G){2}	r(A<->T){2}	r(C<->G){2}	r(C<->T){2}	r(G<->T){2}	r(A<->C){4}	r(A<->G){4}	r(A<->T){4}	r(C<->G){4}	r(C<->T){4}	r(G<->T){4}	pi(A){1}	pi(C){1}	pi(G){1}	pi(T){1}	pi(A){2}	pi(C){2}	pi(G){2}	pi(T){2}	pi(A){3}	pi(C){3}	pi(G){3}	pi(T){3}	pi(A){4}	pi(C){4}	pi(G){4}	pi(T){4}	alpha{1}	alpha{2}	alpha{3}	alpha{4}	pinvar{2}	m{1}	m{2}	m{3}	m{4}
0	-2.068558e+04	6.785573e+01	6.600000e-01	1.000000e+00	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	1.666667e-01	2.500000e-01	2.500000e-01	2.500000e-01	2.500000e-01	2.500000e-01	2.500000e-01	2.500000e-01	2.500000e-01	2.500000e-01	2.500000e-01	2.500000e-01	2.500000e-01	2.500000e-01	2.500000e-01	2.500000e-01	2.500000e-01	1.000000e+00	1.000000e+00	1.000000e+00	1.000000e+00	0.000000e+00	1.000000e+00	1.000000e+00	1.000000e+00	1.000000e+00
```

And this is how a ```.t``` file looks like:
```
#NEXUS
[ID: 8851921581]
[Param: tree{all}]
begin trees;
   translate
       1 Mantis_religiosa,
       2 Carasisus_morosus,
       3 Phyllium_siccifollium,
       4 Peruphasma_schultei,
       5 Aretaon_asperrimus,
       6 Dryococelus_australis,
       7 Spinotectarchus_acornutus,
       8 Micrarchus_hystriculeus,
       9 Micrarchus_sp,
      10 Tectarchus_salebrosus,
      11 Asteliaphasma_jucundum,
      12 Tectarchus_ovobsessus,
      13 Clitarchus_hookeri,
      14 Niveaphasma_anulata,
      15 Acanthoxyla_sp,
      16 Argosarchus_horridus,
      17 Rumulus_artemis,
      18 Medauroidea_extradentata;
   tree gen.0 = [&U] (8:2.000000e-02,(14:2.000000e-02,((12:2.000000e-02,9:2.000000e-02):2.000000e-02,((5:2.000000e-02,(18:2.000000e-02,(6:2.000000e-02,(4:2.000000e-02,(17:2.000000e-02,(7:2.000000e-02,(15:2.000000e-02,(10:2.000000e-02,(11:2.000000e-02,(13:2.000000e-02,(16:2.000000e-02,2:2.000000e-02):2.000000e-02):2.000000e-02):2.000000e-02):2.000000e-02):2.000000e-02):2.000000e-02):2.000000e-02):2.000000e-02):2.000000e-02):2.000000e-02):2.000000e-02,3:2.000000e-02):2.000000e-02):2.000000e-02):2.000000e-02,1:2.000000e-02);
   tree gen.5000 = [&U] (((6:1.463107e-01,(((((14:8.464470e-02,11:8.928425e-02):1.157484e-02,((9:5.397960e-02,8:3.060077e-02):7.490231e-02,10:6.865646e-02):1.468076e-02):2.089574e-03,((16:6.856667e-02,15:3.019294e-02):1.355062e-02,13:7.246397e-02):1.581026e-02):8.312167e-03,12:6.764211e-02):3.917540e-02,7:1.129546e-01):6.890183e-02):7.965970e-02,(4:3.256303e-01,5:4.916759e-01):5.531344e-02):9.363812e-03,(((18:1.629848e-01,17:2.132906e-01):7.415411e-02,2:4.567366e-01):3.808575e-02,3:4.995555e-01):2.006223e-02,1:1.900453e+00);
```

As we have stopped our analysis before its end, we can use a small trick to carry out the post-inference analyses:

1. open the ```.p``` and ```.t``` file: if they have reached the same number of generations go to (3), if not go (2)
2. trim both files at the same number of generations
3. append a line with ```end;``` to both ```.t``` file
4. delete the line which specifies the parameters of the runs / chains.

If the analysis is finished, you don't need to modify anything and you can just proceed.


---

<br/>
<br/>

## Post-Inference diagnostics: 

While there is also plenty of Post-Inference diagnostics approaches which can be applied also to ML analyses, 
in BI framework they are definitively compulsory. 

We can start by recalling MrBayes with ```mb```; btw if you forgot to 
cancel the line which specifies the parameters of the runs / chains the analysis will start again!
We can then type ```sump``` which stands for Summary of Parameters. 
The summary statistics will be calculate using a relative burnin of 0.25 and a written to a ```.nxs.pstat``` file.
Nonetheless some interesting information are also printed to the standard output, such as this nice plot of
the generation (x-axis) versus the lnL (y-axis). This plot is useful to both decide a sensible burn in for the analysis
and to assess the stationarity of the inference. The numbers refer to the two chain and when their value overlap,
an asterisk get printed.
   
```
   +------------------------------------------------------------+ -14372.88
   |                                     1                      |
   |                                                            |
   |                 2                       1              2 2 |
   |       1             2          2   2 2    1         2      |
   |        1                  2     22              1 1       1|
   |  21121     11 2   21   1                 1   2        11  2|
   |1               2     2* 22   21 1  1  1    1          2    |
   |  1   2   21 21  12      1   2     1     2       2  2       |
   |2    1   21    1    21  2  1211    2  1   2 21 11        *1 |
   | 1  2    1            1     1   1 1     2  2  1     11      |
   |        2       1              2                  1         |
   | 2     2   22 2                      2  1       2  2  1     |
   |   2              1       1            2     2    2   2     |
   |                   1                                        |
   |                                               2            |
   +------+-----+-----+-----+-----+-----+-----+-----+-----+-----+ -14405.80
   ^                                                            ^
   100000                                                       410000
```

Then we have a table which include several values for all the parameters, including the Mean and the Estimated Sample Size (ESS) values.
Generally a good run should yield values > 200, but in certain situation (_e.g._ huge trees) values > 150 can be accepted.
An additional convergence diagnostic are the Potential scale Reduction Factor (PSRF) values, which should approach 1 as runs converge.

```
 Parameter       Mean      Variance     Lower       Upper       Median    min ESS*  avg ESS    PSRF+ 
   ----------------------------------------------------------------------------------------------------
   TL{all}      17.329780   16.308377    9.428004   25.028520   16.199610     37.51     49.04    0.993
   kappa{3}      2.169040    0.088139    1.676334    2.797372    2.155969     43.48     53.24    1.006
   r(A<->C){1}   0.096316    0.000232    0.065865    0.122570    0.096614     63.00     63.00    0.995
   r(A<->G){1}   0.224646    0.001095    0.149828    0.273643    0.227926     53.05     58.03    1.013
   r(A<->T){1}   0.201288    0.000403    0.167491    0.242217    0.201533     63.00     63.00    0.998
   r(C<->G){1}   0.029233    0.000385    0.005629    0.063191    0.026151     63.00     63.00    1.005
   r(C<->T){1}   0.427759    0.001132    0.356569    0.485578    0.430849     35.73     43.83    1.020
   r(G<->T){1}   0.020758    0.000119    0.001193    0.039786    0.021054     40.09     51.55    0.994
   r(A<->C){2}   0.110726    0.000550    0.068749    0.150560    0.108394     27.30     45.15    0.993
   r(A<->G){2}   0.139749    0.000942    0.083305    0.186575    0.140128     21.48     28.07    1.004
   r(A<->T){2}   0.066208    0.000217    0.039490    0.097473    0.064489     27.46     45.23    0.992
   r(C<->G){2}   0.042505    0.001237    0.000174    0.124189    0.033859     63.00     63.00    0.992
   r(C<->T){2}   0.569192    0.003012    0.439866    0.641182    0.575364     29.18     29.33    0.996
   r(G<->T){2}   0.071621    0.000433    0.032831    0.113118    0.070081     56.69     57.79    0.994
   r(A<->C){4}   0.010926    0.000056    0.000408    0.024243    0.009126     31.36     43.43    1.026
   r(A<->G){4}   0.143080    0.001426    0.072592    0.202886    0.139407     19.81     22.71    1.071
   r(A<->T){4}   0.005539    0.000007    0.000032    0.009430    0.005334     38.00     42.84    1.005
   r(C<->G){4}   0.040334    0.001535    0.000112    0.120180    0.028911     38.81     50.90    0.992
   r(C<->T){4}   0.753204    0.003654    0.612118    0.834710    0.768088     20.69     41.84    1.052
   r(G<->T){4}   0.046918    0.000592    0.002713    0.090417    0.043384     56.21     59.61    1.014
   pi(A){1}      0.452514    0.000168    0.429403    0.476068    0.452247     39.24     45.88    0.993
   pi(C){1}      0.141084    0.000086    0.124990    0.159765    0.140110     34.01     45.21    0.992
   pi(G){1}      0.074410    0.000073    0.055780    0.090690    0.073533     38.03     50.52    0.993
   pi(T){1}      0.331992    0.000154    0.309703    0.353117    0.333104     63.00     63.00    0.993
   pi(A){2}      0.541550    0.000246    0.515471    0.573513    0.541913     63.00     63.00    0.992
   pi(C){2}      0.070879    0.000070    0.055913    0.090216    0.069623     39.30     45.05    1.011
   pi(G){2}      0.097756    0.000095    0.081971    0.116750    0.096201     46.26     54.63    1.010
   pi(T){2}      0.289815    0.000175    0.260316    0.310865    0.290739     63.00     63.00    0.992
   pi(A){3}      0.218243    0.000313    0.184188    0.248344    0.218887     63.00     63.00    1.005
   pi(C){3}      0.217568    0.000237    0.186255    0.243109    0.218640     59.57     61.29    0.999
   pi(G){3}      0.087879    0.000092    0.069893    0.104906    0.087637     39.35     45.84    0.993
   pi(T){3}      0.476311    0.000472    0.439024    0.521354    0.473956     47.27     55.13    1.013
   pi(A){4}      0.565348    0.000232    0.533703    0.591135    0.565104     27.04     45.02    0.998
   pi(C){4}      0.074129    0.000019    0.066574    0.083400    0.073675     63.00     63.00    0.994
   pi(G){4}      0.026022    0.000007    0.021118    0.030520    0.025841     30.79     46.90    1.013
   pi(T){4}      0.334500    0.000197    0.308913    0.361712    0.333996     22.12     42.56    0.994
   alpha{1}      0.432695    0.001293    0.370616    0.515680    0.433740     56.79     59.90    0.994
   alpha{2}      1.212800    0.056251    0.731057    1.594491    1.199389     36.99     40.10    0.998
   alpha{3}      0.289073    0.001293    0.213586    0.348934    0.281397     40.34     51.67    0.999
   alpha{4}      0.978703    0.011472    0.754036    1.144843    0.971499     52.04     52.08    0.994
   pinvar{2}     0.238929    0.000954    0.194290    0.296642    0.238191     44.30     48.65    0.992
   m{1}          0.176358    0.001653    0.102047    0.248825    0.175473     40.67     43.33    1.002
   m{2}          0.371345    0.007796    0.229560    0.541795    0.371347     33.90     43.42    1.014
   m{3}          0.168411    0.001839    0.089810    0.230713    0.166950     38.69     42.91    1.005
   m{4}          4.259696    0.044155    3.808170    4.583160    4.274428     35.15     40.82    1.008
   ---------------------------------------------------------------------------------------------------- 
```

We can then type the ```sumt``` command which - in a similar way to the ```sump``` command - summarises 
the tree and branch length informations. Here are the files which are created:

* ```.trprobs``` trees that were found during the MCMC search, sorted by posterior probability
* ```.con.tre```  tree file with clade credibility / posterior probability values
* ```.vstat```    summary statistics for all branch and node parameters
* ```.tstat``` summary statistics for informative (_i.d._ non terminals) taxon bipartitions
* ```.parts``` bipartitions

We can also observe the more important summary statistics from the standard output.

Here there is a depiction of bipartitions: a taxon bipartition is defined by removing a branch on the tree, 
dividing the tree into those taxa to the left and right of the removed branch.
In this representation, dots are the taxa that are on one side of the partition and stars are the taxa on the other side. 
As you can observe bipartitions 1 to 18 represent terminal branches (they divide the dataset into a subset with one species
and another subset with all the other species), while the total number of branches (46) and of internal branches (32).
As you can see there are more bipartitions than the 36 branches which we expect ( recall the ```2* n. species -3``` formula).
This is because different trees are sampled throughout the analysis, which have different topologies and a number of 
possible bipartitions grater than the number of bipartitions in a single tree.

```
   ID -- Partition
   ------------------------
    1 -- .*****************
    2 -- .*................
    3 -- ..*...............
    4 -- ...*..............
    5 -- ....*.............
    6 -- .....*............
    7 -- ......*...........
    8 -- .......*..........
    9 -- ........*.........
   10 -- .........*........
   11 -- ..........*.......
   12 -- ...........*......
   13 -- ............*.....
   14 -- .............*....
   15 -- ..............*...
   16 -- ...............*..
   17 -- ................*.
   18 -- .................*
   19 -- ............*.**..
   20 -- ................**
   21 -- ..............**..
   22 -- .......**.........
   23 -- ......**********..
   24 -- ..........*..*....
   25 -- .....***********..
   26 -- .......*********..
   27 -- .......****.****..
   28 -- .***.*************
   29 -- .**.............**
   30 -- .......***........
   31 -- .**..*************
   32 -- ..*.............**
   33 -- .**...............
   34 -- .......****..*....
   35 -- .......**...*.**..
   36 -- .*..............**
   37 -- .......***..*.**..
   38 -- .***............**
   39 -- ..........*.****..
   40 -- .......**.*.****..
   41 -- ...**.............
   42 -- .........**..*....
   43 -- ...*************..
   44 -- ..*..*************
   45 -- ...*.***********..
   46 -- .......**.*..*....
   ------------------------
```

This interesting representation can be better interpreted using the conversion table associated:

```
   List of taxa in bipartitions:                                                 
                                                                                   
      1 -- Mantis_religiosa
      2 -- Carasisus_morosus
      3 -- Phyllium_siccifollium
      4 -- Peruphasma_schultei
      5 -- Aretaon_asperrimus
      6 -- Dryococelus_australis
      7 -- Spinotectarchus_acornutus
      8 -- Micrarchus_hystriculeus
      9 -- Micrarchus_sp
     10 -- Tectarchus_salebrosus
     11 -- Asteliaphasma_jucundum
     12 -- Tectarchus_ovobsessus
     13 -- Clitarchus_hookeri
     14 -- Niveaphasma_anulata
     15 -- Acanthoxyla_sp
     16 -- Argosarchus_horridus
     17 -- Rumulus_artemis
     18 -- Medauroidea_extradentata
```

You can also find some statistics for informative taxon bipartitions, where ```#obs``` stands for
how many times a bipartition has been observed. The asterisk denote bipartition which were just found in 
a single run out of the two.

```
   ID   #obs    Probab.     Sd(s)+      Min(s)      Max(s)   Nruns 
   ----------------------------------------------------------------
   19    78    1.000000    0.000000    1.000000    1.000000    2
   20    78    1.000000    0.000000    1.000000    1.000000    2
   21    78    1.000000    0.000000    1.000000    1.000000    2
   22    78    1.000000    0.000000    1.000000    1.000000    2
   23    78    1.000000    0.000000    1.000000    1.000000    2
   24    78    1.000000    0.000000    1.000000    1.000000    2
   25    78    1.000000    0.000000    1.000000    1.000000    2
   26    78    1.000000    0.000000    1.000000    1.000000    2
   27    64    0.820513    0.036262    0.794872    0.846154    2
   28    64    0.820513    0.036262    0.794872    0.846154    2
   29    59    0.756410    0.090655    0.692308    0.820513    2
   30    46    0.589744    0.072524    0.538462    0.641026    2
   31    43    0.551282    0.126917    0.461538    0.641026    2
   32    28    0.358974    0.036262    0.333333    0.384615    2
   33    24    0.307692    0.000000    0.307692    0.307692    2
   34    22    0.282051    0.072524    0.230769    0.333333    2
   35    22    0.282051    0.036262    0.256410    0.307692    2
   36    20    0.256410    0.000000    0.256410    0.256410    2
   37    16    0.205128    0.072524    0.153846    0.256410    2
   38    14    0.179487    0.000000    0.179487    0.179487    2
   39    12    0.153846    0.072524    0.102564    0.205128    2
   40    11    0.141026    0.054393    0.102564    0.179487    2
   41     7    0.089744    0.018131    0.076923    0.102564    2
   42     7    0.089744    0.054393    0.051282    0.128205    2
   43     6    0.076923    0.036262    0.051282    0.102564    2
   44     5    0.064103    0.054393    0.025641    0.102564    2
   45     4    0.051282    0.072524    0.000000    0.102564    1 *
   46     4    0.051282    0.072524    0.000000    0.102564    1 *
   ----------------------------------------------------------------
```

And here are the summary statistics for branch and node parameters:

```
                                                95% HPD Interval
                                              --------------------
   Parameter           Mean       Variance     Lower       Upper       Median     PSRF+  Nruns
   -------------------------------------------------------------------------------------------
   length{all}[1]     5.025501    1.327553    2.619254    6.749607    5.014459    1.011    2
   length{all}[2]     1.024853    0.053413    0.606395    1.327899    1.043001    0.998    2
   length{all}[3]     1.019522    0.043750    0.652923    1.376403    1.015712    1.048    2
   length{all}[4]     0.789487    0.029700    0.472181    1.102621    0.791617    1.030    2
   length{all}[5]     0.682053    0.046966    0.288199    1.034160    0.656724    1.060    2
   length{all}[6]     0.358092    0.007500    0.204460    0.509097    0.364765    1.035    2
   length{all}[7]     0.227501    0.002687    0.130121    0.327046    0.228461    1.001    2
   length{all}[8]     0.072925    0.000324    0.040592    0.101965    0.071365    1.033    2
   length{all}[9]     0.083527    0.000396    0.046631    0.118552    0.082402    0.995    2
   length{all}[10]    0.135941    0.001019    0.071498    0.175128    0.138679    1.074    2
   length{all}[11]    0.185187    0.001602    0.122649    0.267603    0.186121    1.023    2
   length{all}[12]    0.100245    0.000737    0.054195    0.149384    0.094959    0.988    2
   length{all}[13]    0.118501    0.000849    0.058337    0.164487    0.115537    1.066    2
   length{all}[14]    0.179818    0.001594    0.111610    0.249414    0.177530    1.037    2
   length{all}[15]    0.088897    0.000619    0.047591    0.128974    0.086298    1.037    2
   length{all}[16]    0.132697    0.000825    0.081025    0.186584    0.133709    0.988    2
   length{all}[17]    0.537564    0.011663    0.320523    0.732439    0.533183    1.063    2
   length{all}[18]    0.383215    0.007115    0.212846    0.529269    0.381975    1.001    2
   length{all}[19]    0.029015    0.000225    0.006395    0.052787    0.027513    1.005    2
   length{all}[20]    0.309560    0.006009    0.167760    0.449037    0.301577    1.019    2
   length{all}[21]    0.033070    0.000115    0.017839    0.053328    0.031480    0.987    2
   length{all}[22]    0.172255    0.001205    0.115284    0.243705    0.175538    1.002    2
   length{all}[23]    0.157918    0.002174    0.093674    0.253700    0.155782    0.995    2
   length{all}[24]    0.027919    0.000154    0.008089    0.050635    0.025864    1.015    2
   length{all}[25]    0.238083    0.005005    0.131511    0.364093    0.239854    1.073    2
   length{all}[26]    0.118525    0.001017    0.061124    0.166572    0.118054    1.065    2
   length{all}[27]    0.025352    0.000253    0.003968    0.058422    0.021665    0.989    2
   length{all}[28]    0.303126    0.023794    0.003661    0.604934    0.276902    0.993    2
   length{all}[29]    0.087030    0.001194    0.025831    0.148802    0.084464    0.998    2
   length{all}[30]    0.022545    0.000133    0.003588    0.046518    0.021492    1.003    2
   length{all}[31]    0.055789    0.001557    0.002338    0.129540    0.048530    0.983    2
   length{all}[32]    0.053285    0.001650    0.002934    0.143063    0.048980    0.976    2
   length{all}[33]    0.048459    0.001004    0.000941    0.097250    0.047258    0.959    2
   length{all}[34]    0.010545    0.000028    0.003047    0.021410    0.009932    0.996    2
   length{all}[35]    0.016982    0.000094    0.004665    0.032884    0.017381    0.965    2
   length{all}[36]    0.044277    0.001206    0.003212    0.114249    0.036272    1.026    2
   length{all}[37]    0.008014    0.000061    0.000202    0.014316    0.007694    0.989    2
   length{all}[38]    0.029960    0.001162    0.003229    0.098225    0.012157    0.926    2
   length{all}[39]    0.009688    0.000029    0.001456    0.016265    0.009671    0.915    2
   length{all}[40]    0.014688    0.000169    0.000390    0.031141    0.008365    1.042    2
   length{all}[41]    0.102426    0.009390    0.003782    0.305491    0.073976    1.297    2
   length{all}[42]    0.008689    0.000090    0.001916    0.029206    0.004243    0.915    2
   length{all}[43]    0.065419    0.002263    0.003368    0.131051    0.070553    1.091    2
   length{all}[44]    0.072598    0.000653    0.039144    0.099471    0.079721    0.775    2
   length{all}[45]    0.014515    0.000017    0.008626    0.017673    0.016896     NA      1 *
   length{all}[46]    0.007429    0.000024    0.002922    0.014292    0.007063     NA      1 *
   -------------------------------------------------------------------------------------------
```
At last it is also telling us how may trees sum up to a cumulated probability of 50%, 90% and so on:

```

   Credible sets of trees (63 trees sampled):
      50 % credible set contains 24 trees
      90 % credible set contains 56 trees
      95 % credible set contains 60 trees
      99 % credible set contains 63 trees
```

For both command the burnin can be modified in absolute terms ```burnin = 4000``` (will remove the first 4k generations) or relative terms ```burnin = 0.4``` (will remove the first 40% generations)
by adding this arguments to either  ```sump``` or  ```sumt```.


A nicer way to carry out this post-inference diagnostics (which sometimes are also done during the inference) is Tracer,
which offers a much easier way to explore what's going on in our inference. Nonetheless sometimes using the ```sump```and ```sumt```
commands can be faster when working on remote acces servers, as to use Tracer we have to move very large files. 
In Tracer we want to see the famous "fuzzy caterpillar" - by far the most loved animals by phylogeneticists -
which implies that there's no autocorrelation in our analysis and that it has reached stationarity. You can 
take a look  in Tracer at an [awful](https://github.com/for-giobbe/phy/raw/master/examples/bad_trace_example.p) 
and [cool](https://github.com/for-giobbe/phy/raw/master/examples/good_trace_example.p) sampling coming fromt two
distinct runs of MrBayes, to get a sense of how much they can be different.

---

<br/>
<br/>

## Posterior Probabilities & nodal support overestimation: 

As we have seen before, different metrics of nodal support are associated with different approaches to the phylogenetic inference.
BI analyses are associated to Posterior Probabilities - also called clade credibility values - which are calculated from the frequencies 
with which they occur in the trees sampled at the stationary state, ranging from 0 to 1. PPs are substantially the proportion of the time 
that a bipartition is found and should be considered as an approximation of the posterior probability of the bipartition.

Nonetheless, as we said before, 
different metrics have different pitfalls and PPs habe been often associated with an overestimation 
of support especially when data sets with a lot of characters are used and in the presence of small branchlenghts.
Generally speaking, standard measures of clade support, such as posterior probabilities 
and bootstrap proportions, can support several conflicting hypotheses with high support.
But this apparent confidence can hinder phylogenetic uncertainty, which can be observed when different metrics are compared
and trying to investigate the degree of reproducibility of ou results (different analytical frameworks and inference parameters).
You'll find a couple of interesting papers down below!

---

<br/>
<br/>

## conclusions: 

In this lesson we did quite a lot:

* explored a bit the command-line
* revised the different MSA formats
* used the more established and comprehensive online platform for phylogenetics
* did some hands on with BI and the relative post-inference diagnostics
* focused on how we should treat carefully nodal support metrics, such as PPs

---

<br/>
<br/>

## further reading: 

[Here](http://mrbayes.sourceforge.net/wiki/index.php/Manual) you will find the manual of MrBayes and some tutorials as well.

[Here](http://mrbayes.sourceforge.net/Help/) a very usefull help page for MrBayes.

A couple of nice papers on how the nodal support can be overestimated by [BP](https://doi.org/10.1016/j.ympev.2014.01.018) and [PP](https://doi.org/10.1080/10635150590924208). 