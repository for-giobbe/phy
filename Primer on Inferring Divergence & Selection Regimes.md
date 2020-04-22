# After Tree Building & Selection Regimes Inference



 
---




## intro: 

Phylogenetic systematics aims is to reconstruct the patterns of relationships among different organisms. 
And while phylogenetic trees can be considered by themselves the results of such investigations, they can also be
used to study evolutionary processes with a myriad of possibilities. This is not a comprensive list by any mean, 
just what I stumbled upon during my research:

* molecular evolution / selection regimes
* divergence time analyses
* mapping traits on phylogenies _i.e._ comparative methods

Due to the severe time constrains of the course, we just have time to deal with the analysis of selection regimes,
which may be the topic which is more fit for this course. But remember that the interaction of
molecular phylogenetics and paleontology & morphological studies is extremely exciting and
relies on very serious informatics skills as well. I have selected two papers (respectively on
[molecular dating](https://doi.org/10.1111/brv.12390) and [comparative methods]()), 
along with two excellent tutorials which you will find in the last section of this class.




---




## Inferring selection regimes: 

A powerful approach to molecular evolution comes from the comparison of 
the relative rates of synonymous and nonsynonymous substitutions (dN & dS).
Synonymous mutations do not change the amino acid sequence; hence their substitution rate (dS) is neutral 
with respect to selective pressure on the protein product of a gene.
Nonsynonymous mutations do change the amino acid sequence, so their substitution rate (dN) is a
function of selective pressure on the protein. 

The ratio of these rates (ω = dN/dS) is a measure of selective regime which a gene underwent. 
It's generally accepted that nonsynonymous mutations which are deleterious will be counterselected for
(purifying selection) and dN/dS will be less than 1, whereas if nonsynonymous mutations
are advantageous they will be fixed at a higher rate than synonymous mutations, and dN/dS will
be greater than 1. A dN/dS ratio equal to one is consistent with neutral evolution.


 maximum
likelihood estimation of the dN/dS ratio (Goldman and Yang 1994; Muse and Gaut 1994) and the
likelihood ratio test for positively selected genes (e.g., Nielsen and Yang 1998; Yang et al. 2000)
SCALED TO SYNONIMOUS SUBST.  It is calculated as the ratio of the number of nonsynonymous substitutions per non-synonymous site (Ka), in a given period of time, to the number of synonymous substitutions per synonymous site (Ks), in the same period. The latter are assumed to be neutral, so that the ratio indicates the net balance between deleterious and beneficial mutations. 



PAML [(Ziheng Yang, 2007)](https://academic.oup.com/mbe/article/24/8/1586/1103731) 
is an extremely versatile piece of software whose features include:

* estimating synonymous and nonsynonymous rates
* testing hypotheses concerning dN/dS rate ratios
* ancestral sequence reconstruction (DNA, codon, or AAs)
* various clock models
* simulating nucleotide, codon, or AA sequence data sets
* ...

Just an extensive reding of the [manual](http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf) will disclose all
of its possibilities! Several great and newer alternative exist, the most notable being HyPhy which is also
implemented in the Datamonkey [web server](https://www.datamonkey.org/).

The server is currently down because it's being fully used to study COVID-19: you can follow this effort [here](http://covid19.datamonkey.org/).
It includes an almost real-time update on which sites are found to be undergoing postive selection, which is quite [impressive](https://observablehq.com/@spond/natural-selection-analysis-of-sars-cov-2-covid-19).
This is not just some super-cool bioinformatics but it can prove very useful to identify new strains of the virus.


Anyway, let's check what we need in order to run CODEML:

* a ```.fasta``` alignment
* a ```.nwk``` tree
* a ```.ctl``` file 

Since the last tree we calculated for our PCG contains bootstrap we can quickly recalculate it using
```iqtree -s ND2_p_aligned.n.gb.fasta```. Anyway my[alignment](https://github.com/for-giobbe/phy/blob/master/examples/ND2_p_aligned.n.gb.fasta) and 
[tree]() are available.  Another possible option could be to to use something similar to ```sed "s/\)[0-9]/)/g"``` on the command-line.


Now we just miss the ```.ctl```, which stads for control and has a similar role to the ```.ctl``` we've seen earlier.
Here is how it looks like:

```
seqfile = 								* sequence data filename
treefile = 								* tree structure file name
outfile = 								* main result file name
noisy = 								* 0,1,2,3,9: how much rubbish on the screen
verbose = 								* 1:detailed output
runmode = 								* 0:user defined tree
seqtype = 								* 1:codons
CodonFreq = 								* 0:equal, 1:F1X4, 2:F3X4, 3:F61
model = 								* 0:one omega ratio for all branches
NSsites = 								* 0:one omega ratio (M0 in Tables 2 and 4)
icode = 								* 0:universal code
fix_kappa = 								* 1:kappa fixed, 0:kappa to be estimated
kappa = 								* initial or fixed kappa
fix_omega = 								* 1:omega fixed, 0:omega to be estimated
omega =  								* initial omega
 ```


In the ```.ctl``` found right up here just few parameters are present and it's possible to customize analyses
more deeply. In our first analysis we are going to compare two models - one with omega = 0.1 and one with omega 2.0 - 
using a Likelihood Ratio Test (LRT). 


First we need to modify a ```.ctl``` files to have ```omega = 0.1```, like this:

```
seqfile = ND2_p_aligned.n.gb.fasta					* sequence data filename
treefile = ND2_p_aligned.n.gb.fasta.treefile				* tree structure file name
outfile = ND2_omega_0.1.out						* main result file name
noisy = 9								* 0,1,2,3,9: how much rubbish on the screen
verbose = 1								* 1:detailed output
runmode = 0								* 0:user defined tree
seqtype = 1								* 1:codons
CodonFreq = 2								* 0:equal, 1:F1X4, 2:F3X4, 3:F61
model = 0								* 0:one omega ratio for all branches
NSsites = 0								* 0:one omega ratio (M0 in Tables 2 and 4)
icode = 4								* 0:universal code
fix_kappa = 0								* 1:kappa fixed, 0:kappa to be estimated
kappa = 2								* initial or fixed kappa
fix_omega = 1								* 1:omega fixed, 0:omega to be estimated
omega = 0.1 								* initial omega
```

which we can then save as ```ND2_omega_0.1.ctl```.


Then we need to have another ```.ctl``` with ```omega = 2```     
 
```
seqfile = ND2_p_aligned.n.gb.fasta					* sequence data filename
treefile = ND2_p_aligned.n.gb.fasta.treefile				* tree structure file name
outfile = ND2_omega_2.0.out						* main result file name
noisy = 9								* 0,1,2,3,9: how much rubbish on the screen
verbose = 1								* 1:detailed output
runmode = 0								* 0:user defined tree
seqtype = 1								* 1:codons
CodonFreq = 2								* 0:equal, 1:F1X4, 2:F3X4, 3:F61
model = 0								* 0:one omega ratio for all branches
NSsites = 0								* 0:one omega ratio (M0 in Tables 2 and 4)
icode = 4								* 0:universal code
fix_kappa = 0								* 1:kappa fixed, 0:kappa to be estimated
kappa = 2								* initial or fixed kappa
fix_omega = 1								* 1:omega fixed, 0:omega to be estimated
omega = 2 								* initial omega
```
 
and we can save it as ```ND2_omega_2.0.ctl```. 


Remember to give different names to the ```outfile``` parameter so that the outputs do not get overwritten. 
Also take care in selecting the correct gencode, otherwise it's likely that stop codons will be found and CODEML will give an error 
- I used 4 which is the correct one for the mitogenomes of arthropods but here is the full list:

* 0 for the universal code
* 1 for the mammalian mitochondrial code
* 3 for mold mt., 
* 4 for invertebrate mt.;
* 5 for ciliate nuclear code; 
* 6 for echinoderm mt.; 
* 7 for euplotid mt.;
* 8 for alternative yeast nuclear; 
* 9 for ascidian mt.; 
* 10 for blepharisma nuclear. 

You can also notice how we set ```fix_omega = 1``` in order to force an omega value; 
if you put it to 0 omega will be estimated with a ML search.

To run codeml all we need to do is type 'codeml' in the same folder that the codeml.ctl file is in;
all of the options are already in the codeml.ctl file, including all the other inputs required for the analysis. 

Let's take a look at the outputs:

* the ```.out``` file which we specified is what we are really interested in.
* ```2NG.t```, ```2NG.dS```, ```2NG.dN``` subs rates, dN & dS distances
* ```4fold.nuc``` 4-fold [degenerate](https://en.wikipedia.org/wiki/Codon_degeneracy) sites
* ```rub```, ```rst```, ```rst1``` and ```lnf``` are some rather misterios intermediate files

Here is the line we are interested in from the ```.out``` file:






Subsequently we can compare which was the best model. Of course this is just an example and
if we are really interested in the dNdS value which better describes our data, we should definitively
try more values. We can also make hypotheses different than omega, for example if the 
transition rates k or NSsites.

```
LRT <- -2*(likelihood.summary$Model1LnL-likelihood.summary$Model2LnL)
degrees.of.freedom <- likelihood.summary$Model2np-likelihood.summary$Model1np
p.value <- 1-pchisq(LRT,df=degrees.of.freedom)
adj.p.value <- p.adjust(p.value, method = "hochberg", n = length(p.value))
```







As we can see 0.1 is a much better fit for our data compared to 2, implying that the gene selection regime 
is better described by a strong purifying selection, instead of a positive selection regime. I wanna make you
focus ono a couple more things: in this kind of analyses dNdS values are averaged all throughout a sequence -
tipically a gene - but we know that selection can differ widely between different subsets, such as
functional domains or precise codons. In this perspective a value of 2 is highly unrealistic, as it would 
imply that every codon of the alignment is accumulating 2X coding mutation than non coding. I think
this values are very rarely associated to real positive selection events but can be used to detect
fenomena as pseudogenization and misalignment. The latter brings us back to the first lesson, where 
we learned that among the many ways of filtering out entire alignments dNdS was a possibility.







---




## further reading: 

The free [book](https://lukejharmon.github.io/pcm/) "Learning from trees" By Luke J Harmon, on comparative methods.

The "Taming the BEAST" [platform](https://taming-the-beast.org/) which includes many tutorials on divergence times estimation using BEAST 2.

[Here](http://evomicsorg.wpengine.netdna-cdn.com/wp-content/uploads/2011/08/bielawski_paml_review.pdf) you'll find more on the theoretical background of dNdS calculations. 