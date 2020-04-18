After Tree Building & Selection Regimes Inference



 
---




## intro: 

Phylogenetic systematics aims is to reconstruct the patterns of relationships among different organisms. 
And while phylogenetic trees can be considered by themselves the reults of such investigations they can also be
used to study evolutionary processes with a myriad of possibilities. This is not a comprensive list by any mean, 
just what I stumbled upon during my research:

* molecular evolution / selection regimes
* divergence time analyses
* mapping traits on phylogenies (_i.e._ comparative methods)

Due to the severe time constrains of the course, we just have time to deal with the analysis of selection regimes,
which may be the topic which is more fit for this course. But remember that the interaction of
molecular phylogenetics and palaeontology & morphological studies is extremely exciting and
relies on very serious informatics skills as well. I have selected two papers (respectively on
[molecular dating](https://doi.org/10.1111/brv.12390) and [comparative methods]())
, along with two excellent tutorials which you will find in the last section of this class.




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

The server is currently down because it's being fully used to study : you can follow this effort [here](http://covid19.datamonkey.org/)
It includes an almost real-time update on which sites are found to be undergoing postive selection, which is quite [amazing](http://covid19.datamonkey.org/2020/04/01/covid19-analysis/)


SCALED TO SYNONIMOUS SUBST.  It is calculated as the ratio of the number of nonsynonymous substitutions per non-synonymous site (Ka), in a given period of time, to the number of synonymous substitutions per synonymous site (Ks), in the same period. The latter are assumed to be neutral, so that the ratio indicates the net balance between deleterious and beneficial mutations. 


Running PAML programs
* Sequence data file
* Tree file
* Control file (*.ctl) 


We should already have an alignement file and a tree file from our previous lessons.

To run codeml all we need to do is type 'codeml' in the same folder that the codeml.ctl file is in and all 
of the options are in the codeml.ctl file. Here how a codeml .ctl file looks like:

```
seqfile = seqfile.txt * sequence data filename
treefile = tree.txt * tree structure file name
outfile = results.txt * main result file name
noisy = 9 * 0,1,2,3,9: how much rubbish on the screen
verbose = 1 * 1:detailed output
runmode = 0 * 0:user defined tree
seqtype = 1 * 1:codons
CodonFreq = 2 * 0:equal, 1:F1X4, 2:F3X4, 3:F61
model = 0 * 0:one omega ratio for all branches
NSsites = 0 * 0:one omega ratio (M0 in Tables 2 and 4)
* 1:neutral (M1 in Tables 2 and 4)
* 2:selection (M2 in Tables 2 and 4)
* 3:discrete (M3 in Tables 2 and 4)
* 7:beta (M7 in Tables 2 and 4)
* 8:beta&w; (M8 in Tables 2 and 4)
icode = 0 * 0:universal code
fix_kappa = 0 * 1:kappa fixed, 0:kappa to be estimated
kappa = 2 * initial or fixed kappa
fix_omega = 0 * 1:omega fixed, 0:omega to be estimated
omega = 5 * initial omega
*set ncatG for models M3, M7, and M8!!!
*ncatG = 3 * # of site categories for M3 in Table 4
*ncatG = 10 * # of site categories for M7 and M8 in Table 4
 ```








---




## further reading: 

[Here](http://evomicsorg.wpengine.netdna-cdn.com/wp-content/uploads/2011/08/bielawski_paml_review.pdf) you'll find more on the theoretical background of dNdS calculations. 