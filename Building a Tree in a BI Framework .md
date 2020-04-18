# Building a Tree in a BI Framework



 
---




## intro: 

In this lesson we will use MrBayes (cit) to understande the underlying basisi of Bayesian Inference in phylogenetics.

you should have allready a nexus (link), which is the file format which MrBayes requires.  
Otherwise you may use SeaView or Aliview to open your FASTA files and convert them into Phylip files.

But we have to carry out model selection again, mainly for two reasons:

* MrBayes supports slightly different models than IQ-TREE

* It's really painfull to translate the model from IQ-TREE language to MrBayes one.

So we will take advantage of this situation to leverage PartitionFinder2 which can carry out the model selection for MrBayes.

This is a lesson in usability and how sometimes there is a straight forward process (ModelFinder & IQ-TREE) otherwise is more fragmented
with multiple tools.

The infput of PF2 is a ```.cfg``` file . You can find the one I will use [here]()

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

12S = 1-506;
CO2nd = 2714-3384/3;

## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##
[schemes]

search = greedy;
```

Run .cfg on Cipress





---




Then we need to put the MrBayes block inside the nexus




---




Then we need to put the parameters for the search and run it.




---




While running, MrBayes will print the (estimated) remaining time at the end of each line
showing chain parameters. The MC3
data are printed to the .mcmc file in the same folder of the
NEXUS file. At the end of the run, you can check for convergence by plotting the last column of
that file: the burnin must be selected when the standard deviation of average split frequencies starts
to be (i) relatively very low and (ii) stable.


Go back to the NEXUS file and modify the last lines:




---




Tracer




---

Uncertainty can generally be observed either from scarce nodal support, polytomies  and/or  
sensitivity to the use of independent sampling of species and analytical frameworks (Yuan et al., 2016). 
Moreover, standard measures of clade support, such as posterior probabilities (Lewis et al., 2005) 
and bootstrap proportions (Simmons and Norton, 2014), can support several conflicting hypotheses with high apparent confidence.




---




## further reading: 

[Here](http://mrbayes.sourceforge.net/wiki/index.php/Manual) you will find the manual of MrBayes and some tutorials as well.