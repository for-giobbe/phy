# Building a Tree in a BI Framework


---


## intro: 

In this class, you will use MrBayes to infer the Bayesian phylogeny of your dataset. MrBayes is
really powerful and many options are provided: please refer to the software manual for further
details. At first, you need to convert your dataset into the NEXUS format; again, several options are
available, and three are shown here.

You may use SeaView or Aliview to open your FASTA files and convert them into Phylip files.


While running, MrBayes will print the (estimated) remaining time at the end of each line
showing chain parameters. The MC3
data are printed to the .mcmc file in the same folder of the
NEXUS file. At the end of the run, you can check for convergence by plotting the last column of
that file: the burnin must be selected when the standard deviation of average split frequencies starts
to be (i) relatively very low and (ii) stable.


Go back to the NEXUS file and modify the last lines: