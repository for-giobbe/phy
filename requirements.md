## general requirements


---


[conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/#)

BASH - included natively in all Mac OS X and Linux distributions


---


# conda-installable packages (conda install -c bioconda package):

T-coffe

MAFFT

gblocks

IQ-TREE

MrBayes

EMBOSS

pal2nal

PAML


---


# executables:

[Aliview](https://github.com/AliView)

[PhyUtility](https://code.google.com/p/phyutility/downloads/list)

[PartitionFinder](http://www.robertlanfear.com/partitionfinder/)

[BEAST](tree.bio.ed.ac.uk/software/beast/)

[Tracer](http://tree.bio.ed.ac.uk/software/tracer/)

[FigTree](http://tree.bio.ed.ac.uk/software/figtree/)


---


# dataset building:

In the beginning we need to choose the group whose phylogeny will be inferred and the loci we will use for inferring their phylogeny. 
I suggest a lineage of animals and some mitochondrial markers, for which a lot of data are available. 
It is possible that you will need to refine your choice after checking about data availability.
On [NCBI](https://www.ncbi.nlm.nih.gov/) select “Nucleotide” from the drop-down menu at the top and use the group and the gene names as queries.
For a set of fifteen/twenty species of the group of interest you should download two /three gene sequences which are available for the most of them. 
Possibly, it would be interesting to have a PCG and a ncRNA. When you decide to download a given sequence, you have to place it in a [FASTA-formatted](https://en.wikipedia.org/wiki/FASTA_format) file.
Please ensure that the sequence names are the same in the two files, otherwise the subsequent step of concatenation will not work.
Finally, you have to select outgroup sequences. You can find an example of a starting file [here](https://raw.githubusercontent.com/for-giobbe/phy/master/unaligned_genes/CO1_total.fasta).