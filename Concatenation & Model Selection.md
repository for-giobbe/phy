# Concatenation & Model Selection

## concatenation: 

The more frequent approach for inferring species trees is the use of the so-called "supermatrixes", in which several loci are "concatenated" in a single string, for each OTU.
Other approach are possible, such as inferring single gene trees and then reconciling them into a species tree using tools such as ASTRAL, but in this tutorial we will focus on the supermatrix approach.

After having aligned our genes, we need to concatenate them and we will do it using phyutility

## model of evolution & partitioning scheme selection : 

To carry out the model selection we will use ModelFinder [Kalyaanamoorthy et al., 2017](https://www.nature.com/articles/nmeth.4285).
As for most of the steps which are carried out to build a tree, there is a staggering choice of different tools and is sometimes difficult to decide which to use.
The more widespread tool for model selection is PartitionFinder2 but we are using Model finder for two key reason:

* usability
* speed

Also I think it's a good idea to keep constantly using new and shiny tools.

* -q partition_file: all partitions share the same set of branch lengths.
* -p partition_file: like above but allowing each partition to have its own evolution rate.
* -Q partition_file: each partition has its own set of branch lengths.

```iqtree -s example.phy -m MFP```


