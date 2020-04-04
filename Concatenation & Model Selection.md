# Concatenation & Model Selection

## concatenation: 

Currently, two divergent systematic methods are commonly applied for inferring species trees: the supermatrix approach (concatenation) and the coalescent approach (gene trees and the species tree are not co-estimated). 
You will find an interesting paper on the topic at the end of this tutorial, but here we will focus on the supermatrix approach. 
After having aligned our genes, we will concatenate them using phyutility: the software has a wide array of functions and can be considere a swiss-knife for phylogeneticists.
As you always should, when using a new software, take a look at its [manual](https://github.com/blackrim/phyutility/blob/master/manual.pdf). 

After adjusting the path to the java executable of phyutility, try this string:

```java -jar /Applications/bio/phyutility/phyutility.jar -concat -in *.fasta -out concatenation.nxs```

Then take a look at the outputs:

Get familiar with the nexus [format](http://informatics.nescent.org/wiki/NEXUS_Specification).

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

---

## further reading: 

[Here](http://www.iqtree.org/doc/Tutorial) you will find great tutorials from the authors themselves on ModelFinder and IQ-Tree.

[Resolution of a concatenation/coalescence kerfuffle: partitioned coalescence support and a robust family‐level tree for Mammalia](https://onlinelibrary.wiley.com/doi/full/10.1111/cla.12170?casa_token=X0ctrSm4S1AAAAAA%3AgiB9v0MtJDO6vMWOigdvW9JrgYuJTebMen6zYxg9S0nP8MWIi2zA2fwWfi-lJlMCD9Ir1MDCzkBeyVwg)