
# MP_bioinfo_2021
*Martelossi Jacopo*.  
Teaching assistantship for the molecular phylogenetics classes of the Master in Bioinformatics at the University of Bologna.


## REQUIREMENTS

All packages can be installed using the conda package manager - look [here](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf) for a nice cheet-sheet. What you have to do is just create a new conda environmet and install all dependencies inside it (google the name of the package + conda). 

#### 1st environment: python 3 (deafult conda settings)

```
conda create --name bio_info-env
conda activate bio_info-env
conda install -c bioconda <packagename>
```

  * **Orthofinder**: Conda installation of Orthofinder already provide a lot of usefull packages ([Orthofinder](https://github.com/davidemms/OrthoFinder), [MCL](https://orthomcl.org/orthomcl/app), [raxml](https://cme.h-its.org/exelixis/web/software/raxml/), [fasttree](http://www.microbesonline.org/fasttree/), diamond, [iqtree](http://www.iqtree.org/), [mafft](https://mafft.cbrc.jp/alignment/server/)).
  * **[gblocks](https://mafft.cbrc.jp/alignment/server/)**.
  * **[T-Coffee](http://www.tcoffee.org/Projects/tcoffee/index.html#DOWNLOAD)**. NB: Use [this](https://anaconda.org/bioconda/t-coffee) link.
  * **[AMAS](https://pypi.org/project/amas/)**.
  * **[CAFE](https://hahnlab.github.io/CAFE/manual.html)**
  * **[kalign2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647288/)**
 
**NB** Conda installation of CAFE can fails. In that case you should download and compile it manually. Follow instructions [here](https://github.com/hahnlab/CAFE/blob/master/docs/cafe_manual.pdf). If compilation fails try to add

```
#include <cmath>
```
to cafe/cafe_commands.cpp

If is still failing try to find a solution in their [ghithub page](https://github.com/hahnlab/CAFE/issues?q=is%3Aissue+) - at the end you are bioinformaticians ;). 

#### 2nd environment: python 2.7 (necessary for PartitionFinder)

Creation of an environment with all dependencies necessary for PartitionFinder:

```
conda create --name PF-env python=2.7
conda activate PF-env
conda install numpy pandas pytables pyparsing scipy scikit-learn
```

Download and install the latest version of [PF](https://apolo-docs.readthedocs.io/en/latest/software/applications/partitionFinder/2.1.1/):

```
wget https://github.com/brettc/partitionfinder/archive/v2.1.1.tar.gz
tar xfz v2.1.1.tar.gz
chmod +x partitionfinder-2.1.1/PartitionFinder.py
```
If you encounter any problem you can try to install all dependecies for each envoirnment while creating it. Download [this](https://raw.githubusercontent.com/for-giobbe/phy/master/2021/environments/spec-file_bio_info-env.txt) and [this](https://raw.githubusercontent.com/for-giobbe/phy/master/2021/environments/spec-file_PF-envs.txt) files for respectively the python 3.6 and 2.7 envoirnments, then run the following codes:

```
conda create --name bio_info-env --file spec-file_bio_info-env.txt
conda create --name PF-env --file spec-file_PF-envs.txt
```

>Remember to install manually PartitionFinder and CAFE if you go directly for this solution (the latest one is not present in the spec-file).

#### Super light executables

 * [Aliview](https://github.com/AliView/AliView)
 * [FigTree](http://tree.bio.ed.ac.uk/software/figtree/)

**If you get stuck feel free to ask for help both through Teams and e-mail (jacopo.martelossi2@unibo.it).**

> Small note: As bioinformatics maybe some of you are going to develop softwares during their scientific carrier. For not-bioinformatician people (like me), installation of softwares is one of the biggest obstacles in using a specific, maybe interesting, tool. So please, try to implement the *usability* of conda.

## Table of content

 * **1. Orthology inference, MSA and masking of the alignment.**
 <br/><br/>
 * **2. Model selection, Maximum-likelihood inference and divergence time estimation in a Maximum-likelihood framework.**
 <br/><br/>
 * **3. Bayesian inference.**
 <br/><br/>
 * **4. Gene family evolutionary analyses +  hints of natural selection regimes detection.**
