# MP_bioinfo_2020
*Forni Giobbe*


---



Syllabus:




0. [requirements](https://github.com/for-giobbe/phy/blob/master/requirements.md)

1. [Multiple Sequence Alignement & Filtering](https://github.com/for-giobbe/phy/blob/master/Multiple%20Sequence%20Alignement%20%26%20filtering.md)

2. [Concatenation & Model Selection](https://github.com/for-giobbe/phy/blob/master/Concatenation%20%26%20Model%20Selection.md)

3. [Building a Tree in a ML Framework](https://github.com/for-giobbe/phy/blob/master/Building%20a%20Tree%20in%20a%20ML%20Framework.md)

4. [Building a Tree in a BI Framework](https://github.com/for-giobbe/phy/blob/master/Building%20a%20Tree%20in%20a%20BI%20Framework%20.md)

5. [Post Phylogenetic Inference & Selection Regimes Inference](https://github.com/for-giobbe/phy/blob/master/Primer%20on%20Inferring%20Divergence%20%26%20Selection%20Regimes.md)



---



_arguments which we ignored due to time constrains:_

* [Orthology inference](https://github.com/davidemms/OrthoFinder)

---

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
 
**NB** Conda installation of CAFE can fails. In that case you should download and compile it manually. Follow instructions [here](https://github.com/hahnlab/CAFE/blob/master/docs/cafe_manual.pdf). If compilation fails try to add

```
#include <cmath>
```
to cafe/cafe_commands.cpp

If is still failing try to search a solution in their [ghithub page](https://github.com/hahnlab/CAFE/issues?q=is%3Aissue+) - at the end you are bioinformaticians ;). 

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
If you encounter any problem you can try to install all dependecies for each envoirnment while creating it. Download [this] and [this] files for respectively the python 3.6 and 2.7 envoirnments, then run the following codes:

```
conda create --name bio_info-env --file spec-file_bio_info-env.txt
conda create --name PF-env --file spec-file_PF-envs.txt
```

>Remember to install manually PartitionFinder and CAFE if you go directly for this solution (the latest one is not present in the spec-file).

 if you get stuck feel free to ask for help both through Teams or e-mail (jacopo.martelossi2@unibo.it).

## Table of content

 * 1. Orthology inference, MSA and masking of the alignment.
 * 2. Model selection, Maximum-likelihood inference and divergence time estimation in a Maximum-likelihood framework.
 * 3. Bayesian inference.
 * 4. Gene family evolutionary analyses +  hints of natural selection regimes detection.
