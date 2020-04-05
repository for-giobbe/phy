# Building a Tree in a ML Framework


---


## inferring gene trees

## ingerring species tree

## inferring nodal support


To overcome the computational burden required by the nonparametric bootstrap, IQ-TREE introduces an ultrafast bootstrap approximation (UFBoot) (Minh et al., 2013; Hoang et al., 2018) that is orders of magnitude faster than the standard procedure and provides relatively unbiased branch support values. 

Q-TREE provides an implementation of the SH-like approximate likelihood ratio test (Guindon et al., 2010). To perform this test, run:

iqtree -s example.phy -m TIM2+I+G -alrt 1000