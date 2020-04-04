# Multiple Sequence Alignment & Filtering

---

## intro: 

All phylogenetic methods require a set of homologous characters thus the first step is to infer which nucleotides / codons / amminocid residues are homologous to each other, so that differences among these nucleotides result only from changes that convey descent information. This process is called  Multiple Sequence Alignment (MSA) and is often followed by an additional step of detection and exclusion of those alignment regions whose homology we are uncertain of.
Due to the large amount of data which we process nowadays, this step is often overlooked and as a result it's quite easy to find misaligned loci in modern datasets. Nonetheless this is a key step, which may affect downstream analyses (see further reading for more).

In this tutorial, we will use the most popular tools for multiple sequence alignment:

MAFFT [(Katoh and Standley 2013)](https://academic.oup.com/mbe/article/30/4/772/1073398)

T-Coffee[(Notredame et al., 2000)](https://www.ncbi.nlm.nih.gov/pubmed/10964570)

And subsequently filter the alignment using: 

Gblocks [(Talavera & Castresana, 2007)](https://academic.oup.com/sysbio/article/56/4/564/1682121)

---

## aligning PCGs with MAFFT

At this stage you should have assembled your dataset into a multi-fasta file, which should look like [this](https://raw.githubusercontent.com/for-giobbe/phy/master/unaligned_genes/12S_total.fasta). It's very important that each one of you has it's own dastaset as this may provide natural examples of problems in analyses, but If you still do not have a dataset use [this](https://github.com/for-giobbe/phy/tree/master/unaligned_genes).

As we have a handfull of loci, the first step in MSA should be to think about their characteristics: some genes will be protein-coding genes (PCGs) and other will be non-coding (as rRNAs).

MAFFT is a very fast tool, which implements several different algorithms which makes it also very versatile.

Let's start with a PCG, using the code:

```mafft --auto CO1.fasta > CO1_aligned.fasta```

we can then inspect the outcome using Aliview: Let's take some time to see if there is any stop codon.
Sometimes aligners which are not aware of the underlying structure of codons in PCGs can introduce some artificial stop codons, which are more often than not alignment errors. 
To avoid this pitfall we can translate our PCG from nucleotides into aminoacid residues, align aminoacids and then retrotranslate the alignment to nucleotide:

1. translation from nt to aa
```transeq -sequence $n.n.fasta -outseq $n.p.fasta```
2. aa alignment
``` mafft --auto CO1.fasta > CO1_aligned.fasta ```
3. retrotranslation of the alignment from aa to nt
``` pal2nal.pl $n.mafft.p.aln $n.n.fasta -output fasta >> tmp1.aln ```

---

## aligning non-coding genes with MAFFT


---

## one aligner to rule them all T-coffe: 

T-Coffe has an incredibly nice server[http://tcoffee.crg.cat/] which I strongly encourage to explore.

---

## the tradeoff between speed and accuracy: 

PSIcoffe
---

## automation 

T-Coffe has an incredibly nice server[http://tcoffee.crg.cat/] which I strongly encourage to explore.

---

```mkdir aligned_genes;

for i in unaligned_genes/*fasta; do 
	gene_name=$(echo $i | awk -F “.” ‘{print $1}’); 
	mafft --auto $i > aligned_genes/$name_aligned.fasta;
	done
```

---

## further reading: 

[Evidence of Statistical Inconsistency of Phylogenetic Methods in the Presence of Multiple Sequence Alignment Uncertainty](https://academic.oup.com/gbe/article/7/8/2102/556628)
[Alignment Uncertainty and Genomic Analysis](https://science.sciencemag.org/content/319/5862/473?casa_token=t07ptffISm4AAAAA:j5l4US_y_GHOMduYw6R-MhyM7YUpa__08lw45l455DAU3tGFNKYlV40ZH0Si5w48Xl1gTEqsocLVvaE)

