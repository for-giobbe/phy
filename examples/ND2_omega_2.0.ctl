seqfile = ND2_p_aligned.n.gb.fasta					* sequence data filename
treefile = ND2_p_aligned.n.gb.fasta.treefile				* tree structure file name
outfile = ND2_omega_2.0.out						* main result file name
noisy = 9								* 0,1,2,3,9: how much rubbish on the screen
verbose = 1								* 1:detailed output
runmode = 0								* 0:user defined tree
seqtype = 1								* 1:codons
CodonFreq = 2								* 0:equal, 1:F1X4, 2:F3X4, 3:F61
model = 0								* 0:one omega ratio for all branches
NSsites = 0								* 0:one omega ratio (M0 in Tables 2 and 4)
icode = 4								* 0:universal code
fix_kappa = 0								* 1:kappa fixed, 0:kappa to be estimated
kappa = 2								* initial or fixed kappa
fix_omega = 1								* 1:omega fixed, 0:omega to be estimated
omega = 2 								* initial omega
