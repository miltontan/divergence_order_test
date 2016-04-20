### The DOT Test ###

Original script by D.D. Ackerly and C.O. Webb
Ecological Archives E087-110-S1
Available here: http://esapubs.org/archive/ecol/E087/110/suppl-1.htm

Please read the divergence-order test paper before conducting the test:
D. D. Ackerly, D. S. Schwilk, and C. O. Webb. 2006. Niche evolution and adaptive radiation: testing the order of trait divergence. Ecology 87:S50–S61.

### Modifications from the original script ###
The main modification to this script is the removal of the dependency on ANCML (Schluter et al. 1997), replaced by ace implemented in ape. Other modifications include storing the DOT test as a function that accepts settings (as opposed to requiring modification of the script directly), edits for how tip.label and nodes are specified, and removing the BATCH option and printing all 22 variables (last 9 omitted if zero branch lengths are present) by default.

Note, you may receive warnings when running the script due to NaNs produced during the maximum-likelihood ancestral state reconstruction.

### Instructions ###

1. The R libraries 'ape' and 'msm' are required. Running dot.R and picfixed.R will store dot and picfixed as functions within R.
2. The dot function takes the following commands: 

dot(tree, x, y, nsim = 0, replace=FALSE)

Where "tree" is a phylo object, "x" and "y" are named vectors for the traits of interest (as accepted by ace), "nsim" is the number of bootstrap replicates to simulate (eg. 1000), and "replace" is a logical for whether to sample with replacement in null models.

### Output ###
The dot function returns a matrix of a single row and with named columns for each of the variables. Below README text is copied from the original README of the DOT Test, describing the output variables of the dot function.

Nreps: Number of bootstrap and randomization replicates

Wx.as: Mean node age for trait X, calculated directly from ancestral states

Wy.as: Same for trait Y

D.as: Age difference (Wx.as-Wy.as)

Wx.bs: Mean node age for trait X, averaged across all bootstrap reps

Wx.sd: Standard deviation of node age for X, across the bootstrap reps

Wy.bs: Same for trait Y

Wy.sd: Standard deviation of node age for X, across the bootstrap reps

D.bs: Mean age difference between traits, averaged across all bootstrap reps (Wx.bs-Wy.bs)

D.sd: Standard deviation of the difference in node age, averaged across all bootstrap reps

P.bs: Significance of one-tailed test that D.bs > 0. For D < 0, p = 1-P.bs. For two-tailed test, take lesser of P.bs and (1-P.bs) and divide by 2.

betaX, betaY: Brownian motion rate parameters calculated by ancml for each trait

Wx.pic: Mean node age for trait X, calculated from Felsenstein's independent contrast algorithm

Wy.pic: Same for trait Y

D.pic: Difference in node age, based on Felsenstein contrasts

P.tips: Significance of one-tailed test that D.pic > 0, based on randomization of tip values (see above for D.pic < 0 and two-tailed tests)

P.contrasts: Significance of one-tailed test that D.pic > 0, based on randomization of contrasts across nodes (see above for D.pic < 0 and two-tailed tests)

BMtestX, BMtestY: Correlation of the absolute value of standardized independent contrasts for trait X with the standard deviation of the contrast (square root of subtending branch lengths). This should not be significantly differen from zero (see Garland et al. 1992)

NormTestX, NormTestY: Correlation coefficient for truncated normal probability plot of standardized contrasts of traits X and Y. If trait fits Brownian motion model, the absolute value of the contrasts should exhibit a truncated normal distribution. The normal probability plot would then be a straight line, and this correlation should be close to 1. If it is much lower than 1, the plot should be examined graphically to determine the departures from normality. After running dot.R, this can be done in R with the command: plot(qtnorm(qq,lower=0),sort(abs(Fx[,1]))) (substitute Fy for trait Y).
