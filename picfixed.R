# This function calculates contrasts based on a vector of all tip and
# internal node trait values. Internals may be derived from any source
# (e.g., simulations, ancestral state estimates).
#
# This is a modification of the picfixed.R file written by David Ackerly,
# which was originally modified from the pic function in the R ape library
# written by Emmanuel Paradis paradis@isem.univ-montp2.fr
#
# This script is also extended to calculate multivariate contrasts
# based on McPeek et al. 2008 (DOI: 10.1086/587076).
#
# Original script by David Ackerly, 2005, dackerly@berkeley.edu 
# Modified by Milton Tan, 2016

picfixed<-function (x, phy, scaled = TRUE, var.contrasts = FALSE) 
{
    if (class(phy) != "phylo") 
        stop("object \"phy\" is not of class \"phylo\"")
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    if (nb.node != nb.tip - 1) 
        stop("\"phy\" is not fully dichotomous")
    if (length(as.matrix(x)[,1]) != nb.tip+nb.node) 
        stop("length of phenotypic and of phylogenetic data do not match")
    p <- as.matrix(x)
    rownames(p) <- as.character(c(1:(nb.tip + nb.node)))
    bl <- phy$edge.length
    contr <- as.numeric(rep(NA, nb.node))
    names(contr) <- as.character(unique(phy$edge[, 1]))
    if (var.contrasts) 
        var.con <- contr
    int <- as.numeric(unique(phy$edge[, 1]))
    count <- 0
    for (nod in int) {
        count <- count + 1
        pair.ind <- which(phy$edge[, 1] == nod)
        i <- pair.ind[1]
        j <- pair.ind[2]
        pair <- phy$edge[pair.ind, 2]
        a <- pair[1]
        b <- pair[2]
        if (var.contrasts) 
            var.con[count] <- bl[i] + bl[j]
        if (ncol(p) == 1) {
            pa <- p[rownames(p) == a]
            pb <- p[rownames(p) == b]
            if (scaled) 
                contr[count] <- (pa - pb)/sqrt(bl[i] + bl[j])
            else contr[count] <- pa - pb
        } else {
            z <- NULL
            for (i in 1:ncol(p)) {
                pa <- p[rownames(p) == a,i]
                pb <- p[rownames(p) == b,i]
                z <- c(z,pa-pb)
                }
            if (scaled) 
                contr[count] <- sqrt(sum(z^2))/sqrt(bl[i] + bl[j])
            else contr[count] <- sqrt(sum(z^2))
        }
    }
    if (var.contrasts) 
        return(cbind(contr, var.con))
    else return(contr)
}
