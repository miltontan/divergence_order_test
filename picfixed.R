# This function calculates contrasts based on a vector of all tip and
# internal node trait values. Internals may be derived from any source
# (e.g., simulations, ancestral state estimates).
#
# This is a modification of the picfixed.R file written by David Ackerly,
# which was originally modified from the pic function in the R ape library
# written by Emmanuel Paradis paradis@isem.univ-montp2.fr
#
# David Ackerly, 2005, dackerly@berkeley.edu 


picfixed<-function (x, phy, scaled = TRUE, var.contrasts = FALSE) 
{
    if (class(phy) != "phylo") 
        stop("object \"phy\" is not of class \"phylo\"")
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    if (nb.node != nb.tip - 1) 
        stop("\"phy\" is not fully dichotomous")
    if (length(x) != nb.tip+nb.node) 
        stop("length of phenotypic and of phylogenetic data do not match")
    p <- as.numeric(rep(NA, nb.tip + nb.node))
    p <- x
    names(p) <- as.character(c(1:(nb.tip + nb.node)))
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
        pa <- p[names(p) == a]
        pb <- p[names(p) == b]
        if (scaled) 
            contr[count] <- (pa - pb)/sqrt(bl[i] + bl[j])
        else contr[count] <- pa - pb
        if (var.contrasts) 
            var.con[count] <- bl[i] + bl[j]
    }
    if (var.contrasts) 
        return(cbind(contr, var.con))
    else return(contr)
}