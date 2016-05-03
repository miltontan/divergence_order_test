# This function is modified from the sigma.d function provided by Adams 2014 
# (DOI: 10.1093/sysbio/syt105). This modified sigma.d does not allow specifying
# separate groups of taxa (as originally scripted), and simply calculates sigma.d.all
# for all taxa.
sigma.d<-function(phy,x){
p<-ncol(x)
x<-as.matrix(x)
x<-prcomp(x)$x
N<-length(phy$tip.label)
ones<-array(1,N)
C<-vcv.phylo(phy)
C<-C[rownames(x),rownames(x)]
a.obs<-colSums(solve(C))%*%x/sum(solve(C))
eigC<-eigen(C)
D.mat<-solve(eigC$vectors %*% diag(sqrt(eigC$values))%*%
t(eigC$vectors))
dist.adj<-as.matrix(dist(rbind((D.mat%*%
(x-(ones%*%a.obs))),0)))
vec.d2<-dist.adj[N+1,1:N]^2
sigma.d.all<-sum(vec.d2)/N/p
return(sigma.d.all)
}
