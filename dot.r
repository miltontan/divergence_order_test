# This script reads in a phylogeny and trait values, calculates ancestral 
# trait values with CIs for two traits, and node ages. Then it bootstraps 
# the distribution of ancestral values and calculates the weighted divergence 
# age for each trait, and the difference between these values. This provides a 
# p-value on the DOT statistic.
#
# Modified from script dot.R by David Ackerly 2005
#
# Milton Tan
# Unlike Ackerly's original script, this function only requires ape.
# Ancestral state reconstruction is herein performed in ace
# 
# This script is also extended to calculate the D-statistic given multivariate  
# data (eg. multiple axes of shape). Multivariate contrasts are calculated based 
# on McPeek et al. 2008 (DOI: 10.1086/587076), multivariate ancestral states are 
# calculated based on geomorph scripts (Adams & Otarola-Castillo 2013), and multvariate 
# rates are calculated using a modified script based on sigma.d by Adams 2014
# (DOI: 10.1093/sysbio/syt105). The modified sigma.d does not allow specifying
# separate groups, and calculates sigma.d for the entire dataset.
#
# Note, this script requires ape and msm to be loaded
# It also requires picfixed and sigma.d loaded as functions in the environment

dot <- function (tree, x, y, nsim = 0, replace=FALSE)
{
# Reformat objects if necessary
	x <- as.matrix(as.matrix(x)[tree$tip.label, ])
	y <- as.matrix(as.matrix(y)[tree$tip.label, ])
# Set flag BLzero = 1 if tree has zero length branches; = 0 if not
	ifelse(min(tree$edge.length) == 0,
	{
	BLzero <- 1
	print('WARNING: tree has zero length branches')
	},
	BLzero <- 0)
	
# tree dimensions
	NTips <- length(tree$tip.label)
	NTot <- 2*NTips-1
	Ncont <- NTips-1
	qq<-((1:Ncont)-0.5)/(Ncont) #quantiles for truncated normal prob plots

# Run ace and sigma.d
	if (ncol(x) == 1) {
	ace_x <- ace(as.numeric(x),tree)
	ace_x_ace <- as.matrix(ace_x$ace)
	ace_x_se <- as.matrix((ace_x$CI95[,2]-ace_x$CI95[,1])/3.92)
	beta1 <- ace_x$sigma2[1]
	} else {
	ace_x_ace <- NULL
	ace_x_se <- NULL
	for (i in 1:ncol(x)) {
        x1 <- x[, i]
        ace_tmp<-ace(x1, tree, method = "ML")
        tmp <- ace_tmp$ace
        tmp2 <- (ace_tmp$CI95[,2]-ace_tmp$CI95[,1])/3.92
        ace_x_ace <- cbind(ace_x_ace, tmp)
        ace_x_se <- cbind(ace_x_se, tmp2)
    }
    beta1 <- sigma.d(tree,x)
    }
    if (ncol(y) == 1) {
	ace_y <- ace(as.numeric(y),tree)
	ace_y_ace <- as.matrix(ace_y$ace)
	ace_y_se <- as.matrix((ace_y$CI95[,2]-ace_y$CI95[,1])/3.92)
	beta2 <- ace_y$sigma2[1]
	} else {
	ace_y_ace <- NULL
	ace_y_se <- NULL
	ace_y <- NULL
	for (i in 1:ncol(y)) {
        y1 <- y[, i]
        ace_tmp<-ace(y1, tree, method = "ML")
        tmp <- ace_tmp$ace
        tmp2 <- (ace_tmp$CI95[,2]-ace_tmp$CI95[,1])/3.92
        ace_y_ace <- cbind(ace_y_ace, tmp)
        ace_y_se <- cbind(ace_y_se, tmp2)
    }
    beta2 <- sigma.d(tree,y)
    }

# Parse ages, vars, ses
	ages <- as.numeric(branching.times(tree))
	X <- rbind(x,ace_x_ace)
	Xse <- rbind(matrix(0,nrow=nrow(x),ncol=ncol(x)),ace_x_se)
	Y <- rbind(y,ace_y_ace)
	Yse <- rbind(matrix(0,nrow=nrow(y),ncol=ncol(y)),ace_y_se)

#calculation of mean ages based on ace contrasts 
# ACx is vector of absolute contrasts
# ACageX = Wx in text of paper
	ACx <- abs(picfixed(X,tree,scaled=FALSE))
	ACageX <- weighted.mean(ages,ACx)

	ACy <- abs(picfixed(Y,tree,scaled=FALSE))
	ACageY <- weighted.mean(ages,ACy)

# Save observed contrasts, before running bootstrap reps
	ACxObs <- ACx
	ACyObs <- ACy

#calculation of mean ages based on Felsenstein contrasts
ifelse (BLzero == 0, 
	{
	# Contrasts based on Felsenstein's algorithm are used
	# to check assumptions of Brownian motion, and for significance
	# testing by alternative null models discussed in online archives
	if (ncol(x) == 1) {
	Cx <- abs(pic(x,tree,scaled=FALSE))
	} else {
	Cx <- sqrt(apply(apply(x,2,pic,phy=tree,scaled=FALSE)^2,1,sum))
	}
	if (ncol(y) == 1) {
	Cy <- abs(pic(y,tree,scaled=FALSE))
	} else {
	Cy <- sqrt(apply(apply(y,2,pic,phy=tree,scaled=FALSE)^2,1,sum))
	}
	CageX <- weighted.mean(ages,Cx)
 	CageY <- weighted.mean(ages,Cy)
	CxObs <- Cx
	CyObs <- Cy
	if (ncol(x) == 1) {
	Fx <- pic(x,tree,var.contrasts=TRUE)
	} else {
	Fx <- pic(x[,1],tree,var.contrasts=TRUE)
	Fx[,1] <- sqrt(apply(apply(x,2,pic,phy=tree,scaled=FALSE)^2,1,sum))/sqrt(Fx[,2])
	}
	if (ncol(y) == 1) {
	Fy <- pic(y,tree,var.contrasts=TRUE)
	} else {
	Fy <- pic(y[,1],tree,var.contrasts=TRUE)
	Fy[,1] <- sqrt(apply(apply(y,2,pic,phy=tree,scaled=FALSE)^2,1,sum))/sqrt(Fy[,2])
	}
	# Check if standardized contrasts are correlated with their SD
	BMx <- cor(sqrt(Fx[,2]),abs(Fx[,1]))
	BMy <- cor(sqrt(Fy[,2]),abs(Fy[,1]))
	# Check R value of truncated normal probability plot
	# low values indicate deviations from straight line
	PPx <- cor(qtnorm(qq,lower=0),sort(abs(Fx[,1])))
	PPy <- cor(qtnorm(qq,lower=0),sort(abs(Fy[,1])))
	# Store results
	result <- c(ACageX,ACageY,CageX,CageY,CageX,CageY,BMx,BMy,PPx,PPy)
	},
	# If zero branch lengths present, only store results from ACx
	result <- c(ACageX,ACageY))

# Store cumulative sums and sums of squares to calculate average 
# contrast values, and sd, over bootstrap samples
	cumACx <- ACx
	cumACy <- ACy
	cumACx2 <- ACx*ACx
	cumACy2 <- ACy*ACy

#create variables for bootstrap samples
	Xr<-X
	Yr<-Y

#run bootstrap samples and calculate agediff
	Nreps<-nsim
	for (r in c(1:Nreps)) 
	{
		# bootstrap of ace contrasts
		for (n in c((NTips+1):NTot)) 
		{
			for (i in 1:ncol(x)) {
			Xr[n,i] <- rnorm(1,X[n,i],Xse[n,i])
			Yr[n,i] <- rnorm(1,Y[n,i],Yse[n,i])
			}
		}
		
		ACx <- abs(picfixed(Xr,tree,scaled=FALSE))
		RageX <- cor(ACx,ages)
		mnageX <- weighted.mean(ages,ACx)
	
		ACy <- abs(picfixed(Yr,tree,scaled=FALSE))
		RageY <- cor(ACy,ages)
		mnageY <- weighted.mean(ages,ACy)
		
		ifelse (BLzero == 0, 
		{
		# tip randomization and contrast randomization here
		if (ncol(x) == 1) {
		Cx <- abs(pic(sample(x,replace=replace),tree,scaled=FALSE))
		} else {
		Xrs<-apply(x,2,sample,replace=replace)
		rownames(Xrs)<-rownames(x)
		Cx <- pic(Xrs[,1],tree,var.contrasts=TRUE)
		Cx[,1] <- sqrt(apply(apply(Xrs,2,pic,phy=tree,scaled=FALSE)^2,1,sum))/sqrt(Cx[,2])
		}
		if (ncol(y) == 1) {
		Cy <- abs(pic(sample(y,replace=replace),tree,scaled=FALSE))
		} else {
		Yrs<-apply(y,2,sample,replace=replace)
		rownames(Yrs)<-rownames(y)
		Cy <- pic(Yrs[,1],tree,var.contrasts=TRUE)
		Cy[,1] <- sqrt(apply(apply(Yrs,2,pic,phy=tree,scaled=FALSE)^2,1,sum))/sqrt(Cy[,2])
		}
		CageX1 <- sum(Cx*ages)/sum(Cx)
	 	CageY1 <- sum(Cy*ages)/sum(Cy)
		# randomize contrasts
		Cx <- sample(CxObs,replace=replace)
		Cy <- sample(CyObs,replace=replace)
		CageX2 <- sum(Cx*ages)/sum(Cx)
	 	CageY2 <- sum(Cy*ages)/sum(Cy)
		rand <- c(mnageX,mnageY,CageX1,CageY1,CageX2,CageY2,BMx,BMy,PPx,PPy)
		}, #end BLzero==0
		rand <- c(mnageX,mnageY))
	
		result <- rbind(result,rand)
		
		cumACx <- cumACx + ACx
		cumACy <- cumACy + ACy
		cumACx2 <- cumACx2 + ACx*ACx
		cumACy2 <- cumACy2 + ACy*ACy
	}

# Add one for observed data
	Nreps<-Nreps+1

#calculate differences in mean age under bootstrap
	agediff <- as.array(result[,1]-result[,2])
	dim(agediff) <- c(length(agediff),1)
	colnames(agediff) <- 'bootstrap'

#calculate one-tailed probability that ace agediff is greater than zero
	pAnc <- length(agediff[agediff<=0])

#calculate differences in age under null models
	if (BLzero == 0) {
		pRand <- NULL
		for (i in 2:3) { 
			agediff <- cbind(agediff,result[,2*i-1]-result[,2*i]) 
			pRand <- cbind(pRand,rank(agediff[,i])[1])  
		}
		colnames(pRand) <- c('pTips','pContrasts')
		colnames(agediff)[2:3] <- c('tip.null','contrast.null')
	}

# Generate output
	out <-  cbind(Nreps,result[1,1],result[1,2],
		agediff[1,1],mean(result[,1]),sd(result[,1]),mean(result[,2]),
		sd(result[,2]),mean(agediff[,1]),sd(agediff[,1]),pAnc,
		beta1,beta2)
	colnames(out) <- c("Nreps","Wx.as","Wy.as","D.as","Wx.bs",
			 "Wx.sd","Wy.bs","Wy.sd","D.bs","D.sd","P.bs",
			 "betaX","betaY") 
	if (BLzero == 0) 
		{
			out <-  cbind(out,result[1,3],result[1,4],agediff[1,2],
				pRand,result[1,7],result[1,8],result[1,9],
				result[1,10])
			colnames(out) <- c("Nreps","Wx.as","Wy.as","D.as","Wx.bs",
				 	 "Wx.sd","Wy.bs","Wy.sd","D.bs","D.sd","P.bs",
				 	 "betaX","betaY","Wx.pic","Wy.pic",
					 "D.pic","P.tips","P.contrasts",
					 "BMtestX","BMtestY","NormTestX","NormTestY")
		}
	return(out)
}
