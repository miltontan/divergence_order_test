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
# Note, this script requires ape and msm to be loaded

dot <- function (tree, x, y, nsim = 0, replace=FALSE)
{
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

# Run ace and generate Ttips and Tint directly
	ace_x <- ace(x,tree)
	ace_y <- ace(y,tree)
	beta1 <- ace_x$sigma2[1]
	beta2 <- ace_y$sigma2[1]
	Ttips <- data.frame(node=1:length(tree$tip.label),age=rep(0,length(tree$tip.label)),taxon=tree$tip.label,x=x[tree$tip.label],x_se=rep(0,length(tree$tip.label)),y=y[tree$tip.label],y_se=rep(0,length(tree$tip.label)))
	rownames(Ttips) <- Ttips[,'taxon']
	Tint <- data.frame(node=unique(tree$edge[,1]),age=branching.times(tree),taxon=rep("",length(unique(tree$edge[,1]))),x=ace_x$ace,x_se=(ace_x$CI95[,2]-ace_x$CI95[,1])/3.92,y=ace_y$ace,y_se=(ace_y$CI95[,2]-ace_y$CI95[,1])/3.92)
	rownames(Tint) <- unique(tree$edge[, 1])

# Recombine traits data, tips first then internals
	T2 <- rbind(Ttips,Tint)
	colnames(T2) <- colnames(Ttips)

#parse vars and ses
	ages <- as.numeric(Tint[,'age'])
	X <- as.numeric(T2[,'x'])
	Xse <- as.numeric(T2[,'x_se'])
	Y <- as.numeric(T2[,'y'])
	Yse <- as.numeric(T2[,'y_se'])

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
	Cx <- abs(pic(X[1:NTips],tree,scaled=FALSE))
	Cy <- abs(pic(Y[1:NTips],tree,scaled=FALSE))
	CageX <- weighted.mean(ages,Cx)
 	CageY <- weighted.mean(ages,Cy)
	CxObs <- Cx
	CyObs <- Cy
	Fx <- pic(X[1:NTips],tree,var.contrasts=TRUE)
	Fy <- pic(Y[1:NTips],tree,var.contrasts=TRUE)
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
			#print(n)
			Xr[n] <- rnorm(1,X[n],Xse[n])
			Yr[n] <- rnorm(1,Y[n],Yse[n])
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
		Cx <- abs(pic(sample(X[1:NTips],replace=replace),tree,scaled=FALSE))
		Cy <- abs(pic(sample(Y[1:NTips],replace=replace),tree,scaled=FALSE))
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

#mean and standard deviation of ace contrasts at each node
#ACx and ACy are now averages across all reps
	meanACx <- cumACx/Nreps
	SDx <- (cumACx2-(cumACx^2/Nreps))/(Nreps-1)
	SDx[SDx<0] <- 0
	SDx <- sqrt(SDx)
	
	meanACy <- cumACy/Nreps
	SDy <- (cumACy2-(cumACy^2/Nreps))/(Nreps-1)
	SDy[SDy<0] <- 0
	SDy <- sqrt(SDy)
	
	#standardize contrasts relative to largest
	#SDx <- SDx/max(meanACx)
	#meanACx <- ACx/max(meanACx) 
	
	#SDy <- SDy/max(meanACy)
	#meanACy <- ACy/max(meanACy)

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
