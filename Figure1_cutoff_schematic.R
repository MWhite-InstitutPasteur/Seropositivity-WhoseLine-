#####################################################################
#####################################################################
##                                                                 ##
## CODE FOR FITTING GAUSSIAN MIXTURE MODELS TO ANTIBODY DATA       ##
##                                                                 ##
## Please feel free to share modify the code as you see fit        ##   
## (but please maintain appropriate accreditation)                 ##
##                                                                 ##   
## Michael White, michael.white@pasteur.fr                         ##
## Gaelle Baudemont, gaelle.baudemont@pasteur.fr                   ##
##                                                                 ##
#####################################################################
#####################################################################

rm(list=ls())
library("rstan")
library(MASS)
library(compiler)

set.seed(123)


#############################################
#############################################
##        ##                               ## 
##   ##   ##  ####    ####  ######  ####   ##
##  ###   ##  ## ##  ##  ##   ##   ##  ##  ##
##   ##   ##  ##  ## ######   ##   ######  ##
##   ##   ##  ## ##  ##  ##   ##   ##  ##  ##
##  ####  ##  ####   ##  ##   ##   ##  ##  ##
##        ##                               ##
#############################################
#############################################


#############################################
## 1.1. Sample numbers and data limits

N_pos <- 500
N_neg <- 500

N_sam <- N_pos + N_neg

N_bins <- 30
Ab_min <- log(3)
Ab_max <- log(30000)

Ab_bin_edges <- seq(from=Ab_min, to=Ab_max, length=N_bins+1)

Ab_bin_width <- Ab_bin_edges[2] - Ab_bin_edges[1]


#############################################
## 1.2. Data set A: Clear bimodal pattern

#############################################
## Simulate data

mu_pos_A <- log(3000)
mu_neg_A <- log(50)

sig_pos_A <- 0.8
sig_neg_A <- 1.2

Ab_pos_A <- rnorm( N_pos, mean=mu_pos_A, sd=sig_pos_A )
Ab_neg_A <- rnorm( N_neg, mean=mu_neg_A, sd=sig_neg_A )


#############################################
## Put positive samples in bins

pos_bins_A <- rep(NA, length=N_bins)

for(i in 1:N_bins)
{
	index <- which( Ab_pos_A >  Ab_bin_edges[i] &
                      Ab_pos_A <= Ab_bin_edges[i+1] )

	pos_bins_A[i] <- length(index)	
}

pos_bins_A[1]      <- pos_bins_A[1] + length(which( Ab_pos_A < Ab_bin_edges[1] ))
pos_bins_A[N_bins] <- pos_bins_A[N_bins] + length(which( Ab_pos_A > Ab_bin_edges[N_bins+1] ))


#############################################
## Put negative samples in bins

neg_bins_A <- rep(NA, length=N_bins)

for(i in 1:N_bins)
{
	index <- which( Ab_neg_A >  Ab_bin_edges[i] &
                      Ab_neg_A <= Ab_bin_edges[i+1] )

	neg_bins_A[i] <- length(index)	
}

neg_bins_A[1]      <- neg_bins_A[1] + length(which( Ab_neg_A < Ab_bin_edges[1] ))
neg_bins_A[N_bins] <- neg_bins_A[N_bins] + length(which( Ab_neg_A > Ab_bin_edges[N_bins+1] ))


#############################################
## 1.3. Data set B: Noisier data without bimodal pattern

#############################################
## Simulate data

mu_pos_B <- log(5000)
mu_neg_B <- log(50)

sig_pos_B <- 1.6
sig_neg_B <- 2.5


Ab_pos_B <- rnorm( N_pos, mean=mu_pos_B, sd=sig_pos_B )
Ab_neg_B <- rnorm( N_neg, mean=mu_neg_B, sd=sig_neg_B )

#############################################
## Put positive samples in bins

pos_bins_B <- rep(NA, length=N_bins)

for(i in 1:N_bins)
{
	index <- which( Ab_pos_B >  Ab_bin_edges[i] &
                      Ab_pos_B <= Ab_bin_edges[i+1] )

	pos_bins_B[i] <- length(index)	
}

pos_bins_B[1]      <- pos_bins_B[1] + length(which( Ab_pos_B < Ab_bin_edges[1] ))
pos_bins_B[N_bins] <- pos_bins_B[N_bins] + length(which( Ab_pos_B > Ab_bin_edges[N_bins+1] ))


#############################################
## Put negative samples in bins

neg_bins_B <- rep(NA, length=N_bins)

for(i in 1:N_bins)
{
	index <- which( Ab_neg_B >  Ab_bin_edges[i] &
                      Ab_neg_B <= Ab_bin_edges[i+1] )

	neg_bins_B[i] <- length(index)	
}

neg_bins_B[1]      <- neg_bins_B[1] + length(which( Ab_neg_B < Ab_bin_edges[1] ))
neg_bins_B[N_bins] <- neg_bins_B[N_bins] + length(which( Ab_neg_B > Ab_bin_edges[N_bins+1] ))






########################################
########################################
##          ##                        ## 
##   ####   ##  #####   ####   ####   ##
##  ##  ##  ##  ##  ## ##  ## ##  ##  ##
##     ##   ##  #####  ##  ## ##      ## 
##    ##    ##  ## ##  ##  ## ##  ##  ##
##   #####  ##  ##  ##  ####   ####   ##
##          ##                        ##
########################################
########################################

#####################################
## 2.1 Data set A

Ab_cuts_A <- c( Ab_neg_A, Ab_pos_A)
Ab_cuts_A <- sort( Ab_cuts_A, decreasing=FALSE )

Ab_cuts_A <- c( -100, Ab_cuts_A, 100)

N_cuts_A <- length(Ab_cuts_A)

sens_A <- rep(NA, N_cuts_A)
spec_A <- rep(NA, N_cuts_A)

for(i in 1:N_cuts_A)
{
	spec_A[i] <- length(which(Ab_neg_A < Ab_cuts_A[i]))/length(Ab_neg_A)
	sens_A[i] <- length(which(Ab_pos_A > Ab_cuts_A[i]))/length(Ab_pos_A)	
}

Youden_index_A <- which.max( sens_A + spec_A)

spec98_index_A <- which.min(abs(Ab_cuts_A - sort(Ab_neg_A)[ round(0.977*length(Ab_neg_A)) ]) )


#####################################
## 2.2 Data set B

Ab_cuts_B <- c( Ab_neg_B, Ab_pos_B)
Ab_cuts_B <- sort( Ab_cuts_B, decreasing=FALSE )

Ab_cuts_B <- c( -100, Ab_cuts_B, 100)

N_cuts_B <- length(Ab_cuts_B)

sens_B <- rep(NA, N_cuts_B)
spec_B <- rep(NA, N_cuts_B)

for(i in 1:N_cuts_B)
{
	spec_B[i] <- length(which(Ab_neg_B < Ab_cuts_B[i]))/length(Ab_neg_B)
	sens_B[i] <- length(which(Ab_pos_B > Ab_cuts_B[i]))/length(Ab_pos_B)	
}

Youden_index_B <- which.max( sens_B + spec_B)

spec98_index_B <- which.min(abs(Ab_cuts_B - sort(Ab_neg_B)[ round(0.977*length(Ab_neg_B)) ]) )



##########################################
##########################################
##          ##                          ## 
##   ####   ##   ####  #     # #     #  ##
##  ##  ##  ##  ##     ##   ## ##   ##  ##
##     ##   ##  ## ### ####### #######  ##
##  ##  ##  ##  ##  ## ## # ## ## # ##  ##
##   ####   ##   ####  ##   ## ##   ##  ##
##          ##                          ##
##########################################
##########################################

N_par <- 5

N_chains <- 4
chain_cols <- c("#0033FF", "#0066FF", "#0099FF", "#0099CC")

N_iter <- 5000
N_warmup <- 500

N_witer <-  N_iter - N_warmup

N_subset <- 500


################################################### 
################################################### 
## 3.1 Fit GMM to data set A

print("Creating data for model")  
list_data_model <- list(N_unknown  = N_sam,                    ## Number of unknown samples
                        Ab_unknown = c(Ab_neg_A, Ab_pos_A),    ## Vector of unknown Ab values 
                        Ab_min     = Ab_min,                   ## Maximum Ab value
                        Ab_max     = Ab_max )                  ## Maximum Ab value



GMM_model <- stan_model(file="GMM_unknown_Ab_data.stan")

    
###################################################     
# Defining initial parameters 

init_fun <- function() 
{
	list( theta     = 0.5, 
            mu_Neg    = Ab_min + 0.1*(Ab_max - Ab_min), 
            sigma_Neg = 0.5,
            mu_delta  = 3,
            sigma_Pos = 0.5
          )
}
    
###################################################     
# Fit the model 
    
print("Fitting model")  
GMM_fit <- sampling(GMM_model,                      # Stan program
                    data = list_data_model,         # named list of data
                    chains = N_chains,               # number of Markov chains
                    warmup = N_warmup,              # number of warmup iterations per chain (2500)
                    iter = N_iter,                  # total number of iterations per chain (10000)
                    refresh = 200,                  # show progress every 'refresh' iterations
                    init = init_fun,                # Initial values set at random
                    verbose = TRUE                  # Gives more explicit error messages
)
    
    
print("Model fitting completed")  

###################################################         
# Data Management of the results

GMM_chains <- extract(GMM_fit)[1:(N_par+1)] 
GMM_chains <- matrix(unlist(GMM_chains), ncol = N_par+1) 
colnames(GMM_chains) <- c("theta", "mu_N", "sigma_N", "mu_delta", "sigma_Pos",  "LogL")


###################################################         
# Chain analytics

par(mfcol=c(2,N_par+1))
 	
for(j in 1:(N_par+1))
{
	plot(x=-1, y=-1, 
	xlim=c(0,N_witer), ylim=range(GMM_chains[,j]),
	xlab="MCMC iteration", ylab=colnames(GMM_chains)[j], 
	main=colnames(GMM_chains)[j] )

	for(m in 1:N_chains)
	{
		points(x = 1:N_witer, 
                   y = GMM_chains[1:N_witer + (m-1)*N_witer,j], 
			 type='l', col=chain_cols[m], cex=0.25 )
	}


	DEN_1 <- density( GMM_chains[1:N_witer ,j] )
	
	plot(x=DEN_1$x, y=DEN_1$y, 
           type='l', col=chain_cols[1],
	     xlim=range(GMM_chains[,j]),
	     xlab=colnames(GMM_chains)[j], ylab="density",
	     main=colnames(GMM_chains)[j] )

	for(m in 2:N_chains)
	{
		DEN_m <- density( GMM_chains[1:N_witer + (m-1)*N_witer,j] )
	
		points(x=DEN_m$x, y=DEN_m$y, 
            type='l', col=chain_cols[m])			
	}
}

################################################### 
###################################################         
# Compare model to data set A

GMM_chains_subset <- GMM_chains[seq(from=1, to=nrow(GMM_chains), length=N_subset),]

Ab_seq <- seq(from=Ab_min, to=Ab_max, length=100)

Neg_dist_A <- matrix(NA, nrow=N_subset, ncol=length(Ab_seq))
Pos_dist_A <- matrix(NA, nrow=N_subset, ncol=length(Ab_seq))

	

for(i in 1:N_subset)
{
	Neg_dist_A[i,] <- (1-GMM_chains_subset[i,1])*
                         dnorm(Ab_seq, mean=GMM_chains_subset[i,2], 
                                     sd=GMM_chains_subset[i,3])

	Pos_dist_A[i,] <- GMM_chains_subset[i,1]*
                        dnorm(Ab_seq, mean=GMM_chains_subset[i,2]+GMM_chains_subset[i,4], 
                                    sd=GMM_chains_subset[i,5])
}

Neg_dist_quant_A <- matrix(NA, nrow=3, ncol=100)
Pos_dist_quant_A <- matrix(NA, nrow=3, ncol=100)

for(j in 1:length(Ab_seq))
{
	Neg_dist_quant_A[,j] <- quantile( Neg_dist_A[,j], prob=c(0.025, 0.5, 0.975) )
	Pos_dist_quant_A[,j] <- quantile( Pos_dist_A[,j], prob=c(0.025, 0.5, 0.975) )
}


#####################################
## Calulate cut-off

GMM_cut_A <- GMM_chains_subset[,2] + 2*GMM_chains_subset[,3]

GMM_cut_median_A <- median( GMM_cut_A )


GMM_cut_A <- density( GMM_cut_A )



################################################### 
###################################################         
# Plot fit of model to data set A

par(mfrow=c(1,1))


plot(x=1e10, y=1e10, 
pch=15, cex=1, col="orangered",
xlim=log(c(2,50000)), ylim=c(0,120), 
#xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n',
xlab="antibody levels (MFI)", ylab="",  
main="" )



for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=c(0, 0, neg_bins_A[i], neg_bins_A[i]), 
		   col="royalblue" )
}

for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=neg_bins_A[i] + c(0, 0, pos_bins_A[i], pos_bins_A[i]), 
		   col="orangered" )
}


mtext(side = 2, line = 2.5, 
cex=1.75, 
text="sample set A")



points(x=Ab_seq, y=Ab_bin_width*N_sam*Neg_dist_quant_A[2,], 
       type='l', lwd=2, col="blue")

polygon( x=c( Ab_seq, rev(Ab_seq) ),
         y=Ab_bin_width*c( N_sam*Neg_dist_quant_A[1,], rev(N_sam*Neg_dist_quant_A[3,]) ), 
	   col=rgb(0,0,1,0.5), border=NA)

points(x=Ab_seq, y=Ab_bin_width*N_sam*Pos_dist_quant_A[2,], 
       type='l', lwd=2, col="red")

polygon( x=c( Ab_seq, rev(Ab_seq) ),
         y=Ab_bin_width*c( N_sam*Pos_dist_quant_A[1,], rev(N_sam*Pos_dist_quant_A[3,]) ), 
	   col=rgb(1,0,0,0.5), border=NA)


polygon( x=c( GMM_cut_A$x, rev(GMM_cut_A$x) ),
         y=(120/max(GMM_cut_A$y))*c( rep(0,length(GMM_cut_A$y)), rev(GMM_cut_A$y) ), 
	   col=rgb(0,1,0,0.5), border=NA)

points( x = rep( GMM_cut_median_A, 2),
        y = (120/max(GMM_cut_A$y))*c(0, GMM_cut_A$y[which.min(abs(GMM_cut_A$x - GMM_cut_median_A))] ),
        type='l', col="green", lwd=2 )



################################################### 
################################################### 
## 3.2 Fit GMM to data set B

print("Creating data for model")  
list_data_model <- list(N_unknown  = N_sam,                    ## Number of unknown samples
                        Ab_unknown = c(Ab_neg_B, Ab_pos_B),    ## Vector of unknown Ab values 
                        Ab_min     = Ab_min,                   ## Maximum Ab value
                        Ab_max     = Ab_max )                  ## Maximum Ab value



GMM_model <- stan_model(file="GMM_unknown_Ab_data.stan")

    
###################################################     
# Defining initial parameters 

init_fun <- function() 
{
	list( theta     = 0.5, 
            mu_Neg    = Ab_min + 0.1*(Ab_max - Ab_min), 
            sigma_Neg = 0.5,
            mu_delta  = 3,
            sigma_Pos = 0.5
          )
}
    
###################################################     
# Fit the model 
    
print("Fitting model")  
GMM_fit <- sampling(GMM_model,                      # Stan program
                    data = list_data_model,         # named list of data
                    chains = N_chains,               # number of Markov chains
                    warmup = N_warmup,              # number of warmup iterations per chain (2500)
                    iter = N_iter,                  # total number of iterations per chain (10000)
                    refresh = 200,                  # show progress every 'refresh' iterations
                    init = init_fun,                # Initial values set at random
                    verbose = TRUE                  # Gives more explicit error messages
)
    
    
print("Model fitting completed")  

###################################################         
# Data Management of the results

GMM_chains <- extract(GMM_fit)[1:(N_par+1)] 
GMM_chains <- matrix(unlist(GMM_chains), ncol = N_par+1) 
colnames(GMM_chains) <- c("theta", "mu_N", "sigma_N", "mu_delta", "sigma_Pos",  "LogL")


###################################################         
# Chain analytics

par(mfcol=c(2,N_par+1))
 	
for(j in 1:(N_par+1))
{
	plot(x=-1, y=-1, 
	xlim=c(0,N_witer), ylim=range(GMM_chains[,j]),
	xlab="MCMC iteration", ylab=colnames(GMM_chains)[j], 
	main=colnames(GMM_chains)[j] )

	for(m in 1:N_chains)
	{
		points(x = 1:N_witer, 
                   y = GMM_chains[1:N_witer + (m-1)*N_witer,j], 
			 type='l', col=chain_cols[m], cex=0.25 )
	}


	DEN_1 <- density( GMM_chains[1:N_witer,j] )
	
	plot(x=DEN_1$x, y=DEN_1$y, 
           type='l', col=chain_cols[1],
	     xlim=range(GMM_chains[,j]),
	     xlab=colnames(GMM_chains)[j], ylab="density",
	     main=colnames(GMM_chains)[j] )

	for(m in 2:N_chains)
	{
		DEN_m <- density( GMM_chains[1:N_witer + (m-1)*N_witer,j] )
	
		points(x=DEN_m$x, y=DEN_m$y, 
            type='l', col=chain_cols[m])			
	}
}

################################################### 
###################################################         
# 2.2 Compare model to data set B

GMM_chains_subset <- GMM_chains[seq(from=1, to=nrow(GMM_chains), length=N_subset),]

Ab_seq <- seq(from=Ab_min, to=Ab_max, length=100)

Neg_dist_B <- matrix(NA, nrow=N_subset, ncol=length(Ab_seq))
Pos_dist_B <- matrix(NA, nrow=N_subset, ncol=length(Ab_seq))

	

for(i in 1:N_subset)
{
	Neg_dist_B[i,] <- (1-GMM_chains_subset[i,1])*
                       dnorm(Ab_seq, mean=GMM_chains_subset[i,2], 
                                     sd=GMM_chains_subset[i,3])

	Pos_dist_B[i,] <- GMM_chains_subset[i,1]*
                      dnorm(Ab_seq, mean=GMM_chains_subset[i,2]+GMM_chains_subset[i,4], 
                                    sd=GMM_chains_subset[i,5])
}

Neg_dist_quant_B <- matrix(NA, nrow=3, ncol=100)
Pos_dist_quant_B <- matrix(NA, nrow=3, ncol=100)

for(j in 1:length(Ab_seq))
{
	Neg_dist_quant_B[,j] <- quantile( Neg_dist_B[,j], prob=c(0.025, 0.5, 0.975) )
	Pos_dist_quant_B[,j] <- quantile( Pos_dist_B[,j], prob=c(0.025, 0.5, 0.975) )
}


#####################################
## Calulate cut-off

GMM_cut_B <- GMM_chains_subset[,2] + 2*GMM_chains_subset[,3]

GMM_cut_median_B <- median( GMM_cut_B )


GMM_cut_B <- density( GMM_cut_B )



################################################### 
###################################################         
# Plot fit of model to data set B

par(mfrow=c(1,1))


plot(x=1e10, y=1e10, 
pch=15, cex=1, col="orangered",
xlim=log(c(2,50000)), ylim=c(0,120), 
#xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n',
xlab="antibody levels", ylab="",  
main="" )



for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=c(0, 0, neg_bins_B[i], neg_bins_B[i]), 
		   col="royalblue" )
}

for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=neg_bins_B[i] + c(0, 0, pos_bins_B[i], pos_bins_B[i]), 
		   col="orangered" )
}


mtext(side = 2, line = 2.5, 
cex=1.75, 
text="sample set B")



points(x=Ab_seq, y=Ab_bin_width*N_sam*Neg_dist_quant_B[2,], 
       type='l', lwd=2, col="blue")

polygon( x=c( Ab_seq, rev(Ab_seq) ),
         y=Ab_bin_width*c( N_sam*Neg_dist_quant_B[1,], rev(N_sam*Neg_dist_quant_B[3,]) ), 
	   col=rgb(0,0,1,0.5), border=NA)

points(x=Ab_seq, y=Ab_bin_width*N_sam*Pos_dist_quant_B[2,], 
       type='l', lwd=2, col="red")

polygon( x=c( Ab_seq, rev(Ab_seq) ),
         y=Ab_bin_width*c( N_sam*Pos_dist_quant_B[1,], rev(N_sam*Pos_dist_quant_B[3,]) ), 
	   col=rgb(1,0,0,0.5), border=NA)


polygon( x=c( GMM_cut_B$x, rev(GMM_cut_B$x) ),
         y=(120/max(GMM_cut_B$y))*c( rep(0,length(GMM_cut_B$y)), rev(GMM_cut_B$y) ), 
	   col=rgb(0,1,0,0.5), border=NA)

points( x = rep( GMM_cut_median_B, 2),
        y = (120/max(GMM_cut_B$y))*c(0, GMM_cut_B$y[which.min(abs(GMM_cut_B$x - GMM_cut_median_B))] ),
        type='l', col="green", lwd=2 )




##############################################
##############################################
##          ##                              ## 
##  ##      ##  #####  ##     ####  ######  ##
##  ## ##   ##  ##  ## ##    ##  ##   ##    ##
##  ######  ##  #####  ##    ##  ##   ##    ## 
##     ##   ##  ##     ##    ##  ##   ##    ##
##     ##   ##  ##     #####  ####    ##    ##
##          ##                              ##
##############################################
##############################################

line_seq_y = log(c(3, 10, 30, 100, 300, 1000, 3000, 10000, 30000))
line_seq_x = c(0, 20, 40, 60, 80, 100, 120)



main.size = 2.5
lab.size  = 2
axis.size = 1.0
line.size = 3
text.size = 1


tiff( file="Figure1_cutoff_schematic.tif", width=30, height=18, units="cm", res=500)



lay.mat <- rbind( c( 1, 1, 2, 3), 
                  c( 4, 5, 6, 7), 
                  c( 8, 9,10,11), 
                  c(12,12,12,12) )
layout(lay.mat, widths=c(7,5,6,6), heights=c(1,6,6,1))
layout.show(12)



############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))


plot.new()
title( "(1) ROC Curves", cex.main=2.0, line=-2)

plot.new()
title( "(2) Negative Controls", cex.main=2.0, line=-2)

plot.new()
title( "(3) Mixture Model", cex.main=2.0, line=-2)



#####################################
#####################################
##                                 ##
## PANEL 1: nice data, ROC curves  ##
##                                 ##
#####################################
#####################################

par(mar=c(5,5,1,1))
par(mgp=c(2, 0.65, 0))

plot(x=1e10, y=1e10, 
pch=15, cex=1, col="orangered",
xlim=log(c(2,50000)), ylim=c(0,120), 
xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n',
xlab="antibody levels (MFI)", ylab="",  
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size )

for(i in 1:length(line_seq_x))
{
	points(x=c(0.001,100000), y=rep(line_seq_x[i],2), type='l', lwd=1, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(0.1,100000), type='l', lwd=1, col="grey", lty="dashed")
}


for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=c(0, 0, neg_bins_A[i], neg_bins_A[i]), 
		   col="royalblue" )
}

for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=neg_bins_A[i] + c(0, 0, pos_bins_A[i], pos_bins_A[i]), 
		   col="orangered" )
}


points( x=rep(Ab_cuts_A[spec98_index_A],2), y=c(0,100000), 
        type='l', lwd=line.size, col="green3")

mtext(side = 2, line = 2.5, 
cex=1.75, 
text="sample set A")


axis(1, at=log(c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000)), 
        label=c("", "", "10", "", "100", "", "1000", "", "10000", "", ""), cex.axis=axis.size )


axis(2, at=c(0, 20, 40, 60, 80, 100, 120), 
        label=c("0", "20", "40", "60", "80", "100", "120"), cex.axis=axis.size, las=2 )


#####################################
#####################################
##                                 ##
## PANEL 2: Data set A, ROC curves ##
##                                 ##
#####################################
#####################################

par(mar=c(5,4,1,1))

line_seq <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)

plot(x=c(0,1), y=c(0,1), type='l', lty="dashed", lwd=1,
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="")

mtext(side = 1, line = 2.25, 
cex=1.25, 
text="1 - specificity")

mtext(side = 2, line = 2.5, 
cex=1.25, 
text="sensitivity")

for(i in 1:length(line_seq))
{
	points(x=c(-1,1), y=rep(line_seq[i],2), type='l', lwd=1, col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(-1,1), type='l', lwd=1, col="grey", lty="dashed")
}



points( x= 1 - spec_A, y=sens_A, 
    	    type='S', lwd=4, col="black" )


points( x= 1 - spec_A[spec98_index_A], y=sens_A[spec98_index_A], 
    	    pch=19, cex=2, col="green3" )


axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        cex.axis=0.9, lwd = 1, lwd.ticks = 1) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"),
        las=2, cex.axis=1, lwd = 1, lwd.ticks = 1 ) 



############################################
############################################
##                                        ##
## PANEL 3: data set A, Negative controls ##
##                                        ##
############################################
############################################

par(mar=c(5,3,1,1))


plot(x=1e10, y=1e10, 
pch=15, cex=1, col="orangered",
xlim=log(c(2,50000)), ylim=c(0,120), 
xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n',
xlab="antibody levels (MFI)", ylab="",  
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size )

for(i in 1:length(line_seq_x))
{
	points(x=c(0.001,100000), y=rep(line_seq_x[i],2), type='l', lwd=1, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(0.1,100000), type='l', lwd=1, col="grey", lty="dashed")
}


for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=c(0, 0, neg_bins_A[i], neg_bins_A[i]), 
		   col="royalblue" )
}

	
points( x=rep(mean(Ab_neg_A) + 2*sd(Ab_neg_A),2), y=c(0,100000), 
        type='l', lwd=line.size, col="green3")

points( x=rep(mean(Ab_neg_A) + 3*sd(Ab_neg_A),2), y=c(0,100000), 
        type='l', lwd=line.size, col="green3", lty="dashed")




axis(1, at=log(c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000)), 
        label=c("", "", "10", "", "100", "", "1000", "", "10000", "", ""), cex.axis=axis.size )


axis(2, at=c(0, 20, 40, 60, 80, 100, 120), 
        label=c("0", "20", "40", "60", "80", "100", "120"), cex.axis=axis.size, las=2 )



###########################################
###########################################
##                                       ##
## PANEL 4: data set A, Gaussian mixture ##
##                                       ##
###########################################
###########################################



plot(x=1e10, y=1e10, 
pch=15, cex=1, col="orangered",
xlim=log(c(2,50000)), ylim=c(0,120), 
xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n',
xlab="antibody levels (MFI)", ylab="",  
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size )

for(i in 1:length(line_seq_x))
{
	points(x=c(0.001,100000), y=rep(line_seq_x[i],2), type='l', lwd=1, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(0.1,100000), type='l', lwd=1, col="grey", lty="dashed")
}




for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=c(0, 0, neg_bins_A[i], neg_bins_A[i]), 
		   col="grey" )
}

for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=neg_bins_A[i] + c(0, 0, pos_bins_A[i], pos_bins_A[i]), 
		   col="grey" )
}


##mtext(side = 2, line = 2.5, 
##cex=1.75, 
##text="sample set A")



points(x=Ab_seq, y=Ab_bin_width*N_sam*Neg_dist_quant_A[2,], 
       type='l', lwd=2, col="blue")

polygon( x=c( Ab_seq, rev(Ab_seq) ),
         y=Ab_bin_width*c( N_sam*Neg_dist_quant_A[1,], rev(N_sam*Neg_dist_quant_A[3,]) ), 
	   col=rgb(0,0,1,0.33), border=NA)

points(x=Ab_seq, y=Ab_bin_width*N_sam*Pos_dist_quant_A[2,], 
       type='l', lwd=2, col="red")

polygon( x=c( Ab_seq, rev(Ab_seq) ),
         y=Ab_bin_width*c( N_sam*Pos_dist_quant_A[1,], rev(N_sam*Pos_dist_quant_A[3,]) ), 
	   col=rgb(1,0,0,0.33), border=NA)


polygon( x=c( GMM_cut_A$x, rev(GMM_cut_A$x) ),
         y=(120/max(GMM_cut_A$y))*c( rep(0,length(GMM_cut_A$y)), rev(GMM_cut_A$y) ), 
	   col=rgb(0,1,0,0.33), border=NA)

points( x = rep( GMM_cut_median_A, 2),
        y = (120/max(GMM_cut_A$y))*c(0, GMM_cut_A$y[which.min(abs(GMM_cut_A$x - GMM_cut_median_A))] ),
        type='l', col="green", lwd=2 )



	
axis(1, at=log(c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000)), 
        label=c("", "", "10", "", "100", "", "1000", "", "10000", "", ""), cex.axis=axis.size )


axis(2, at=c(0, 20, 40, 60, 80, 100, 120), 
        label=c("0", "20", "40", "60", "80", "100", "120"), cex.axis=axis.size, las=2 )




#####################################
#####################################
##                                 ##
## PANEL 5: data set B, ROC curves ##
##                                 ##
#####################################
#####################################


par(mar = c(5,5,1,1))


plot(x=1e10, y=1e10, 
pch=15, cex=1, col="orangered",
xlim=log(c(2,50000)), ylim=c(0,120), 
xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n',
xlab="antibody levels (MFI)", ylab="",  
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size )

for(i in 1:length(line_seq_x))
{
	points(x=c(0.001,100000), y=rep(line_seq_x[i],2), type='l', lwd=1, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(0.1,100000), type='l', lwd=1, col="grey", lty="dashed")
}


for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=c(0, 0, neg_bins_B[i], neg_bins_B[i]), 
		   col="royalblue" )
}

for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=neg_bins_B[i] + c(0, 0, pos_bins_B[i], pos_bins_B[i]), 
		   col="orangered" )
}

	
points( x=rep(Ab_cuts_B[spec98_index_B],2), y=c(0,100000), 
        type='l', lwd=line.size, col="green3")


mtext(side = 2, line = 2.5, 
cex=1.75, 
text="sample set B")


axis(1, at=log(c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000)), 
        label=c("", "", "10", "", "100", "", "1000", "", "10000", "", ""), cex.axis=axis.size )


axis(2, at=c(0, 20, 40, 60, 80, 100, 120), 
        label=c("0", "20", "40", "60", "80", "100", "120"), cex.axis=axis.size, las=2 )


#####################################
#####################################
##                                 ##
## PANEL 6: data set B, ROC curves ##
##                                 ##
#####################################
#####################################

par(mar=c(5,4,1,1))


line_seq <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed", lwd=1,
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="")

mtext(side = 1, line = 2.25, 
cex=1.25, 
text="1 - specificity")

mtext(side = 2, line = 2.5, 
cex=1.25, 
text="sensitivity")

for(i in 1:length(line_seq))
{
	points(x=c(-1,1), y=rep(line_seq[i],2), type='l', lwd=1, col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(-1,1), type='l', lwd=1, col="grey", lty="dashed")
}



points( x= 1 - spec_B, y=sens_B, 
    	    type='S', lwd=4, col="black" )


points( x= 1 - spec_B[spec98_index_B], y=sens_B[spec98_index_B], 
    	    pch=19, cex=2, col="green3" )


axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        cex.axis=0.9, lwd = 1, lwd.ticks = 1) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"),
        las=2, cex.axis=1, lwd = 1, lwd.ticks = 1 ) 




#####################################
#####################################
##                                 ##
## PANEL 7: data set B, Neg con    ##
##                                 ##
#####################################
#####################################

par(mar=c(5,3,1,1))


plot(x=1e10, y=1e10, 
pch=15, cex=1, col="orangered",
xlim=log(c(2,50000)), ylim=c(0,120), 
xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n',
xlab="antibody levels (MFI)", ylab="",  
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size )

for(i in 1:length(line_seq_x))
{
	points(x=c(0.001,100000), y=rep(line_seq_x[i],2), type='l', lwd=1, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(0.1,100000), type='l', lwd=1, col="grey", lty="dashed")
}


for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=c(0, 0, neg_bins_B[i], neg_bins_B[i]), 
		   col="royalblue" )
}

	
points( x=rep(mean(Ab_neg_B) + 2*sd(Ab_neg_B),2), y=c(0,100000), 
        type='l', lwd=line.size, col="green3")


points( x=rep(mean(Ab_neg_B) + 3*sd(Ab_neg_B),2), y=c(0,100000), 
        type='l', lwd=line.size, col="green3", lty="dashed")

axis(1, at=log(c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000)), 
        label=c("", "", "10", "", "100", "", "1000", "", "10000", "", ""), cex.axis=axis.size )


axis(2, at=c(0, 20, 40, 60, 80, 100, 120), 
        label=c("0", "20", "40", "60", "80", "100", "120"), cex.axis=axis.size, las=2 )


#####################################
#####################################
##                                 ##
## PANEL 8: data set B, GMM        ##
##                                 ##
#####################################
#####################################


plot(x=1e10, y=1e10, 
pch=15, cex=1, col="orangered",
xlim=log(c(2,50000)), ylim=c(0,120), 
xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n',
xlab="antibody levels (MFI)", ylab="",  
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size )

for(i in 1:length(line_seq_x))
{
	points(x=c(0.001,100000), y=rep(line_seq_x[i],2), type='l', lwd=1, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(0.1,100000), type='l', lwd=1, col="grey", lty="dashed")
}



for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=c(0, 0, neg_bins_B[i], neg_bins_B[i]), 
		   col="grey" )
}

for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=neg_bins_B[i] + c(0, 0, pos_bins_B[i], pos_bins_B[i]), 
		   col="grey" )
}


##mtext(side = 2, line = 2.5, 
##cex=1.75, 
##text="sample set B")



points(x=Ab_seq, y=Ab_bin_width*N_sam*Neg_dist_quant_B[2,], 
       type='l', lwd=2, col="blue")

polygon( x=c( Ab_seq, rev(Ab_seq) ),
         y=Ab_bin_width*c( N_sam*Neg_dist_quant_B[1,], rev(N_sam*Neg_dist_quant_B[3,]) ), 
	   col=rgb(0,0,1,0.33), border=NA)

points(x=Ab_seq, y=Ab_bin_width*N_sam*Pos_dist_quant_B[2,], 
       type='l', lwd=2, col="red")

polygon( x=c( Ab_seq, rev(Ab_seq) ),
         y=Ab_bin_width*c( N_sam*Pos_dist_quant_B[1,], rev(N_sam*Pos_dist_quant_B[3,]) ), 
	   col=rgb(1,0,0,0.33), border=NA)


polygon( x=c( GMM_cut_B$x, rev(GMM_cut_B$x) ),
         y=(120/max(GMM_cut_B$y))*c( rep(0,length(GMM_cut_B$y)), rev(GMM_cut_B$y) ), 
	   col=rgb(0,1,0,0.33), border=NA)

points( x = rep( GMM_cut_median_B, 2),
        y = (120/max(GMM_cut_B$y))*c(0, GMM_cut_B$y[which.min(abs(GMM_cut_B$x - GMM_cut_median_B))] ),
        type='l', col="green", lwd=2 )



	
axis(1, at=log(c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000)), 
        label=c("", "", "10", "", "100", "", "1000", "", "10000", "", ""), cex.axis=axis.size )


axis(2, at=c(0, 20, 40, 60, 80, 100, 120), 
        label=c("0", "20", "40", "60", "80", "100", "120"), cex.axis=axis.size, las=2 )



##############################################
## Add the legend to the bottom of the plot ##
##############################################

par(mar = c(0,0,0,0))
plot.new()
legend(x='center', 
       legend = c( "negative samples", "positive samples", "unknown status" ), 
       fill = c("royalblue", "orangered", "grey"),
       border = c("royalblue", "orangered", "grey"),
	 ncol=3, bty="n", cex=2)





dev.off()







