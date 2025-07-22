
par(ask=FALSE)

NTD_data <- read.csv("20plexNTD_all_data_20220419.csv")

NTD_short_names <- colnames( NTD_data )[9:28]
N_NTD <- length(NTD_short_names)



NTD_long_names <- c("LF: Bm33", "Crypto: Cp23", "Schisto: SEA", "Oncho: Ov16", "Dengue: NS1-3",
                    "Dengue: NS1-4", "Entamoeba: LecA", "Zika: NS1", "Dengue: NS1-1", "Giardia: VSP3",
                    "LF: Bm14", "LF: Wb123", "Trachoma: CT694", "Giardia: VSP5", "Chikungunya: E1",
                    "Crypto: Cp17", "Trachoma: pgp3", "Dengue: NS1-2", "Strongy: NIE", "Schisto: Sm25")

NTD_reorder <- c(20,  3, 11,  1, 12,
                  4, 16,  2, 10, 14,
                  7, 19, 17, 13, 15,
                  8,  9, 18,  5,  6)

NTD_short_names <- NTD_short_names[NTD_reorder]                 

NTD_data[,9:28] <- NTD_data[,(8 + NTD_reorder)]

NTD_long_names <- NTD_long_names[NTD_reorder]

colnames( NTD_data )[9:28] <- NTD_short_names 



NTD_pos <- NTD_data[which( NTD_data$sample %in% c("POS1", "POS2", "POS3") ),]



NTD_neg <- read.csv("20plexNTD_Covidoise_NEGcontrols_plate2_session2.csv")

NTD_neg <- NTD_neg[,c(1, 21, 4, 12, 2, 13, 
                          5,17,  3,11, 15,
                          8,20, 18,14, 16,
                          9,10, 19, 6,  7) ]

colnames(NTD_neg)[2:21] <- colnames(NTD_data)[9:28] 
                         

m2sd_cuts <- apply(X=log(NTD_neg[,-1]), MARGIN=2, FUN=mean) +
             2*apply(X=log(NTD_neg[,-1]), MARGIN=2, FUN=sd)  

m2sd_cuts <- exp(m2sd_cuts)


m3sd_cuts <- apply(X=log(NTD_neg[,-1]), MARGIN=2, FUN=mean) +
             3*apply(X=log(NTD_neg[,-1]), MARGIN=2, FUN=sd)  

m3sd_cuts <- exp(m3sd_cuts)






##########################################
##########################################
##          ##                          ## 
##   ####   ##   ####  #     # #     #  ##
##  ##  ##  ##  ##     ##   ## ##   ##  ##
##     ##   ##  ## ### ####### #######  ##
##    ##    ##  ##  ## ## # ## ## # ##  ##
##   #####  ##   ####  ##   ## ##   ##  ##
##          ##                          ##
##########################################
##########################################

library(MASS)
library(compiler)
library("rstan")

N_bins <- 30
Ab_min <- log(3)
Ab_max <- log(30000)

Ab_bin_edges <- seq(from=Ab_min, to=Ab_max, length=N_bins+1)

Ab_bin_width <- Ab_bin_edges[2] - Ab_bin_edges[1]


N_par <- 5

N_chains <- 4
chain_cols <- c("#0033FF", "#0066FF", "#0099FF", "#0099CC")

N_iter <- 5000
N_warmup <- 500

N_witer <-  N_iter - N_warmup

N_subset <- 500



################################################### 
## 2.1 Put Senegal and France data in bins (anti-pgp)


Sen_bins_pgp <- rep(NA, length=N_bins)

for(i in 1:N_bins)
{
	index <- which( log(NTD_data$pgp) >  Ab_bin_edges[i] &
                      log(NTD_data$pgp) <= Ab_bin_edges[i+1] )

	Sen_bins_pgp[i] <- length(index)	
}


Fra_bins_pgp <- rep(NA, length=N_bins)

for(i in 1:N_bins)
{
	index <- which( log(NTD_neg$pgp) >  Ab_bin_edges[i] &
                      log(NTD_neg$pgp) <= Ab_bin_edges[i+1] )

	Fra_bins_pgp[i] <- length(index)	
}





################################################### 
## 2.2 Fit GMM to trachoma data (anti-pgp)

N_sam <- nrow(NTD_neg) + nrow(NTD_data)

print("Creating data for model")  
list_data_model <- list(N_unknown  = N_sam,     ## Number of unknown samples
                        Ab_unknown = log(c(NTD_neg$pgp, NTD_data$pgp)),  ## Vector of unknown Ab values 
                        Ab_min     = Ab_min,                             ## Maximum Ab value
                        Ab_max     = Ab_max )                            ## Maximum Ab value



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
                   y = GMM_chains[1:N_witer+ (m-1)*N_witer,j], 
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
# Compare model to anti-pgp data

GMM_chains_subset <- GMM_chains[seq(from=1, to=nrow(GMM_chains), length=N_subset),]

Ab_seq <- seq(from=Ab_min, to=Ab_max, length=100)

Neg_dist_pgp <- matrix(NA, nrow=N_subset, ncol=length(Ab_seq))
Pos_dist_pgp <- matrix(NA, nrow=N_subset, ncol=length(Ab_seq))

	

for(i in 1:N_subset)
{
	Neg_dist_pgp[i,] <- (1-GMM_chains_subset[i,1])*
                          dnorm(Ab_seq, mean=GMM_chains_subset[i,2], 
                                     sd=GMM_chains_subset[i,3])
 
	Pos_dist_pgp[i,] <- GMM_chains_subset[i,1]*
                         dnorm(Ab_seq, mean=GMM_chains_subset[i,2]+GMM_chains_subset[i,4], 
                                    sd=GMM_chains_subset[i,5])
}

Neg_dist_quant_pgp <- matrix(NA, nrow=3, ncol=100)
Pos_dist_quant_pgp <- matrix(NA, nrow=3, ncol=100)

for(j in 1:length(Ab_seq))
{
	Neg_dist_quant_pgp[,j] <- quantile( Neg_dist_pgp[,j], prob=c(0.025, 0.5, 0.975) )
	Pos_dist_quant_pgp[,j] <- quantile( Pos_dist_pgp[,j], prob=c(0.025, 0.5, 0.975) )
}


#####################################
## Calulate cut-off

GMM_cut_pgp <- GMM_chains_subset[,2] + 2*GMM_chains_subset[,3]

GMM_cut_median_pgp <- median( GMM_cut_pgp )


GMM_cut_pgp <- density( GMM_cut_pgp )



################################################### 
###################################################         
# Plot fit of model to anti-pgp data

par(mfrow=c(1,1))


plot(x=1e10, y=1e10, 
pch=15, cex=1, col="orangered",
xlim=log(c(2,50000)), ylim=c(0,200), 
#xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n',
xlab="antibody levels (MFI)", ylab="",  
main="" )



for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=c(0, 0, Fra_bins_pgp[i], Fra_bins_pgp[i]), 
		   col="orangered" )
}

for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=Fra_bins_pgp[i] + c(0, 0, Sen_bins_pgp[i], Sen_bins_pgp[i]), 
		   col="royalblue" )
}


mtext(side = 2, line = 2.5, 
cex=1.75, 
text="anti-pgp antibody levels")



points(x=Ab_seq, y=Ab_bin_width*N_sam*Neg_dist_quant_pgp[2,], 
       type='l', lwd=2, col="blue")

polygon( x=c( Ab_seq, rev(Ab_seq) ),
         y=Ab_bin_width*c( N_sam*Neg_dist_quant_pgp[1,], rev(N_sam*Neg_dist_quant_pgp[3,]) ), 
	   col=rgb(0,0,1,0.5), border=NA)

points(x=Ab_seq, y=Ab_bin_width*N_sam*Pos_dist_quant_pgp[2,], 
       type='l', lwd=2, col="red")

polygon( x=c( Ab_seq, rev(Ab_seq) ),
         y=Ab_bin_width*c( N_sam*Pos_dist_quant_pgp[1,], rev(N_sam*Pos_dist_quant_pgp[3,]) ), 
	   col=rgb(1,0,0,0.5), border=NA)


polygon( x=c( GMM_cut_pgp$x, rev(GMM_cut_pgp$x) ),
         y=(200/max(GMM_cut_pgp$y))*c( rep(0,length(GMM_cut_pgp$y)), rev(GMM_cut_pgp$y) ), 
	   col=rgb(0,1,0,0.5), border=NA)

points( x = rep( GMM_cut_median_pgp, 2),
        y = (200/max(GMM_cut_pgp$y))*c(0, GMM_cut_pgp$y[which.min(abs(GMM_cut_pgp$x - GMM_cut_median_pgp))] ),
        type='l', col="green", lwd=2 )




################################################### 
## 2.3 Put Senegal and France data in bins (anti-CT694)


Sen_bins_CT694 <- rep(NA, length=N_bins)

for(i in 1:N_bins)
{
	index <- which( log(NTD_data$CT694) >  Ab_bin_edges[i] &
                      log(NTD_data$CT694) <= Ab_bin_edges[i+1] )

	Sen_bins_CT694[i] <- length(index)	
}


Fra_bins_CT694 <- rep(NA, length=N_bins)

for(i in 1:N_bins)
{
	index <- which( log(NTD_neg$CT694) >  Ab_bin_edges[i] &
                      log(NTD_neg$CT694) <= Ab_bin_edges[i+1] )

	Fra_bins_CT694[i] <- length(index)	
}





################################################### 
## 2.2 Fit GMM to trachoma data (anti-pgp)

N_sam <- nrow(NTD_neg) + nrow(NTD_data)

print("Creating data for model")  
list_data_model <- list(N_unknown  = N_sam,     ## Number of unknown samples
                        Ab_unknown = log(c(NTD_neg$CT694, NTD_data$CT694)),  ## Vector of unknown Ab values 
                        Ab_min     = Ab_min,                             ## Maximum Ab value
                        Ab_max     = Ab_max )                            ## Maximum Ab value



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
                   y = GMM_chains[1:N_witer+ (m-1)*N_witer,j], 
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
# Compare model to anti-pgp data

GMM_chains_subset <- GMM_chains[seq(from=1, to=nrow(GMM_chains), length=N_subset),]

Ab_seq <- seq(from=Ab_min, to=Ab_max, length=100)

Neg_dist_CT694 <- matrix(NA, nrow=N_subset, ncol=length(Ab_seq))
Pos_dist_CT694 <- matrix(NA, nrow=N_subset, ncol=length(Ab_seq))

	

for(i in 1:N_subset)
{
	Neg_dist_CT694[i,] <- (1-GMM_chains_subset[i,1])*
                          dnorm(Ab_seq, mean=GMM_chains_subset[i,2], 
                                     sd=GMM_chains_subset[i,3])
 
	Pos_dist_CT694[i,] <- GMM_chains_subset[i,1]*
                         dnorm(Ab_seq, mean=GMM_chains_subset[i,2]+GMM_chains_subset[i,4], 
                                    sd=GMM_chains_subset[i,5])
}

Neg_dist_quant_CT694 <- matrix(NA, nrow=3, ncol=100)
Pos_dist_quant_CT694 <- matrix(NA, nrow=3, ncol=100)

for(j in 1:length(Ab_seq))
{
	Neg_dist_quant_CT694[,j] <- quantile( Neg_dist_CT694[,j], prob=c(0.025, 0.5, 0.975) )
	Pos_dist_quant_CT694[,j] <- quantile( Pos_dist_CT694[,j], prob=c(0.025, 0.5, 0.975) )
}


#####################################
## Calulate cut-off

GMM_cut_CT694 <- GMM_chains_subset[,2] + 2*GMM_chains_subset[,3]

GMM_cut_median_CT694 <- median( GMM_cut_CT694 )


GMM_cut_CT694 <- density( GMM_cut_CT694 )



################################################### 
###################################################         
# Plot fit of model to anti-CT694 data

par(mfrow=c(1,1))


plot(x=1e10, y=1e10, 
pch=15, cex=1, col="orangered",
xlim=log(c(2,50000)), ylim=c(0,200), 
#xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n',
xlab="antibody levels (MFI)", ylab="",  
main="" )



for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=c(0, 0, Fra_bins_CT694[i], Fra_bins_CT694[i]), 
		   col="orangered" )
}

for(i in 1:N_bins)
{
	polygon( x=c(Ab_bin_edges[i], Ab_bin_edges[i+1], Ab_bin_edges[i+1], Ab_bin_edges[i]),
               y=Fra_bins_CT694[i] + c(0, 0, Sen_bins_CT694[i], Sen_bins_CT694[i]), 
		   col="royalblue" )
}


mtext(side = 2, line = 2.5, 
cex=1.75, 
text="anti-CT694 antibody levels")



points(x=Ab_seq, y=Ab_bin_width*N_sam*Neg_dist_quant_CT694[2,], 
       type='l', lwd=2, col="blue")

polygon( x=c( Ab_seq, rev(Ab_seq) ),
         y=Ab_bin_width*c( N_sam*Neg_dist_quant_CT694[1,], rev(N_sam*Neg_dist_quant_CT694[3,]) ), 
	   col=rgb(0,0,1,0.5), border=NA)

points(x=Ab_seq, y=Ab_bin_width*N_sam*Pos_dist_quant_CT694[2,], 
       type='l', lwd=2, col="red")

polygon( x=c( Ab_seq, rev(Ab_seq) ),
         y=Ab_bin_width*c( N_sam*Pos_dist_quant_CT694[1,], rev(N_sam*Pos_dist_quant_CT694[3,]) ), 
	   col=rgb(1,0,0,0.5), border=NA)


polygon( x=c( GMM_cut_CT694$x, rev(GMM_cut_CT694$x) ),
         y=(200/max(GMM_cut_CT694$y))*c( rep(0,length(GMM_cut_CT694$y)), rev(GMM_cut_CT694$y) ), 
	   col=rgb(0,1,0,0.5), border=NA)

points( x = rep( GMM_cut_median_CT694, 2),
        y = (200/max(GMM_cut_CT694$y))*c(0, GMM_cut_CT694$y[which.min(abs(GMM_cut_CT694$x - GMM_cut_median_CT694))] ),
        type='l', col="green", lwd=2 )




#################################
#################################
##                             ## 
##  ####   ####  ####  ######  ##
##  ## ##   ##  ##       ##    ##
##  ##  ##  ##   ####    ##    ## 
##  ## ##   ##      ##   ##    ##
##  ####   ####  ####    ##    ##
##                             ##
#################################
#################################


SP_cuts <- m2sd_cuts
SP_cuts[c(3,4)] <- m3sd_cuts[c(3,4)]
SP_cuts[13] <- exp(GMM_cut_median_pgp)
SP_cuts[14] <- exp(GMM_cut_median_CT694)

SP_cuts[c(7,8)] <- NA


#################################
## Calculate binned data

N_bins <- 20

MFI_bin_edges <- exp(seq(from=log(3), to=log(30000), length=N_bins+1))

Dielmo_MFI <- matrix(NA, nrow=N_bins, ncol=N_NTD)
colnames(Dielmo_MFI) <- NTD_short_names

for(j in 1:N_NTD)
{
	for(i in 1:N_bins)
	{
		index <- which( NTD_data$village == "Dielmo" &
                            NTD_data[,8+j] > MFI_bin_edges[i] &
                            NTD_data[,8+j] <= MFI_bin_edges[i+1] )

		Dielmo_MFI[i,j] <- length(index)	
	}
}


Ndiop_MFI <- matrix(NA, nrow=N_bins, ncol=N_NTD)
colnames(Ndiop_MFI) <- NTD_short_names

for(j in 1:N_NTD)
{
	for(i in 1:N_bins)
	{
		index <- which( NTD_data$village == "Ndiop" &
                            NTD_data[,8+j] > MFI_bin_edges[i] &
                            NTD_data[,8+j] <= MFI_bin_edges[i+1] )

		Ndiop_MFI[i,j] <- length(index)	
	}
}





Neg_MFI <- matrix(NA, nrow=N_bins, ncol=N_NTD)
colnames(Neg_MFI) <- NTD_short_names

for(j in 1:N_NTD)
{
	for(i in 1:N_bins)
	{
		index <- which( NTD_neg[,1+j] > MFI_bin_edges[i] &
                            NTD_neg[,1+j] <= MFI_bin_edges[i+1] )

		Neg_MFI[i,j] <- length(index)	
	}
}

Neg_MFI <- 3*Neg_MFI


#################################
## Plot data distribution

line_seq_y = c(3, 10, 30, 100, 300, 1000, 3000, 10000, 30000)
line_seq_x = c(0, 100, 200, 300, 400, 500, 600, 700)


main.size = 1.5
lab.size  = 1.25
axis.size = 0.7
line.size = 1
text.size = 1



tiff( file="Figure2_AB_dist_20panel.tif", width=30, height=20, units="cm", res=500)

lay.mat <- rbind( c( 1, 2, 3, 4, 5), 
                  c( 6, 7, 8, 9,10), 
                  c(11,12,13,14,15),
                  c(16,17,18,19,20), 
                  c(21,21,21,21,21) )
layout(lay.mat, heights=c(6,6,6,6,1))
layout.show(21)

par(mar=c(4,3,2,0.5))
par(mgp=c(2, 0.5,0))

for(j in 1:N_NTD)
{
	plot(x=1e10, y=1e10, 
	pch=15, cex=1, col="orangered",
	xlim=c(10,30000), ylim=c(0,500), log="x",
	xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n',
	xlab="antibody levels (MFI)", ylab="samples",  
	main=NTD_long_names[j],
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
		polygon( x=c(MFI_bin_edges[i], MFI_bin_edges[i+1], MFI_bin_edges[i+1], MFI_bin_edges[i]),
                     y=c(0, 0, Neg_MFI[i,j], Neg_MFI[i,j]), 
		         col="royalblue" )
	}

	for(i in 1:N_bins)
	{
		polygon( x=c(MFI_bin_edges[i], MFI_bin_edges[i+1], MFI_bin_edges[i+1], MFI_bin_edges[i]),
                     y=Neg_MFI[i,j] + c(0, 0, Dielmo_MFI[i,j], Dielmo_MFI[i,j]), 
		         col="tan3" )
	}

	for(i in 1:N_bins)
	{
		polygon( x=c(MFI_bin_edges[i], MFI_bin_edges[i+1], MFI_bin_edges[i+1], MFI_bin_edges[i]),
                     y=Neg_MFI[i,j] + Dielmo_MFI[i,j] + c(0, 0, Ndiop_MFI[i,j], Ndiop_MFI[i,j]), 
		         col="orangered" )
	}


	points( x=c(SP_cuts[j], SP_cuts[j]), y=c(0,1000), 
           	  type='l', lwd=3, col="green3" )
	
	
	axis(1, at=c(10, 30, 100, 300, 1000, 3000, 10000, 30000), 
              label=c(10, 30, 100, 300, 1000, 3000, 10000, 30000), cex.axis=axis.size )

	axis(2, at=c(0, 100, 200, 300, 400, 500), label=c("0", "100", "200", "300", "400", "500"), 
           cex.axis=1.4*axis.size, las=2)

}


##############################################
## Add the legend to the bottom of the plot ##
##############################################

par(mar = c(0,0,0,0))
plot.new()
legend(x='center', 
       legend = c( "France", "Dielmo", "Ndiop"  ), 
       fill = c("royalblue", "tan3", "orangered"),
       border = c("royalblue", "tan3", "orangered"),
	 ncol=3, bty="n", cex=2)


dev.off()



save.image("Figure2_AB_dist.RData")






