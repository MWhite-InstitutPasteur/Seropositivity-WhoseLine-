

load("C:\\U\\MultiAssay\\WhoseLine\\Figures_v2\\Figure2_AB_dist.RData")

###########################################
###########################################
##                                       ##
##   ####  #####     ####   ####  #####  ##
##  ##     ##  ##   ##  ## ##     ##     ##
##   ####  #####    ###### ## ### ####   ##
##      ## ##       ##  ## ##  ## ##     ##
##   ####  ##       ##  ##  ####  #####  ##
##                                       ##
###########################################
###########################################

library(binom)

###############################################
## 0.1 Prepare data for plotting

age_bins     <- c(0, 5, 10, 15, 20, 30, 40, 50, 60, 100)
age_bins_mid <- 0.5*( age_bins[-1] + age_bins[-length(age_bins)] )
 
N_bins <- length(age_bins) - 1 


SP_Dielmo_bins_med <- matrix(NA, nrow=N_bins, ncol=N_NTD)
colnames(SP_Dielmo_bins_med) <- NTD_short_names

SP_Dielmo_bins_lwr <- matrix(NA, nrow=N_bins, ncol=N_NTD)
colnames(SP_Dielmo_bins_lwr) <- NTD_short_names

SP_Dielmo_bins_upr <- matrix(NA, nrow=N_bins, ncol=N_NTD)
colnames(SP_Dielmo_bins_upr) <- NTD_short_names


for(j in 1:N_NTD)
{
	for(i in 1:N_bins)
	{
		index <- which( NTD_data[,5] > age_bins[i] & 
                            NTD_data[,5] <= age_bins[i+1] &
                            NTD_data[,7] == "Dielmo"    ) 
	
		temp  <- NTD_data[index,8+j]


		SP_est <- as.numeric(as.vector(
            	                	 binom.confint( length(which(temp > SP_cut_select[j])), length(temp), method="wilson")[1,4:6]
            	                ))

		SP_Dielmo_bins_med[i,j] <- SP_est[1]
		SP_Dielmo_bins_lwr[i,j] <- SP_est[2]
		SP_Dielmo_bins_upr[i,j] <- SP_est[3]
	}
}




SP_Ndiop_bins_med <- matrix(NA, nrow=N_bins, ncol=N_NTD)
colnames(SP_Ndiop_bins_med) <- NTD_short_names

SP_Ndiop_bins_lwr <- matrix(NA, nrow=N_bins, ncol=N_NTD)
colnames(SP_Ndiop_bins_lwr) <- NTD_short_names

SP_Ndiop_bins_upr <- matrix(NA, nrow=N_bins, ncol=N_NTD)
colnames(SP_Ndiop_bins_upr) <- NTD_short_names


for(j in 1:N_NTD)
{
	for(i in 1:N_bins)
	{
		index <- which( NTD_data[,5] > age_bins[i] & 
                            NTD_data[,5] <= age_bins[i+1] &
                            NTD_data[,7] == "Ndiop"    ) 
	
		temp  <- NTD_data[index,8+j]


		SP_est <- as.numeric(as.vector(
            	                	 binom.confint( length(which(temp > SP_cut_select[j])), length(temp), method="wilson")[1,4:6]
            	                ))

		SP_Ndiop_bins_med[i,j] <- SP_est[1]
		SP_Ndiop_bins_lwr[i,j] <- SP_est[2]
		SP_Ndiop_bins_upr[i,j] <- SP_est[3]
	}
}


#################################
## Plot data distribution

line_seq_x = c(0, 0.2, 0.4, 0.6, 0.8, 1)
line_seq_y = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)


main.size = 1.2
lab.size  = 1.5
axis.size = 0.9
line.size = 1.2
text.size = 1




tiff( file="Figure4_SP_curves_20panel.tif", width=30, height=20, units="cm", res=500)

lay.mat <- rbind( c( 1, 2, 3, 4, 5), 
                  c( 6, 7, 8, 9,10), 
                  c(11,12,13,14,15),
                  c(16,17,18,19,20), 
                  c(21,21,21,21,21) )
layout(lay.mat, heights=c(6,6,6,6,1))
layout.show(21)

par(mar=c(4,4,2,0.5))
par(mgp=c(2.25, 0.6,0))

for(j in 1:N_NTD)
{
	par(mgp=c(2.0, 0.6,0))

	plot(x=1e10, y=1e10, 
	pch=15, cex=1, col="orangered",
	xlim=c(0,82), ylim=c(0,1.02), 
	xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n',
	xlab="age (years)", ylab="",  
	main=NTD_long_names[j],
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size )

	par(mgp=c(2.5, 0.6,0))
	title( ylab="seroprevalence", cex.lab=lab.size )  


	for(i in 1:length(line_seq_x))
	{
		points(x=c(0.001,100000), y=rep(line_seq_x[i],2), type='l', lwd=1, col="grey", lty="dashed")
	}

	for(i in 1:length(line_seq_y))
	{
		points(x=rep(line_seq_y[i],2), y=c(0,100000), type='l', lwd=1, col="grey", lty="dashed")
	}

	if( j %in% c( 1:6, 9:20) )
	{

		points(x=age_bins_mid - 0.5, y=SP_Dielmo_bins_med[,j], 
		pch=15, cex=1.5, col="tan3")

		points(x=age_bins_mid - 0.5, y=SP_Dielmo_bins_med[,j], 
		type='l', lwd=2,  col="tan3")

		for(i in 1:N_bins)
		{
			arrows(x0=age_bins_mid[i] - 0.5, y0=SP_Dielmo_bins_lwr[i,j], 
		             x1=age_bins_mid[i] - 0.5, y1=SP_Dielmo_bins_upr[i,j], 
		             length=0.03, angle=90, code=3, col="tan3", lwd=1)	
		}

		points(x=age_bins_mid + 0.5, y=SP_Ndiop_bins_med[,j], 
		pch=15, cex=1.5,  col="orangered")

		points(x=age_bins_mid + 0.5, y=SP_Ndiop_bins_med[,j], 
		type='l', lwd=2,  col="orangered")

		for(i in 1:N_bins)
		{
			arrows(x0=age_bins_mid[i] + 0.5, y0=SP_Ndiop_bins_lwr[i,j], 
		             x1=age_bins_mid[i] + 0.5, y1=SP_Ndiop_bins_upr[i,j], 
		             length=0.03, angle=90, code=3, col="orangered", lwd=1)	
		}
	}


	axis(1, at=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), 
              label=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), cex.axis=axis.size )

	axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), 
              label=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=axis.size, las=2 )

}


##############################################
## Add the legend to the bottom of the plot ##
##############################################

par(mar = c(0,0,0,0))
plot.new()
legend(x='center', 
       legend = c( "Dielmo", "Ndiop" ), 
       fill = c("tan3", "orangered"),
       border = c("tan3", "orangered"),
	 ncol=2, bty="n", cex=2)


dev.off()






