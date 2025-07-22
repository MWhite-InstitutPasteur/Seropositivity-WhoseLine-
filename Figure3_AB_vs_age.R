


load("C:\\U\\MultiAssay\\WhoseLine\\Figures_v2\\Figure2_AB_dist.RData")



###########################################
###########################################
##                                       ##
##   ####  #####     ####   ####  #####  ##
##  ##  ## ##  ##   ##  ## ##     ##     ##
##  ###### #####    ###### ## ### ####   ##
##  ##  ## ##  ##   ##  ## ##  ## ##     ##
##  ##  ## #####    ##  ##  ####  #####  ##
##                                       ##
###########################################
###########################################


NTD_data$age[which(NTD_data$age > 80)] <- 80


###############################################
## 0.1 Prepare data for plotting

#################################
## Plot data distribution

line_seq_x = c(3, 10, 30, 100, 300, 1000, 3000, 10000, 30000)
line_seq_y = c(-15, -5, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)


main.size = 1.2
lab.size  = 1.0
axis.size = 0.75
line.size = 1.2
text.size = 1
point.size = 0.4


tiff( file="Figure3_AB_vs_age_20panel.tif", width=30, height=20, units="cm", res=500)


lay.mat <- rbind( c( 1, 2, 3, 4, 5), 
                  c( 6, 7, 8, 9,10), 
                  c(11,12,13,14,15),
                  c(16,17,18,19,20), 
                  c(21,21,21,21,21) )
layout(lay.mat, heights=c(6,6,6,6,1))
layout.show(21)

par(mar=c(2.75,3.25,2,0.5))
par(mgp=c(2.25, 0.6,0))

for(j in 1:N_NTD)
{
	par(mgp=c(2.25, 0.6,0))

	plot(x=1e10, y=1e10, 
	pch=15, cex=1, col="orangered",
	xlim=c(-20,80), ylim=c(10,30000), log="y",
	xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n',
	xlab="", ylab="antibody levels (MFI)",  
	main=NTD_long_names[j],
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size )


	par(mgp=c(1.5, 0.6,0))
	title(xlab="     age (years)" )

	for(i in 1:length(line_seq_x))
	{
		points(x=c(-100,100000), y=rep(line_seq_x[i],2), type='l', lwd=1, col="grey", lty="dashed")
	}

	for(i in 1:length(line_seq_y))
	{
		points(x=rep(line_seq_y[i],2), y=c(0.1,100000), type='l', lwd=1, col="grey", lty="dashed")
	}


	points( x = NTD_data[which(NTD_data[,7]=="Dielmo"),5], 
              y = NTD_data[which(NTD_data[,7]=="Dielmo"),8+j], 
	pch=19, cex=point.size, col="tan3")

	points( x = NTD_data[which(NTD_data[,7]=="Ndiop"),5], 
              y = NTD_data[which(NTD_data[,7]=="Ndiop"),8+j], 
	pch=19, cex=point.size, col="orangered")


	points( x = runif(n=nrow(NTD_neg), min=-15, max=-5), 
              y = NTD_neg[,1+j], 
	pch=19, cex=point.size, col="royalblue")

	

	points( x=c(-100,1000), y=c(SP_cut_select[j], SP_cut_select[j]), 
              type='l', lwd=2, col="green3" )

	axis(1, at=c(0, 10, 20, 30, 40, 50, 60, 70, 80), 
              label=c(0, 10, 20, 30, 40, 50, 60, 70, 80), cex.axis=axis.size )

	axis(1, at=c(-15, -5), 
              label=c("       neg", ""), cex.axis=axis.size )

	axis(2, at=c(10, 30, 100, 300, 1000, 3000, 10000, 30000), 
              label=c(10, 30, 100, 300, 1000, 3000, 10000, 30000), 
              cex.axis=axis.size, las=2 )

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




