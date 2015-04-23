# [gwennm@bloom Cell_Division] 
# code to run on cluster from this folder, use the following in the command line
# for i in $(seq 0 1 24); do echo "Rscript ~/DeepDOM/ssPopModel/run_model.R $i prochloro DeepDOM ~/DeepDOM/Cell_Division" | qsub -lwalltime=78:00:00,nodes=1:ppn=1 -N proGR$i -d.; done

#  library(rgl)
library(DEoptim)
library(zoo)

#home <- "/Volumes/gwennm/DeepDOM/Cell_division/"; folder <- NULL; cruise <- "DeepDOM"
home <- '~/DeepDOM/Cell_Division/'; folder <- NULL

#source('/Volumes/gwennm/DeepDOM/ssPopModel/functions_model.R', chdir = TRUE)
source('~/DeepDOM/ssPopModel/functions_model.R', chdir = TRUE)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",	"#FF7F00", "red", "#7F0000"))


args <- commandArgs(TRUE)
t <- as.numeric(args[1])
phyto <- as.character(args[2])
cruise <- as.character(args[3])
outdir <- as.character(args[4])



###############################
m= 57
time.interval <- 60*60*24 #number of seconds in 1 day
#m <- 2^6 # number of size class
###############################




###############
## REFERENCE ##
###############
# library(R.matlab); mat <- readMat("/Users/francois/Documents/DATA/SeaFlow/Cell_Division/Matlab/day733320data.mat"); res <- readMat("/Users/francois/Documents/DATA/SeaFlow/Cell_Division/Matlab/results.mat")
# volbins <- mat$volbins[1,]
# Edata <- mat$Edata
# V.hists <- mat$Vhists
# N.dist <- mat$N.dist
# Vproj <- res$Vproj
# para <- V.hists; percentile <- cut(para, 100); plot3d(log(rep(volbins , breaks)), rep(1:ncol(para), each=nrow(para)), z=matrix(para), col=jet.colors(100)[percentile], type='l', lwd=6)


	##############
	## PAR DATA ##
	##############
	
	Par.path <- paste(home, folder,"PAR_",cruise,sep="")
	Par <- read.csv(Par.path, sep=",")
	Par$time <- as.POSIXct(Par$time, tz="GMT")
	Par$num.time <- as.numeric(Par$time)


	#######################	
	## SIZE DISTRIBUTION ##
	#######################

	# t <- 0	
	# phyto <- "prochloro"
   
    print(paste("time delay:", t))
	print(paste("phytoplankton population:",phyto))
	
	load(paste(home,folder, phyto,"_dist_Ncat",m,"_",cruise,sep=""))
	Vhists <- distribution[[1]]
	Vhists <- sweep(Vhists, 2, colSums(Vhists), '/') # Normalize each column of VHists to 1
	N_dist <- distribution[[2]]

	volbins <- as.numeric(row.names(Vhists))
			sizebins <- 2*(volbins*3/(pi*4))^(1/3)# to check the actual diameter
	volbins <- volbins/max(volbins) # to make sure values are never > 1, for compatibility issue with the Delta function

	time.numc <- as.numeric(colnames(Vhists))	
	time <- as.POSIXct(time.numc, origin="1970-01-01" ,tz="GMT")
	time.day <- seq(range(time)[1], range(time)[2], by=time.interval)	
	n.day <- round(diff(range(time))); print(paste("Number of days in the dataset:",n.day))

	# para <- Vhists; percentile <- cut(unlist(para), 100); plot3d(log(rep(as.numeric(row.names(para)), dim(para)[2])), rep(as.numeric(colnames(para)), each=dim(para)[1]) , Vhists , col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")
	
	

	##############################
	## RUN size.model.functions ##
	##############################

	resol <-  60 # number of minutes per interval
	breaks <- 25*60/resol

	model <- array(NA, dim=c(4,1))

# t <- 1

	for(i in time.day){
		#i <- time.day[1]
		start <- i + t*60*60
		end <- start + time.interval
		print(paste("calculating growth projection from ",start , "to",end))
	
	
	#plot(Par$time, Par$par, type='o'); points(c(start, end),c(0,0), col='red',pch=16, cex=2)

		### SELECT SIZE DISTRIBUTION for DAY i
		start.i <- findInterval(start, time)
		end.i <- findInterval(end, time)
		
		print(paste("the time series has", end.i-start.i +1, "/24 data points"))
		if(end.i-start.i+1 < 18){
			print(paste("Note enough data points, skipping to next 24hr-period"))
			next
		}
		V.hists <- Vhists[,c(start.i:end.i)]
		N.dist <- N_dist[,c(start.i:end.i)]
	    # para <- V.hists; percentile <- cut(unlist(para), 100)
	    #plot3d(log(rep(as.numeric(row.names(para))), dim(para)[2]), rep(as.numeric(colnames(para), each=dim(para)[1])), unlist(para) , col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")


		### SELECT PAR corresponding to each sample
		interp.par.mean <- approx(Par[,"time"], Par[,"par"], xout=time[start.i:end.i], ties=mean)
		Edata <- matrix(cbind(interp.par.mean$x, interp.par.mean$y), ncol=2)
        
		
		### RUN size.class.model_functions
		proj <- try(determine.opt.para(V.hists=V.hists,N.dist=N.dist,Edata=Edata,volbins=volbins))
		
		#para <- proj$Vproj; percentile <- cut(unlist(para), 100); plot3d(log(rep(volbins, 24)), rep(1:ncol(para), each=nrow(para)), z=matrix(para), col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")
		
		if(class(proj) !='try-error'){
		model <- matrix(cbind(as.array(model), as.array(proj)), nrow=4,ncol=ncol(model)+1)
	    save(model, file=paste(outdir, "/",phyto,"_modelHD_growth_",cruise,"_Ncat",m,"_t",t, sep=""))

	  }else{print("error during optimization")}
}
