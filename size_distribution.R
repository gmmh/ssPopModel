#GMH 04.19.15
# local machine version

library(popcycle)

cruise <- "DeepDOM"
home <- "~/Desktop/Cruises/DeepDOM_2013/seaflow/Cell_Division"
set.project.location("/Volumes/gwennm/popcycle")
flag.file <- "/Volumes/gwennm/popcycle/sqlite/flag_file.txt"
phyto <- "synecho"
para <- "fsc_small"
n.breaks <- 57
time.interval <- 60 #minutes

#################################
### Get the range of 'para' for 'phyto' ###
#################################

# Get the time range
stat <- get.stat.table()
stat$time <- as.POSIXct(stat$time,format="%FT%T",tz='GMT')
phyto.stat <- subset(stat, pop == phyto)
flag <- read.table(flag.file, header =T) #flag info to take out bad files
flag.good <- subset(flag, flag==0)
phyto.stat <- subset(phyto.stat, file %in% flag.good$file)
time.range <- range(phyto.stat$time)
time <- seq(time.range[1],time.range[2] , by=60*time.interval) # cut the time series according to time interval

# Get the range of 'para' for 'phyto' 
#changed this to use the gate ranges saved in the log rather than dowloading opp files
log.gate <- list.files(log.gate.location, phyto)
para.phyto <- as.data.frame(matrix(nrow=length(log.gate), ncol=2))
for(i in 1:length(log.gate)){
	log.table <- read.table(paste0(log.gate.location,"/", log.gate[i]), sep=",", header=T)
	para.phyto[i,] <- range(log.table[,para])
}
#para.phyto <- get.opp.by.date(time.range[1], time.range[2], pop=phyto, channel=para)[,para]
para.range <- c(min(para.phyto[,1]), max(para.phyto[,2]))


#########################
### SMOOTH bead signal  ###
#########################

# Get the beads data
#m.beads <- median(subset(stat, pop =='beads' | time > time.range[1] & time < time.range[2])[,para])
# Plot the light scattering of beads over time
plot.time(stat, popname='beads',param='fsc_small')

# BInned data according to 'time'interval'
beads <- subset(stat, pop=='beads' & time > time.range[1] & time < time.range[2])
time.binned <- cut(beads$time, time, labels=F)
para.beads.binned <- as.vector(tapply(beads$fsc_small, time.binned, median))
time.beads.binned <- as.POSIXct(as.vector(tapply(beads$time, time.binned, mean)), origin="1970-01-01", tz='GMT')
points(time.beads.binned , para.beads.binned, type='o', col=2, pch=16)

# Smooth the data
spar <- 0.45 # smooothing parameter, the higher the more smoothing is applied.
smooth <- smooth.spline(time.beads.binned, para.beads.binned,spar=spar)
smooth.para.beads.binned <- spline(as.POSIXct(smooth$x,origin="1970-01-01",tz="GMT"), smooth$y, xout=as.POSIXct(smooth$x,origin="1970-01-01",tz="GMT"))
lines(smooth.para.beads.binned ,col=3,lwd=3) # visualize the smooth data and re-adjust the spar parameter if necessary.


###############################################################
### Generate SIZE distribution, binned into 'n.breaks' for each time interval ###
###############################################################

Vhist <- Ndist  <- Time <- NULL
for( t in time){

    #get the opp for phyto
    t <- as.POSIXct(t, origin="1970-01-01", tz='GMT'); print(paste(t))
    pop <- try(get.opp.by.date(t, t+60*time.interval, pop=phyto, channel=para))

    if(class(pop) == "try-error" | nrow(pop) < 10){
        print("not enough opp, skipping hour")
        next
        }

    # get Beads signal
    m.beads1 <- smooth.para.beads.binned$y[which(smooth.para.beads.binned$x > as.numeric(t))[1]]
    m.beads2 <- smooth.para.beads.binned$y[max(which(smooth.para.beads.binned$x < as.numeric(t)+60*time.interval))]
    m.beads <- mean(c(m.beads1, m.beads2), na.rm=T)
    # get opp/evt ratio (used to calculate Ndist)
  ## opp.evt.ratio <- median(get.opp.evt.ratio.by.date(t, t+60*time.interval)$ratio)

    # create the size distribution of normalized forward scatter, using a Gaussian filter
    dens <- density(log2(pop[,para]/m.beads), n=n.breaks,from=log2(para.range[1]/m.beads) , to=log2(para.range[2]/m.beads), kernel='gaussian')
    freq.dist <-  dens$y*diff(dens$x)[1] # convert density to frequency
    freq.dist <- freq.dist/sum(freq.dist) # normailize the frequency to 1
        Vhist <- data.frame(cbind(Vhist, freq.dist))
     size.dist <- round(freq.dist * nrow(pop))
        Ndist <- data.frame(cbind(Ndist, size.dist))
        Time <- c(Time, t)

       } 

#################################################################
### CONVERT normalized forward SCATTER by 1 micron beads to VOLUME ###
#################################################################       
norm.fsc <- 2^dens$x
   if(phyto == "synecho" | phyto == "pico" | phyto == "prochloro"){
        volbins <- round(10^(0.524*log10(norm.fsc) + 0.283),3)
        # Size$volume <- 10^(0.5*log10(Size$stages/Size$fsc_beads))# MIE THEORY
        }
    
  
    if(phyto == "crypto" | phyto =='nano'){
      #Size$volume <- 10^(0.75*log10(Size$stages/Size$fsc_beads)) # MIE THEORY
      volbins <- round(10^(1.2384*log10(norm.fsc) + 1.003),3)
    }
  


####################
### SAVE the output ###
####################

colnames(Vhist) <- colnames(Ndist) <- as.character(Time)
rownames(Vhist) <- rownames(Ndist) <- volbins
   
distribution <- list()
    distribution[[1]] <- Vhist
    distribution[[2]] <- Ndist

save(distribution, file=paste(home,"/",phyto,"_dist_Ncat",n.breaks,"_",cruise,sep=""))
print(paste("saving ", home,"/",phyto,"_dist_Ncat",n.breaks,"_",cruise,sep=""))



#########################
### SHOW size distribution ###
#########################
# home <- "/Volumes/gwennm/DeepDOM/Cell_Division"
# load(file=paste(home,"/",phyto,"_dist_Ncat",n.breaks,"_",cruise,sep=""))

library(rgl)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

para <- distribution[[1]]
percentile <- cut(unlist(para), 100)

# in linear scale
plot3d(rep(as.numeric(row.names(para)), dim(para)[2]), 
            rep(as.numeric(colnames(para)), each=dim(para)[1]) , 
            unlist(para), 
            col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")

# in log scale
plot3d(log2(rep(as.numeric(row.names(para)), dim(para)[2])), 
            rep(as.numeric(colnames(para)), each=dim(para)[1]) , 
            unlist(para), 
            col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")



