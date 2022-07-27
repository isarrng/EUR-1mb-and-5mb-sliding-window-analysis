

library(vegan)
library(snpStats)

### Inputs for analysis, ALSO INPUT SAME PLINK FILES BELOW ###

plink.gen=read.plink("individuals_m45_g50_18.bed",
                     "individuals_m45_g50_18.bim",
                     "individuals_m45_g50_18.fam")
GEO.DIST.MAT <- read.table("ind_m45_g50_Geo.txt", head = T) # Loading the Geographic distance matrix
data <- read.table("ind_m45_g50_dates.txt", head = T) # Loading sample date information
nindividuals <- length(data[,1])

total.gen.length.from=min(plink.gen$map$position)
total.gen.length.to=max(plink.gen$map$position)
block.length=5000000
min.snps.in.block=100



############################### Genetic Blocks #################################


plink.snps=plink.gen$map$position
snp.order=sort(plink.snps)

blocks=split(snp.order, cut(plink.gen$map$position, seq(total.gen.length.from, 
                                                        total.gen.length.to, 
                                                        block.length)))
bsize=sapply(blocks, length)
blocks=blocks[bsize>=min.snps.in.block]
block.nr=length(blocks)
b.names=names(blocks)
block.array=array(data=NA, dim=c(nindividuals,nindividuals,length(blocks)),
                  dimnames=list(NULL,NULL,b.names)) #Block array for Loog analysis#


#Genetic distance of blocks#

for (a in 1:block.nr){
  print(paste("block nr ", a))
  block.x=blocks[a]
  names(block.x)[names(block.x)==b.names[a]]="positions"
  block.x$names=plink.gen$map$snp.name[match(block.x$positions,plink.gen$map$position)]
  b.df=data.frame(block.x$names,block.x$positions)
  b.df=b.df[!duplicated(b.df),]
  
  plink.block.x=read.plink("individuals_m45_g50_18.bed",
                           "individuals_m45_g50_18.bim",
                           "individuals_m45_g50_18.fam",select.snps=b.df$block.x.names)
  
  block.ibs=ibsCount(plink.block.x$genotypes,
                     uncertain=FALSE)
  
  block.ibs.dist=ibsDist(block.ibs)
  
  block.ibs.dist.m=as.matrix(block.ibs.dist)
  
  block.array[,,a]=block.ibs.dist.m
  
}

# Drop individuals with NA values in genetic distance matrix #

test.GEO.DIST.MAT=GEO.DIST.MAT
for (s in 1:block.nr){
  print(paste("run nr ", s))
  if (sum(is.na(block.array[,,s]))>0){
    ind.positions.r=which(rowSums(is.na(block.array[,,s]))<1)
    ind.positions.c=which(colSums(is.na(block.array[,,s]))<1)
    ind.positions=intersect(ind.positions.r,ind.positions.c)
    block.array.1=array(data=NA, dim=c(length(ind.positions),
                                       length(ind.positions),block.nr),
                        dimnames=list(NULL,NULL,b.names))
    for (d in 1:block.nr){
      print(paste("!: ", d))
      block.array.1[,,d]=block.array[ind.positions,ind.positions,d]
    }
    data=data[ind.positions,]
    GEO.DIST.MAT=GEO.DIST.MAT[ind.positions,ind.positions]
    block.array=block.array.1
  }
}

diff.test=setdiff(colnames(test.GEO.DIST.MAT),colnames(GEO.DIST.MAT))

nindividuals <- length(data[,1])

############################## Loog Analysis ###################################
tstart=date()
############################## START ##############################
nr.run <- 1 # Number of resampling runs
window.moves = 65 # Number of windows

results.array = array(data=NA, dim=c(window.moves,9,block.nr),
                      dimnames=list(NULL,NULL,b.names))
corvec.array = array(data=NA, dim=c(window.moves,90,block.nr))

for (m in seq(1,nr.run)){
  print(paste("run nr ", m))
  ############################## Draw Dates From Sample Date Ranges ##############################
  
  datefrom <- as.numeric(data[,3])
  dateto <- as.numeric(data[,2])
  data[,4]<-((datefrom+dateto)/2)
  
  ##############################  CALCULATE TEMPORAL DISTANCES ##############################
  
  TIME.DIST.MAT <- matrix(nrow= nindividuals,ncol= nindividuals)
  
  for(i in seq(1,nindividuals)){
    
    for(j in seq(1,nindividuals)){
      TIME.DIST.MAT[i,j] <- abs(data[i,4] - data[j,4])
    }
  }
  ############################## SLIDING WINDOW PARAMETERS ##############################
  
  start.windows.at = 9000 # Starting time of the Sliding Window analyses (Years Ago)
  gap.window = 3000 # Width of the sliding window (in Years)
  move.window.by = 100 # Stride (in Years)
  
  for(w in seq(1: window.moves)) {
    ############################## Define the Window ##############################
    
    time.from = start.windows.at
    time.to = start.windows.at - gap.window
    start.windows.at = time.from - move.window.by
    
    data.positions <- which((data[,4] <= time.from) & (data[,4] >= time.to))
    nsamples <- length(data.positions)
    
    #		data.positions <- sample(data.positions,1, replace = F, ) # <- Use this line for jackknifing
    
    
    ############################## Geographic Distances ##############################
    
    geo.dist.mat <- matrix(NA, nsamples, nsamples)
    geo.dist.mat <- GEO.DIST.MAT[data.positions,data.positions]
    
    ###################################### Temporal Distances ###################################
    
    time.dist.mat <- matrix(NA, nsamples, nsamples)
    time.dist.mat <- TIME.DIST.MAT[data.positions,data.positions]
    
    ############################## Genetic Distances ##############################
    
    for (b in 1:block.nr){
      print(paste("run nr ", b))
      gen.dist.mat <- matrix(NA, nsamples, nsamples)
      gen.dist.mat <- block.array[data.positions, data.positions, b]
      
      ############################## MANTEL TESTING GEOG/TIME EUCLUD DERRIV ##############################
      
      nsim=90 # number of angles where correlation is calculated
      mcor_resvec = rep(0,nsim)
      
      scalevec <- tan(seq(0, ((pi/2) * (1-1e-8)), length.out= nsim))
      ang <- seq(0, (90 * (1-1e-8)), length.out= nsim)
      
      for(i in seq(1,nsim)){
        geotimedist <- (((geo.dist.mat )^2) + ((time.dist.mat * scalevec[i])^2))^0.5
        mcor_resvec[i] <- mantel(gen.dist.mat, geotimedist, permutations=0, method="pearson",
                                 parallel=getOption("mc.cores"))$statistic
      }
      mcor_resvec[is.na(mcor_resvec)] <- 0
      maxratio = ang[mcor_resvec==max(mcor_resvec)]
      maxcor = max(mcor_resvec)
      mincor = min(mcor_resvec)
      ends = c(mcor_resvec[1],  mcor_resvec[nsim])
      peak = maxcor - max(ends)
      if(length(maxratio) > 1){
        maxratio = 0
      }
      
      
      ############################## RECORD RESULTS  ##############################
      results.array[w,1,b] <- w # window number
      results.array[w,2,b] <- time.from # Starting time of the window (Years Ago)
      results.array[w,3,b] <- time.to # End time of the window (Years Ago)
      results.array[w,4,b] <- nsamples # Number of samples in a window
      results.array[w,5,b] <- maxratio # Angle value at the point of highest correlation
      results.array[w,6,b] <- maxcor # Highest correlation
      results.array[w,7,b] <- peak # Difference between the highest correlation and the correlation between the highest of the ends (correlation eith either geographic or temporal distance matrix)
      results.array[w,8,b] <- mcor_resvec[1] # Correlation with the geographic distance matrix
      results.array[w,9,b] <- mcor_resvec[nsim] # Correlation with temporal distance matrix
      
      corvec.array[w,,b] <-mcor_resvec # Vector of correlation values given each angle (n = nsim)
      
    }
  }
}

tend=date()
tstart
tend


saveRDS(results.array, file="ind_m45_g50_18_5mb_res.RData")

if (length(diff.test)>0){
  saveRDS(diff.test, file="ind_m45_g50_18_5mb_diff.RData")
}



