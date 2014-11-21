##########################
# Submission FlowCAP IV  #
# Step 1 - Preprocessing #
##########################

FCSloc <- "/group/irc/shared/data/FlowCAPIV"
dir.create(paste(FCSloc,"dambi",sep="/"),showWarnings=FALSE)

#############################################
# Load the required libraries and functions #
#############################################

library(flowCore)
library(flowDensity)

# Quality control function
medianTime_select <- function(flowFrame, parameter, nIntervals, maxDiff){
  # Split the data in equal size intervals
  cuts <- cut(exprs(flowFrame)[,"Time"],breaks=nIntervals)
  splitted <- split(as.data.frame(exprs(flowFrame)),cuts)
  
  # Count the number of cells in each interval
  counts <- sapply(splitted,function(flowDataFrame){dim(flowDataFrame)[1]})
  enoughCells <- counts > median(counts)-(2*sd(counts))
  
  # Calculate the median values of the cells in each interval
  medians <- sapply(splitted,function(flowDataFrame){median(flowDataFrame[,parameter])})
  stableMedian <- abs(diff(medians))<maxDiff
  stableMedian <- c(T,stableMedian) & c(stableMedian,T)
  
  # Return only the cells which meet the conditions
  selection <- which(enoughCells & stableMedian)
  return(flowFrame[which(as.numeric(cuts) %in% selection),])
}

# Function to remove events on the margins
removeMargins <- function(flowFrame,dimensions){
  # Look op the accepted ranges for the dimensions
  meta <- pData(flowFrame@parameters)
  rownames(meta) <- meta[,"name"]
  
  # Initialize variables
  selection <- rep(TRUE,times=dim(flowFrame)[1])
  e <- exprs(flowFrame)
  
  # Make selection
  for(d in dimensions){
    selection <- selection & 
      e[,d] > max(meta[d,"minRange"],min(e[,d])) &
      e[,d] < min(meta[d,"maxRange"],max(e[,d]))
  }
  return(flowFrame[selection,])
}

# Function to select single cells
removeDoublets <- function(flowFrame, d1="FSC-A", d2="FSC-H", w=NULL,silent=TRUE){
  # Calculate the ratios
  ratio <- exprs(flowFrame)[,d1] / (1+ exprs(flowFrame)[,d2])
  
  # Define the region that is accepted
  r <- median(ratio)
  if(is.null(w)){ w <- 2*sd(ratio) }
  if(!silent){
    print(r)
    print(w)
  }
  
  # Make selection
  selection <- which(ratio < r+w)
  return(flowFrame[selection,])
}

######################
# Read the meta data #
######################

meta_patients <- read.csv(paste(FCSloc,"MetaDataTrain.csv",sep="/"),as.is=TRUE)

filenames <- c(rbind(meta_patients[,"Stim"],meta_patients[,"Unstim"]))

#############################
# Initialize some variables #
#############################

colnames_meta <- c("File","Original","Bad quality","On margins","doublets","Alive T cells",
                   "% Good quality","% Not on margins","% Single cells","% Alive T cells")
meta <- matrix(0,nrow=length(filenames),ncol=10,dimnames=list(filenames,colnames_meta))

#######################
# Preprocess the data #
#######################

for(file in filenames){ 
  flowFrame <- read.FCS(paste(FCSloc,file,sep="/"))
  
  # a. Quality Control
  flowFrame_q <- medianTime_select(flowFrame,"FSC-A",100,10000)
  
  # b. Remove margin events
  flowFrame_m <- removeMargins(flowFrame_q,dimensions=colnames(flowFrame))
  
  # c. Select single cells
  flowFrame_d <- removeDoublets(flowFrame_m,d1="FSC-A",d2="FSC-H")
  
  # d. Compensate
  flowFrame_c <- compensate(flowFrame_d,flowFrame@description$SPILL)
  
  # e. Transform
  flowFrame_t <- transform(flowFrame_c,transformList(c("SSC-A",colnames(flowFrame@description$SPILL)), logicleTransform()))
  
  # f. Select alive T-cells
  flowFrame_a <- tryCatch({ # lots of try catch situations to make sure there are cells selected
    print("flowFrame_tmp")
    flowFrame_tmp <- flowFrame_t[(flowDensity(flowFrame_t,channels=c("V450-A","R780-A"),position=c(FALSE,TRUE)))@index,]
    if(dim(flowFrame_tmp)[1]<100) stop("No cells selected")
    flowFrame_tmp
  }, error = function(e){
    print("Failed on combined directions")
    return(NULL)
  })
    
  if(is.null(flowFrame_a)){
    flowFrame_tmp <- tryCatch({
      print("flowFrame_tmp2")
      flowFrame_tmp2 <- flowFrame_t[(flowDensity(flowFrame_t,channels=c("V450-A","R780-A"),position=c(FALSE,NA)))@index,]
      if(dim(flowFrame_tmp2)[1]<1) stop("No cells selected")  
      flowFrame_tmp2
    }, error = function(e){
      print("Failed on V450-A")
      flowFrame_tmp2 <- flowFrame_t[(flowDensity(flowFrame_t,channels=c("V450-A","R780-A"),position=c(FALSE,NA),upper=c(FALSE,NA)))@index,]  
      return(flowFrame_tmp2)
    })
    
    flowFrame_a <- tryCatch({
      print("flowFrame_tmp3")
      flowFrame_tmp3 <- flowFrame_tmp[(flowDensity(flowFrame_tmp,channels=c("V450-A","R780-A"),position=c(NA,TRUE)))@index,]
      if(dim(flowFrame_tmp3)[1]<1) stop("No cells selected")  
      flowFrame_tmp3
    }, error = function(e){
      print("Failed on R780-A")
      flowFrame_tmp3 <- flowFrame_tmp[(flowDensity(flowFrame_tmp,channels=c("V450-A","R780-A"),position=c(NA,TRUE),percentile=c(NA,0.5)))@index,]  
      return(flowFrame_tmp3)
    })  
  }
  
  
  # Save the result
  write.FCS(flowFrame_a,paste(FCSloc,"/dambi/",gsub(".fcs","",file),"_preprocessed.fcs",sep=""))
  # Save some metadata
  meta[file,] <- c(file,dim(flowFrame)[1],
                   dim(flowFrame)[1]-dim(flowFrame_q)[1],dim(flowFrame_q)[1]-dim(flowFrame_m)[1],
                   dim(flowFrame_m)[1]-dim(flowFrame_d)[1],dim(flowFrame_a)[1],
                   dim(flowFrame_q)[1]/dim(flowFrame)[1],dim(flowFrame_m)[1]/dim(flowFrame_q)[1],
                   dim(flowFrame_d)[1]/dim(flowFrame_m)[1],dim(flowFrame_a)[1]/dim(flowFrame_t)[1])
  print(meta[file,])
}
write.csv(meta,file=paste(FCSloc,"/dambi/MetaData_preprocessing.csv", sep=""))
