###############################
# Submission FlowCAP IV       #
# Step 2 - Feature extraction #
###############################

FCSloc <- "/group/irc/shared/data/FlowCAPIV"

#############################################
# Load the required libraries and functions #
#############################################

library(flowCore)
library(flowDensity)
library(flowType)

# Function to extract features from one file
getFeaturesFromFile <- function(file){
  
  # Intitialize variables
  markers <- c("FSC-A","FSC-H","SSC-A","B710-A","IFNg","TNFa","CD4","CD27","CD107a","CD154","CD3","CCR7","IL2","CD8","V705-A","V655-A","V605-A","CD57","V565-A","CD45RO","Vivid/CD14")
  prop <- c(1,3,7,8,9,10,12,14,18,20)
  mfi <- c(1,3,5:10,12:14,18,20)

  # Read the preprocessed file
  print(file)
  file_preprocessed <- paste(gsub(".fcs","",file),"_preprocessed.fcs",sep="")
  flowFrame <- read.FCS(paste(FCSloc,"dambi",file_preprocessed,sep="/"))
  
  # Calculate the splits with the flowDensity package
  tresholds <- NULL
  for(d in colnames(flowFrame)[prop]){
    tresholds[d] <- deGate(flowFrame,d)
  }
  
  # Calculate the features with the flowType package
  res <- flowType(flowFrame, PropMarkers = prop, MFIMarkers = mfi, Methods = 'thresholds', Thresholds= as.list(tresholds), MarkerNames=colnames(flowFrame))
  MFIs <- res@MFIs
  proportions <- res@CellFreqs / max(res@CellFreqs)
  
  # Define the feature names
  resultNames <- res@PhenoCodes
  for(i in mfi){
    resultNames <- c(resultNames,paste(res@PhenoCodes,markers[i]))
  }
  
  # Put everything in a vector
  result <- rep(0,length=(1+length(mfi))*length(proportions))
  names(result) <- resultNames
  result[1:length(proportions)] <- proportions
  for(i in 1:length(mfi)){
    result[(i*length(proportions)+1):((i+1)*length(proportions))] <- MFIs[,i]
  }
  return(result)
}

######################
# Read the meta data #
######################

meta_patients <- read.csv(paste(FCSloc,"MetaDataTrain.csv",sep="/"),as.is=TRUE)
nPatients<-nrow(meta_patients)

filenames <- sort(c(meta_patients[,"Stim"],meta_patients[,"Unstim"]))

########################
# Extract the features #
########################

# Initialize variables
nFeatures <- 3^10 * 14 
allfeatures <- matrix(ncol=nPatients,nrow=2+(3*nFeatures)) # rows: survival time + censor + (stim + unstim + diff)*nFeatures
tmp <- getFeaturesFromFile("001.fcs")
featureNames <- names(tmp)
rownames(allfeatures) <- c("Survival.Time","Status",paste("Stim",featureNames), paste("Unstim",featureNames), paste("diff",featureNames))

# Extract features for each patient
for (i in 1:nPatients) {
  print(i)
  file_stim = meta_patients[i,"Stim"]
  features_stim <- getFeaturesFromFile(file_stim)
  file_unstim = meta_patients[i,"Unstim"]
  features_unstim <- getFeaturesFromFile(file_unstim)
  allfeatures[,i]<-unlist(c(meta_patients[i,"Survival.Time"],
                            meta_patients[i,"Status"],
                            features_stim,
                            features_unstim,
                            features_stim - features_unstim))
}

# Save features in separate files
for(i in seq(1,(3*nFeatures),by=10000)){
  lastindex<-min(i+10001,3*nFeatures)
  write.csv(allfeatures[(i+2):lastindex,],file=paste(FCSloc,"/dambi/features_allpatients_",i,"-",i+9999,".csv",sep=""))
}
