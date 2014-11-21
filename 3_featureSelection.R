##############################
# Submission FlowCAP IV      #
# Step 3 - Feature selection #
##############################

FCSloc <- "/group/irc/shared/data/FlowCAPIV"

#############################################
# Load the required libraries and functions #
#############################################

library(survival)

scoreFeature <- function(survivalTime,status,feature,na.action="impute"){
  if(sum(is.na(feature))==length(feature)){
    return(1)
  }
  if(na.action == "impute"){
    feature[is.na(feature)] <- mean(feature,na.rm=TRUE)  
  }
  data <- list(time=survivalTime,status=status,x=feature)
  fit <- coxph(Surv(time,status)~x, data)
  pvalue <- summary(fit)$logtest[3]
  return(pvalue)
}

######################
# Read the meta data #
######################

meta_patients <- read.csv(paste(FCSloc,"MetaDataTrain.csv",sep="/"),as.is=TRUE)
nPatients<-nrow(meta_patients)

train_ids <- 1:191
train_ids <- 1:2 ############## TODO
train_survivalTime <- meta_patients[train_ids,"Survival.Time"]
train_survivalTime <- train_survivalTime - min(train_survivalTime) + 1 # To ensure all survival times are positive
train_status <- meta_patients[train_ids,"Status"]

########################################
# Calculate the score for the features #
########################################

bestScores <- rep(1,500)
bestFeatures <- matrix(0,nrow=500,ncol=383)
nFeatures <- 3^10 * 14 

# For all the separate feature files
for(i in seq(1,(3*nFeatures),by=10000)){
  # Read the features
  allfeatures <- read.csv(paste(FCSloc,"/dambi/features_allpatients_",i,"-",i+9999,".csv",sep=""),header=T,stringsAsFactors=FALSE)
  rownames(allfeatures)=allfeatures[,1]
  allfeatures=as.matrix(allfeatures[,-1])
  
  # Initialize variables
  scores <- rep(1,nrow(allfeatures))
  names(scores) <- rownames(allfeatures)
  
  # Score all these features
  for (f in 1:nrow(allfeatures)) {
    scores[f] <- scoreFeature(train_survivalTime,train_status,allfeatures[f,train_ids])
    #print(paste(f,scores[f]))
  }
  
  scores <- sort(scores)
  
  # Save all the scores
  write.csv(scores,file=paste(FCSloc,"/dambi/featureScores_",i,"-",i+9999,".csv",sep=""))
  
  # Only keep those features without NA values for the event-instances
  nrNAs <- rowSums(is.na(allfeatures[names(scores),which(train_status==1)]))
  scores <- scores[nrNAs == 0]
  
  # Select the best 1000 features
  if(i == 1){
    bestScores <- scores[1:1000]
    bestFeatures <- t(allfeatures[names(bestScores),])
  } else {
    bestScores <- sort(c(bestScores,scores[1:1000]))[1:1000]
    newFeatures <- bestScores[names(bestScores) %in% rownames(allfeatures)]
    oldFeatures <- bestScores[names(bestScores) %in% colnames(bestFeatures)]
    bestFeatures <- cbind(bestFeatures[,names(oldFeatures)],t(allfeatures[names(newFeatures),]))
    bestFeatures <- bestFeatures[,names(bestScores)]
  }
  
  
}

write.csv(bestFeatures,file=paste(FCSloc,"dambi","best1000Features.csv",sep="/"))
write.csv(bestScores,file=paste(FCSloc,"dambi","best1000Scores.csv",sep="/"))

# Select uncorrelated features
selected_features <- c(1,2)
for(i in 3:1000){
  print(i)
  max <- 0
  for(f in selected_features){
    if(sd(bestFeatures[,f]) > 0){
      tmp <- abs(cor(bestFeatures[train_ids,i],bestFeatures[train_ids,f]))
      if(tmp > max) max <- tmp
    } else {
      max <- 1
    }
  }
  if(max < 0.2){
    selected_features <- c(selected_features,i)
    print(selected_features)
  } 
}

# Save this selection of features
write.csv(bestFeatures[,selected_features],file=paste(FCSloc,"dambi","bestUncorrelatedFeatures.csv",sep="/"))