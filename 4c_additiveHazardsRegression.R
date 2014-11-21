#########################################
# Submission FlowCAP IV                 #
# Step 4c - Additive Hazards regression #
#########################################

FCSloc <- "/group/irc/shared/data/FlowCAPIV"

#############################################
# Load the required libraries and functions #
#############################################

library(survival)
library(ahaz)

fitModel <- function(survivalTime,status,features){
  data <- list(time=survivalTime,status=status,x=features)
  fit <- coxph(Surv(time,status)~x, data)
  return(fit)
}

fitAdditiveModel <- function(survivalTime,status,features){
  fit <- tune.ahazpen(Surv(survivalTime, status), features)
  return(fit)
}

######################
# Read the meta data #
######################

meta_patients <- read.csv(paste(FCSloc,"MetaDataTrain.csv",sep="/"),as.is=TRUE)
nPatients<-nrow(meta_patients)

train_ids <- 1:191
train_survivalTime <- meta_patients[train_ids,"Survival.Time"]
train_survivalTime <- train_survivalTime - min(train_survivalTime) + 1 # To ensure all survival times are positive
train_status <- meta_patients[train_ids,"Status"]

##############################
# Read the selected features #
##############################

features <- read.csv(paste(FCSloc,"dambi","best1000Features.csv",sep="/"),as.is=TRUE)
rownames(features) <- features[,1]
features <- as.matrix(features[,-1])

selected_features <- 1:100
######################################
# Try leave one out cross-validation #
######################################

# Make a prediction for each sample in the training set
prediction <- NULL
for (i in train_ids) {
  fold <- setdiff(train_ids,i)
  
  fit <- fitAdditiveModel(train_survivalTime[fold],train_status[fold],features[fold,selected_features])
  
  toPredict <- list(time=train_survivalTime[i],status=train_status[i],x=features[i,selected_features,drop=F])
  prediction <- c(prediction,predict(fit,newX=rbind(toPredict$x),lambda=fit$lambda.min,type="lp"))
  
  print(paste(i,train_survivalTime[i],prediction[i]))
}

# Evaluate the result
fit <- fitModel(train_survivalTime,train_status,prediction)
pvalue <- summary(fit)$logtest[3]
cindex <- summary(fit)$concordance
correlation <- cor(prediction[train_status==1], train_survivalTime[which(train_status==1)])

# Print evaluation measures
print(paste("p value: ",pvalue," c-index: ",cindex[1]," correlation: ",correlation))

#################################
# Generate results for test set #
#################################

fit <- fitModel(train_survivalTime,train_status,features[train_ids,selected_features])
toPredict <- list(time=rep(NA,383),status=rep(NA,383),x=features[,selected_features])
prediction <- predict(fit,toPredict)
meta_patients[,"Prediction"] <- prediction
write.csv(meta_patients,file=paste(FCSloc,"dambi","predictionAdditiveHazards.csv",sep="/"))