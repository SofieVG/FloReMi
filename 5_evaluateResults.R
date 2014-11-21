coxph <- read.csv("../Dropbox/Phd_Sofie/FlowCAP_IV/predictions/predictionCoxPH.csv")
addh <- read.csv("../Dropbox/Phd_Sofie/FlowCAP_IV/predictions/predictionAdditiveHazards.csv")
rf <- read.csv("../Dropbox/Phd_Sofie/FlowCAP_IV/predictions/predictionRandomSurvivalForest.csv")
rf1000 <- read.csv("../Dropbox/Phd_Sofie/FlowCAP_IV/predictions/predictionRandomSurvivalForest1000.csv")
rf1000_m <- read.csv("../Dropbox/Phd_Sofie/FlowCAP_IV/predictions/predictionRandomSurvivalForest1000_moreFeatures.csv")

glmnet <- read.csv("../Dropbox/Phd_Sofie/FlowCAP_IV/predictions/predictionGlmnet.csv")

gold <- read.csv("../Dropbox/Phd_Sofie/FlowCAP_IV/predictions/MetaDataFull.csv")


library(survival)

stime <- gold[,"Survival.Time"]
sstatus <- gold[,"Status"]

features <- list(coxph[,"Prediction"], addh[,"Prediction"], rf[,"Prediction"])
features <- list(rf1000[,"Prediction"],rf1000_m[,"Prediction"])
features <- list(glmnet[,"Prediction"])
features_transformed = lapply(features, function(f){ 1 - ((f - (min(f))) / max(f))})

for(feature in features_transformed){
  
  data_all <- list(time=stime, status=sstatus, x=feature)
  data_train <- list(time=stime[1:191], status=sstatus[1:191], x=feature[1:191])
  data_test <- list(time=stime[192:383], status=sstatus[192:383], x=feature[192:383])
  for(data in list(data_all,data_train,data_test)){
    fit<-coxph(Surv(time, status) ~ x, data)
    pvalue=summary(fit)$logtest[3]
    print(paste(length(data[[1]]),pvalue, summary(fit)$concordance[1]))
  }
}


coxph[,"Prediction"] <- features_transformed[[1]]
addh[,"Prediction"] <- features_transformed[[2]]
rf[,"Prediction"] <- features_transformed[[3]]


write.csv(coxph, file="../Dropbox/Phd_Sofie/FlowCAP_IV/predictions/predictionCoxPH_t.csv")
write.csv(addh, file="../Dropbox/Phd_Sofie/FlowCAP_IV/predictions/predictionAdditiveHazards_t.csv")
write.csv(rf, file="../Dropbox/Phd_Sofie/FlowCAP_IV/predictions/predictionRandomSurvivalForest_t.csv")
