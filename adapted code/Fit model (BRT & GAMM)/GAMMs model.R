library(mgcv) #GAMMs
library(lubridate) #datetimes


# 1. ------------------- > Load and prep data ####
# Tracking data provided is a subset of 10 randomly selected individuals
# Response variable name is 'presabs' for Presence-Absence
tracks <- read.csv("2_Data/blwh_subset_dataset_ps_env.csv") 
head(tracks)
tracks$track_hame <- as.factor(tracks$track_hame)
tracks$month <- lubridate::month(tracks$dt)
tracks<-na.omit(tracks)

# 2. ------------------- > Fit candidate GAMMs ####
# Individuals (variable name 'ptt') are nested as a random effect

# sst,z,z_sd,ssh_sd,ild,eke (Model 1)
gam.mod1 <- mgcv::gam(PresAbs~s(sst, bs="ts")+s(z, bs="ts")+s(z_sd, bs="ts")+s(ssh_sd, bs="ts")+s(ild, bs="ts")+ s(EKE, bs="ts")+s(ptt, bs = "re"),family=binomial, data=tracks, method = "REML", select = T)

# sst,z,z_sd,ssh_sd,ild,eke,te(lon,lat) (Model 2)
gam.mod2 <- mgcv::gam(PresAbs ~ s(sst, bs="ts")+s(z, bs="ts")+s(z_sd, bs="ts")+s(ssh_sd, bs="ts")+s(ild, bs="ts")+ s(EKE, bs="ts")+te(lon,lat,bs="ts")+s(ptt, bs = "re"),family=binomial, data=tracks, method = "REML", select = T)

# z,z_sd,ssh_sd,ild,eke,te(sst,lat) (Model 3)
gam.mod3 <- mgcv::gam(PresAbs ~ s(z, bs="ts")+s(z_sd, bs="ts")+s(ssh_sd, bs="ts")+s(ild, bs="ts")+s(EKE, bs="ts")+te(sst, lat,bs="ts")+s(ptt, bs = "re"),family=binomial, data=tracks, method = "REML", select = T)

# sst, z,z_sd
gam.mod4 <- mgcv::gam(PresAbs ~ s(sst, bs="ts")+s(z, bs="ts")+s(z_sd, bs="ts")+s(ptt, bs = "re"),family=binomial, data=tracks, method = "REML", select = T)

# curl,sst,ssh,ssh_sd, sst_sd,z,z_sd,ild,eke
gam.mod5 <- mgcv::gam(PresAbs ~ s(curl, bs="ts")+s(sst, bs="ts")+s(ssh, bs="ts")+s(ssh_sd, bs="ts")+s(sst_sd, bs="ts")+s(z, bs="ts")+s(z_sd, bs="ts")+ s(ild, bs="ts")+s(EKE, bs="ts")+s(ptt, bs = "re"),family=binomial, data=tracks, method = "REML", select = T)

# sst, ssh_sd, z, z_sd, ild, eke, slope
gam.mod6 <- mgcv::gam(PresAbs~s(sst, bs="ts")+s(ssh_sd, bs="ts")+s(z, bs="ts")+s(z_sd, bs="ts")+s(ild, bs="ts")+ s(EKE, bs="ts")+ s(slope, bs="ts")+s(ptt, bs = "re"),family=binomial, data=tracks, method = "REML", select = T)

# sst, ssh_sd, z, z_sd, ild, eke, aspect
gam.mod7 <- mgcv::gam(PresAbs~s(sst, bs="ts")+s(ssh_sd, bs="ts")+s(z, bs="ts")+s(z_sd, bs="ts")+s(ild, bs="ts")+ s(EKE, bs="ts")+s(aspect, bs="ts")+s(ptt, bs = "re"),family=binomial, data=tracks, method = "REML", select = T)

# sst, ssh_sd, z*slope, z_sd, ild, eke
gam.mod8 <- mgcv::gam(PresAbs~s(sst, bs="ts")+s(ssh_sd, bs="ts")+s(z, slope,bs="ts")+ s(z_sd, bs="ts")+s(ild, bs="ts")+ s(EKE, bs="ts")+s(ptt, bs = "re"),family=binomial, data=tracks, method = "REML", select = T)

# bv, ssh_sd, z*slope, z_sd, ild, eke
gam.mod9 <- mgcv::gam(PresAbs~s(BV, bs="ts")+s(ssh_sd, bs="ts")+s(z, slope,bs="ts")+ s(z_sd, bs="ts")+s(ild, bs="ts")+ s(EKE, bs="ts")+s(ptt, bs = "re"),family=binomial, data=tracks, method = "REML", select = T)

# sst, ssh_sd, z, z_sd, ild, eke, month
gam.mod10<-mgcv::gam(PresAbs~s(sst,bs="ts")+s(ssh_sd,bs="ts")+s(z,bs="ts")+s(z_sd,bs="ts")+s(ild,bs="ts")+ s(EKE, bs="ts")+s(month, bs="cc") + s(ptt, bs = "re"),family=binomial, data=tracks, method = "REML", select = T)


#by season
library(dplyr)
library(lubridate)

# Assuming your date column is named 'date' and is already of Date or POSIXct class
your_dataframe <- tracks %>%
  mutate(month = month(date),
         season = case_when(
           month %in% c(12, 1, 2, 3, 4, 5, 6) ~ "winter_spring",
           month %in% c(7, 8, 9, 10, 11) ~ "summer_fall"
         ))

# Split into two datasets
winter_spring_df <- your_dataframe %>% filter(season == "winter_spring")
summer_fall_df   <- your_dataframe %>% filter(season == "summer_fall")


#winter-spring

summer_fall_df <- dismo::gbm.step(data=summer_fall_df, 
                                  gbm.x=c("curl","ild", "ssh", "sst","sst_sd", "ssh_sd", "z", "z_sd", "EKE","bv","slope", "aspect"), 
                                  gbm.y="PresAbs", ### response variable
                                  family = "bernoulli",
                                  tree.complexity = 3, ### complexity of the interactions that the model will fit
                                  learning.rate = 0.05,  ### optimized to end up with >1000 trees
                                  bag.fraction = 0.6) ### recommended by Elith, amount of input data used each time


dev_eval(summer_fall_df)

summary(summer_fall_df) #This function give us the % importance of each env. variable to explain the distribution of our species
summary_res <- summary(summer_fall_df)
setwd("5_Fit_model")
write.csv(summary_res, file = "summary_results_summer_fall_df.csv", row.names = FALSE)

dev_eval(summer_fall_df)# deviance explained by the model
summary_dev <- dev_eval(summer_fall_df)
write.csv(summary_dev, file = "summary_dev_summer_fall_df.csv", row.names = FALSE)

windows(20,20)
gbm.plot(summer_fall_df, smooth=TRUE, plot.layout = c(4,4), write.title=T) #response curves/ fitted functions
dev.print(png, file = "partial_curves_summer_fall_df.png", width = 12, height = 8, units = "in", res = 300,bg = "white")
dev.off()


saveRDS(summer_fall_df,"BLWH.res1.tc3.lr05.summer_fall_df.rds")



#winter_spring_df

winter_spring_df <- dismo::gbm.step(data=winter_spring_df, 
                                    gbm.x=c("curl","ild", "ssh", "sst","sst_sd", "ssh_sd", "z", "z_sd", "EKE","bv","slope", "aspect"), 
                                    gbm.y="PresAbs", ### response variable
                                    family = "bernoulli",
                                    tree.complexity = 3, ### complexity of the interactions that the model will fit
                                    learning.rate = 0.05,  ### optimized to end up with >1000 trees
                                    bag.fraction = 0.6) ### recommended by Elith, amount of input data used each time


summary(winter_spring_df)

dev_eval=function(model_object){
  null <- model_object$self.statistics$mean.null
  res <- model_object$self.statistics$mean.resid
  dev=((null - res)/null)*100 
  return(dev)
}

dev_eval(winter_spring_df)

summary(winter_spring_df) #This function give us the % importance of each env. variable to explain the distribution of our species
summary_res <- summary(winter_spring_df)
setwd("5_Fit_model")
write.csv(summary_res, file = "summary_results_winter_spring_df.csv", row.names = FALSE)

dev_eval(winter_spring_df)# deviance explained by the model
summary_dev <- dev_eval(winter_spring_df)
write.csv(summary_dev, file = "summary_dev_winter_spring_df.csv", row.names = FALSE)

windows(20,20)
gbm.plot(winter_spring_df, smooth=TRUE, plot.layout = c(4,4), write.title=T) #response curves/ fitted functions
dev.print(png, file = "partial_curves_winter_spring_df.png", width = 12, height = 8, units = "in", res = 300,bg = "white")
dev.off()


saveRDS(winter_spring_df,"BLWH.res1.tc3.lr05.winter_spring_df.rds")




#evaluation
#Leave One year Out analysis
LOO_eval <- function(DataInput, gbm.x, gbm.y, lr=lr, tc){
  DataInput <- tracks
  DataInput$date = as.POSIXct(DataInput$date, format = '%Y-%m-%d')
  DataInput$Year <- format(DataInput$date, "%Y")
  Evaluations_LOO <- as.data.frame(matrix(data=0,nrow=1,ncol=4))
  colnames(Evaluations_LOO) <- c("k","Deviance","AUC","TSS")
  counter=1
  for (y in min(DataInput$Year):max(DataInput$Year)){
    print(y)
    DataInput_train <- DataInput[DataInput$Year!=y,]
    DataInput_test <- DataInput[DataInput$Year==y,]
    DataInput.loo <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = c("PresAbs"), 
                              family="bernoulli", tree.complexity=tc,
                              learning.rate = lr, bag.fraction = 0.6)
    preds <- predict.gbm(DataInput.loo, DataInput_test,
                         n.trees=DataInput.loo$gbm.call$best.trees, type="response")
    dev <- calc.deviance(obs=DataInput_test$PresAbs, pred=preds, calc.mean=TRUE)
    d <- cbind(DataInput_test$PresAbs, preds)
    pres <- d[d[,1]==1,2]
    abs <- d[d[,1]==0,2]
    if(length(pres)>0 & length(abs)>0){
      e <- evaluate(p=pres, a=abs)
      
      Evaluations_LOO[counter,1] <- y
      Evaluations_LOO[counter,2] <- dev
      Evaluations_LOO[counter,3] <- e@auc
      Evaluations_LOO[counter,4] <- max(e@TPR + e@TNR-1)
      counter=counter+1 
    }
  }
  return(Evaluations_LOO)}


#######
#Make function to 75/25 split AUC test. 
#This is to 'repeat' what Elliot thinks they did for EcoCast (Kylie disagrees, see next function)
eval_7525 <- function(DataInput, gbm.x, gbm.y, lr, tc=tc, family){
  DataInput <- tracks
  Evaluations_7525 <- as.data.frame(matrix(data=0,nrow=1,ncol=5))
  colnames(Evaluations_7525) <- c("Deviance","AUC","TSS","Sensitivity", "Specificity")
  DataInput_bound <- floor((nrow(DataInput)/4)*3)         #define % of training and test set
  DataInput_train<- DataInput[sample(nrow(DataInput),DataInput_bound),]
  DataInput_test<- sqldf('SELECT * FROM DataInput EXCEPT SELECT * FROM DataInput_train')
  dim(DataInput_test)
  dim(DataInput_train)
  DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = "PresAbs", 
                               family="bernoulli", tree.complexity=tc,
                               learning.rate = lr, bag.fraction = 0.6)
  preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                       n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
  dev <- calc.deviance(obs=DataInput_test$PresAbs, pred=preds, calc.mean=TRUE)
  d <- cbind(DataInput_test$PresAbs, preds)
  pres <- d[d[,1]==1,2]
  abs <- d[d[,1]==0,2]
  e <- evaluate(p=pres, a=abs)
  Evaluations_7525[1,1] <- dev
  Evaluations_7525[1,2] <- e@auc
  Evaluations_7525[1,3] <- max(e@TPR + e@TNR-1)
  Evaluations_7525[1,4] <- mean(e@TPR)
  Evaluations_7525[1,5] <- mean(e@TNR)
  return(Evaluations_7525)}


#Make function to 100/100  AUC test
eval_100_percent <- function(dataInput, gbm.x, gbm.y, lr=lr, tc){
  DataInput <- tracks
  Evaluations_100_percent <- as.data.frame(matrix(data=0,nrow=1,ncol=3))
  colnames(Evaluations_100_percent) <- c("Deviance","AUC","TSS")
  DataInput_train<- DataInput
  DataInput_test<- DataInput
  DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                               family="bernoulli", tree.complexity=tc,
                               learning.rate = lr, bag.fraction = 0.6)
  preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                       n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
  dev <- calc.deviance(obs=DataInput_test$PresAbs, pred=preds, calc.mean=TRUE)
  d <- cbind(DataInput_test$PresAbs, preds)
  pres <- d[d[,1]==1,2]
  abs <- d[d[,1]==0,2]
  e <- evaluate(p=pres, a=abs)
  Evaluations_100_percent[1,1] <- dev
  Evaluations_100_percent[1,2] <- e@auc
  Evaluations_100_percent[1,3] <- max(e@TPR + e@TNR-1)
  
  return(Evaluations_100_percent)}

library(sqldf)
BLWH.loo.eval_full <-LOO_eval(tracks,gbm.x=gbm.x, gbm.y="PresAbs",lr=0.05, tc=3)
BLWH.7525.eval_full <- eval_7525(tracks,gbm.x=gbm.x, gbm.y="PresAbs",lr=0.05, tc=3)
BLWH.100.eval_full <- eval_100_percent(tracks,gbm.x=gbm.x, gbm.y="PresAbs",lr=0.05, tc=3)


saveRDS(BLWH.loo.eval_full,paste("BLWH.loo.eval_full.rds"))
saveRDS(BLWH.7525.eval_full,paste("BLWH.7525.eval_full.rds"))
saveRDS(BLWH.100.eval_full,paste("BLWH.100.eval_full.rds"))


BLWH.loo.eval_fullBLWH.loo.eval_fullSDR.7525.eval
write.csv(BLWH.loo.eval_full, file = "BLWH.loo.eval_full.csv", row.names = FALSE)

BLWH.7525.eval_full <- BLWH.7525.eval_full
write.csv(BLWH.7525.eval_full, file = "BLWH.7525.eval_full.csv", row.names = FALSE)


BLWH.100.eval_full <- BLWH.100.eval_full
write.csv(BLWH.100.eval_full, file = "BLWH.100.eval_full.csv", row.names = FALSE)


