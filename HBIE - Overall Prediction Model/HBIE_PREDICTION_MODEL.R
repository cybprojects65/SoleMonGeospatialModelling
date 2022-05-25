library(sqldf)

#structure of the file, a CSV with the following columns:
#station	lat	lon	biomIndex2019	diva	maxent	ssa	reconstructed	divaerror	ssauncert
#NOTE: change the maxent omission rate threshold

fileIdxes<-c()
maxEntOR<-c()

fileIdx = "./All Inputs/Sepia_officinalis_index_reconstructed_diva_maxent_ssa2019.csv"
maxEnt_1Percent_Omission_Rate<-0.021
fileIdxes<-c(fileIdxes,fileIdx)
maxEntOR<-c(maxEntOR,maxEnt_1Percent_Omission_Rate)

fileIdx = "./All Inputs/Sepia_offcinalis_index_reconstructed_diva_maxent_ssa_fixed2020.csv"
maxEnt_1Percent_Omission_Rate<-0.022
fileIdxes<-c(fileIdxes,fileIdx)
maxEntOR<-c(maxEntOR,maxEnt_1Percent_Omission_Rate)

fileIdx = "./All Inputs/Solea_solea_index_reconstructed_diva_maxent_ssa_optimal2019.csv"
maxEnt_1Percent_Omission_Rate<-0.036
fileIdxes<-c(fileIdxes,fileIdx)
maxEntOR<-c(maxEntOR,maxEnt_1Percent_Omission_Rate)

fileIdx = "./All Inputs/Solea_solea_index_reconstructed_diva_maxent_ssa_2020_conservative_2020ref.csv"
maxEnt_1Percent_Omission_Rate<-0.037
fileIdxes<-c(fileIdxes,fileIdx)
maxEntOR<-c(maxEntOR,maxEnt_1Percent_Omission_Rate)

fileIdx = "./All Inputs/Squilla_mantis__index_reconstructed_diva_maxent_ssa2019_fixed.csv"
maxEnt_1Percent_Omission_Rate<-0.031
fileIdxes<-c(fileIdxes,fileIdx)
maxEntOR<-c(maxEntOR,maxEnt_1Percent_Omission_Rate)

fileIdx = "./All Inputs/Squilla_mantis__index_reconstructed_diva_maxent_ssa2020_fixed.csv"
maxEnt_1Percent_Omission_Rate<-0.016
fileIdxes<-c(fileIdxes,fileIdx)
maxEntOR<-c(maxEntOR,maxEnt_1Percent_Omission_Rate)

fileIdx = "./All Inputs/Pecten_jacobaeus__index_reconstructed_diva_maxent_ssa2019.csv"
maxEnt_1Percent_Omission_Rate<-0.033
fileIdxes<-c(fileIdxes,fileIdx)
maxEntOR<-c(maxEntOR,maxEnt_1Percent_Omission_Rate)

fileIdx = "./All Inputs/Pecten_jacobaeus__index_reconstructed_diva_maxent_ssa2020.csv"
maxEnt_1Percent_Omission_Rate<-0.037
fileIdxes<-c(fileIdxes,fileIdx)
maxEntOR<-c(maxEntOR,maxEnt_1Percent_Omission_Rate)


for (i in 1:length(fileIdxes)){

  fileIdx<- fileIdxes[i]
  maxEnt_1Percent_Omission_Rate<- maxEntOR[i]
  
fileOut<-paste0(gsub(".csv", "_predictionmodel_val",fileIdx),".csv")
data<-read.csv(file = fileIdx)
reconstructed<-data[data$reconstructed=="Y",]
truevalues<-data[data$reconstructed=="N",]

lowthreshold<<-maxEnt_1Percent_Omission_Rate##1% omission rate
midthreshold<<-0.21	#10% omission rate
minweight<<-0.4#0.5
midweight<<-1
highweight<<-1

thresholdMaxEnt<-function(x){
  if (x<lowthreshold)
    return (minweight)
  else if (x<midthreshold)  
    return (midweight)
  else
    return (highweight)
}

Wpresence =  (lapply(reconstructed$maxent,thresholdMaxEnt))
WpresenceVector<-unlist(Wpresence)
reconstructed$wpresence<-WpresenceVector


alphas = c(0.6) #seq(0.1,1,by=0.1)
betas = c(1) #seq(0.1,1,by=0.1)

minerror = 100
bestalpha <-0
minrerrorrelative<-100
minaveragerelativeerror<-100

for (alpha in alphas){
  for (beta in betas){
    trueval<-reconstructed$biomIndex2019
    
    diva <- reconstructed$diva
    #pag 44 DIVA user guide When relative errors are demanded, one could divide by σ^2 Bii; https://github.com/gher-ulg/Diva-User-Guide/blob/master/DivaUserGuide.pdf
    divarelativeerror<-sqrt(reconstructed$divaerror)/(sd(truevalues$biomIndex2019))
    divaabsoluteerror<-divarelativeerror*reconstructed$diva
    sigmamultiplier=1.96
    #sigmamultiplier=sigmamultiplier*sigmamultiplier
    confid<-sqrt(reconstructed$divaerror)*sigmamultiplier*sd(truevalues$biomIndex2019)
    divaMin <- reconstructed$diva-confid
    divaMax <- reconstructed$diva+confid
    bonus=1.2#1.24
    malus=0.8#0.7
    divaerror <- reconstructed$divaerror
    for (d in 1:length(reconstructed$diva)){
      if (divaerror[d] == 0){
        diva[d]=0
        divaMin[d]=0
      }
      else if (divaerror[d] < 0.1){
        diva[d]=bonus*diva[d]
        divaMin[d]=bonus*divaMin[d]
        divaMax[d]=bonus*divaMax[d]
      }else if (divaerror[d] < 0.15){
        diva[d]=1*diva[d]
        divaMin[d]=1*divaMin[d]
        divaMax[d]=1*divaMax[d]
      }else{
        diva[d]=malus*diva[d]
        divaMin[d]=malus*divaMin[d]
        divaMax[d]=malus*divaMax[d]
      }
      
    }
    
    #diva<-alpha*reconstructed$diva
    ssa<-beta*reconstructed$ssa
    
    #reconstructed$averagebio<-(alpha*reconstructed$diva+beta*reconstructed$ssa)/2
    reconstructed$averagebio<-(diva+ssa)/2#2
    #THIS IS ERROR FIELD NOT CONFIDENCE INTERVAL
    averagebioMin<-(divaMin+ssa)/2
    averagebioMax<-(divaMax+ssa)/2
    
    predictionModel <- reconstructed$wpresence*reconstructed$averagebio
    reconstructed$prediction<-predictionModel
    reconstructed$predictionMin <- reconstructed$wpresence*averagebioMin
    reconstructed$predictionMin[which(reconstructed$predictionMin<0)]<-0
    reconstructed$predictionMax <- reconstructed$wpresence*averagebioMax
    
    reconstructed$discrepancy2019<-(reconstructed$biomIndex2019-reconstructed$prediction)
    reconstructed$discrepancy2019abs<-abs(reconstructed$biomIndex2019-reconstructed$prediction)
    reconstructed$discrepancy2019rel<-reconstructed$discrepancy2019abs*100/reconstructed$biomIndex2019
    reconstructed$discrepancy2019rel[which(is.nan(reconstructed$discrepancy2019rel))]<-0
    reconstructed$discrepancy2019rel[which(is.infinite(reconstructed$discrepancy2019rel))]<-NA
    averageRelativeError<-mean(reconstructed$discrepancy2019rel,na.rm = T)
    
    npoints<-dim(reconstructed)[1]
    
    error2019<-sqrt(sum((reconstructed$discrepancy2019)*(reconstructed$discrepancy2019))/npoints)
    
    merror2019<-mean(reconstructed$discrepancy2019abs)
    mean2019<-mean(reconstructed$biomIndex2019)
    merror2019relative <- merror2019/mean2019
    sigma<-sd(reconstructed$discrepancy2019)
    
    #cat("Alpha:",alpha,"\n")
    #cat("Mean Squared Error:",error2019,"\n")
    #cat("Mean absolute error:",merror2019,"\n")
    #cat("Mean relative error:",merror2019relative,"\n")
    #cat("Average relative error:",averageRelativeError,"\n")
    #cat("Sigma Error:",sigma,"\n")
    
    #complete<-rbind(reconstructed,truevalues)
    write.csv(file=fileOut,x = reconstructed,row.names = F)
    #if (error2019<minerror)
    if (averageRelativeError<minaveragerelativeerror)
    {
      minerror<-error2019
      bestalpha<-alpha
      bestbeta<-beta
      minrerrorrelative<-merror2019relative
      minaveragerelativeerror<-averageRelativeError
    }
  }
}

nestimationscaught<-0
nestimations<-length(reconstructed$biomIndex2019)

for (g in 1:nestimations){
  
  check<-(reconstructed$biomIndex2019[g]>=reconstructed$predictionMin[g] && 
            reconstructed$biomIndex2019[g]<=reconstructed$predictionMax[g])
  if (check){
    cat("OK:",reconstructed$biomIndex2019[g],"vs",reconstructed$prediction[g],"(",reconstructed$predictionMin[g],";",reconstructed$predictionMax[g],")",";",reconstructed$divaerror[g],"\n")
    nestimationscaught=nestimationscaught+1
  }
  else
    cat("KO:",reconstructed$biomIndex2019[g],"vs",reconstructed$prediction[g],"(",reconstructed$predictionMin[g],";",reconstructed$predictionMax[g],")",";",reconstructed$divaerror[g],"\n")
}


cat("Correct estimations:",nestimationscaught,"over",nestimations,"\n")
accuracy = nestimationscaught*100/nestimations

cat("ACCURACY:",accuracy,"%\n")

cat("\n")

#cat("Best Beta:",bestbeta,"\n")
#cat("Optimal MSE:",minerror,"\n")

#cat("Optimal mean Relative error/mean biomass :",minrerrorrelative*100,"%\n")
#cat("Optimal Mean pointwise relative error :",minaveragerelativeerror,"%\n")
#cat("Relative error range:",min(reconstructed$discrepancy2019rel),";",max(reconstructed$discrepancy2019rel),"\n")

discrepancy2019<-reconstructed$biomIndex2019-reconstructed$prediction
mse<-sqrt(sum((discrepancy2019)*(discrepancy2019))/npoints)
mse_disc<-mse/mean(reconstructed$biomIndex2019)

#cat("MSE:",mse,"\n")
#cat("MSE/MeanBio2019:",mse_disc,"\n")

averageRelativeErrorOftheModel<-(reconstructed$predictionMax-reconstructed$predictionMin)/2
averageRelativeErrorOftheModel<-averageRelativeErrorOftheModel/reconstructed$prediction
averageRelativeErrorOftheModel[is.infinite(averageRelativeErrorOftheModel)]<-0
averageRelativeErrorOftheModelm<-mean(averageRelativeErrorOftheModel)
#cat("AVG Relative Error of the estimates:",averageRelativeErrorOftheModelm,"\n")


d0<-reconstructed$discrepancy2019rel

de<-reconstructed$divaerror/max(reconstructed$divaerror)
di<-reconstructed$discrepancy2019abs/max(reconstructed$discrepancy2019abs)
ds<-reconstructed$ssauncert/max(reconstructed$ssauncert)

#plot (y=de,x=seq(1:length(de)),type='l')
#lines(y=di,x=seq(1:length(de)),col='red')
#lines(y=ds,x=seq(1:length(de)),col='green')


totalbiomasstrue<-sum(data$biomIndex2019)
totalbiomassestimated<-sum(c(truevalues$biomIndex2019,reconstructed$prediction))
totalMINbiomassestimated<-sum(c(truevalues$biomIndex2019,reconstructed$predictionMin))
totalMAXbiomassestimated<-sum(c(truevalues$biomIndex2019,reconstructed$predictionMax))
percbiomasstrue<-sum(reconstructed$biomIndex2019)*100/totalbiomasstrue

#cat("Total biomass in missing locations in 2019 :",sum(reconstructed$biomIndex2019),"vs",sum(truevalues$biomIndex2019)," over total ",totalbiomasstrue, percbiomasstrue,"%","\n")
#cat("Total biomass 2019 :",totalbiomasstrue," vs ",totalbiomassestimated,"(",totalMINbiomassestimated,";",totalMAXbiomassestimated,")","\n")


############################
# Import Data
#HaulBiomass <- read.csv("HaulBiomass.csv") # Biomass Index by Haul

calcBiomass<-function(HaulBiomass){
  StrataWeight <- read.csv("StrataWeight.csv") # Area of strata
  HaulData <- read.csv("HaulData.csv") # Swept Areas
  
  # Compute Index
  Index=merge(HaulBiomass, HaulData, by=c("year", "Station"))
  
  Index$BiomRaw=Index$BiomIndex*Index$SweptArea # return from Biomass Index to raw weight of species in the haul
  
  Index=Index[!is.na(Index$BiomRaw),]
  
  Index=aggregate(list(Biomass=Index$BiomRaw, SweptArea=Index$SweptArea), 
                  by = list(year=Index$year, Stratum=Index$Stratum), 
                  FUN=sum) # summarize weight and swept area by year and stratum
  
  Index$BiomassStratum=Index$Biomass/Index$SweptArea # Stratum Index 
  
  Index=merge(Index, StrataWeight, by="Stratum") 
  
  Index$BiomassStratumWeighted=Index$BiomassStratum*Index$StratumWeight # Weight Stratum Index by the relative area of the stratum
  
  Index=aggregate(list(Index=Index$BiomassStratumWeighted), by=list(year=Index$year), FUN=sum) # Get Index
  return (Index)
}

HaulBiomass<-data.frame(year="2019",Station=truevalues$station,BiomIndex=truevalues$biomIndex2019)
Index2019Known2020<-calcBiomass(HaulBiomass)
HaulBiomass<-data.frame(year="2019",Station=c(truevalues$station,reconstructed$station),BiomIndex=c(truevalues$biomIndex2019,reconstructed$biomIndex2019))
Index2019True<-calcBiomass(HaulBiomass)
HaulBiomass<-data.frame(year="2019",Station=c(reconstructed$station),BiomIndex=c(reconstructed$biomIndex2019))
Index2019Missing2020<-calcBiomass(HaulBiomass)

HaulBiomass<-data.frame(year="2019",Station=c(truevalues$station,reconstructed$station),BiomIndex=c(truevalues$biomIndex2019,reconstructed$prediction))
Index2019Predicted<-calcBiomass(HaulBiomass)
HaulBiomass<-data.frame(year="2019",Station=c(truevalues$station,reconstructed$station),BiomIndex=c(truevalues$biomIndex2019,reconstructed$predictionMin))
Index2019PredictedMin<-calcBiomass(HaulBiomass)
HaulBiomass<-data.frame(year="2019",Station=c(truevalues$station,reconstructed$station),BiomIndex=c(truevalues$biomIndex2019,reconstructed$predictionMax))
Index2019PredictedMax<-calcBiomass(HaulBiomass)

HaulBiomass<-data.frame(year="2019",Station=c(reconstructed$station),BiomIndex=c(reconstructed$prediction))
Index2019PredictedMissing2020<-calcBiomass(HaulBiomass)
HaulBiomass<-data.frame(year="2019",Station=c(reconstructed$station),BiomIndex=c(reconstructed$predictionMin))
Index2019PredictedMissing2020Min<-calcBiomass(HaulBiomass)
HaulBiomass<-data.frame(year="2019",Station=c(reconstructed$station),BiomIndex=c(reconstructed$predictionMax))
Index2019PredictedMissing2020Max<-calcBiomass(HaulBiomass)

#cat("Biomass Index in 2019 (true value):",Index2019True[2][[1]],"\n")
#cat("Biomass Index in 2019 (predicted value):",Index2019Predicted[2][[1]],"(",Index2019PredictedMin[2][[1]],";",Index2019PredictedMax[2][[1]],")","\n")

#cat("Biomass Index in 2019 in known 2020 locations (true value):",Index2019Known2020[2][[1]],"\n")

#cat("Biomass Index in 2019 in missing 2020 locations (true value):",Index2019Missing2020[2][[1]],"\n")
#cat("Biomass Index in 2019 in missing 2020 locations (predicted value):",Index2019PredictedMissing2020[2][[1]],"(",Index2019PredictedMissing2020Min[2][[1]],";",Index2019PredictedMissing2020Max[2][[1]],")","\n")

#format(round(x, 2), nsmall = 2)

ff<-function(x){
  return (round(x*100)/100)
}


Total_biomass_true = ff(totalbiomasstrue)
Total_biomass_predicted = ff(totalbiomassestimated)
Total_biomass_predicted_min = ff(totalMINbiomassestimated)
Total_biomass_predicted_max = ff(totalMAXbiomassestimated)

Total_biomass_in_missing_hauls_true = ff(sum(reconstructed$biomIndex2019))
Total_biomass_in_missing_hauls_perc = ff(percbiomasstrue)
Total_biomass_in_missing_hauls_predicted<-ff(sum(c(reconstructed$prediction)))
Total_biomass_in_missing_hauls_predicted_min<-ff(sum(c(reconstructed$predictionMin)))
Total_biomass_in_missing_hauls_predicted_max<-ff(sum(c(reconstructed$predictionMax)))

Total_biomass_in_known_hauls = ff(sum(truevalues$biomIndex2019))
Biomass_index_in_known_hauls = ff(Index2019Known2020[2][[1]])

Biomass_index_in_all_hauls_true = ff(Index2019True[2][[1]])
Biomass_index_in_all_hauls_predicted = ff(Index2019Predicted[2][[1]])
Biomass_index_in_all_hauls_predicted_min = ff(Index2019PredictedMin[2][[1]])
Biomass_index_in_all_hauls_predicted_max = ff(Index2019PredictedMax[2][[1]])

Biomass_index_in_missing_hauls_true = ff(Index2019Missing2020[2][[1]])
Biomass_index_in_missing_hauls_predicted = ff(Index2019PredictedMissing2020[2][[1]])
Biomass_index_in_missing_hauls_predicted_min = ff(Index2019PredictedMissing2020Min[2][[1]])
Biomass_index_in_missing_hauls_predicted_max = ff(Index2019PredictedMissing2020Max[2][[1]])

Accuracy_on_hauls_prediction = ff(accuracy)

cat("Scenario :",fileIdx,"\n")
cat("Accuracy at predicting the true values in the hauls:",Accuracy_on_hauls_prediction,"%","\n")
cat("Total biomass in known hauls (measured)    :",Total_biomass_in_known_hauls,"\n")
cat("Biomass Index in known hauls (measured)    :",Biomass_index_in_known_hauls,"\n")

cat("Total Biomass (measured)                   :",Total_biomass_true,"\n")
cat("Total Biomass (predicted)                  :",Total_biomass_predicted,"(",ff(totalMINbiomassestimated),";",ff(totalMAXbiomassestimated),")","\n")
cat("Biomass Index (measured)                   :",Biomass_index_in_all_hauls_true,"\n")
cat("Biomass Index (predicted)                  :",Biomass_index_in_all_hauls_predicted,"(",Biomass_index_in_all_hauls_predicted_min,";",Biomass_index_in_all_hauls_predicted_max,")","\n")

cat("Total biomass in missing hauls (measured)  :",Total_biomass_in_missing_hauls_true,"(",Total_biomass_in_missing_hauls_perc,"% of total)","\n")
cat("Total biomass in missing hauls (predicted) :",Total_biomass_in_missing_hauls_predicted,"(",Total_biomass_in_missing_hauls_predicted_min,";",Total_biomass_in_missing_hauls_predicted_max,")","\n")
cat("Biomass Index in missing hauls (measured)  :",Biomass_index_in_missing_hauls_true,"\n")
cat("Biomass Index in missing hauls (predicted) :",Biomass_index_in_missing_hauls_predicted,"(",Biomass_index_in_missing_hauls_predicted_min,";",Biomass_index_in_missing_hauls_predicted_max,")","\n")

cat("\n")

################################
}
