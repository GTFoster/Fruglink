library(tidyverse)
library(magrittr)
library(igraph)
library(pROC)
library(BIEN)
library(np)
library(PVR)
library(ape)
library("phytools")
library(tictoc)
library(corrplot)

load("Data.nosync/FullData/fullReplicateRF_rfobs_2.Rda")

#' Summarize Gini impurity score across one-hot encoded variables
#' Does this by created a weighted mean impurity score, where weight is based on frequency of the variable
#'
#' @citation I DID NOT WRITE THIS CODE! Credit to Max Ghenis on Stack Overflow 
#' https://stats.stackexchange.com/questions/92419/relative-importance-of-a-set-of-predictors-in-a-random-forests-classification-in
#' 
#' @param rf.obj Random forest object you want to compute variable importance for
#' @param groups a list, where each element is a character vector of rownames you want to combine
#' 
#' @return A 1 dimensional array, where rownames are the groups of variables and associated values are weighted average importance values (accoring to Gini coefficient)
#' 
group.importance <- function(rf.obj, groups) {
  var.share <- function(rf.obj, members) {
    count <- table(rf.obj$forest$bestvar)[-1]
    names(count) <- names(rf.obj$forest$ncat)
    share <- count[members] / sum(count[members])
    return(share)
  }
  var.imp <- as.matrix(sapply(groups, function(g) {
    sum(randomForest::importance(rf.obj, 2)[g, ]*var.share(rf.obj, g))
  }))
  colnames(var.imp) <- "MeanDecreaseGini"
  return(var.imp)
}

output <- rf_output

temp <- NULL
for(i in 1:2){
  # Phy
  temp <- rbind(temp, data.frame(giniDec=output[[i]]$Phy$importance, model="Phy", run=i, trait=c("FrugPhy1", "FrugPhy2", "FrugPhy3", "Plantphy1", "Plantphy2", "Plantphy3","Plantphy4")))
  
  
  # Traits
  traitstemp <- group.importance(rf.obj=output[[i]]$Traits, groups=list(
    growthform=c("tree","liana", "palm", "scrub"), 
    fruitcolor=c("yellow", "red", "black", "brown", "green"),
    LipidScore=c("Lipid_Score"),
    fruit_diameter=c("fruit_diameter"),
    fruit_length=c("fruit_length"),
    Frug_Body_Mass="Frug_Body_Mass",
    Frug_Mean_Gape_Size="Frug_Mean_Gape_Size",
    FrugScore="Frugivory_score")
  )
  temp <- rbind(temp, data.frame(giniDec=traitstemp, model="Traits", run=i, trait=rownames(traitstemp)))
  rm(traitstemp)
  #Latent
  temp <- rbind(temp, data.frame(giniDec=output[[i]]$Latent$importance, model="Latent", run=i, trait=c("p_latent1", "p_latent2","p_latent3","f_latent1","f_latent2","f_latent3")))
  
  #PhyTraits
  phytraitstemp <- group.importance(rf.obj=output[[i]]$PhyTraits, groups=list(
    growthform=c("tree","liana", "palm", "scrub"), 
    fruitcolor=c("yellow", "red", "black", "brown", "green"),
    LipidScore=c("Lipid_Score"),
    fruit_diameter=c("fruit_diameter"),
    fruit_length=c("fruit_length"),
    Frug_Body_Mass="Frug_Body_Mass",
    Frug_Mean_Gape_Size="Frug_Mean_Gape_Size",
    FrugScore="Frugivory_score",
    FrugPhy1="fc1",
    FrugPhy2="fc2",
    FrugPhy3="fc3",
    Plantphy1="c1PlDecomp",
    Plantphy2="c2PlDecomp",
    Plantphy3="c3PlDecomp",
    Plantphy4="c4PlDecomp"))
  temp <- rbind(temp, data.frame(giniDec=phytraitstemp, model="PhyTraits", run=i, trait=rownames(phytraitstemp)))
  rm(phytraitstemp)
  #TraitsLatent
  traitslatenttemp <- group.importance(rf.obj=output[[i]]$TraitsLatent, groups=list(
    growthform=c("tree","liana", "palm", "scrub"), 
    fruitcolor=c("yellow", "red", "black", "brown", "green"),
    LipidScore=c("Lipid_Score"),
    fruit_diameter=c("fruit_diameter"),
    fruit_length=c("fruit_length"),
    Frug_Body_Mass="Frug_Body_Mass",
    Frug_Mean_Gape_Size="Frug_Mean_Gape_Size",
    FrugScore="Frugivory_score",
    p_latent1="Psvd1",
    p_latent2="Psvd2",
    p_latent3="Psvd3",
    f_latent1="Fsvd1",
    f_latent2="Fsvd2",
    f_latent3="Fsvd3"))
  
  temp <- rbind(temp, data.frame(giniDec=traitslatenttemp, model="TraitsLatent", run=i, trait=rownames(traitslatenttemp)))
  rm(traitslatenttemp)
  #PhyLatent
  phylatenttemp <- group.importance(rf.obj=output[[i]]$PhyLatent, groups=list(
    p_latent1="Psvd1",
    p_latent2="Psvd2",
    p_latent3="Psvd3",
    f_latent1="Fsvd1",
    f_latent2="Fsvd2",
    f_latent3="Fsvd3",
    FrugPhy1="fc1",
    FrugPhy2="fc2",
    FrugPhy3="fc3",
    Plantphy1="c1PlDecomp",
    Plantphy2="c2PlDecomp",
    Plantphy3="c3PlDecomp",
    Plantphy4="c4PlDecomp"))
  temp <- rbind(temp, data.frame(giniDec=phylatenttemp, model="PhyLatent", run=i, trait=rownames(phylatenttemp)))
  rm(phylatenttemp)
  #Trio
  triotemp <- group.importance(rf.obj=output[[i]]$Trio, groups=list(
    growthform=c("tree","liana", "palm", "scrub"), 
    fruitcolor=c("yellow", "red", "black", "brown", "green"),
    LipidScore=c("Lipid_Score"),
    fruit_diameter=c("fruit_diameter"),
    fruit_length=c("fruit_length"),
    Frug_Body_Mass="Frug_Body_Mass",
    Frug_Mean_Gape_Size="Frug_Mean_Gape_Size",
    FrugScore="Frugivory_score",
    p_latent1="Psvd1",
    p_latent2="Psvd2",
    p_latent3="Psvd3",
    f_latent1="Fsvd1",
    f_latent2="Fsvd2",
    f_latent3="Fsvd3",
    FrugPhy1="fc1",
    FrugPhy2="fc2",
    FrugPhy3="fc3",
    Plantphy1="c1PlDecomp",
    Plantphy2="c2PlDecomp",
    Plantphy3="c3PlDecomp",
    Plantphy4="c4PlDecomp"))
  temp <- rbind(temp, data.frame(giniDec=triotemp, model="Trio", run=i, trait=rownames(triotemp)))
  rm(traitslatenttemp)
}

varimport <- temp

#varimporSummary <- varimport %>% dplyr::group_by(., model, trait) %>% dplyr::summarise(avg=(mean(MeanDecreaseGini)/max(MeanDecreaseGini)), sd=sd((MeanDecreaseGini)/MeanDecreaseGini))

varimporSummary <- varimport %>% dplyr::group_by(., model, trait) %>% dplyr::summarise(avg=(mean(MeanDecreaseGini)), sd=sd((MeanDecreaseGini)))

varimporSummary$trait %<>% as.factor()

save(varimporSummary, file="Data.nosyncVarImportSummary.Rda")