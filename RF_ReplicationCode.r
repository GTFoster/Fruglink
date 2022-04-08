install.packages("randomForest")
install.packages("pROC")
install.packages("magrittr")
install.packages("plyr")
install.packages("plyr")
install.packages("ape")
install.packages("tidyr")
install.packages("tictoc")
install.packages("PVR")
install.packages("tibble")

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#library(tidyverse)
library(randomForest)
library(pROC)
library(magrittr)
library(plyr)
library(dplyr)
library(ape)
library(tidyr)
library(tictoc)
library(PVR)
library(tibble)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
tictoc::tic()
dat <- read.csv("Data.nosync/ATLANTIC_frugivory.csv")
dat %<>% filter(., Frugivore_Species != "Carollia castanea") #This bat has an incorrect gape size
dat <- dplyr::select(dat, -ID, -Latitude, -Longitude, -Study_Location, -Precision, -Study_Method, -Study.reference, -Doi.Link, -Frug_Population_Trend, -Frug_Migration_status)

dat %<>% unique(.)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("BIEN_subtree.Rda")
Mammal_trees <- ape::read.nexus(file="Data.nosync/VertLife_FrugMam/output.nex")
Bird_trees <- ape::read.nexus(file="Data.nosync/VertLife_FrugBird/output.nex")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
dat <- dplyr::filter_at(dat, vars(Plant_Form, Frugivory_score, Lipid_Score, Fruit_color, Plant_origin), all_vars(!is.na(.))) #make sure we have scores for all our binary variables
dat %<>% dplyr::filter(., Fruit_color !="") #This Remove blank fruit color entries
dat$Plant_origin[dat$Plant_origin !="native"] <- "nonnative" #condense all nonnative plants into one (can revisit later)

#dat$real <- 1 # note that these are real interactions; important when we expand later

dat_old <- dat
dat$Frugivory_score <- paste("FrugScore", dat$Frugivory_score, sep="") #Make value names unambiguous as column names
dat$Lipid_Score <- paste("LipScore", dat$Lipid_Score, sep="")


categories <- list(Forms=unique(dat$Plant_Form), #take the unique entries of our categorical data and save them in a list; (less typing to assign column names based on them below)
                   Frugivory_score=unique(dat$Frugivory_score),
                   Lipid_Score=unique(dat$Lipid_Score),
                   Fruit_color=unique(dat$Fruit_color),
                   Plant_origin=unique(dat$Plant_origin)
                   )
dat_old <- dat

dat <- dat %>% mutate(bin=1) %>% pivot_wider(., names_from = Plant_Form, values_from = bin)
dat <- dat %>% mutate(bin=1) %>% pivot_wider(., names_from = Frugivory_score, values_from = bin)
dat <- dat %>% mutate(bin=1) %>% pivot_wider(., names_from = Lipid_Score, values_from = bin)
dat <- dat %>% mutate(bin=1) %>% pivot_wider(., names_from = Fruit_color, values_from = bin)
dat <- dat %>% mutate(bin=1) %>% pivot_wider(., names_from = Plant_origin, values_from = bin)


dat[, which(colnames(dat)=="liana"):ncol(dat)][is.na(dat[,which(colnames(dat)=="liana"):ncol(dat)])==TRUE] <- 0 #Replae all NA's in our newly created columns with 0

dat <- left_join(dat, dat_old) #add back in our factor columns just in case we want them later


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
plant_litter <- PVR::PVRdecomp(phy=BIEN_subtree, type="newick", scale=TRUE)

sum((round((plant_litter@Eigen$values)/sum((plant_litter@Eigen$values)),3)*100)[1:5]) #first three vectors contain about 42.6% of variation :(

plant_PhyEig <-data.frame(BIEN_subtree$tip.label, plant_litter@Eigen$vectors[,1:5])
colnames(plant_PhyEig) <- c("Plant_Species", paste(colnames(plant_PhyEig)[2:ncol(plant_PhyEig)], "PlDecomp", sep=""))

plant_PhyEig$Plant_Species <- gsub(pattern="_", replace=" ", x=plant_PhyEig$Plant_Species)

dat <- left_join(dat, plant_PhyEig, by="Plant_Species")
birds <- dplyr::filter(dat, Frug_Class=="Aves")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
vegstuff <- function(dat){
mat_assym <- as.data.frame.matrix(table(dat$Frugivore_Species, dat$Plant_Species)) #Rows = Frugivores, Columns = Plants
frugSim <- rowMeans(as.matrix(vegan::vegdist(mat_assym)))	 #What's my average distance from all other frugivores (higher means I'm more unique, lower means I'm more similiar)
plantSim <- rowMeans(as.matrix(vegan::vegdist(t(mat_assym)))) #Same thing for plants
simMat <- frugSim %*% t(plantSim)	#uniqueness of the plant * uniqueness of the frugivore
rownames(simMat) <- rownames(mat_assym)
data.frame(frugSim) %>% rownames_to_column("Frugivore_Species")

dat <- data.frame(frugSim) %>% rownames_to_column("Frugivore_Species") %>% left_join(dat, by="Frugivore_Species")
dat_new <- data.frame(plantSim) %>% rownames_to_column("Plant_Species") %>% left_join(dat, by="Plant_Species")
return(dat_new)
}

birds <- vegstuff(birds)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
svData <- function(dat){
mat_assym <- as.data.frame.matrix(table(dat$Frugivore_Species, dat$Plant_Species)) #Rows = Frugivores, Columns = Plants

decomp <- svd(mat_assym)#run svd decomposition
u <- decomp$u #extract first 3 axis
v <- data.frame(decomp$v) #I think this is the right dimension?

plantsSVD <-data.frame(colnames(mat_assym), v[,1], v[,2], v[,3])
colnames(plantsSVD) <- c("Plant_Species", "Psvd1", "Psvd2", "Psvd3")
frugSVD <- data.frame(colnames(t(mat_assym)), u[,1], u[,2], u[,3])
colnames(frugSVD) <- c("Frugivore_Species", "Fsvd1", "Fsvd2", "Fsvd3")

dat_wP <- left_join(dat, plantsSVD, by="Plant_Species")
dat_new <- left_join(dat_wP, frugSVD, by="Frugivore_Species")

return(dat_new)
}

birds <- svData(birds)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
p_phylo <- c("c1PlDecomp", "c2PlDecomp", "c3PlDecomp", "c4PlDecomp")
p_traits <- c("fruit_diameter", "fruit_length", "tree","liana", "palm", "scrub", "yellow", "red", "black", "brown", "green", "LipScore1", "LipScore2", "LipScore3")
p_latent <- c("Psvd1", "Psvd2", "Psvd3")

f_phylo <- c("fc1", "fc2", "fc3")
f_traits <- c("Frug_Body_Mass","Frug_Mean_Gape_Size", "FrugScore1", "FrugScore2", "FrugScore3")
f_latent <- c("Fsvd1", "Fsvd2", "Fsvd3")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
woodedWalk <- function(dat, FrugTraits, PlantTraits, class_balancing=FALSE, balance_ratio=3, returnTestTrain=FALSE){
  #Set up our data into a fully expanded edgelist
  require(tidyr)
  dat$real <- 1 #make a new column denoting all of these edges are real; important when we expand out to all pairwise connections later
  dat <- dplyr::filter_at(dat, vars(c(FrugTraits, PlantTraits)), all_vars(!is.na(.))) #make sure we have data for all our predictors
  dat <- dplyr::select_at(dat, vars(c(Frugivore_Species, Plant_Species, real, FrugTraits, PlantTraits))) #Select our relevant predictors
  dat %<>% unique(.) #Make sure we only have unique entries
  full_L <- tidyr::expand(dat, Frugivore_Species, Plant_Species) #Expand to include all possible pairwise plant-frugivore interactions.
  full_real <- dplyr::select(dat, Frugivore_Species, Plant_Species, real) %>% left_join(full_L, ., by = c("Frugivore_Species", "Plant_Species")) #notating which of our edges in this expanded data are real 
  full_real$real[is.na(full_real$real)==TRUE] <- 0

  full_real_frugs <- dat[,c("Frugivore_Species", FrugTraits)] %>% unique() %>% left_join(full_real, ., by="Frugivore_Species") #Add in our frugivore traits

  full_real_both <- dat[,c("Plant_Species", PlantTraits)] %>% unique() %>% left_join(full_real_frugs, ., by="Plant_Species") #add in our plant traits
  
  full_real_both$real <- as.factor(full_real_both$real)
  
  if(class_balancing==TRUE){
    # Subsample our edgelist
    ones <- full_real_both %>% dplyr::filter(., real==1)
    #one_IDs <- sample(round(nrow(ones)*0.8), replace=F)
    #testOnes <- ones[-one_IDs,]
    #trainOnes <- ones[one_IDs,]
    
    zeroes <- full_real_both %>% dplyr::filter(., real==0)
    zero_IDs <- sample(x=1:nrow(zeroes), size=nrow(ones)*balance_ratio, replace=F) #Ratio of train to test 1:3
    keptZeroes <- zeroes[zero_IDs,]
    
    #nrow(rbind(trainZeroes, testZeroes, trainOnes, testOnes)) == nrow(full_real_both) #Making sure we're not losing anything
    balanced <- rbind(ones, keptZeroes)
    
    balanced_IDs <- sample(x=1:nrow(balanced),size=round(nrow(balanced)*0.8), replace=F) #subsample out 80% of our edges. 
    test <- balanced[-c(balanced_IDs),] #assign 0% train set
    train <- balanced[balanced_IDs,] #assign 80% train set
  }
    
  if(class_balancing==FALSE){
       # Just Subsample our edgelist as is-don't do any balancing
    sample_IDs <- sample(x=1:nrow(full_real_both), size=round(nrow(full_real_both)*0.8), replace=F) #subsample out 80% of our edges. 
    test <- full_real_both[-sample_IDs,] #assign 80% train set
    train <- full_real_both[sample_IDs,]
  }

  rf<-randomForest::randomForest(real ~. -Frugivore_Species - Plant_Species, data=train, ntree=100)
  trainROC <- roc(train$real,rf$votes[,2])
  predictions <- stats::predict(rf, newdata=test, type="prob")
  test$S <- predictions[,2]
  
  ##########################
  #Part where we throw a bunch out for ease of replication
  ##########################
  roc(test$real, test$S)
  AUC <- as.numeric(roc(test$real, test$S)$auc)
  varimports <- rf$importance
  
  if(returnTestTrain==TRUE){
    model=list(rfModel=rf,test=test, train=train)
  }
  if(returnTestTrain==FALSE){
    #model<- rf
    model <- list(AUC, varimports)
  }
  return(model)
}


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
output <- list()
for(i in 1:100){
  #Set up Phylogeny
  frug_litter <- PVR::PVRdecomp(phy=Bird_trees[[round(i)]], type="newick", scale=TRUE)
  frug_PhyEig <-data.frame(Bird_trees[[5]]$tip.label, frug_litter@Eigen$vectors[,1:3]) #take first three vectors and make them a dataframe
  colnames(frug_PhyEig) <- c("Frugivore_Species", "fc1", "fc2", "fc3")
  frug_PhyEig$Frugivore_Species <- gsub(pattern="_", replace=" ", x=frug_PhyEig$Frugivore_Species) 
  birds_temp <- left_join(birds, frug_PhyEig, by="Frugivore_Species") #Make a new object that includes our phylo information
  
  rebalancedRF_Phy <- woodedWalk(birds_temp, FrugTraits = c(f_phylo), PlantTraits = c(p_phylo), class_balancing = TRUE)
  rebalancedRF_Traits <- woodedWalk(birds_temp, FrugTraits = c(f_traits), PlantTraits = c(p_traits), class_balancing = TRUE)
  rebalancedRF_Latent <- woodedWalk(birds_temp, FrugTraits = c(f_latent), PlantTraits = c(p_latent), class_balancing = TRUE)
  rebalancedRF_PhyTraits <- woodedWalk(birds_temp, FrugTraits = c(f_phylo, f_traits), PlantTraits = c(p_phylo, p_traits), class_balancing = TRUE)
  rebalancedRF_TraitsLatent <- woodedWalk(birds_temp, FrugTraits = c(f_latent, f_traits), PlantTraits = c(p_latent, p_traits), class_balancing = TRUE)
  rebalancedRF_PhyLatent <- woodedWalk(birds_temp, FrugTraits = c(f_latent, f_phylo), PlantTraits = c(p_latent, p_phylo), class_balancing = TRUE)
  rebalancedRF_Trio <- woodedWalk(birds_temp, FrugTraits = c(f_phylo, f_latent, f_traits), PlantTraits = c(p_phylo, p_latent, p_traits), class_balancing = TRUE)
  
  temp <- list(Phy=rebalancedRF_Phy, 
               Traits=rebalancedRF_Traits,
               Latent=rebalancedRF_Latent,
               PhyTraits=rebalancedRF_PhyTraits,
               TraitsLatent=rebalancedRF_TraitsLatent,
               PhyLatent=rebalancedRF_PhyLatent,
               Trio=rebalancedRF_Trio,
               run=i
               )
  output[[i]] <- temp
}


save(output, file="ReplicateRFsmaller.Rda")

time <- tictoc::toc()

save(time, file="replicateRuntimeob.Rda")
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
output <- data.frame(models[[7]]$importance)
output$model <- names(models)[7]
output <- rownames_to_column(output, var = "Variable")
  

for(i in 1:6){
  temp <- data.frame(models[[i]]$importance)
  temp$model <- names(models)[i]
  temp <- rownames_to_column(temp, var = "Variable")
  output <- rbind(output, temp)
}

output$MeanDecreaseGini <- as.numeric(output$MeanDecreaseGini)

