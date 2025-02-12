load(file="rereqs.rda")

library(tidyverse)
library(PVR)
library(ape)
library("phytools")
library(picante)
library(magrittr)

suitout <- NULL
for(ratioval in c(1, 3, 10)){
  for(i in 1:500){
    nas <- birds[gsub(" ", "_", birds$Frugivore_Species) %in%  Bird_trees[[round(i)]]$tip.label==FALSE,] %>% dplyr::select(., Frugivore_Species)  %>% unique() #Extract Bird taxa that don't occur in our tree
    genera <- sub("_[^ ]+$", "", Bird_trees[[round(i)]]$tip.label) %>% unique() #Genera list present in our bird tree
    nas$genus <- sub(" [^ ]+$", "", nas$Frugivore_Species) #Extract only genera of the nas
    nas <- nas[nas$genus %in% genera,] #Remove species without congeneric already in the tree
    nas <- nas[grepl(" ", nas$Frugivore_Species)==TRUE,] #Remove interactions only assigned to the genus level-this leaves just 1 warbler species
    ult_birdTree <- force.ultrametric(Bird_trees[[round(i)]]) #Our tree should already be ultrametric, but just in case
    polytree <- ult_birdTree #new tree object
    for(j in 1:nrow(nas)){ #add in our missing genera as polytomies
      polytree <- phytools::add.species.to.genus(tree=polytree, species=nas$Frugivore_Species[j], where="root") #Create polytomy at the genera root. 
    }
    #length(polytree$tip.label) - length(ult_birdTree$tip.label) #Number of net gained tips
    
    #Set up Phylogeny
    frug_litter <- PVR::PVRdecomp(phy=polytree, type="newick", scale=TRUE)
    frug_PhyEig <-data.frame(polytree$tip.label, frug_litter@Eigen$vectors[,1:3]) #take first three vectors and make them a dataframe
    colnames(frug_PhyEig) <- c("Frugivore_Species", "fc1", "fc2", "fc3")
    frug_PhyEig$Frugivore_Species <- gsub(pattern="_", replace=" ", x=frug_PhyEig$Frugivore_Species) 
    birds_temp <- left_join(birds, frug_PhyEig, by="Frugivore_Species") #Make a new object that includes our phylo information
    
    rebalancedRF_Phy <- woodedWalk(birds_temp, FrugTraits = c(f_phylo), PlantTraits = c(p_phylo), class_balancing = TRUE, returnTestTrain=TRUE, balance_ratio = ratioval, output_type = "rfobject")
    rebalancedRF_Phy_Suit <- rebalancedRF_Phy$test
    rebalancedRF_Phy_Suit$S <- predict(rebalancedRF_Phy$rfModel, rebalancedRF_Phy_Suit, type="prob")[,1]
    rebalancedRF_Phy_Suit$model <- "Phy"
    rebalancedRF_Phy_Suit %<>% dplyr::select(., Frugivore_Species, Plant_Species, model, real, S)
    
    rebalancedRF_Traits <- woodedWalk(birds_temp, FrugTraits = c(f_traits), PlantTraits = c(p_traits), class_balancing = TRUE,  returnTestTrain=TRUE, balance_ratio = ratioval, output_type = "rfobject")
    rebalancedRF_Traits_Suit <- rebalancedRF_Traits$test
    rebalancedRF_Traits_Suit$S <- predict(rebalancedRF_Traits$rfModel, rebalancedRF_Traits_Suit, type="prob")[,1]
    rebalancedRF_Traits_Suit$model <- "Traits"
    rebalancedRF_Traits_Suit %<>% dplyr::select(., Frugivore_Species, Plant_Species, model, real, S)
    
    rebalancedRF_Latent <- woodedWalk(birds_temp, FrugTraits = c(f_latent), PlantTraits = c(p_latent), class_balancing = TRUE,  returnTestTrain=TRUE, balance_ratio = ratioval, output_type = "rfobject")
    rebalancedRF_Latent_Suit <- rebalancedRF_Latent$test
    rebalancedRF_Latent_Suit$S <- predict(rebalancedRF_Latent$rfModel, rebalancedRF_Latent_Suit, type="prob")[,1]
    rebalancedRF_Latent_Suit$model <- "Latent"
    rebalancedRF_Latent_Suit %<>% dplyr::select(., Frugivore_Species, Plant_Species, model, real, S)
    
    rebalancedRF_PhyTraits <- woodedWalk(birds_temp, FrugTraits = c(f_phylo, f_traits), PlantTraits = c(p_phylo, p_traits), class_balancing = TRUE,  returnTestTrain=TRUE, balance_ratio = ratioval, output_type = "rfobject")
    rebalancedRF_PhyTraits_Suit <- rebalancedRF_PhyTraits$test
    rebalancedRF_PhyTraits_Suit$S <- predict(rebalancedRF_PhyTraits$rfModel, rebalancedRF_PhyTraits_Suit, type="prob")[,1]
    rebalancedRF_PhyTraits_Suit$model <- "PhyTraits"
    rebalancedRF_PhyTraits_Suit %<>% dplyr::select(., Frugivore_Species, Plant_Species, model, real, S)
    
    rebalancedRF_TraitsLatent <- woodedWalk(birds_temp, FrugTraits = c(f_latent, f_traits), PlantTraits = c(p_latent, p_traits), class_balancing = TRUE,  returnTestTrain=TRUE, balance_ratio = ratioval, output_type = "rfobject")
    rebalancedRF_TraitsLatent_Suit <- rebalancedRF_TraitsLatent$test
    rebalancedRF_TraitsLatent_Suit$S <- predict(rebalancedRF_TraitsLatent$rfModel, rebalancedRF_TraitsLatent_Suit, type="prob")[,1]
    rebalancedRF_TraitsLatent_Suit$model <- "TraitsLatent"
    rebalancedRF_TraitsLatent_Suit %<>% dplyr::select(., Frugivore_Species, Plant_Species, model, real, S)
    
    rebalancedRF_PhyLatent <- woodedWalk(birds_temp, FrugTraits = c(f_latent, f_phylo), PlantTraits = c(p_latent, p_phylo), class_balancing = TRUE,  returnTestTrain=TRUE, balance_ratio = ratioval, output_type = "rfobject")
    rebalancedRF_PhyLatent_Suit <- rebalancedRF_PhyLatent$test
    rebalancedRF_PhyLatent_Suit$S <- predict(rebalancedRF_PhyLatent$rfModel, rebalancedRF_PhyLatent_Suit, type="prob")
    rebalancedRF_PhyLatent_Suit$model <- "PhyLatent"
    rebalancedRF_PhyLatent_Suit %<>% dplyr::select(., Frugivore_Species, Plant_Species, model, real, S)
    
    rebalancedRF_Trio <- woodedWalk(birds_temp, FrugTraits = c(f_phylo, f_latent, f_traits), PlantTraits = c(p_phylo, p_latent, p_traits), class_balancing = TRUE,  returnTestTrain=TRUE, balance_ratio = ratioval, output_type = "rfobject")
    rebalancedRF_Trio_Suit <- rebalancedRF_Trio$test
    rebalancedRF_Trio_Suit$S <- predict(rebalancedRF_Trio$rfModel, rebalancedRF_Trio_Suit, type="prob")[,1]
    rebalancedRF_Trio_Suit$model <- "Trio"
    rebalancedRF_Trio_Suit %<>% dplyr::select(., Frugivore_Species, Plant_Species, model, real, S)
    
    temp <- rbind(rebalancedRF_Traits_Suit,
                  rebalancedRF_Phy_Suit,
                  rebalancedRF_Latent_Suit,
                  rebalancedRF_PhyTraits_Suit,
                  rebalancedRF_PhyLatent_Suit,
                  rebalancedRF_TraitsLatent_Suit,
                  rebalancedRF_Trio_Suit)
    temp$run <- i
    temp$balanceratio <- ratioval
    suitout <- rbind(temp, suitout)
    print(paste(i, "complete of 100", sep=" "))
  }
}
save(suitout, file="suitout.rda")