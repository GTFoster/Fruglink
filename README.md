# Research Questions:
### 1. How well are we able to predict missing frugivory links?
### 2. How do link prediction methods using trait-based, phylogenetic, or latent graph information compare? Do they overlap in the kinds of links they tend to get right, or can we combine them to be more than the sum of their parts?


### Abstract:

Due to the constraints of limited effort and sampling error, observed species interaction networks will always be an imperfect representation of the ``true'' underlying community. Link prediction methods allow us to construct a potentially more complete representation of a given empirical network by guiding targeted sampling of predicted links, as well as offer insight into potential interactions that may occur as species' ranges shift. Various data types can predict interactions; understanding how different kinds of information compare in their ability to predict links between different types of nodes is important. To this end, we compare random-forest regression models informed by combinations of phylogenetic structure, species traits, and latent network structure in their ability to predict interactions in a diverse network of fruiting plants and frugivorous birds in Brazil's Atlantic forest. We found that for our dataset, latent structure was the most important determinant of model predictive performance. While incorporating trait or phylogenetic information alongside latent features had little effect on discriminatory power, they did meaningfully increase overall model accuracy. Our results highlight the potential importance of latent structural features for predicting mutualistic interactions, and encourage a clear conceptual link between prediction performance metrics and the overall goal of predicting cryptic links.


### Repository Structure
```{bash}
├── README.md. 
├── Full_Analysis.Rmd. Main analysis file; should recreate the entire analysis start to finish provide you install the appropriate packages (easiest way to do this is to set initial chunk eval=TRUE)
├── Data.nosync
│   ├── Results. These files are all are created as outputs of my analyses. These are saved as intermediate .Rda or .csv values since the analysis can be somewhat computationally intensive
│   ├── ├── ttsplitReplicateRF_withJ.Rda Nested list; results of 100 iterations of model training of all 7 model variaitions after test-train split. Top level is 1:100 replites, lower level is the 7 model structures. W
│   ├── ├── fullReplicateRF_predicts.Rda Nested list; same approach as above, but saving the output randomforest objects. This is quite large, but necessary for calculating variable importance
│   ├── ├── fullNetPredictSummary.csv Model predictions of all 7 models averaged over 100 iterations trained on and predicting all potential interactions. 
│   ├── DataSources. Data products used to perform this analyses. Citations found below. 
│   ├── Old Stuff. The ether. Ive worked on this project for a while, and am not great at using git branches for alternate ideas or directions. This is the idea graveyard. 
├── Figures Figures as seen in the manuscript
```


### Data Citations:
Bello, Carolina, Mauro Galetti, Denise Montan, Marco A Pizo, Tatiane C Mariguela, Laurence Culot, Felipe Bufalo, et al. 2017. “Atlantic Frugivory: A Plant–Frugivore Interaction Data Set for the Atlantic Forest.” Ecology 98 (6): 1729. https://doi.org/https://doi.org/10.1002/ecy.1818.

Jetz, Walter, Gavin H Thomas, Jeffery B Joy, Klaas Hartmann, and Arne O Mooers. 2012. “The Global Diversity of Birds in Space and Time.” Nature 491 (7424): 444–48. https://doi.org/https://doi.org/10.1038/nature11631.

Maitner, Brian S, Brad Boyle, Nathan Casler, Rick Condit, John Donoghue, Sandra M Durán, Daniel Guaderrama, et al. 2018. “The Bien r Package: A Tool to Access the Botanical Information and Ecology Network (BIEN) Database.” Methods in Ecology and Evolution 9 (2): 373–79. https://doi.org/https://doi.org/10.1111/2041-210X.12861.
