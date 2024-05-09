---
title: "README"
author: "Grant Foster"
date: "9/10/2021"
output: html_document
---


#Research Questions:
## 1. How well are we able to predict missing frugivory links across diverse taxa?
## 2. How do link prediction methods using trait-based, phylogenetic, or latent graph information compare? Do they overlap in the kinds of links they tend to get right, or can we combine them to be more than the sum of their parts?
## 3. Compared to the entire network, are we able to predict better or worse on particular taxonomic or functional groups? How about for more generalist or specialist nodes?
## 4. If we train network reconstruction methods on only native interactions, how do they perform when applied to nonnative species?


## Methods:
Using Atlantic Forest data from Bello 2021: (7344 links with native plants, 857 nonnative links). Split into test and training sets, and train a variety of network reconstruction model types including: 
a) methods agnostic to node biology (preferential-attatchment models, hierarchical structure (Clauset 2008)
b )Models incorporating phylogenetic information (in the style of Elmasri 2020)
c) Trait-based methods, and
d) Combinations of the above. 

Comparing the performance of these competing models using network-wide measures of link prediction performance such as AUC will allow us to answer research questions 1&2. We can then hopefully answer question 3 by comparing the accuracy of particular link sets (ex. Can we detect a phylogenetic signal in how well we’re able to predict host sets using a metric like Moran’s I). We can then apply our best fitting link prediction models to the entire set of native hosts, and then test performance on non-native species to evaluate question 4. In addition to the main research questions outlined above, this project also allows us to predict possible unobserved links, as well as how incorporating of those links may influence overall network structure or stability. 

## Model Construction:
In order to make different methods both comparable across models as well as combinable into complex models, we have to establish a common link prediction model structure for calculating suitabilities. So far, I've done a bit of messing around with spectral curve fitting methods that use latent graph information, but these seem hard to bake into a common structure (at least right now).

### Plug & Play
I've also messed with reusing Dallas et al. 2017's Plug & Play method, which uses a bayesian justification for comparing the covariance structure of all real training nodes (where $y_{ij}=1$), compared to all possible node combinations. This approach can easily accept both continuous and binary variables, and phylogenetic decomposition allows us to turn phylogenetic information into a series of continuous eigenvectors, albeit with a significant loss of information. Still, these allows for them to be easily incorporated into a plug & play framework. 

While the plug and play method does allow a way to incorporate phylogenetic information into prediction models in comparable way, it still lives two problems. First, as I currently see it the model does not allow for ordinal data; even its method of allowing binary factorial data is by converting it to 0's and 1's when calculating the covariance structure. Maybe there is an easy workaround for this that I haven't found yet; the problem itself occurs when we're calculating covariacne matrices for the two point sets (Dallas 2017 does this using `stats::cov()`). Second, I'm not sure (yet) how you'd be able to incorporate latent methods. Maybe there's some sort of graph decomposition approach similar to the phylogenetic eigen approach that would allow you to condense that informaiton into a continous variable, or maybe a variable you get out of spectral curve fitting methods that can then be fed into the plug and play mode? I'm not sure yet. 

### Data: Bell et al 2017 (https://doi.org/10.1002/ecy.1818)
As of right now, I've been just pulling trees of Tree-Of-Life, but I'll try to get better trees soon. I'm hoping the vertebrate tree will be easy to find, and Cleber suggested the BIEN tree for the plants which seems like a really good resource. 


### Papers Methods to look at:
Preferrential Attatchment:
Spectral Curve fitting: https://towardsdatascience.com/link-prediction-in-bipartite-graph-ad766e47d75c, which is largely built around Kunegis 2011 (https://github.com/kunegis/phd)
Heirarchical: Clauset 2008 (https://doi.org/10.1038/nature06830)
Phylogenetic: Elmasri 2020 (DOI: 10.1214/19-AOAS1296)
Trait Based: 



