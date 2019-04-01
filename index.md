---
title: "Hierarchical Bayesian Climate Reconstruction Using Tree-Ring Data"
author: "Daniel J Hocking and Laura G Smith"
date: "01 April, 2019"
output: 
  html_document: 
    keep_md: yes
    template: sheds-template.html
---

## Abstract

Tree growth is often limited in part by temperature or precipitation. Scientists have exploited this relationship to reconstruct climate over centuries and even millennia. However, annual tree ring growth is also a result of allometric growth patterns and non-climate related environmental variables. The traditional approach to reconstruction is to remove the influence of aging while maintaining the effect of climate common to the selected trees. The resulting chronologies are then regressed against the climate variable of interest. This has proven to be an effective strategy, but the uncertainty in the standardization and chronology building steps are not propagated to the reconstruction, thereby overestimating reconstruction confidence. Conducting the full process of standardization and reconstruction in one hierarchical Bayesian model allows for full estimation of paleoclimate uncertainty. We explored models using splines, negative exponential, linear, and constant age-related detrending using independent, partially-pooled, and regional-curve standardizations. We also examined models assuming stable climate or time-varying climate using splines. Finally, we created a changepoint model to allow variation in the relationship between tree growth and climate. We applied these models to Scots pine from Tornetrask, Sweden and to a multispecies reconstruction using all ITRDB records from the NOAA New Mexico Northwest Mountain climate division (29-04). We found that reconstructions differed less depending on the standardization method in the model than on the underlying assumption of the climate model. The uncertainty among models was greater than the uncertainty within models, which may not be clearly resolved from the validation data in the instrumental record.

## Introduction

$$
log(rwl_{ikt}) = \alpha_{0ik} + \alpha_{1ik}age_{ikt} + \beta_{k}climate_{t} + \epsilon_{ikt}
$$

where $rwl_{ikt}$ is the ring-width length for tree $i$ of species $k$ in year $t$. The error is additive on the logrithmic scale and therefore multiplicative on the arthimetic scale with $\alpha_{0ik} + \alpha_{1ik}age_{ikt}$ representing the modified negative exponential growth curve varying by tree and species, $\beta_{k}$ the log-linear effect of climate on tree ring growth, and the iid error is $\epsilon_{ikt}$ such that 

$$
\epsilon_{ikt} \sim \mathcal{N}(0, \sigma_{ik}^2)
$$

with the standard deviation $\sigma$ varying by tree and species.

## Methods

### Data

### Code

## Results

## Discussion

## Literature Cited

## More Coming Soon...
