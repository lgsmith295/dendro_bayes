---
title: "Hierarchical Bayesian Climate Reconstruction Using Tree-Ring Data"
author: "Daniel J Hocking and Laura G Smith"
date: "05 April, 2019"
output:
  html_document:
    keep_md: yes
    template: sheds-template.html
  word_document: default
---

## Abstract

Tree growth is often limited in part by temperature or precipitation. Scientists have exploited this relationship to reconstruct climate over centuries and even millennia. However, annual tree ring growth is also a result of allometric growth patterns and non-climate related environmental variables. The traditional approach to reconstruction is to remove the influence of aging while maintaining the effect of climate common to the selected trees. The resulting chronologies are then regressed against the climate variable of interest. This has proven to be an effective strategy, but the uncertainty in the standardization and chronology building steps are not propagated to the reconstruction, thereby overestimating reconstruction confidence. Conducting the full process of standardization and reconstruction in one hierarchical Bayesian model allows for full estimation of paleoclimate uncertainty. We explored models using splines, negative exponential, linear, and constant age-related detrending using independent, partially-pooled, and regional-curve standardizations. We also examined models assuming stable climate or time-varying climate using splines. Finally, we created a changepoint model to allow variation in the relationship between tree growth and climate. We applied these models to Scots pine from Tornetrask, Sweden and to a multispecies reconstruction using all ITRDB records from the NOAA New Mexico Southwest Mountain climate division (29-04). We found that reconstructions differed less depending on the standardization method in the model than on the underlying assumption of the climate model. The uncertainty among models was greater than the uncertainty within models, which may not be clearly resolved from the validation data in the instrumental record.

## Introduction

Tree growth is often limited in part by temperature or precipitation. Scientists have exploited this relationship to reconstruct climate over centuries and even millennia. However, annual tree ring growth is also a result of allometric growth patterns and non-climate related environmental variables. The traditional approach to reconstruction is to remove the influence of aging while maintaining the effect of climate common to the selected trees. The resulting chronologies are then regressed against the climate variable of interest. This has proven to be an effective strategy, but the uncertainty in the standardization and chronology building steps are not propagated to the reconstruction, thereby overestimating reconstruction confidence. Conducting the full process of standardization and reconstruction in one hierarchical Bayesian model allows for full estimation of paleoclimate uncertainty.

### Model

Our approach is based on the classic multiplicative dendroclimatal tree-ring model with `Ring Width Length = biological growth x climate response x error`. We follow the formulation developed by Schofield et al. (2016) and generalize it to incorporate multiple species as described below.

$$
log(rwl_{ikt}) = \alpha_{0ik} + \alpha_{1ik}age_{ikt} + \beta_{k}climate_{t} + \epsilon_{ikt}
$$

where $rwl_{ikt}$ is the ring-width length for tree $i$ of species $k$ in year $t$. The error is additive on the logrithmic scale and therefore multiplicative on the arthimetic scale with $\alpha_{0ik} + \alpha_{1ik}age_{ikt}$ representing the modified negative exponential growth curve varying by tree and species, $\beta_{k}$ the log-linear effect of climate on tree ring growth, and the iid error is $\epsilon_{ikt}$ such that 

$$
\epsilon_{ikt} \sim \mathcal{N}(0, \sigma_{i}^2)
$$

with the standard deviation $\sigma$ varying by tree.

The full model with a negative exponential growth model with partial pooling (one exponential decay curve for each tree but all are drawn from a random normal distribution), a linear climate relationship, and assuming climate varying around a stable mean is represented as:

$$
log(rwl_{ikt}) = \alpha_{0ik} + \alpha_{1ik}age_{ikt} + \eta_{t} + \epsilon_{ikt}
$$
$$
\epsilon_{ikt} \sim \mathcal{N}(0, \sigma_{i}^2)
$$

$$
\eta_{t} \sim \mathcal{N}(\beta_{k}climate_{t}, \sigma_{k \eta}^2)
$$

$$
x_{t} \sim \mathcal{N}(\mu_{x}, \sigma_{x}^2)
$$

$$
\alpha_{0ik} \sim \mathcal{N}(\mu_{\alpha_0k}, \sigma_{\alpha_0k}^2)
$$

$$
\alpha_{1ik} \sim \mathcal{N}(\mu_{\alpha_1k}, \sigma_{\alpha_1k}^2)
$$

### Priors

The priors can be varied depending on prior knowledge. We take the approach of using relatively uninformative priors. Additional information can be used to constrain the model by using more informative priors. This information often comes from previous studies and expert knowledge. 

$$
\mu_{\alpha_0k} \sim \mathcal{N}(0, 10^2)
$$

$$
\mu_{\alpha_1k} \sim \mathcal{N}(0, 10^2) T(-\infty, 0)
$$

$$
\beta_{k} \sim \mathcal{N}(0, 10^2)
$$

The standard deviations from each part of the model follow a half-cauchy distribution with a scale of 5.

$$
\sigma \sim cauchy(0, 5^2) T(0, \infty)
$$

### Autoregressive (AR1)
 
Tree ring growth is often autocorrelated from one year to the next. If this isn't accounted for in the model the residuals may be autocorrelated, which violates assumptions of linear regression models. Steinschneider et al. (2017) added an autoregressive approach with a simple AR1 model and we follow that approach here.

$$
log(rwl_{ikt}) = \alpha_{0ik} + \alpha_{1ik}age_{ikt} + \eta_{t} + \delta_i log(rwl_{ikt-1}) + \epsilon_{ikt}
$$

where $\delta_i$ is the autoregressive correlation between year $t$ and $t + 1$. The prior is $\delta_i \sim unif(0, 1)$. The correlation varies by series so that each tree is allowed to have its own independent autocorrelation.

## Tornestrask

We used the Scots Pine (*Pinus sylvestris*) data from the Tornetrask area of Sweden as a well known example for examining different model options.

The first model was the same as model MB_TS_CON in Schofield et al. (2016). Model 3 in Steinschneider et al. (2017) is similar but with a Box-Cox transformation rather than a log transformation. This model assumes a modified negative exponential biological growth function, a linear relationship with climate, and that climate is stable following a normal distribution over the past 500 years (length of these data).

**This model overcomes the segment length curse.** Estimating the components of the model simultaneously and having known ages with partial pooling overcomes the segment length curse and is described in the supplements of Schofield et al. (2016). Partial pooling is a result of a hierarchical model on $\alpha_{0i}$ and $\alpha_{1i}$ so that each tree gets its own detrending curve but those curves are drawn from a normal distribution with a common mean and variance.

<div class="figure">
<img src="Results/Figures/JAGS/negexp_norm_paper.png" alt="Mean annual temperature reconstruction for Tornestrask, Sweden assuming negative exponential biological growth function, a linear relationship with climate, and that climate is stable following a normal distribution." width="75%" />
<p class="caption">Mean annual temperature reconstruction for Tornestrask, Sweden assuming negative exponential biological growth function, a linear relationship with climate, and that climate is stable following a normal distribution.</p>
</div>

### Spline Detrending

This is a similar model but with 2/3 cubic spline detrending rather than negative exponential. One limitation currently is that there is no restriction to have monotonically decreasing curves and some recently sampled young trees experiencing rapid climate warming will have increasing biological growth curves, thereby obscurring some recent climate warming.

Below is an example of the detrending from three single series that were done simultaneously with the climate reconstruction.


```r
include_graphics("Results/Figures/Detrend/Splines/detrend_tree_spl_4.pdf")
```


```r
include_graphics("Results/Figures/Detrend/Splines/detrend_tree_spl_240.pdf")
```


```r
include_graphics("Results/Figures/Detrend/Splines/detrend_tree_spl_245.pdf")
```

Below is the resulting reconstruction.

<div class="figure">
<img src="Results/Figures/JAGS/spline_norm_paper.pdf" alt="Mean annual temperature reconstruction for Tornestrask, Sweden assuming a smooth biological growth function (spline), a linear relationship with climate, and that climate is stable following a normal distribution." width="75%" />
<p class="caption">Mean annual temperature reconstruction for Tornestrask, Sweden assuming a smooth biological growth function (spline), a linear relationship with climate, and that climate is stable following a normal distribution.</p>
</div>

### Negative Exponential with AR1

<div class="figure">
<img src="Results/Figures/JAGS/negexp_linear_ar_paper.png" alt="Mean annual temperature reconstruction for Tornestrask, Sweden assuming negative exponential biological growth function, a linear relationship with climate, and that climate is stable following a normal distribution. The model also allows for autocorrelation with the previous year (AR1)." width="75%" />
<p class="caption">Mean annual temperature reconstruction for Tornestrask, Sweden assuming negative exponential biological growth function, a linear relationship with climate, and that climate is stable following a normal distribution. The model also allows for autocorrelation with the previous year (AR1).</p>
</div>

### Negative Exponential with unstable climate

To relax the assumption of a stable climate, Schofield et al. (2016) developed a model with a cubic B-spline. We followed that approach here for this reconstruction. 

<div class="figure">
<img src="Results/Figures/JAGS/negexp_spl25_paper.pdf" alt="Mean annual temperature reconstruction for Tornestrask, Sweden assuming negative exponential biological growth function, a linear relationship with climate, and that climate is unstable varying smoothly and estimated with a cubic B-Spline." width="75%" />
<p class="caption">Mean annual temperature reconstruction for Tornestrask, Sweden assuming negative exponential biological growth function, a linear relationship with climate, and that climate is unstable varying smoothly and estimated with a cubic B-Spline.</p>
</div>

A problem with this model is that estimates go into unrealistic space. Trees stop growth at very low temperatures and this changes the relationship between temperature and tree growth. Therefore this extrapolation of a linear model allows for biologically unrealistic relationships. The following models attempt to address this issue.

### Negative Exponential with variable climate relationship

For this model we assume a stable climate but allow the relationship between climate and tree growth to change with temperature. Schofield et al. (2016) used a piecewise linear regression to accomplish this because they had information from other studies to indicate that below 4 C the trees stop growth. So they had a fixed changepoint based on existing knowledge. 

<div class="figure">
<img src="Results/Figures/JAGS/negexp_1change_paper.pdf" alt="Mean annual temperature reconstruction for Tornestrask, Sweden assuming negative exponential biological growth function, a piecewise linear relationship with climate with an unknown change point, and that climate is stable following a normal distribution." width="75%" />
<p class="caption">Mean annual temperature reconstruction for Tornestrask, Sweden assuming negative exponential biological growth function, a piecewise linear relationship with climate with an unknown change point, and that climate is stable following a normal distribution.</p>
</div>

This model in an improvement allowing for more biological realism but keeping the estimated within the temperatures that allow for tree growth and survival (otherwise there would not be living trees producing rings in those years).

### Regional Curve Standardization with a Negative Exponential and Unstable Climate

Our final model using a negative exponential regional curve standardization (NegExp-RCS) where all trees are assumed to have the same growth at a given age following a negative exponential function. For this the parameters $\alpha_{0i}$ and $\alpha_{1i}$ are reduced to $\alpha_{0}$ and $\alpha_{1}$ such that they are the same across all trees of a given age. We used a basis spline with knots every 25 years.

<div class="figure">
<img src="Results/Figures/JAGS/rcs_spline_25_paper.pdf" alt="Mean annual temperature reconstruction for Tornestrask, Sweden assuming NegExp-RCS biological growth function, a linear relationship with climate, and that climate is unstable but varies smoothly over time (estimated with a B-spline)." width="75%" />
<p class="caption">Mean annual temperature reconstruction for Tornestrask, Sweden assuming NegExp-RCS biological growth function, a linear relationship with climate, and that climate is unstable but varies smoothly over time (estimated with a B-spline).</p>
</div>

This model also keeps the reconstructed climate within biologically realistic bounds but relaxes the assumption of a stable climate. **There is a large amount of variation in reconstructions depending on model assumptions, even among models that validate well and produce biologically realistic estimates.**

In the future it may be worth trying a NegExp-RCS model with a piecewise linear climate relationship, and a climate that can vary smoothly over time and adding an AR1 term (maybe unnecessary with a spline).

## Methods

After extending some previous models and testing them on the well-known Tornetrask data, we wanted to further expand this framework to incorporate multispecies reconstructions. 

We downloaded ITRDB data from New Mexico and clipped it to the NOAA Southwest Climate Division area (29-04).

<div class="figure">
<img src="Results/Figures/NM/nm_study_sites.png" alt="Tree ring data from the New Mexico Southwest Climate Division (29-04)." width="75%" />
<p class="caption">Tree ring data from the New Mexico Southwest Climate Division (29-04).</p>
</div>

We used the [`dplR`](https://cran.r-project.org/web/packages/dplR/index.html) package in R to temporarily detrend each series with a negative exponential function. We then created a chronology from these detrended RWI.

We downloaded the associated [monthly precipitation data](https://www.esrl.noaa.gov/psd/cgi-bin/data/timeseries/timeseries.pl?ntype=2&typediv=2&state=+29&averaged=11&division=4&year1=1895&year2=2019&anom=0&iseas=0&mon1=0&mon2=0&typeout=1&y1=&y2=&plotstyle=0&Submit=Create+Timeseries).  

We then used the `treeclim` package to examine the relationship of the chronology to the monthly temperatures and found that January - July correlated signficantly with the chronology. Therefore, we use the mean precipitation in centimeters for these months in our Bayesian Hierarchical model.

We removed any series that did not correlate with these climate data (spearman $\rho < 0.4$) and also any series that did not correlate with each other ($\rho < 0.4$). We then used the raw ring width lengths from these series in our model.

### Data

### Code

## Results

## Discussion

## Literature Cited

## More Coming Soon...
