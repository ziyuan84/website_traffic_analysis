---
title: "Analysis of Company Website Traffic"
date: "2023-03-14"
output:
  html_document:
    keep_md: true
  extra_dependencies: ['amsmath', 'hyperref']
monofont: "Roboto Mono"
header-includes:
  - \usepackage{fontspec}
  - \newfontfamily\urlfont{Roboto Mono}

---



## Abstract

It has been suggested that dwell times on webpages may be modeled using a Weibull distribution, and that Web browsing exhibits a significant "negative ageing" phenomenon, that is, the rate of Webpage abandonment decreases over time  <span class="citation">(Liu, White, and
Dumais 2010)</span>. The goal of this project is to assess the suitability of the normal and Weibull distributions for modeling the average session duration of visitors to a specific website `(http://www.polygraphis.com)`, and to determine whether average session durations for this website have a negative ageing effect. Furthermore, we classify visitors by the region they come from and use a hierarchical Weibull model to study the average session duration of visitors from each region. We then fit the model using JAGS and R and perform convergence diagnostics as well as residual analysis, and assess predictive performance of the model. These results are then compared to those obtained by fitting a normal distribution to the data. 

<p>Our results are as follows. First, it is observed that in both the Weibull and normal models, MCMC convergence of the monitored parameters is achieved after 5000 iterations, and that the corresponding autocorrelation values decrease sharply to 0. Second, the Deviance Information Criterion (DIC) value for the Weibull model turns out to be smaller than that for the normal model. Third, the residual plots for both models do not appear to follow any pattern, but the mean observational level residuals for the normal model are generally lower than those in the Weibull model. That is to say, it appears that the normal model is a better fit to the data than the Weibull model. Fourth, we confirm that the posterior Weibull distribution for each region is indeed negatively ageing, that is, the hazard function for the pdf in each case is decreasing.</p>


## Introduction

The main questions addressed are: first, is the Weibull distribution appropriate for modeling the average session duration of visitors to the website `http://www.polygraphis.com`; second, do average session durations on this website exhibit a negative ageing effect? We propose to analyze the daily average session duration of visitors, classified by geographic region, that arrived at the website via organic search from 6 March 2021 to 12 March 2023. We will apply a hierarchical model to the data, where each region has its own Weibull distribution with a common normal distribution on the shape parameter and a gamma distribution $\Gamma(\alpha,\beta)$ on the scale parameter such that the mean $\frac{\alpha}{\beta}$ and standard deviation $\frac{\sqrt{\alpha}}{\beta}$ are modeled by gamma distributions themselves. Using JAGS and R, we will then fit the model and address the first question by performing convergence diagnostics as well as assessing predictive performance. The second question will be addressed by determining the hazard function of the distribution for each region using the mean parameters of the model obtained.    

## Data

The data were collected from the company's Google Analytics account. The head of the extracted raw dataset is show below.


```
##   Continent     Date Users Avg..Session.Duration
## 1      Asia 20220429    70              00:01:53
## 2      Asia 20220430    65              00:02:16
## 3      Asia 20210514    50              00:01:30
## 4      Asia 20210512    49              00:01:00
## 5      Asia 20220428    42              00:01:29
## 6      Asia 20210420    41              00:01:36
```
The dataset only included organic search traffic. A *session* is defined to be a period of time during which a user interacts with the website; it is initiated when a user that is currently not in any active session views a page,and it ends after 30 minutes of user inactivity, when the user enters the last page without an engagement hit, or at the moment of the last engagement hit on the last page. The *average session duration* for a given date and continent is calculated by dividing the total duration of all sessions on that date of users from the given continent by the corresponding number of sessions <span class="citation">(<span>“Session Duration, Avg - Analytics Help”</span>
2023)</span>. The average session durations in the dataset were converted into seconds, and only rows with a well-defined continent (one of the 5 continents - Africa, Americas, Asia, Europe or Oceania) and an average session duration of at least two seconds were kept. Each continent was converted into a number which denoted its position in the lexicographical ordering of all 5 continents (so Africa was mapped to 1, Americas to 2, Asia to 3, Europe to 4 and Oceania to 5). The cleaned dataset had 1128 rows and 2 columns, one for the continent and the other for the average session duration (in seconds). As a side note, we remark that the dwell times considered by <span class="citation">(Liu, White,
and Dumais 2010)</span> are calculated for single webpages and not for whole sessions.  


```
##   Continent_number average_session_duration
## 1                3                      113
## 2                3                      136
## 3                3                       90
## 4                3                       60
## 5                3                       89
## 6                3                       96
```
Before taking into account the distribution of the average session durations across different regions, we observed that the density line of the data appeared to increase sharply at first, and after reaching its peak seemed to approximately follow a Weibull distribution with negative ageing. The boxplots showed the presence of quite a few outliers for each continent. 

![](website_traffic_analysis_files/figure-html/unnamed-chunk-3-1.png)<!-- -->![](website_traffic_analysis_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

## Model
The Weibull distribution seems to be quite a natural model for dwell times on webpages. The idea is that web users often quickly scan through a webpage for relevance; they abandon it as soon as they discover that it is irrelevant, but are more likely to stay on after the short window for assessing its relevance is over and when they dwell longer on the webpage. The full hierarchical Weibull model is specified in the following code.

```r
mod_string = " model {
for (i in 1:length(Continent)) {
  average_session_duration[i] ~ dweib(v,lam[Continent[i]])
}

for (j in 1:5) {
  lam[j] ~ dgamma(alpha, beta)
  mu[j] = (exp(loggam(1+1/v)))*lam[j]^(-1/v)
}

alpha = mu0^2 / sig0^2
beta = mu0 / sig0^2
v ~ dnorm(3.0,1.0)

mu0 ~ dgamma(.1, 10.0)
sig0 ~ dgamma(.1,10.0)

} "
```
We postulated that the daily average session duration for visitors from each continent $i$ follows a Weibull distribution with shape parameter $v \sim N(3.0,1.0)$ and scale parameter $\lambda_i \sim \Gamma(\alpha,\beta)$, where the mean $\mu_0 = \frac{\alpha}{\beta}$ and standard deviation $\sigma_0 = \frac{\sqrt{\alpha}}{\beta}$ both follow a gamma distribution, in particular $\frac{\alpha}{\beta} \sim \Gamma(0.1,10.0)$ and $\frac{\sqrt{\alpha}}{\beta} \sim \Gamma(3.0,2.5)$. The parameters $a,b,a',b'$ for the distributions $\Gamma(a,b)$ and $\Gamma(a',b')$ of $\mu_0$ and $\sigma_0$ respectively were chosen after simulating draws from prior distributions for different values of $a,b,a',b'$ between $0$ and $10$. The mean $\mu_i$ of the Weibull distribution with shape parameter $v$ and scale parameter $\lambda_i$ was calculated using the formula $\mu_i = \Gamma(1+\frac{1}{v})\cdot \lambda_i^{-\frac{1}{v}}$; here $\Gamma$ denotes the gamma function. Further, we note that when the shape parameter $v$ is known, the conjugate prior of the likelihood function of $\lambda_i$ conditioned on a set of data points is a gamma distribution <span class="citation">(Fink 1997)</span>; thus the gamma distribution was chosen as the prior distribution for the scale parameter. Another reason for choosing the gamma distribution is to ensure that the values of the scale parameter are always positive; otherwise, the Weibull pdf for the average session duration might be undefined. Based on the mean estimates of the model parameters $\lambda_i$ and $v$ for each continent $i$, the first question of this project could be addressed by calculating the mean residuals for each continent $i$. The second question could be addressed by determining, for each continent $i$, whether or not the hazard function $h(t) = v\lambda_it^{v-1}$ is decreasing; if so, then this would confirm the negative ageing effect of average session durations, and otherwise it would suggest that positive ageing is observed instead.

Having specified the model, we set a seed of 113 and fitted the model using JAGS and R. The monitored parameters were $v$, $\mu$, $\lambda$, $\mu_0$ and $\sigma_0$ ($\mu$ and $\lambda$ are vectors of length 5, where each coordinate corresponds to a continent). The algorithm was run for three chains with a burn-in period of 1000 iterations, and each chain was run for 1e4 iterations. It was observed that convergence for each parameter was achieved after 1e4 iterations. The Gelman-Rubin diagnostic potential reduction factor for every parameter except $\mu_0$ and $\sigma_0$ also converged to $1$; the shrink factor for $\mu_0$ reduced to between 1.03 and 1.09 while that for $\sigma_0$ reduced to between 1.07 and 1.15. The trace, density and autocorrelation plots for $\lambda_1$ and $v$ (for all three chains) are shown below.   


![](website_traffic_analysis_files/figure-html/unnamed-chunk-6-1.png)<!-- -->![](website_traffic_analysis_files/figure-html/unnamed-chunk-6-2.png)<!-- -->![](website_traffic_analysis_files/figure-html/unnamed-chunk-6-3.png)<!-- -->![](website_traffic_analysis_files/figure-html/unnamed-chunk-6-4.png)<!-- -->

We also calculated a value of 12848 for the Deviance Information Criterion (DIC); this value would later be compared to the DIC value obtained for the alternative, normal model. 

In order to carry out residual analysis, we divided all the daily average session durations into 5 sublists, one for each of the 5 continents, and subtracted the mean average session duration $\mu_i$ for each continent $i$ from every element in the corresponding sublist for continent $i$. The residual plots are shown below.   

![](website_traffic_analysis_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

The mean as well as the standard deviation of the residuals for the 5 continents are listed below in ascending order of the continent indices. 


```
## Mean residuals: 7.016164 3.157353 -2.69756 19.45388 -21.53142
```

```
## Standard deviation of residuals: 327.9693 144.0456 127.6912 501.8171 99.26842
```

Apart from the presence of outliers, there did not appear to be any unusual patterns in the residual plots. We also compared the continent average session duration means to the overall mean average session duration. The latter was calculated using the formula $\mu = \Gamma(1+\frac{1}{v})\cdot\mu_0^{-\frac{1}{v}}$; here $\mu_0$ denotes the overall mean of the scale parameter $\lambda$. 
The plot of the differences $\mu_i - \mu$ for $i \in \{1,2,3,4,5\}$ is shown below.

![](website_traffic_analysis_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

We also used the posterior samples to get Monte Carlo estimates of the mean average session duration. This was done by drawing $3e4$ samples of the scale parameter $\lambda$ from the gamma distribution with shape $\frac{\mu_0^2}{\sigma_0^2}$ and rate $\frac{\mu_0}{\sigma_0^2}$ (thus taking into account the uncertainty in the parameters $\alpha$ and $\beta$), then, for the $j$-th sample $\lambda(j)$ drawn, calculating the corresponding mean $\lambda(j)^{-\frac{1}{v_j}}\cdot\Gamma(1+\frac{1}{v_j})$, where $v_j$ is the $j$-th sample drawn from the distribution $N(3,1)$. Quantiles of the estimates, ranging from 0% to 100% in intervals of 5%, are shown below. While the 95% quantile of the simulated distribution was approximately 650.8, there were a number of extreme outliers (the maximum value being 7.238549e14), resulting in the mean being unrealistically large (24128943778). We would thus suggest using the median of 147.4 (correct to 1 decimal place) as a measure of central tendency. In addition, we used the 3e4 samples of $\mu_0$ (the mean of the distribution of the scale parameter $\lambda$) to simulate the average  session duration for a new continent; the histogram of the simulated posterior predictive distribution is shown below.


```
##           0%           5%          10%          15%          20%          25% 
## 4.902418e+00 5.136580e+01 6.668267e+01 7.829630e+01 8.814285e+01 9.750929e+01 
##          30%          35%          40%          45%          50%          55% 
## 1.067941e+02 1.160488e+02 1.253589e+02 1.351033e+02 1.460321e+02 1.584667e+02 
##          60%          65%          70%          75%          80%          85% 
## 1.721120e+02 1.882833e+02 2.061260e+02 2.303778e+02 2.638039e+02 3.138271e+02 
##          90%          95%         100% 
## 4.076356e+02 6.430286e+02 2.037462e+11
```

![](website_traffic_analysis_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

For comparison, we tried fitting an alternative model to the data using the normal distribution. The mean and precision were modeled using gamma distributions. The mean and standard deviation of the parameters of the mean distribution were also modeled using gamma distributions.


```r
mod1_string = " model {
for (i in 1:length(Continent)) {
  average_session_duration[i] ~ dnorm(mu[Continent[i]],prec[Continent[i]])
}

for (j in 1:5) {
  mu[j] ~ dgamma(alpha, beta)
  prec[j] ~ dgamma(50.0, 5*1.0/2.0)

}

alpha = mu0^2 / sig0^2
beta = mu0 / sig0^2
sig = 1/sqrt(prec)

mu0 ~ dgamma(124.0, 1.0)
sig0 ~ dgamma(50.0,1/10.0)


} "
```

As was the case for the Weibull model, convergence was observed after 1e4 iterations. However, it was observed that even though the DIC value for the normal model was 14821 - larger than the DIC value of 12848 for the Weibull model, the mean residuals were smaller in the normal model, as shown below. 




```
## Mean residuals: 2.8081 0.5070523 0.3057574 5.304358 0.8125322
```

```
## Standard deviation of residuals: 327.9693 144.0456 127.6912 501.8171 99.26842
```

The posterior samples were again drawn to get Monte Carlo estimates of the mean average session duration. In comparison to the Weibull model, the range of values of the samples of the mean drawn for the normal model seemed to be more reasonable (lying between 0 and 9.051076e3). Using these samples, we simulated the distribution of the  average session duration (thus taking into account the uncertainty in the mean average session duration).  


```
##           0%           5%          10%          15%          20%          25% 
## 6.900533e-73 3.577576e-12 1.790361e-08 2.088687e-06 4.992267e-05 6.529384e-04 
##          30%          35%          40%          45%          50%          55% 
## 4.763267e-03 2.378462e-02 9.840787e-02 3.099742e-01 9.045118e-01 2.420467e+00 
##          60%          65%          70%          75%          80%          85% 
## 5.625104e+00 1.284010e+01 2.631738e+01 5.246707e+01 9.933563e+01 1.855677e+02 
##          90%          95%         100% 
## 3.537972e+02 7.664404e+02 1.040256e+04
```

![](website_traffic_analysis_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

## Results and Conclusions

The following is the summary of the posterior distribution for the Weibull model.

```
## 
## Iterations = 1:30000
## Thinning interval = 1 
## Number of chains = 1 
## Sample size per chain = 30000 
## 
## 1. Empirical mean and standard deviation for each variable,
##    plus standard error of the mean:
## 
##             Mean        SD  Naive SE Time-series SE
## lam[1]   0.01229  0.002116 1.221e-05      4.185e-05
## lam[2]   0.03150  0.003254 1.878e-05      8.864e-05
## lam[3]   0.02497  0.002498 1.442e-05      7.011e-05
## lam[4]   0.01284  0.001894 1.094e-05      4.573e-05
## lam[5]   0.02454  0.006275 3.623e-05      7.858e-05
## mu[1]  270.76565 45.425168 2.623e-01      2.724e-01
## mu[2]   82.36395  5.081795 2.934e-02      3.117e-02
## mu[3]  109.77293  5.880096 3.395e-02      3.378e-02
## mu[4]  253.39095 29.547687 1.706e-01      1.767e-01
## mu[5]  121.24571 38.251433 2.208e-01      2.467e-01
## mu0      0.02225  0.006285 3.629e-05      1.284e-04
## sig0     0.01131  0.006612 3.818e-05      1.480e-04
## v        0.80713  0.018104 1.045e-04      5.548e-04
## 
## 2. Quantiles for each variable:
## 
##             2.5%       25%       50%       75%     97.5%
## lam[1] 8.599e-03 1.079e-02 1.215e-02   0.01361   0.01685
## lam[2] 2.550e-02 2.926e-02 3.136e-02   0.03358   0.03837
## lam[3] 2.036e-02 2.325e-02 2.486e-02   0.02658   0.03022
## lam[4] 9.508e-03 1.151e-02 1.270e-02   0.01402   0.01692
## lam[5] 1.413e-02 2.008e-02 2.390e-02   0.02834   0.03853
## mu[1]  1.955e+02 2.386e+02 2.658e+02 297.64106 372.99270
## mu[2]  7.295e+01 7.884e+01 8.218e+01  85.70368  92.81407
## mu[3]  9.883e+01 1.057e+02 1.096e+02 113.63038 121.81466
## mu[4]  2.019e+02 2.326e+02 2.513e+02 271.76231 317.07145
## mu[5]  6.604e+01 9.448e+01 1.147e+02 140.71506 213.41181
## mu0    1.300e-02 1.822e-02 2.130e-02   0.02504   0.03765
## sig0   4.811e-03 7.269e-03 9.475e-03   0.01306   0.02939
## v      7.717e-01 7.948e-01 8.072e-01   0.81923   0.84255
```

Here is the summary of the posterior distribution for the normal model.


```
## 
## Iterations = 1:30000
## Thinning interval = 1 
## Number of chains = 1 
## Sample size per chain = 30000 
## 
## 1. Empirical mean and standard deviation for each variable,
##    plus standard error of the mean:
## 
##          Mean     SD Naive SE Time-series SE
## mu[1]  274.97 26.422  0.15255        0.19599
## mu[2]   85.01  6.475  0.03738        0.04921
## mu[3]  106.77  5.115  0.02953        0.03795
## mu[4]  267.54 34.524  0.19932        0.25655
## mu[5]   98.90  9.089  0.05247        0.06795
## mu0    131.98 11.269  0.06506        0.08350
## sig[1] 195.18 11.252  0.06496        0.06496
## sig[2] 128.95  4.108  0.02372        0.02372
## sig[3] 117.50  3.272  0.01889        0.01869
## sig[4] 368.53 17.846  0.10303        0.10625
## sig[5]  33.89  2.273  0.01312        0.01312
## sig0   421.12 63.115  0.36439        0.46906
## 
## 2. Quantiles for each variable:
## 
##          2.5%   25%    50%    75%  97.5%
## mu[1]  223.35 257.2 275.02 292.81 326.42
## mu[2]   72.21  80.7  85.02  89.35  97.56
## mu[3]   96.87 103.3 106.71 110.23 116.96
## mu[4]  199.62 244.2 267.34 290.83 335.00
## mu[5]   81.11  92.8  98.85 105.00 116.83
## mu0    110.71 124.2 131.73 139.31 155.07
## sig[1] 174.63 187.4 194.64 202.48 218.58
## sig[2] 121.21 126.1 128.82 131.63 137.25
## sig[3] 111.28 115.3 117.43 119.66 124.11
## sig[4] 335.63 356.2 367.84 380.11 405.58
## sig[5]  29.79  32.3  33.76  35.34  38.67
## sig0   306.52 376.8 417.60 461.81 553.37
```

The mean values of the posterior average session durations for each of the 5 continents in the Weibull model appeared to be fairly similar to the corresponding values in the normal model. However, while the residual errors for the normal model were smaller than those for the Weibull model, the DIC value of 12484 for the Weibull model  was smaller than the DIC value of 14821 for the normal model, suggesting that the Weibull model may be preferable to the normal model. 

Finally, the second question of this project may be answered by looking at the mean value of the shape parameter $v$ as well as the mean value of the scale parameter $\lambda_i$ ($i \in \{1,2,3,4,5\}$) for each of the 5 continents. The hazard function of the Weibull distribution with shape parameter $v$ and scale parameter $\lambda$ is given by $h(t) = v\lambda_it^{v-1}$, and $h(t)$ is a positive decreasing function iff $v < 1$ and $\lambda_i > 0$. The plots below show that the hazard functions for all 5 distributions are indeed decreasing, thus confirming the phenomenon of negative ageing.

![](website_traffic_analysis_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

As briefly mentioned earlier, one caveat is that the Weibull model we fitted has quite a number of outliers. Furthermore, a limitation of the model is that the shape parameter $v$ follows a normal distribution with fixed parameters. It might be reasonable to allow the parameters of the distribution of $v$ to have distributions themselves, thus reflecting the variability of the mean values of $v$ for the different continents. Another possible improvement is to consider fitting a mixed model, for example a mixture of two Weibull distributions, which might reduce the residual errors.       

## References

<div id="refs" class="references csl-bib-body hanging-indent">
<div id="ref-Fink97" class="csl-entry">
Fink, Daniel. 1997. <span>“A Compendium of Conjugate Priors.”</span> <a href="https://www.johndcook.com/CompendiumOfConjugatePriors.pdf">https://www.johndcook.com/CompendiumOfConjugatePriors.pdf</a>.
</div>
<div id="ref-Liu10" class="csl-entry">
Liu, Chao, Ryen W. White, and Susan Dumais. 2010. <span>“Understanding
Web Browsing Behaviors Through Weibull Analysis of Dwell Time.”</span>
In <em>Proceedings of the 33rd International ACM SIGIR Conference on
Research and Development in Information Retrieval</em>, 379–86. SIGIR
’10. New York, NY, USA: Association for Computing Machinery.
</div>
<div id="ref-Google_analytics_avg_session_duration" class="csl-entry">
<span>“Session Duration, Avg - Analytics Help.”</span> 2023. <a href="https://support.google.com/analytics/answer/1006253" class="uri">https://support.google.com/analytics/answer/1006253</a>.
</div>
</div>
