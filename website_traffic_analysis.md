---
title: "Analysis of Company Website Traffic"
date: "2023-03-14"
output:
  html_document:
    keep_md: true
  pdf_document:
  extra_dependencies: ['amsmath', 'hyperref']
monofont: "Roboto Mono"
header-includes:
  - \usepackage{fontspec}
  - \newfontfamily\urlfont{Roboto Mono}
bibliography: bayesianbib.bib
---



\urlstyle{tt}

## Abstract

<!-- <style>pre {width: 600px; border: 0; display: inline; margin: 0; padding: 0; white-space: normal}</style> -->

It has been suggested that dwell times on webpages may be modeled using a Weibull distribution, and that Web browsing exhibits a significant "negative ageing" phenomenon, that is, the rate of Webpage abandonment decreases over time [@Liu10]. The goal of this project is to assess the suitability of the normal and Weibull distributions for modeling the average session duration of visitors to a specific website `(http://www.polygraphis.com)`, and to determine whether average session durations for this website have a negative ageing effect. Furthermore, we classify visitors by the region they come from and use a hierarchical Weibull model to study the average session duration of visitors from each region. We then fit the model using JAGS and R and perform convergence diagnostics as well as residual analysis, and assess predictive performance of the model. These results are then compared to those obtained by fitting a normal distribution to the data. 

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
The dataset only included organic search traffic. A *session* is defined to be a period of time during which a user interacts with the website; it is initiated when a user that is currently not in any active session views a page,and it ends after 30 minutes of user inactivity, when the user enters the last page without an engagement hit, or at the moment of the last engagement hit on the last page. The *average session duration* for a given date and continent is calculated by dividing the total duration of all sessions on that date of users from the given continent by the corresponding number of sessions [@Google_analytics_avg_session_duration]. The average session durations in the dataset were converted into seconds, and only rows with a well-defined continent (one of the 5 continents - Africa, Americas, Asia, Europe or Oceania) and an average session duration of at least two seconds were kept. Each continent was converted into a number which denoted its position in the lexicographical ordering of all 5 continents (so Africa was mapped to 1, Americas to 2, Asia to 3, Europe to 4 and Oceania to 5). The cleaned dataset had 1128 rows and 2 columns, one for the continent and the other for the average session duration (in seconds). As a side note, we remark that the dwell times considered by [@Liu10] are calculated for single webpages and not for whole sessions.  


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
We postulated that the daily average session duration for visitors from each continent $i$ follows a Weibull distribution with shape parameter $v \sim N(3.0,1.0)$ and scale parameter $\lambda_i \sim \Gamma(\alpha,\beta)$, where the mean $\mu_0 = \frac{\alpha}{\beta}$ and standard deviation $\sigma_0 = \frac{\sqrt{\alpha}}{\beta}$ both follow a gamma distribution, in particular $\frac{\alpha}{\beta} \sim \Gamma(0.1,10.0)$ and $\frac{\sqrt{\alpha}}{\beta} \sim \Gamma(3.0,2.5)$. The parameters $a,b,a',b'$ for the distributions $\Gamma(a,b)$ and $\Gamma(a',b')$ of $\mu_0$ and $\sigma_0$ respectively were chosen after simulating draws from prior distributions for different values of $a,b,a',b'$ between $0$ and $10$. The mean $\mu_i$ of the Weibull distribution with shape parameter $v$ and scale parameter $\lambda_i$ was calculated using the formula $\mu_i = \Gamma(1+\frac{1}{v})\cdot \lambda_i^{-\frac{1}{v}}$; here $\Gamma$ denotes the gamma function. Further, we note that when the shape parameter $v$ is known, the conjugate prior of the likelihood function of $\lambda_i$ conditioned on a set of data points is a gamma distribution [@Fink97]; thus the gamma distribution was chosen as the prior distribution for the scale parameter. Another reason for choosing the gamma distribution is to ensure that the values of the scale parameter are always positive; otherwise, the Weibull pdf for the average session duration might be undefined. Based on the mean estimates of the model parameters $\lambda_i$ and $v$ for each continent $i$, the first question of this project could be addressed by calculating the mean residuals for each continent $i$. The second question could be addressed by determining, for each continent $i$, whether or not the hazard function $h(t) = v\lambda_it^{v-1}$ is decreasing; if so, then this would confirm the negative ageing effect of average session durations, and otherwise it would suggest that positive ageing is observed instead.

Having specified the model, we set a seed of 113 and fitted the model using JAGS and R. The monitored parameters were $v$, $\mu$, $\lambda$, $\mu_0$ and $\sigma_0$ ($\mu$ and $\lambda$ are vectors of length 5, where each coordinate corresponds to a continent). The algorithm was run for three chains with a burn-in period of 1000 iterations, and each chain was run for 1e4 iterations. It was observed that convergence for each parameter was achieved after 1e4 iterations. The Gelman-Rubin diagnostic potential reduction factor for every parameter except $\mu_0$ and $\sigma_0$ also converged to $1$; the shrink factor for $\mu_0$ reduced to between 1.03 and 1.09 while that for $\sigma_0$ reduced to between 1.07 and 1.15. The trace, density and autocorrelation plots for $\lambda_1$ and $v$ (for all three chains) are shown below.   


![](website_traffic_analysis_files/figure-html/unnamed-chunk-6-1.png)<!-- -->![](website_traffic_analysis_files/figure-html/unnamed-chunk-6-2.png)<!-- -->![](website_traffic_analysis_files/figure-html/unnamed-chunk-6-3.png)<!-- -->![](website_traffic_analysis_files/figure-html/unnamed-chunk-6-4.png)<!-- -->

We also calculated a value of 12848 for the Deviance Information Criterion (DIC); this value would later be compared to the DIC value obtained for the alternative, normal model. 

In order to carry out residual analysis, we divided all the daily average session durations into 5 sublists, one for each of the 5 continents, and subtracted the mean average session duration $\mu_i$ for each continent $i$ from every element in the corresponding sublist for continent $i$. The residual plots are shown below.   

![](website_traffic_analysis_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

The mean as well as the standard deviation of the residuals for the 5 continents are listed below in ascending order of the continent indices. 


```
## Mean residuals: 6.968559 3.203226 -2.686607 19.38497 -21.31251
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
## 1.936160e+00 5.110010e+01 6.671620e+01 7.812275e+01 8.794898e+01 9.759376e+01 
##          30%          35%          40%          45%          50%          55% 
## 1.069973e+02 1.163594e+02 1.255360e+02 1.355689e+02 1.466733e+02 1.586251e+02 
##          60%          65%          70%          75%          80%          85% 
## 1.725716e+02 1.884885e+02 2.087711e+02 2.328303e+02 2.668738e+02 3.177312e+02 
##          90%          95%         100% 
## 4.110857e+02 6.616254e+02 5.381240e+12
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
## Mean residuals: 3.178891 0.4537799 0.2232925 4.830155 0.8981627
```

```
## Standard deviation of residuals: 327.9693 144.0456 127.6912 501.8171 99.26842
```

The posterior samples were again drawn to get Monte Carlo estimates of the mean average session duration. In comparison to the Weibull model, the range of values of the samples of the mean drawn for the normal model seemed to be more reasonable (lying between 0 and 9.051076e3). Using these samples, we simulated the distribution of the  average session duration (thus taking into account the uncertainty in the mean average session duration).  


```
##           0%           5%          10%          15%          20%          25% 
## 1.820272e-73 4.634810e-12 1.829416e-08 1.896221e-06 6.093720e-05 7.103372e-04 
##          30%          35%          40%          45%          50%          55% 
## 5.000559e-03 2.582620e-02 1.082508e-01 3.585969e-01 1.027618e+00 2.595177e+00 
##          60%          65%          70%          75%          80%          85% 
## 6.077616e+00 1.311816e+01 2.682163e+01 5.131927e+01 9.743540e+01 1.845555e+02 
##          90%          95%         100% 
## 3.590833e+02 7.751596e+02 8.471048e+03
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
## lam[1]   0.01225  0.002107 1.216e-05      4.221e-05
## lam[2]   0.03145  0.003259 1.882e-05      9.049e-05
## lam[3]   0.02491  0.002498 1.442e-05      7.267e-05
## lam[4]   0.01279  0.001873 1.081e-05      4.713e-05
## lam[5]   0.02451  0.006236 3.600e-05      7.761e-05
## mu[1]  270.81326 45.156914 2.607e-01      2.746e-01
## mu[2]   82.31808  5.142321 2.969e-02      2.995e-02
## mu[3]  109.76197  5.838832 3.371e-02      3.371e-02
## mu[4]  253.45985 29.266929 1.690e-01      1.728e-01
## mu[5]  121.02680 37.880893 2.187e-01      2.464e-01
## mu0      0.02225  0.006835 3.946e-05      1.668e-04
## sig0     0.01150  0.008509 4.913e-05      2.415e-04
## v        0.80760  0.018075 1.044e-04      5.536e-04
## 
## 2. Quantiles for each variable:
## 
##             2.5%       25%       50%       75%     97.5%
## lam[1] 8.583e-03 1.077e-02 1.209e-02   0.01356   0.01681
## lam[2] 2.555e-02 2.916e-02 3.128e-02   0.03357   0.03828
## lam[3] 2.035e-02 2.317e-02 2.477e-02   0.02655   0.03011
## lam[4] 9.513e-03 1.147e-02 1.267e-02   0.01398   0.01678
## lam[5] 1.408e-02 2.009e-02 2.387e-02   0.02829   0.03848
## mu[1]  1.959e+02 2.390e+02 2.660e+02 297.51580 372.27709
## mu[2]  7.276e+01 7.879e+01 8.214e+01  85.64489  92.85728
## mu[3]  9.893e+01 1.057e+02 1.096e+02 113.59632 121.68102
## mu[4]  2.024e+02 2.328e+02 2.514e+02 271.67450 317.29917
## mu[5]  6.584e+01 9.446e+01 1.148e+02 141.01399 213.95250
## mu0    1.280e-02 1.814e-02 2.122e-02   0.02491   0.03758
## sig0   4.750e-03 7.314e-03 9.511e-03   0.01298   0.02974
## v      7.719e-01 7.951e-01 8.078e-01   0.81973   0.84285
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
## mu[1]  274.60 26.693  0.15411        0.19305
## mu[2]   85.07  6.457  0.03728        0.04765
## mu[3]  106.85  5.052  0.02917        0.03716
## mu[4]  268.01 34.462  0.19897        0.25104
## mu[5]   98.82  9.146  0.05281        0.06942
## mu0    132.10 11.444  0.06607        0.08427
## sig[1] 195.27 11.286  0.06516        0.06516
## sig[2] 129.00  4.118  0.02377        0.02377
## sig[3] 117.50  3.275  0.01891        0.01891
## sig[4] 368.32 18.081  0.10439        0.10439
## sig[5]  33.90  2.272  0.01312        0.01340
## sig0   420.49 63.531  0.36680        0.47844
## 
## 2. Quantiles for each variable:
## 
##          2.5%    25%    50%    75%  97.5%
## mu[1]  222.09 256.52 274.79 292.65 327.02
## mu[2]   72.37  80.76  85.11  89.41  97.62
## mu[3]   96.94 103.48 106.84 110.24 116.73
## mu[4]  200.58 244.96 267.70 291.26 335.32
## mu[5]   80.72  92.71  98.83 105.00 116.62
## mu0    110.29 124.34 131.76 139.69 155.32
## sig[1] 174.55 187.42 194.81 202.50 218.95
## sig[2] 121.24 126.18 128.92 131.68 137.37
## sig[3] 111.31 115.25 117.42 119.68 124.12
## sig[4] 334.91 355.77 367.62 380.08 405.67
## sig[5]  29.85  32.29  33.78  35.36  38.69
## sig0   305.75 376.09 416.93 461.39 554.64
```

The mean values of the posterior average session durations for each of the 5 continents in the Weibull model appeared to be fairly similar to the corresponding values in the normal model. However, while the residual errors for the normal model were smaller than those for the Weibull model, the DIC value of 12484 for the Weibull model  was smaller than the DIC value of 14821 for the normal model, suggesting that the Weibull model may be preferable to the normal model. 

Finally, the second question of this project may be answered by looking at the mean value of the shape parameter $v$ as well as the mean value of the scale parameter $\lambda_i$ ($i \in \{1,2,3,4,5\}$) for each of the 5 continents. The hazard function of the Weibull distribution with shape parameter $v$ and scale parameter $\lambda$ is given by $h(t) = v\lambda_it^{v-1}$, and $h(t)$ is a positive decreasing function iff $v < 1$ and $\lambda_i > 0$. The plots below show that the hazard functions for all 5 distributions are indeed decreasing, thus confirming the phenomenon of negative ageing.

![](website_traffic_analysis_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

As briefly mentioned earlier, one caveat is that the Weibull model we fitted has quite a number of outliers. Furthermore, a limitation of the model is that the shape parameter $v$ follows a normal distribution with fixed parameters. It might be reasonable to allow the parameters of the distribution of $v$ to have distributions themselves, thus reflecting the variability of the mean values of $v$ for the different continents. Another possible improvement is to consider fitting a mixed model, for example a mixture of two Weibull distributions, which might reduce the residual errors.       

## References
