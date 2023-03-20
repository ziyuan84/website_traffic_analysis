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
bibliography: bayesianbib.bib
---



## Abstract

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
## Mean residuals: 6.783184 3.125954 -2.701038 20.14128 -21.52738
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
## 6.144106e+00 5.183820e+01 6.691371e+01 7.833263e+01 8.835362e+01 9.810520e+01 
##          30%          35%          40%          45%          50%          55% 
## 1.073620e+02 1.164208e+02 1.263255e+02 1.361756e+02 1.469466e+02 1.591016e+02 
##          60%          65%          70%          75%          80%          85% 
## 1.727186e+02 1.883973e+02 2.079350e+02 2.328437e+02 2.661226e+02 3.180218e+02 
##          90%          95%         100% 
## 4.073108e+02 6.685608e+02 6.755625e+11
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
## Mean residuals: 3.011816 0.5432441 0.2039191 5.201533 0.7201194
```

```
## Standard deviation of residuals: 327.9693 144.0456 127.6912 501.8171 99.26842
```

The posterior samples were again drawn to get Monte Carlo estimates of the mean average session duration. In comparison to the Weibull model, the range of values of the samples of the mean drawn for the normal model seemed to be more reasonable (lying between 0 and 9.051076e3). Using these samples, we simulated the distribution of the  average session duration (thus taking into account the uncertainty in the mean average session duration).  


```
##           0%           5%          10%          15%          20%          25% 
## 5.632670e-66 4.007364e-12 1.996959e-08 2.181921e-06 5.230344e-05 6.644559e-04 
##          30%          35%          40%          45%          50%          55% 
## 4.626948e-03 2.327692e-02 1.029172e-01 3.602710e-01 1.060681e+00 2.578360e+00 
##          60%          65%          70%          75%          80%          85% 
## 6.049204e+00 1.306816e+01 2.691566e+01 5.197088e+01 9.829807e+01 1.849997e+02 
##          90%          95%         100% 
## 3.571570e+02 7.690441e+02 9.018937e+03
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
## lam[1]   0.01231  0.002136 1.233e-05      4.581e-05
## lam[2]   0.03156  0.003275 1.891e-05      9.141e-05
## lam[3]   0.02502  0.002514 1.451e-05      7.262e-05
## lam[4]   0.01289  0.001878 1.084e-05      4.639e-05
## lam[5]   0.02459  0.006247 3.607e-05      7.992e-05
## mu[1]  270.99863 45.646220 2.635e-01      2.795e-01
## mu[2]   82.39535  5.135042 2.965e-02      3.065e-02
## mu[3]  109.77641  5.856358 3.381e-02      3.381e-02
## mu[4]  252.70355 29.166230 1.684e-01      1.746e-01
## mu[5]  121.24166 38.303982 2.211e-01      2.515e-01
## mu0      0.02216  0.006035 3.484e-05      1.182e-04
## sig0     0.01122  0.006369 3.677e-05      1.425e-04
## v        0.80673  0.018155 1.048e-04      5.548e-04
## 
## 2. Quantiles for each variable:
## 
##             2.5%       25%       50%       75%     97.5%
## lam[1] 8.549e-03   0.01081 1.217e-02   0.01365   0.01688
## lam[2] 2.557e-02   0.02928 3.142e-02   0.03367   0.03844
## lam[3] 2.043e-02   0.02326 2.493e-02   0.02665   0.03023
## lam[4] 9.516e-03   0.01157 1.280e-02   0.01408   0.01690
## lam[5] 1.421e-02   0.02013 2.397e-02   0.02833   0.03846
## mu[1]  1.956e+02 238.63535 2.663e+02 298.54108 373.25200
## mu[2]  7.284e+01  78.85702 8.222e+01  85.73504  92.90414
## mu[3]  9.900e+01 105.70178 1.095e+02 113.58189 121.85616
## mu[4]  2.016e+02 232.37749 2.505e+02 270.46702 316.03191
## mu[5]  6.617e+01  94.44466 1.148e+02 140.37619 214.11537
## mu0    1.292e-02   0.01819 2.136e-02   0.02507   0.03630
## sig0   4.798e-03   0.00732 9.511e-03   0.01302   0.02787
## v      7.715e-01   0.79453 8.065e-01   0.81894   0.84254
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
## mu[1]  274.77 26.543  0.15325        0.19994
## mu[2]   84.98  6.497  0.03751        0.04734
## mu[3]  106.87  5.041  0.02910        0.03670
## mu[4]  267.64 34.461  0.19896        0.25557
## mu[5]   98.99  9.119  0.05265        0.06824
## mu0    132.12 11.441  0.06606        0.08586
## sig[1] 195.18 11.256  0.06499        0.06398
## sig[2] 128.98  4.076  0.02353        0.02310
## sig[3] 117.47  3.248  0.01875        0.01861
## sig[4] 368.36 17.971  0.10376        0.10376
## sig[5]  33.91  2.274  0.01313        0.01332
## sig0   421.02 63.251  0.36518        0.48166
## 
## 2. Quantiles for each variable:
## 
##          2.5%    25%    50%    75%  97.5%
## mu[1]  222.50 257.08 274.90 292.52 326.82
## mu[2]   72.10  80.65  85.05  89.35  97.51
## mu[3]   97.09 103.47 106.83 110.21 116.92
## mu[4]  199.53 244.52 267.65 290.95 335.16
## mu[5]   81.07  92.90  99.01 105.13 116.92
## mu0    110.54 124.30 131.80 139.61 155.63
## sig[1] 174.82 187.38 194.61 202.28 218.93
## sig[2] 121.22 126.20 128.90 131.65 137.28
## sig[3] 111.31 115.26 117.37 119.58 124.08
## sig[4] 334.84 355.98 367.77 379.87 405.71
## sig[5]  29.82  32.32  33.78  35.37  38.73
## sig0   306.68 376.93 417.83 462.19 553.36
```

The mean values of the posterior average session durations for each of the 5 continents in the Weibull model appeared to be fairly similar to the corresponding values in the normal model. However, while the residual errors for the normal model were smaller than those for the Weibull model, the DIC value of 12484 for the Weibull model  was smaller than the DIC value of 14821 for the normal model, suggesting that the Weibull model may be preferable to the normal model. 

Finally, the second question of this project may be answered by looking at the mean value of the shape parameter $v$ as well as the mean value of the scale parameter $\lambda_i$ ($i \in \{1,2,3,4,5\}$) for each of the 5 continents. The hazard function of the Weibull distribution with shape parameter $v$ and scale parameter $\lambda$ is given by $h(t) = v\lambda_it^{v-1}$, and $h(t)$ is a positive decreasing function iff $v < 1$ and $\lambda_i > 0$. The plots below show that the hazard functions for all 5 distributions are indeed decreasing, thus confirming the phenomenon of negative ageing.

![](website_traffic_analysis_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

As briefly mentioned earlier, one caveat is that the Weibull model we fitted has quite a number of outliers. Furthermore, a limitation of the model is that the shape parameter $v$ follows a normal distribution with fixed parameters. It might be reasonable to allow the parameters of the distribution of $v$ to have distributions themselves, thus reflecting the variability of the mean values of $v$ for the different continents. Another possible improvement is to consider fitting a mixed model, for example a mixture of two Weibull distributions, which might reduce the residual errors.       

## References
