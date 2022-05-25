# Bayes factors for Comparing Variances
[![Build Status](https://travis-ci.org/fdabl/bfvartest.svg?branch=master)](https://travis-ci.org/fdabl/bfvartest)

This R package allows one to test equality, inequality, and mixed hypotheses on independent population variances. The method is described in Dablander, F.<sup>&#11089;</sup>, van den Bergh, D.<sup>&#11089;</sup>, Wagenmakers, E.-J., & Ly, A. ([2022](https://arxiv.org/abs/2003.06278)) Default Bayes Factors for Testing the (In)equality of Several Population Variances.

Specifically, the package allows testing hypotheses of the form:

<p align='center' style='padding-top: 1em;'>
  <img src='Variances-Math.png' width=300/>
</p>

## Overview
Below you find the code to reproduce the analyses in the paper.

```r
devtools::install_github('fdabl/bfvartest', build_vignettes = TRUE)
library('bfvartest')

# 5.1 Sex Differences in Personality
twosd_test(n1 = 969, n2 = 716, sd1 = sqrt(15.6), sd2 = sqrt(19.9), u = 0.50)


# 5.2 Testing Against a Single Value
x <- c(6.2, 5.8, 5.7, 6.3, 5.9, 5.8, 6.0)
n <- length(x)
sd_x <- sd(x) # use rounded 0.22 in the paper

## (i) BF_{+0}
onesd_test(
    n = n, s = sd_x, popsd = sqrt(0.10),
    u = 0.50, alternative_interval = c(1, Inf), log = FALSE
)

## (ii) BF_{10}
onesd_test(
    n = n, s = sd_x, popsd = sqrt(0.10),
    u = 0.50, alternative_interval = c(0, Inf), log = FALSE
)

## (iii) BF_{+0} informed
onesd_test(
    n = n, s = sd_x, popsd = sqrt(0.10),
    u = 2.16, alternative_interval = c(1, Inf), log = FALSE
)


# 5.3 Comparing Measurement Precision
n <- 990
sdigit <- 0.98
slaser <- 0.89

## (i) BF_{+0}
twosd_test(
    n1 = n, n2 = n, sd1 = slaser, sd2 = sdigit,
    u = 0.50, alternative_interval = c(1, Inf), log = FALSE
)


## (ii) BF'_{0+} non-overlapping interval
1 / twosd_test(
    n1 = n, n2 = n, sd1 = slaser, sd2 = sdigit, u = 0.50, log = FALSE,
    null_interval = c(0.90, 1.10), alternative_interval = c(1.10, Inf)
)



# 5.4 The "Standardization" Hypothesis in Archeology
ns <- c(117, 171, 55)
sds <- c(12.74, 8.13, 5.83)
hyp <- c('1=2=3', '1,2,3', '1>2>3')
res <- ksd_test(hyp = hyp, ns = ns, sds = sds, u = 0.50, iter = 6000)
res$BF

## (i) log BF_{10} ~ 20
res$BF[2, 1]

## (ii) log BF_{r0} ~ 21.80
res$BF[3, 1]

## (iii) log BF_{r1} ~ 1.79
res$BF[3, 2]


# 5.5 Increased Variability in Mathematical Ability
ns <- c(3280, 6007, 7549, 9160, 9395, 6410)
sds <- c(5.99, 5.39, 4.97, 4.62, 3.69, 3.08)
hyp <- c('1=2=3=4=5=6', '1,2,3,4,5,6', '1>2>3>4>5>6')
res <- ksd_test(hyp = hyp, ns = ns, sds = sds, u = 0.50, iter = 6000)
res$BF

## (i) log BF_{10} ~ 1660.53
res$BF[2, 1]

## (ii) log BF_{r0} ~ 1667.11
res$BF[3, 1]

## (iii) log BF_{r1} ~ 6.58
res$BF[3, 2]


# For more, browse the vignette
browseVignettes('bfvartest')
```
