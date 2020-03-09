# Bayes factors for Testing Variances
This R package allows one to test equality, inequality, and mixed hypotheses on independent population variances. The method is described in Dablander, F.<sup>&#8902;</sup>, van den Bergh, D.<sup>&#8902;</sup>, Ly, A., Wagenmakers, E.-J. (2020) Default Bayes Factors for Testing the (In)equality of Several Population Variances.


## Overview
```{r}
devtools::install_github('fdabl/bfvartest')
library('bfvartest')

# 2.6.1 Testing Against a Single Value
x <- c(6.2, 5.8, 5.7, 6.3, 5.9, 5.8, 6.0)
onesd_test(
    n = length(x), s = sd(x), popsd = sqrt(0.10),
    alpha = 2.16, alternative_interval = c(1, Inf)
)


# 2.6.2 Comparing Measurement Precision
n <- 990
sdigit <- 0.98
slaser <- 0.89
twosd_test(
    n1 = n, n2 = n, sd1 = sdigit, sd2 = slaser,
    alpha = 0.50, alternative_interval = c(1, Inf)
) # H_+ vs H_0

twosd_test(
    n1 = n, n2 = n, sd1 = sdigit, sd2 = slaser, alpha = 0.50,
    null_interval = c(0.90, 1.10), alternative_interval = c(1.10, Inf)
) # H'_1 vs H'_0


# 2.6.3 Sex Differences in Personality
twosd_test(n1 = 969, n2 = 716, sd1 = 15.6, sd2 = 20, alpha = 4.50)


# 3.1.1 The "Standardization" Hypothesis in Archeology
ns <- c(55, 171, 117)
sds <- c(5.83, 8.13, 12.74)
hyp <- c('1=2=3', '1>2>3')
res <- ksd_test(hyp = hyp, ns = ns, sds = sds, alpha = 0.50)
res$BF


# 3.1.2 Increased Variability in Mathematical Ability
ns <- c(3280, 6007, 7749, 9160, 9395, 6410)
sds <- c(5.983, 5.385, 4.970, 4.626, 3.688, 3.082)
hyp <- c('1=2=3=4=5=6', '1,2,3,4,5,6', '1>2>3>4>5>6')
res <- ksd_test(hyp = hyp, ns = ns, sds = sds, alpha = 0.50)
res$BF
```
