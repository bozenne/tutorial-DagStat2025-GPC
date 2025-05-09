library(BuyseTest)
library(survival)
library(prodlim)
library(ggplot2) ## graphical illustration

## * load data
data("prodige", package = "BuyseTest")
head(prodige)


## * Marginal analyses
## ** Traditional analyses
## efficacy
e.HR <- coxph(Surv(OS, statusOS) ~ treatment, data = prodige)
cbind(estimate = exp(coef(e.HR)), exp(confint(e.HR)), p.value = summary(e.HR)$coef[,"Pr(>|z|)"])
##             estimate     2.5 %    97.5 %     p.value
## treatmentT 0.7783213 0.6582702 0.9202666 0.003366709

e.glmTox <- glm(toxicity>=3 ~ treatment, family = binomial(link = "identity"),  data = prodige)
cbind(estimate = coef(e.glmTox), confint(e.glmTox), p.value = summary(e.glmTox)$coef[,"Pr(>|z|)"])
## Waiting for profiling to be done...
##              estimate       2.5 %    97.5 %       p.value
## (Intercept) 0.5323383 0.483474459 0.5807906 1.603621e-101
## treatmentT  0.0733624 0.005721244 0.1406032  3.319334e-02

## * Benefit risk analysis
## ** main
e.BR <- BuyseTest(treatment ~ tte(OS, statusOS, threshold = 6) + cont(toxicity, operator = "<0", threshold = 2) + tte(OS, statusOS, threshold = 2) + cont(toxicity, operator = "<0"),
                  data = prodige)
plot(e.BR)
## plot(e.BR, type = "racetrack")
summary(e.BR)

table.BR <- model.tables(e.BR)
data.frame(favorable = round(table.BR$favorable,1),
           unfavorable = round(table.BR$unfavorable,1),
           neutral = round(table.BR$neutral + table.BR$uninf,1),
           delta = round(100*table.BR$delta,1),
           Delta = paste0(round(100*table.BR$Delta,1)," (p=",format.pval(table.BR$p.value, digits = 3, eps = 1e-5),")"))
##   favorable unfavorable neutral delta             Delta
## 1      28.2        15.6    56.2  12.7 12.7 (p=0.000401)
## 2      10.3        13.0    32.9  -2.7   10 (p=0.011969)
## 3       5.1         5.0    22.8   0.1 10.1 (p=0.015047)
## 4       5.7         4.7    12.3   1.0 11.1 (p=0.010297)
sum(table.BR$favorable)       ## [1] 49.3886
sum(table.BR$unfavorable)      ## [1] 38.26632

## ** decomposition
e.BRfav <- coef(e.BR, statistic = "favorable", cumulative = FALSE)
e.BRfav
##      OS_t6 toxicity_t2       OS_t2    toxicity 
## 0.28238768  0.10317657  0.05095573  0.05736603 
e.BRunfav <- coef(e.BR, statistic = "unfavorable", cumulative = FALSE)
e.BRunfav
##      OS_t6 toxicity_t2       OS_t2    toxicity 
## 0.15552191  0.12995841  0.05029263  0.04689025 
cumsum(e.BRfav - e.BRunfav) ## match the Net Treatment Benefit
##     OS_t6 toxicity_t2       OS_t2    toxicity 
## 0.1268658   0.1000839   0.1007470   0.1112228 

coef(e.BR, statistic = "winRatio")
##    OS_t6 toxicity_t2       OS_t2    toxicity 
## 1.815742    1.350581    1.300045    1.290655 
cumprod(e.BRfav/e.BRunfav) ## do not match the Win Ratio
##    OS_t6 toxicity_t2       OS_t2    toxicity 
## 1.815742    1.441554    1.460560    1.786865 
cumsum(e.BRfav/e.BRunfav)  ## do not match the Win Ratio
##    OS_t6 toxicity_t2       OS_t2    toxicity 
## 1.815742    2.609662    3.622847    4.846257 


## ** sensitivity analysis
## requires the riskRegression package to be install for multiple testing adjustment
M.threshold <- cbind(OS_t6 = c(3:5,3:5),
                     toxicity_t2 = c(2,2,2,3,3,3),
                     OS_t2 = 1,
                     toxicity = 0)
M.threshold
eBR.Se <- sensitivity(e.BR, band = TRUE,
                      threshold = M.threshold)
plot(eBR.Se)
eBR.Se

