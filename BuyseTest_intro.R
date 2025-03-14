library(BuyseTest)

## * Software for GPC, part 1: measures of treatment effect
## ** create a dataset to showcase the BuyseTest function
set.seed(10)
data <- simBuyseTest(100, n.strata = 2)
head(data)
##       id treatment  eventtime status toxicity      score strata
##    <num>    <fctr>      <num>  <num>   <fctr>      <num> <fctr>
## 1:     1         C 0.17392093      1      yes -2.1250686      a
## 2:     2         C 0.16255166      0      yes  0.5211787      a
## 3:     3         C 0.08302502      1      yes -0.0464229      b
## 4:     4         C 0.22204972      0       no -1.1494717      b
## 5:     5         C 0.11669726      1       no  0.6293383      a
## 6:     6         C 0.11885540      1      yes -0.7264715      a

## ** Formula interface + Output of the GPC procedure
e.BT <- BuyseTest(treatment ~ tte(eventtime, status),
                  data = data)
summary(e.BT)
##       Generalized pairwise comparisons with 1 endpoint
##
## - statistic       : net treatment benefit  (delta: endpoint specific, Delta: global) 
## - null hypothesis : Delta == 0 
## - confidence level: 0.95 
## - inference       : H-projection of order 1 after atanh transformation 
## - treatment groups: T (treatment) vs. C (control) 
## - censored pairs  : probabilistic score based on the survival curves
## - results
##  endpoint total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  Delta CI [2.5% ; 97.5%] p.value 
## eventtime      100        57.39          42.61          0        0 0.1479  [-0.0293;0.3161] 0.10151 

summary(e.BT, percentage = FALSE) ## number of pairs
## [...]
##  endpoint total favorable unfavorable neutral uninf  Delta CI [2.5% ; 97.5%] p.value 
## eventtime 10000   5739.39     4260.61       0     0 0.1479  [-0.0293;0.3161] 0.10151 

## ** Other measures of treatment effect
confint(e.BT, statistic = "winRatio") ## win ratio
##           estimate        se  lower.ci upper.ci null   p.value
## eventtime 1.347081 0.2450411 0.9430953 1.924118    1 0.1014458

e.BThalf <- BuyseTest(treatment ~ tte(eventtime, status),
                      data = data, add.halfNeutral = TRUE, trace = FALSE)
model.tables(e.BThalf, statistic = "favorable") ## probabilistic index
##    endpoint total favorable unfavorable neutral uninf     Delta  lower.ci  upper.ci   p.value
## 1 eventtime   100  57.39388    42.60612       0     0 0.5739388 0.4852354 0.6581263 0.1019135
coef(e.BThalf, statistic = "winRatio") ## win odds
## [1] 1.347081

## ** More options
## example of use with a prioritized bivariate outcome
e.MBT <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 1) + bin(toxicity, operator = "<0"),
                   data = data, trace = 0)
model.tables(e.MBT)
print(model.tables(e.MBT),digits = 3)
##    endpoint threshold total favorable unfavorable neutral uninf   delta  Delta lower.ci upper.ci p.value
## 1 eventtime     1e+00 100.0      10.2        2.55    87.2     0  0.0768 0.0768 -0.00928    0.162  0.0803
## 3  toxicity     1e-12  87.2      18.8       24.72    43.7     0 -0.0590 0.0178 -0.13396    0.169  0.8192

## graphical display
plot(e.MBT)
plot(e.MBT, type = "racetrack")

## * Software for GPC, part 2: statistical inference
## ** Inference based on U-statistic theory
## toy example
df10 <- rbind(data.frame(Y = c(-1.2,-0.5,-0.8,0.3,1.1,1.2,0.7,-0.5,0.6,-1.2), group = "C"),
              data.frame(Y = c(-0.6,-2.2,-0.7,-2.1,-1.3,-0.4,-0.7,-0.9,0.1,-0.3), group = "E"))

## first order H-projection
e.Ustat <- BuyseTest(group ~ cont(Y), data = df10, trace = 0)
coef(e.Ustat, statistic = "favorable") ## [1] 0.26
coef(e.Ustat, statistic = "unfavorable") ## [1] 0.74
coef(e.Ustat, statistic = "netBenefit") ## [1] -0.48

eIid.Ustat <- getIid(e.Ustat)*10
(eIid.Ustat[1:10])^2 ## [1] 0.7744 0.0064 0.4624 0.2704 0.2704 0.2704 0.2704 0.0064 0.2704 0.7744
(eIid.Ustat[11:20])^2 ## [1] 0.0064 0.2704 0.0064 0.2704 0.2704 0.2304 0.0064 0.0144 0.2304 0.2304

crossprod(eIid.Ustat) ## [1] 4.912
confint(e.Ustat)[,"se"]^2 ## [1] 0.04912

## second order H-projection
BuyseTest.options(order.Hprojection = 2)
e.Ustat2 <- BuyseTest(group ~ cont(Y), data = df10, trace = 0)
confint(e.Ustat2)[,"se"]^2 ## [1] 0.051904
BuyseTest.options(order.Hprojection = 1)

## ** Inference based on U-statistic theory (code)
rbind(confint(e.BT, transformation = TRUE),
      confint(e.BT, transformation = FALSE))
##             estimate         se    lower.ci  upper.ci null    p.value
## eventtime  0.1478776 0.08897931 -0.02931684 0.3160612    0 0.10150573
## eventtime1 0.1478776 0.08897931 -0.02651861 0.3222739    0 0.09652625

NTB <- coef(e.BT)
sigma.NTB <- sqrt(crossprod(getIid(e.BT)))
sigmaTrans.NTB <- sigma.NTB/(1-NTB^2)
c(estimate = NTB, se = sigmaTrans.NTB,
  p.value = 2*(1-pnorm(NTB/sigma.NTB)),
  pTrans.value = 2*(1-pnorm(atanh(NTB)/sigmaTrans.NTB)))
##   estimate           se      p.value pTrans.value 
## 0.14787764   0.09096860   0.09652625   0.10150573 

## ** Equivalence GPC and Wilcoxon rank sum test
eBT.perm <- BuyseTest(treatment ~ cont(score), data = data,
                      method.inference = "varexact permutation")
model.tables(eBT.perm)
##   endpoint total favorable unfavorable neutral uninf  Delta   p.value
## 1    score   100     53.67       46.33       0     0 0.0734 0.3698664

wilcox.test(score ~ treatment, data = data, correct = FALSE)$p.value
## [1] 0.3698664

## ** Argument method.inference in Buyse Test

## change default behavior
BuyseTest.options(method.inference = "permutation", n.resampling = 1000, statistic = "winRatio")
BuyseTest.options(reinitialise = TRUE)


## * Examples revisited
## see files revisit_CHARM.R
##           revisit_EB.R
##           revisit_prodige.R

## * Advanced Topics, part 1: censoring

## ** Peron scoring rule (example)
library(prodlim)
df.exPeron <- rbind(data.frame(time = c(3,5,7), event = 1, group = "T"),
                    data.frame(time = c(1,2,3,4,5), event = 1, group = "C"))
KM.exPeron <- prodlim(Hist(time,event) ~ group, data = df.exPeron)
plot(KM.exPeron)

e.exPeron <- BuyseTest(group ~ tte(time, event),
                       model.tte = KM.exPeron, method.inference = "none", 
                       data = data.frame(time = c(4.5,6,1.5), event = c(1,1,0), group = c("T","T","C")),
                       trace = FALSE, keep.pairScore = TRUE)
getPairScore(e.exPeron)
##    index.C index.T favorable unfavorable neutral uninf weight
##      <num>   <num>     <num>       <num>   <num> <num>  <num>
## 1:       3       1      0.75        0.25       0     0      1
## 2:       3       2      1.00        0.00       0     0      1

## ** Using build-in imputation approach
data(prodige, package = "BuyseTest")
e.NTB_Gehan <- BuyseTest(treatment ~ tte(OS, statusOS), scoring.rule = "Gehan", data = prodige,
                         keep.pairScore = TRUE, trace = FALSE)
getPairScore(e.NTB_Gehan)[1:2,]
##    index.C index.T favorable unfavorable neutral uninf weight
##      <num>   <num>     <num>       <num>   <num> <num>  <num>
## 1:       1     403         1           0       0     0      1
## 2:       2     403         0           0       0     1      1

e.NTB_Peron <- BuyseTest(treatment ~ tte(OS, statusOS), scoring.rule = "Peron", data = prodige,
                         keep.pairScore = TRUE, trace = FALSE)
getPairScore(e.NTB_Peron)[1:2,]
##    index.C index.T favorable unfavorable neutral        uninf weight
##      <num>   <num>     <num>       <num>   <num>        <num>  <num>
## 1:       1     403 1.0000000     0.00000       0 0.0000000000      1
## 2:       2     403 0.5286551     0.47068       0 0.0006648516      1

## ** Using your own imputation approach
e.NTB_Latta <- BuyseTest(treatment ~ tte(OS, statusOS), scoring.rule = "Peron",
                         data = prodige, trace = FALSE,
                         model.tte = prodlim(Hist(OS, statusOS) ~ 1, data = prodige))

e.NTB_Efron <- BuyseTest(treatment ~ tte(OS, statusOS), scoring.rule = "Efron", data = prodige, trace = FALSE)

library(survival)
e.NTB_Weibull <- BuyseTest(treatment ~ tte(OS, statusOS),
                           data = prodige, trace = FALSE,
                           model.tte = survreg(Surv(OS, statusOS) ~ treatment, data = prodige, dist = "weibull"))

M.results <- rbind(Gehan = model.tables(e.NTB_Gehan),
                   Peron = model.tables(e.NTB_Peron),
                   Efron = model.tables(e.NTB_Efron),
                   Latta = model.tables(e.NTB_Latta),
                   Weibull = model.tables(e.NTB_Weibull))
print(M.results, digits = 3)
##         endpoint total favorable unfavorable  neutral    uninf  Delta lower.ci upper.ci p.value
## Gehan         OS   100      35.2        31.8 0.000591 33.09049 0.0339 -0.03057   0.0981  0.3026
## Peron         OS   100      54.7        45.2 0.000967  0.02500 0.0951  0.00838   0.1804  0.0317
## Efron         OS   100      54.7        45.2 0.025971  0.00000 0.0951  0.00834   0.1805  0.0317
## Latta         OS   100      52.8        47.0 0.099400  0.04895 0.0582 -0.01539   0.1313  0.1210
## Weibull       OS   100      50.9        43.3 5.855714  0.00413 0.0758 -0.00830   0.1589  0.0772

## ** What about administrative censoring?

e.NTB_restricted <- BuyseTest(treatment ~ tte(OS, statusOS, restriction = 24), scoring.rule = "Peron", data = prodige, trace = FALSE)
model.tables(e.NTB_restricted)
##   endpoint restriction total favorable unfavorable   neutral uninf      Delta    lower.ci  upper.ci    p.value
## 1       OS          24   100  54.50976    45.04014 0.4501017     0 0.09469619 0.008047478 0.1799335 0.03224151


## * Advanced Topics, part 2: covariate adjustment

## see file revisit_prodige.R


## * Designing trial

## ** Power/sample size calculation with BuyseTest
simFCT <- function(n.C, n.T){
    df.C <- data.frame(id = paste0("C",1:n.C), group = 0,
                       tox = sample(1:6, n.C, replace=TRUE, prob= c(16.09, 15.42, 33.26, 26.18, 8.38, 0.67)/100),                       
                       time = rweibull(n.C, scale = 9.995655, shape = 1.28993),
                       event = 1)
    df.T <- data.frame(id = paste0("T",1:n.T), group = 1,
                       tox = sample(1:6, n.T, replace=TRUE, prob= c(8.21, 13.09, 31.29, 30.87, 12.05, 4.49)/100),                       
                       time = rweibull(n.T, scale = 13.16543, shape = 1.575269),
                       event = 1)
    return(rbind(df.C,df.T))
}
set.seed(10)
simFCT(2,2)
##   id group tox      time event
## 1 C1     0   4  8.821945     1
## 2 C2     0   3  4.591318     1
## 3 T1     1   3 15.495787     1
## 4 T2     1   3 15.557655     1

## ** Power calculation with BuyseTest
e.power <- powerBuyseTest(formula = group ~ tte(time, event, threshold = 1) + cont(tox, operator = "<0"),
                          sim = simFCT, sample.size = c(10,25,50),
                          n.rep = 100, seed = 10, cpus = 1)
summary(e.power)
##         Simulation study with Generalized pairwise comparison
##         with 100 samples

##  - net benefit statistic (null hypothesis Delta=0)
##  endpoint threshold n.T n.C mean.estimate sd.estimate mean.se rejection.rate
##       tox     1e-12  10  10        0.2156      0.2656  0.2468           0.13
##                      25  25        0.2032      0.1677  0.1582            0.2
##                      50  50        0.2015      0.1228  0.1121           0.43

##  n.T          : number of observations in the treatment group
##  n.C          : number of observations in the control group
##  mean.estimate: average estimate over simulations
##  sd.estimate  : standard deviation of the estimate over simulations
##  mean.se      : average estimated standard error of the estimate over simulations
##  rejection    : frequency of the rejection of the null hypothesis over simulations
## (standard error: H-projection of order 1| p-value: after transformation) 

## more precise estimation of the power by simulating more datasets
## use parallel calculation to reduce the of waiting time
e.power2 <- powerBuyseTest(formula = group ~ tte(time, event, threshold = 1) + cont(tox, operator = "<0"),
                          sim = simFCT, sample.size = c(10,25,50),
                          n.rep = 1000, seed = 10, cpus = 4)
summary(e.power2)
##         Simulation study with Generalized pairwise comparison
##         with 1000 samples

##  - net benefit statistic (null hypothesis Delta=0)
##  endpoint threshold n.T n.C mean.estimate sd.estimate mean.se rejection.rate
##       tox     1e-12  10  10        0.2031       0.253  0.2486          0.103
##                      25  25         0.202      0.1568  0.1586          0.219
##                      50  50         0.202      0.1132  0.1122          0.413

##  n.T          : number of observations in the treatment group
##  n.C          : number of observations in the control group
##  mean.estimate: average estimate over simulations
##  sd.estimate  : standard deviation of the estimate over simulations
##  mean.se      : average estimated standard error of the estimate over simulations
##  rejection    : frequency of the rejection of the null hypothesis over simulations
## (standard error: H-projection of order 1| p-value: after transformation) 


## ** Sample size calculation with BuyseTest
e.n <- powerBuyseTest(formula = group ~ tte(time,event, threshold = 1) + cont(tox, operator = "<0"),
                      sim = simFCT, power = 0.8,
                      n.rep = c(1000,10), seed = 10, trace = 2, cpus = 1)
summary(e.n)
##          Sample size calculation with Generalized pairwise comparison
##         for a power of 0.8 and type 1 error rate of 0.05 

##  - estimated sample size (mean [min;max]): 126 [91;155] controls
##                                            126 [91;155] treated

##  - net benefit statistic (null hypothesis Delta=0)
##  endpoint threshold n.T n.C mean.estimate sd.estimate mean.se rejection.rate
##       tox     1e-12 126 126        0.2049       0.069  0.0707          0.818

##  n.T          : number of observations in the treatment group
##  n.C          : number of observations in the control group
##  mean.estimate: average estimate over simulations
##  sd.estimate  : standard deviation of the estimate over simulations
##  mean.se      : average estimated standard error of the estimate over simulations
##  rejection    : frequency of the rejection of the null hypothesis over simulations
## (standard error: H-projection of order 1| p-value: after transformation) 


## more accurate sample size calculation using larger datasets to estimate the asymptotic variance
e.n2 <- powerBuyseTest(group ~ tte(time,event, threshold = 1) + cont(tox, operator = "<0"),
                       sim = simFCT, power = 0.8, max.sample.size = 5000,
                       n.rep = c(1000,10), seed = 10, trace = 2, cpus = 1)
summary(e.n2)
##         Sample size calculation with Generalized pairwise comparison
##         for a power of 0.8 and type 1 error rate of 0.05 

##  - estimated sample size (mean [min;max]): 118 [100;138] controls
##                                            118 [100;138] treated

##  - net benefit statistic (null hypothesis Delta=0)
##  endpoint threshold n.T n.C mean.estimate sd.estimate mean.se rejection.rate
##       tox     1e-12 118 118         0.203      0.0743  0.0731          0.776

##  n.T          : number of observations in the treatment group
##  n.C          : number of observations in the control group
##  mean.estimate: average estimate over simulations
##  sd.estimate  : standard deviation of the estimate over simulations
##  mean.se      : average estimated standard error of the estimate over simulations
##  rejection    : frequency of the rejection of the null hypothesis over simulations
## (standard error: H-projection of order 1| p-value: after transformation) 
