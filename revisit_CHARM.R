library(dplyr)
library(survival)
library(survminer)
library(BuyseTest)
BuyseTest.options(order.Hprojection=2) ## Use second order H-project

## * load data
data("CHARM", package = "BuyseTest")
head(CHARM)

## * Tabulate number of events
trt.arm <- 1 * (CHARM[,'treatment'] == 'T') 
cens.ind.all <- 1 * (CHARM[,'statusComposite'] == '1')
cens.mort <- 1 * (CHARM[,'statusMortality'] == '1')
cens.hosp <- 1 * (CHARM[,'statusHospitalization'] == '1')
cens.both <- 1 * (CHARM[,'statusHospitalization'] == '1' & CHARM[,'statusMortality'] == '1')
 
addmargins(table(cens.ind.all, trt.arm))
##             trt.arm
## cens.ind.all    0    1  Sum
##          0   1143 1181 2324
##          1    366  333  699
##          Sum 1509 1514 3023
table(cens.mort, trt.arm)
##          trt.arm
## cens.mort    0    1
##         0 1339 1344
##         1  170  170
table(cens.hosp, trt.arm)
##          trt.arm
## cens.hosp    0    1
##         0 1233 1273
##         1  276  241
table(cens.both, trt.arm)
##          trt.arm
## cens.both    0    1
##         0 1429 1436
##         1   80   78

## * time to first event analysis
km_fit <- survfit(Surv(Composite, statusComposite) ~ treatment, data = CHARM)

km_fit %>%
  ggsurvplot(
    data = CHARM,
    fun = "pct",
    # linetype = "strata", # Change line type by groups
    # pval = TRUE,
    # conf.int = TRUE,
    risk.table = TRUE,
    fontsize = 3, # used in risk table
    surv.median.line = "hv", # median horizontal and vertical ref lines
    ggtheme = theme_light(),
    title = "Kaplan-Meier Survival Function Estimate",
    legend.title = "",
    legend.labs = levels(data$treatment)
  )

survdiff( Surv(Composite, statusComposite) ~ treatment, data = CHARM)


## * Time to worse event analysis
BT_CHARM <- BuyseTest(treatment~tte(Mortality,statusMortality) + tte (Hospitalization,statusHospitalization), 
                      data=CHARM, scoring.rule = "Gehan", trace=0)
summary(BT_CHARM)
##       Generalized pairwise comparisons with 2 prioritized endpoints

## - statistic       : net treatment benefit  (delta: endpoint specific, Delta: global) 
## - null hypothesis : Delta == 0 
## - confidence level: 0.95 
## - inference       : H-projection of order 1 after atanh transformation 
## - treatment groups: T (treatment) vs. C (control) 
## - censored pairs  : deterministic score or uninformative
## - uninformative pairs: no contribution at the current endpoint, analyzed at later endpoints
## - results
##        endpoint total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  delta  Delta CI [2.5% ; 97.5%]  p.value  
##       Mortality   100.00         9.51           9.08          0    81.41 0.0042 0.0042  [-0.0157;0.0241] 0.676327  
## Hospitalization    81.41        10.58           7.94          0    62.90 0.0264 0.0306   [0.0029;0.0582] 0.030108 *

BT14_CHARM <- BuyseTest(treatment~tte(Mortality,statusMortality, threshold=14) +
                            tte(Hospitalization,statusHospitalization, threshold=14), 
                        data=CHARM, scoring.rule = "Gehan", trace=0)
summary(BT14_CHARM)
##       Generalized pairwise comparisons with 2 prioritized endpoints

## - statistic       : net treatment benefit  (delta: endpoint specific, Delta: global) 
## - null hypothesis : Delta == 0 
## - confidence level: 0.95 
## - inference       : H-projection of order 1 after atanh transformation 
## - treatment groups: T (treatment) vs. C (control) 
## - censored pairs  : deterministic score or uninformative
## - uninformative pairs: no contribution at the current endpoint, analyzed at later endpoints
## - results
##        endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  delta  Delta CI [2.5% ; 97.5%]  p.value  
##       Mortality        14   100.00         9.47           8.95       0.03    81.55 0.0052 0.0052   [-0.0147;0.025] 0.609892  
## Hospitalization        14    81.58        10.51           7.94       0.04    63.09 0.0257 0.0308   [0.0033;0.0584] 0.028366 *

## * Exit
BuyseTest.options(order.Hprojection=1) ## re-set default options
