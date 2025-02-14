library(dplyr)
library(survival)
library(survminer)
library(BuyseTest)

## * load data
data("CHARM")
head(CHARM)

## * Tabulate number of events
trt.arm <- 1 * (CHARM[,'treatment'] == 'T') 
cens.ind.all <- 1 * (CHARM[,'statusComposite'] == '1')
cens.mort <- 1 * (CHARM[,'statusMortality'] == '1')
cens.hosp <- 1 * (CHARM[,'statusHospitalization'] == '1')
cens.both <- 1 * (CHARM[,'statusHospitalization'] == '1' & CHARM[,'statusMortality'] == '1')
 
table(cens.ind.all, trt.arm)
table(cens.mort, trt.arm)
table(cens.hosp, trt.arm)
table(cens.both, trt.arm)

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

BT_CHARM <- BuyseTest(treatment~tte(Mortality,statusMortality, threshold=14) +
                        tte (Hospitalization,statusHospitalization, threshold=14), hierarchical=FALSE, 
                      data=CHARM, scoring.rule = "Gehan", trace=0)
summary(BT_CHARM)
