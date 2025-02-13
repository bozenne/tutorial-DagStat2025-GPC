library(dplyr)
library(survival)
library(survminer)


charm <- read.csv("../Data/CHARM_sim.csv")
head(charm)

# Tabulate number of events
trt.arm <- 1 * (charm[,'treatment'] == 'T') 
cens.ind.all <- 1 * (charm[,'statusComposite'] == '1')
cens.mort <- 1 * (charm[,'statusMortality'] == '1')
cens.hosp <- 1 * (charm[,'statusHospitalization'] == '1')
cens.both <- 1 * (charm[,'statusHospitalization'] == '1' & charm[,'statusMortality'] == '1')
 
table(cens.ind.all, trt.arm)
table(cens.mort, trt.arm)
table(cens.hosp, trt.arm)
table(cens.both, trt.arm)

#time to first event analysis
km_fit <- survfit(Surv(Composite, statusComposite) ~ treatment, data = charm)

km_fit %>%
  ggsurvplot(
    data = charm,
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

survdiff( Surv(Composite, statusComposite) ~ treatment, data = charm)


#Time to worse event analysis
BT_charm <- BuyseTest(treatment~tte(Mortality,statusMortality) + tte (Hospitalization,statusHospitalization), 
                      data=charm, scoring.rule = "Gehan", trace=0)
summary(BT_charm)

BT_charm <- BuyseTest(treatment~tte(Mortality,statusMortality, threshold=14) +
                        tte (Hospitalization,statusHospitalization, threshold=14), hierarchical=FALSE, 
                      data=charm, scoring.rule = "Gehan", trace=0)
summary(BT_charm)
