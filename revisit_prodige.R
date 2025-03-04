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

chisq.test(table(prodige$toxicity,prodige$treatment))
## X-squared = 42.015, df = 5, p-value = 5.849e-08

## ** GPC alternative (different estimands)

e.NTB <- BuyseTest(treatment ~ tte(OS, statusOS), data = prodige)
summary(e.NTB)
 ## endpoint total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  Delta CI [2.5% ; 97.5%]  p.value  
 ##       OS      100        54.74          45.23          0     0.03 0.0951   [0.0084;0.1804] 0.031662 *

e.NTH <- BuyseTest(treatment ~ cont(toxicity, operator = "<0"), data = prodige)
summary(e.NTH)
 ## endpoint total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   Delta CI [2.5% ; 97.5%]    p.value    
 ## toxicity      100        33.68          47.96      18.36        0 -0.1427 [-0.2181;-0.0656] 0.00030318 ***

## * Benefit risk analysis
## ** main
e.BR <- BuyseTest(treatment ~ tte(OS, statusOS, threshold = 6) + cont(toxicity, operator = "<0", threshold = 2) + tte(OS, statusOS, threshold = 2) + cont(toxicity, operator = "<0"),
                  data = prodige)
plot(e.BR, type = "hist")
plot(e.BR, type = "racetrack")

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

## other measure of treatment effect
summary(e.BR, statistic = "winRatio")
coef(e.BR, statistic = "winRatio", cumulative = TRUE)
##    OS_t6 toxicity_t2       OS_t2    toxicity 
## 1.815742    1.593925    1.499402    1.392608 
coef(e.BR, statistic = "winRatio", cumulative = FALSE)
##     OS_t6 toxicity_t2       OS_t2    toxicity 
## 1.8157421   1.2595728   1.0131848   0.8173871 

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

## * Censoring
e.NTB_Gehan <- BuyseTest(treatment ~ tte(OS, statusOS), scoring.rule = "Gehan", data = prodige,
                         keep.pairScore = TRUE, trace = FALSE)
getPairScore(e.NTB_Gehan)[1:2,]

e.NTB_Peron <- BuyseTest(treatment ~ tte(OS, statusOS), scoring.rule = "Peron", data = prodige,
                         keep.pairScore = TRUE, trace = FALSE)
getPairScore(e.NTB_Peron)[1:2,]


e.NTB_Efron <- BuyseTest(treatment ~ tte(OS, statusOS), scoring.rule = "Efron", data = prodige, trace = FALSE)
e.NTB_Latta <- BuyseTest(treatment ~ tte(OS, statusOS), scoring.rule = "Peron", data = prodige, trace = FALSE,
                         model.tte = prodlim(Hist(OS, statusOS) ~ 1, data = prodige))

M.results <- rbind(Gehan = model.tables(e.NTB_Gehan),
                   Peron = model.tables(e.NTB_Peron),
                   Efron = model.tables(e.NTB_Efron),
                   Latta = model.tables(e.NTB_Latta))
print(M.results, digits = 3)
##              endpoint total favorable unfavorable  neutral  uninf  Delta lower.ci upper.ci p.value
## Gehan              OS   100      35.2        31.8 0.000591 33.090 0.0339 -0.03057   0.0981  0.3026
## Peron              OS   100      54.7        45.2 0.000967  0.025 0.0951  0.00838   0.1804  0.0317
## Efron              OS   100      54.7        45.2 0.025971  0.000 0.0951  0.00834   0.1805  0.0317
## Latta              OS   100      52.8        47.0 0.099400  0.049 0.0582 -0.01539   0.1313  0.1210

e.NTB_restricted <- BuyseTest(treatment ~ tte(OS, statusOS, restriction = 24), scoring.rule = "Peron", data = prodige, trace = FALSE)
model.tables(e.NTB_restricted)
##   endpoint restriction total favorable unfavorable   neutral uninf      Delta    lower.ci  upper.ci    p.value
## 1       OS          24   100  54.50976    45.04014 0.4501017     0 0.09469619 0.008047478 0.1799335 0.03224151

## * Extra: more precise comparison of proportions
library(exact2x2)

tableTox <- xtabs(cbind(n = 1, tox3 = toxicity>=3) ~ treatment, data = prodige)
binomMeld.test(n1 = tableTox["C","n"], n2 = tableTox["T","n"],
               x1 = tableTox["C","tox3"], x2 = tableTox["T","tox3"],
               conf.int = TRUE, parmtype = "difference")
## proportion 1 = 0.53234, proportion 2 = 0.6057, p-value = 0.03988
## difference (p2-p1) 
##          0.0733624 
## 95 percent confidence interval:
##  0.003272581 0.142782537

## exact test: uncondExact2x2pn

## * Extra: descriptive benefit-risk illustration
library(ggpubr)
library(riskRegression)
library(survival)

e.coxph <- coxph(Surv(OS,statusOS)~strata(treatment), data = prodige, x = TRUE, y = TRUE)
pred.coxph <- predictCox(e.coxph, type = "survival", keep.newdata = TRUE)
ggEx2_KM <- autoplot(pred.coxph, group.by = "strata", plot = FALSE)$plot
ggEx2_KM <- ggEx2_KM + scale_colour_manual(name = "",
                                           values = c("T" = "darkblue",
                                                      "C" = "darkorange"),
                                           labels = c("T" = "arm Folfirinox",
                                                      "C" = "arm Gemcitabine"))
ggEx2_KM <- ggEx2_KM + scale_shape_manual(name = "",
                                          breaks = c(0,1),
                                          values = c(3,18),
                                          labels = c("censoring","death"))
ggEx2_KM <- ggEx2_KM + labs(x = "Months", y = "Probability of survival")
ggEx2_KM <- ggEx2_KM + theme(text = element_text(size=20),
                             axis.line = element_line(size = 1.25),
                             axis.ticks = element_line(size = 2),
                             axis.ticks.length=unit(.25, "cm"),
                             legend.box.margin = margin(c(-40,0,0,0)),
                             legend.position="bottom",
                             legend.direction = "vertical",
                             legend.background = element_rect(fill = NA))
ggEx2_KM

## ** toxiciy: barplot
## proportion of each grade within each treatment arm
prodigeTtox <- as.data.frame(prop.table(table(prodige$treatment, prodige$toxicity), margin = 1))
names(prodigeTtox) <- c("treatment","grade","probability")

## color scale
scale_G2R <- scales::seq_gradient_pal(rgb(r=0,g=0.9,b=0), rgb(r=0.9,g=0,b=0),"Lab")(seq(0,1,length.out=6))

## graphical display
ggEx2_tox <- ggplot(prodigeTtox, aes(fill = grade, y = probability, x = treatment))
ggEx2_tox <- ggEx2_tox + geom_col(position = position_stack(reverse = TRUE))
ggEx2_tox <- ggEx2_tox + scale_fill_manual("Worst adverse event", values = scale_G2R)
ggEx2_tox <- ggEx2_tox + xlab("") + ylab("Relative frequency") 
ggEx2_tox <- ggEx2_tox + scale_y_continuous(labels = scales::percent) 
ggEx2_tox <- ggEx2_tox + scale_x_discrete(breaks = c("C","T"), labels = c("C" = "arm Gemcitabine", "T" = "arm Folfirinox")) 

ggEx2_tox <- ggEx2_tox + theme(text = element_text(size=20),
                               axis.text.x = element_text(colour = c("C" = "orange", "T" = "darkblue")),
                               axis.line = element_line(size = 1.25),
                               axis.ticks = element_line(size = 2),
                               axis.ticks.length=unit(.25, "cm"),
                               legend.box.margin = margin(c(-35,0,0,0)),
                               legend.position="bottom",
                               legend.direction = "vertical")
ggEx2_tox <- ggEx2_tox + guides(fill=guide_legend(nrow=1,byrow=TRUE))
ggEx2_tox

## ** assemble
ggEx2_KMtox <- ggarrange(ggEx2_KM, ggEx2_tox)
ggEx2_KMtox

## ggsave(ggEx2_KMtox, filename = file.path("figures","gg-Ex2-KMtox.png"), width = 10, height = 5, dpi = 150)
