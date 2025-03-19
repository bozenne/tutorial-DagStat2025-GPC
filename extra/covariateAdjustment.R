### covariateAdjustment.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 19 2025 (18:42) 
## Version: 
## Last-Updated: mar 19 2025 (18:42) 
##           By: Brice Ozenne
##     Update #: 1
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(BuyseTest)
library(survival)
library(prodlim)
library(ggplot2) ## graphical illustration

## * load data
data("prodige", package = "BuyseTest")
head(prodige)

## * Advanced Topics, part 2: covariate adjustment

## ** Revisit Prodige trial with sex stratification
s.BR <- BuyseTest(treatment ~ tte(OS, statusOS, threshold = 6) + cont(toxicity, operator = "<0", threshold = 2) + tte(OS, statusOS, threshold = 2) + cont(toxicity, operator = "<0") + strata(sex),
                  pool.strata = "CMH", data = prodige)
summary(s.BR, digit = c(2,3,2)) ## digit: number of decimal for percentage, delta, and p-value
## endpoint threshold strata total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  delta Delta CI [2.5% ; 97.5%] p.value    
##       OS         6 global   100.00        28.22          15.64      56.04     0.10  0.126 0.126     [0.056;0.195] 0.00048 ***
##                         M    49.04        13.79           6.84      28.33     0.08  0.142                                    
##                         F    50.96        14.43           8.80      27.71     0.02  0.111                                    
## toxicity         2 global    56.14        10.30          13.07      32.77     0.00 -0.028 0.098      [0.02;0.176] 0.01410   *
##                         M    28.41         5.52           6.28      16.60     0.00 -0.016                                    
##                         F    27.73         4.77           6.78      16.17     0.00 -0.039                                    
##       OS         2 global    32.77         4.99           5.02      22.72     0.05  0.000 0.098     [0.016;0.178] 0.01856   *
##                         M    16.60         2.22           2.46      11.87     0.05 -0.005                                    
##                         F    16.17         2.77           2.56      10.85     0.00  0.004                                    
## toxicity           global    22.76         5.71           4.71      12.35     0.00  0.010 0.108     [0.023;0.192] 0.01312   *
##                         M    11.92         2.83           2.39       6.70     0.00  0.009                                    
##                         F    10.85         2.88           2.31       5.65     0.00  0.011                                    

nobs(s.BR)
##   C     T pairs 
## 402   421 84456 
nobs(s.BR, strata = TRUE)
##     C   T pairs
## M 190 218 41420
## F 212 203 43036
confint(s.BR, strata = TRUE)
##                 estimate         se     lower.ci  upper.ci null     p.value
## OS_t6.M       0.14163397 0.04738957  0.047772987 0.2330156    0 0.003192326
## OS_t6.F       0.11055413 0.05319843  0.005450657 0.2132417    0 0.039286663
## toxicity_t2.M 0.12608947 0.05325181  0.020703383 0.2287044    0 0.019148860
## toxicity_t2.F 0.07116027 0.05904555 -0.045004977 0.1854268    0 0.229711152
## OS_t2.M       0.12124870 0.05565770  0.011132977 0.2284589    0 0.031001199
## OS_t2.F       0.07528238 0.06110224 -0.044985322 0.1933990    0 0.219667852
## toxicity.M    0.13014476 0.05847593  0.014300729 0.2425413    0 0.027778835
## toxicity.F    0.08636752 0.06352251 -0.038834727 0.2088999    0 0.176098720

## ** Non-parametric adjustment
marginal.BR <- BuyseTest(treatment ~ tte(OS, statusOS, threshold = 6) + cont(toxicity, operator = "<0", threshold = 2) + tte(OS, statusOS, threshold = 2) + cont(toxicity, operator = "<0"),
                         data = prodige, trace = FALSE)
strata.BR <- BuyseTest(treatment ~ bin(sex), data = prodige, trace = FALSE)
## diff(tapply(prodige$sex == "F", prodige$treatment, mean)) - coef(BT.strata)

Vmat2 <- crossprod(cbind(getIid(marginal.BR, endpoint = "toxicity"),getIid(strata.BR)))
## point estimate
coef(marginal.BR, endpoint = "toxicity") -  Vmat2[1,-1,drop=FALSE] %*% solve(Vmat2[-1,-1,drop=FALSE]) %*% coef(strata.BR) ## 0.1122007
## variance
Vmat2[1,1] - Vmat2[1,-1,drop=FALSE] %*% solve(Vmat2[-1,-1,drop=FALSE]) %*% Vmat2[-1,1,drop=FALSE] ## 0.001847724

## ** Alternative: standardization
tableN <- addmargins(table(treatment = prodige$treatment,sex = prodige$sex))
df.tableN <- as.data.frame(tableN[1:2,1:2])
rownames(df.tableN) <- paste(df.tableN$treatment,df.tableN$sex, sep = ":")

prodige$weight <- sqrt(tableN["C","Sum"]*tableN["T","Sum"])/tableN["Sum","Sum"]*tableN["Sum",as.character(prodige$sex)]/df.tableN[paste(prodige$treatment,prodige$sex,sep=":"),"Freq"]
std.BR <- BuyseTest(treatment ~ tte(OS, statusOS, threshold = 6) + cont(toxicity, operator = "<0", threshold = 2) + tte(OS, statusOS, threshold = 2) + cont(toxicity, operator = "<0"),
                    data = prodige, weightObs = prodige$weight, trace = FALSE)
confint(std.BR, endpoint = "toxicity") ## WARNING 1: uncertainty ignore variability of the weights, i.e., assume known covariate distribution
                                       ## WARNING 2: survival probabilities only treatment dependent but not sex dependent
##           estimate        se   lower.ci  upper.ci null   p.value
## toxicity 0.1119014 0.0431152 0.02678982 0.1954023    0 0.0100622

#### EXTRA: equivalence between explicit standardization and standardization using weights            ####
####        require same survival model - illustrated only with survival only depending on treatment  ####

## approach 1 with weights
e.ref <- BuyseTest(treatment ~ tte(OS, statusOS, threshold = 6) + cont(toxicity, operator = "<0", threshold = 2),
                   model.tte = prodlim(Hist(OS,statusOS)~treatment, data = prodige), data = prodige, trace = FALSE)

confint(e.ref, endpoint = "toxicity_t2")
##              estimate         se   lower.ci  upper.ci null   p.value
## toxicity_t2 0.1000839 0.03955895 0.02209791 0.1768593    0 0.0119687

## approach 2 explicit standardization (omitting the argument model.tte is prefered as survival probabilities would be treatment and sex dependent)
prodige$sex2 <- ifelse(prodige$treatment=="T",prodige$sex=="M",prodige$sex=="F")
e.strata1 <- BuyseTest(treatment ~ tte(OS, statusOS, threshold = 6) + cont(toxicity, operator = "<0", threshold = 2) + strata(sex),
                       model.tte = prodlim(Hist(OS,statusOS)~treatment, data = prodige), data = prodige, trace = FALSE)
e.strata2 <- BuyseTest(treatment ~ tte(OS, statusOS, threshold = 6) + cont(toxicity, operator = "<0", threshold = 2) + strata(sex2),
                       model.tte = prodlim(Hist(OS,statusOS)~treatment, data = prodige), data = prodige, trace = FALSE)

eDelta.strata1 <- coef(e.strata1, strata = TRUE, endpoint = "toxicity_t2")
eDelta.strata2 <- coef(e.strata2, strata = TRUE, endpoint = "toxicity_t2")

eTable.std <- rbind(data.frame(sex.C = names(eDelta.strata1),
                               pi.C = tableN["C",names(eDelta.strata1)]/tableN["C","Sum"],
                               sex.T = names(eDelta.strata1),
                               pi.T = tableN["T",names(eDelta.strata1)]/tableN["T","Sum"],
                               Delta = eDelta.strata1),
                    data.frame(sex.C = c("M","F"),
                               pi.C = tableN["C",c("M","F")]/tableN["C","Sum"],
                               sex.T = c("F","M"),
                               pi.T = tableN["T",c("F","M")]/tableN["T","Sum"],
                               Delta = eDelta.strata2)
                    )
sum(eTable.std$pi.C*eTable.std$pi.T*eTable.std$Delta)
## [1] 0.1000839




##----------------------------------------------------------------------
### covariateAdjustment.R ends here
