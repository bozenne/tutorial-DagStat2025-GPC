### covariateAdjustment.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 19 2025 (18:42) 
## Version: 
## Last-Updated: May  8 2025 (11:23) 
##           By: Brice Ozenne
##     Update #: 4
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

## ** standardization

## *** build-in
std.BR <- BuyseTest(treatment ~ tte(OS, statusOS, threshold = 6) + cont(toxicity, operator = "<0", threshold = 2) + tte(OS, statusOS, threshold = 2) + cont(toxicity, operator = "<0") + strata(sex),
                    pool.strata = "standardization", data = prodige)
summary(std.BR, digit = c(2,3,2))
## endpoint threshold strata total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  delta Delta CI [2.5% ; 97.5%] p.value    
##       OS         6 global   100.00        28.16          15.71      56.04     0.09  0.126 0.126     [0.054;0.197] 0.00067 ***
##                         M    24.47         6.88           3.42      14.14     0.04  0.142                                    
##                       F.M    27.31         7.42           4.90      14.97     0.02  0.092                                    
##                       M.F    22.79         6.66           3.00      13.11     0.02  0.160                                    
##                         F    25.43         7.20           4.39      13.83     0.01  0.111                                    
## toxicity         2 global    56.13        10.27          13.01      32.85     0.00 -0.028 0.098     [0.019;0.176] 0.01571   *
##                         M    14.18         2.76           3.14       8.29     0.00 -0.016                                    
##                       F.M    14.99         2.79           3.25       8.94     0.00 -0.017                                    
##                       M.F    13.13         2.34           3.24       7.55     0.00 -0.040                                    
##                         F    13.84         2.38           3.38       8.07     0.00 -0.039                                    
##       OS         2 global    32.85         5.08           5.02      22.71     0.04  0.001 0.099     [0.016;0.181] 0.01977   *
##                         M     8.29         1.11           1.23       5.92     0.02 -0.005                                    
##                       F.M     8.94         1.37           1.39       6.19     0.00 -0.001                                    
##                       M.F     7.55         1.22           1.12       5.18     0.02  0.004                                    
##                         F     8.07         1.38           1.28       5.41     0.00  0.004                                    
## toxicity           global    22.75         5.73           4.69      12.33     0.00  0.010 0.109     [0.022;0.194] 0.01385   *
##                         M     5.95         1.41           1.19       3.34     0.00  0.009                                    
##                       F.M     6.19         1.64           1.24       3.30     0.00  0.015                                    
##                       M.F     5.20         1.24           1.10       2.86     0.00  0.006                                    
##                         F     5.41         1.44           1.15       2.82     0.00  0.011                                    

coef(s.BR) ## marginal
##      OS_t6 toxicity_t2       OS_t2    toxicity 
## 0.12592885  0.09833291  0.09802122  0.10802346 
coef(std.BR)  ## stratified followed by standardization
##      OS_t6 toxicity_t2       OS_t2    toxicity 
## 0.12604879  0.09807780  0.09885171  0.10900857 

## *** manual
## WARNING: does not match the build-in approach due to a different censoring model (survival probabilities only treatment dependent but not sex dependent)
##          uncertainty not properly computed as the software mis-understand the number of independent observations (taken to be the sum of the weights instead of the number of lines in the dataset)
tableN <- addmargins(table(treatment = prodige$treatment,sex = prodige$sex))
df.tableN <- as.data.frame(tableN[1:2,1:2])
rownames(df.tableN) <- paste(df.tableN$treatment,df.tableN$sex, sep = ":")

prodige$weight <- sqrt(tableN["C","Sum"]*tableN["T","Sum"])/tableN["Sum","Sum"]*tableN["Sum",as.character(prodige$sex)]/df.tableN[paste(prodige$treatment,prodige$sex,sep=":"),"Freq"]
std2.BR <- BuyseTest(treatment ~ tte(OS, statusOS, threshold = 6) + cont(toxicity, operator = "<0", threshold = 2) + tte(OS, statusOS, threshold = 2) + cont(toxicity, operator = "<0"),
                    data = prodige, weightObs = prodige$weight, trace = FALSE)
coef(std2.BR, endpoint = "toxicity") 
## [1] 0.1119014

## to be compared to
coef(BuyseTest(treatment ~ tte(OS, statusOS, threshold = 6) + cont(toxicity, operator = "<0", threshold = 2) + tte(OS, statusOS, threshold = 2) + cont(toxicity, operator = "<0") + strata(sex),
               pool.strata = "standardization", data = prodige, model.tte =  prodlim(Hist(OS, statusOS)~treatment, data = prodige), trace = FALSE))
##     OS_t6 toxicity_t2       OS_t2    toxicity 
## 0.1282043   0.1008404   0.1016326   0.1119014 

## another approach would be to fit a GPC for each combination of strata among treatment groups: (active.control)=(M.M, F.M, M.F, F.F) and pool with appropriate weights

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




##----------------------------------------------------------------------
### covariateAdjustment.R ends here
