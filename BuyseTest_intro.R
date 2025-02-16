library(BuyseTest)

## * create a dummy dataset to showcase the BuyseTest function
set.seed(10)
data <- simBuyseTest(100, n.strata = 2)
head(data)

## * example of use with a single outcome
e.BT <- BuyseTest(treatment ~ tte(eventtime, status),
                  data = data)
summary(e.BT)
summary(e.BT, percentage = FALSE) ## number of pairs

## other measures of treatment effect
confint(e.BT, statistic = "winRatio") ## win ratio
e.BThalf <- BuyseTest(treatment ~ tte(eventtime, status),
                      data = data, add.halfNeutral = TRUE, trace = FALSE)
model.tables(e.BThalf, statistic = "favorable") ## probabilistic index
coef(e.BThalf, statistic = "winRatio") ## win odds

## other types of outcomes
BuyseTest(treatment ~ cont(score, threshold = 2), data = data)
BuyseTest(treatment ~ bin(toxicity, operator = "<0")data = data)



## * example of use with a prioritized bivariate outcome
e.MBT <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 1) + bin(toxicity, operator = "<0"),
                   data = data, trace = 0)
model.tables(e.MBT)

print(model.tables(e.MBT),digits = 3)

plot(e.MBT)
plot(e.MBT, type = "racetrack")

## * Obtain score-matrix
e.BTindiv <- BuyseTest(treatment ~ tte(eventtime, status = status),
                   data = data, keep.pairScore = TRUE)
getPairScore(e.BTindiv)


## * Matching wilcoxon test
eBT.perm <- BuyseTest(treatment ~ cont(score), data = data,
                      method.inference = "varexact permutation")
model.tables(eBT.perm)

wilcox.test(score ~ treatment, data = data, correct = FALSE)$p.value

## * No tanh transformation
rbind(confint(e.BTindiv, transformation = TRUE),
      confint(e.BTindiv, transformation = FALSE))

NTB <- coef(e.BTindiv)
sigma.NTB <- sqrt(crossprod(getIid(e.BTindiv)))
sigmaTrans.NTB <- sigma.NTB/(1-NTB^2)
c(estimate = NTB, se = sigmaTrans.NTB,
  p.value = 2*(1-pnorm(NTB/sigma.NTB)),
  pTrans.value = 2*(1-pnorm(atanh(NTB)/sigmaTrans.NTB)))
##   estimate           se      p.value pTrans.value 
## 0.14787764   0.09096860   0.09652625   0.10150573 

## * update default behavior
BuyseTest.options(method.inference = "permutation", n.resampling = 1000, statistic = "winRatio")
