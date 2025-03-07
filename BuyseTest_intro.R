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
BuyseTest.options(reinitialise = TRUE)

## * designing trial
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

e.power <- powerBuyseTest(formula = group ~ tte(time, event, threshold = 1) + cont(tox, operator = "<0"),
                          sim = simFCT, sample.size = c(10,25,50),
                          n.rep = 100, seed = 10, cpus = 1)
summary(e.power)

e.n <- powerBuyseTest(formula = group ~ tte(time,event, threshold = 1) + cont(tox, operator = "<0"),
                      sim = simFCT, power = 0.8,
                      n.rep = c(1000,10), seed = 10, trace = 2, cpus = 1)
summary(e.n)

e.n2 <- powerBuyseTest(group ~ tte(time,event, threshold = 1) + cont(tox, operator = "<0"),
                       sim = simFCT, power = 0.8, max.sample.size = 5000,
                       n.rep = c(1000,10), seed = 10, trace = 2, cpus = 1)
summary(e.n2)
