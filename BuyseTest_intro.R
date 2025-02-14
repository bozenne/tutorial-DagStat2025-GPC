library(BuyseTest)

## * create a dummy dataset to showcase the BuyseTest function
set.seed(10)
data <- simBuyseTest(100, n.strata = 2)
head(data)

## * example of use with a single outcome
e.BT <- BuyseTest(treatment ~ tte(eventtime, status), data = data)
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


## ** Net Treatment Benefit
summary(e.BT)
model.tables(e.BT, percentage = TRUE)
model.tables(e.BT, percentage = FALSE)
confint(e.BT)
## graphical display
plot(e.BT, type = "hist")
plot(e.BT, type = "pie")
plot(e.BT, type = "racetrack")

## ** Win ratio
confint(e.BT, statistic="winRatio")

## ** Probabilistic index (PI) and success odds (SO)
e.BT2 <- BuyseTest(treatment ~ tte(eventtime, status = "status") + bin(toxicity, operator = "<0"),
                   scoring.rule = "Gehan", data = data, trace = 0,
                   add.halfNeutral = TRUE)
confint(BT2, statistic = "favorable") ## PI
confint(BT2, statistic = "winRatio") ## SO

## * Obtain score-matrix
e.BT3 <- BuyseTest(treatment ~ tte(eventtime, status = "status") + bin(toxicity, operator = "<0"),
                   scoring.rule = "Gehan", data = data, trace=0,
                   keep.pairScore = TRUE)
getPairScore(e.BT3)

