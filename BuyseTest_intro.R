## * install BuyseTest package
## install.packages ("BuyseTest", quiet=TRUE)
library(BuyseTest)
## Default settings of the BuyseTest package
BuyseTest.options(order.Hprojection=2) ## request the second order U-statistic inference

## * create a dummy dataset to showcase the BuyseTest function
set.seed(10)
data <- simBuyseTest(100)

## * example of use with a single outcome
BuyseTest(treatment ~ bin(toxicity, operator = "<0"), data = data)

BuyseTest(treatment ~ cont(score, operator = ">0", threshold = 2), data = data)

BuyseTest(treatment ~ tte(eventtime, status = "status", threshold = 7), data = data)

## * example of use with a prioritized bivariate outcome
e.BT <- BuyseTest(treatment ~ tte(eventtime, status = "status") + bin(toxicity, operator = "<0"),
                  scoring.rule = "Gehan", data = data, trace = 0)


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

