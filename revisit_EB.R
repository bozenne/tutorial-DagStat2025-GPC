library(BuyseTest)

## * Default settings of the BuyseTest package
BuyseTest.options(order.Hprojection=2) ## request the second order U-statistic inference

## * load data
data("EB")
head(EB)


## * prioritized
print(BuyseTest(Group~b(Bin)+c(DiffQoL), data=EB,method.inference="permutation",n.resampling=10000))

print(BuyseTest(Group~c(StdDiffCount)+c(DiffQoL), data=EB,method.inference="permutation",n.resampling=10000))

print(BuyseTest(Group~c(StdDiffCount, threshold=0.5)+c(DiffQoL), data=EB,method.inference="studentized permutation",n.resampling=10000))

print(BuyseTest(Group~b(Bin)+c(DiffQoL), data=EB,method.inference="u-statistic"), percentage=FALSE)

