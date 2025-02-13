EB <- read.csv("../Data/EB.csv")
head(EB)

BuyseTest.options(order.Hprojection=2)

#prioritized
print(BuyseTest(Group~b(Bin)+c(DiffQoL), data=EB,method.inference="permutation",n.resampling=10000))

print(BuyseTest(Group~c(StdDiffCount)+c(DiffQoL), data=EB,method.inference="permutation",n.resampling=10000))

print(BuyseTest(Group~c(StdDiffCount, threshold=0.5)+c(DiffQoL), data=EB,method.inference="studentized permutation",n.resampling=10000))

print(BuyseTest(Group~b(Bin)+c(DiffQoL), data=EB,method.inference="u-statistic"), percentage=FALSE)

