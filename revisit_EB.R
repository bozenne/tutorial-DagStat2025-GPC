library(BuyseTest)

## * Default settings of the BuyseTest package
BuyseTest.options(order.Hprojection=2) ## request the second order U-statistic inference

## * load data
data("EB", package = "BuyseTest")
head(EB)


## * prioritized
print(BuyseTest(Group~b(Bin)+c(DiffQoL), data=EB,method.inference="varexact-permutation"))
## endpoint total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  delta  Delta   p.value
##      Bin   100.00           44          10.67      45.33     0.00 0.3333 0.3333 0.0701057
##  DiffQoL    45.33           32           6.22       5.33     1.78 0.2578 0.5911 0.0051302

print(BuyseTest(Group~c(StdDiffCount)+c(DiffQoL), data=EB,method.inference="varexact-permutation"))
##     endpoint total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  delta  Delta  p.value
## StdDiffCount   100.00        64.00          27.11       2.22     6.67 0.3689 0.3689 0.069625
##      DiffQoL     8.89         6.22           0.44       1.78     0.44 0.0578 0.4267 0.040017

print(BuyseTest(Group~c(StdDiffCount, threshold=0.2)+c(DiffQoL), data=EB,method.inference="varexact-permutation"))
##     endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  delta  Delta  p.value
## StdDiffCount       0.2      100        56.44          19.56      17.33     6.67 0.3689 0.3689 0.059927
##      DiffQoL                 24        15.11           2.67       3.56     2.67 0.1244 0.4933 0.016440

print(BuyseTest(Group~b(Bin)+c(DiffQoL), data=EB,method.inference="u-statistic"), percentage=FALSE)
## endpoint total favorable unfavorable neutral uninf  delta  Delta CI [2.5% ; 97.5%]   p.value
##      Bin   225        99          24     102     0 0.3333 0.3333  [-0.0291;0.6183] 0.0706270
##  DiffQoL   102        72          14      12     4 0.2578 0.5911   [0.1931;0.8221] 0.0059238

## * Exit
BuyseTest.options(order.Hprojection=1) ## re-set default options
