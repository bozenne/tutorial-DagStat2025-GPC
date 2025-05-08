### descriptivePlot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 19 2025 (18:42) 
## Version: 
## Last-Updated: May  8 2025 (11:21) 
##           By: Brice Ozenne
##     Update #: 5
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(ggpubr)
library(riskRegression)
library(prodlim)
library(BuyseTest)
library(survival)

## * load data
data("prodige", package = "BuyseTest")
head(prodige)

## * slide 56: descriptive benefit-risk illustration

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
                             axis.line = element_line(linewidth = 1.25),
                             axis.ticks = element_line(linewidth = 2),
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
                               axis.line = element_line(linewidth = 1.25),
                               axis.ticks = element_line(linewidth = 2),
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


## * slide 73: Peron score rule example
## ** Peron scoring rule (example)
library(prodlim)
df.exPeron <- rbind(data.frame(time = c(3,5,7), event = 1, group = "E"),
                    data.frame(time = c(1,2,3,4,5), event = 1, group = "C"))
KM.exPeron <- prodlim(Hist(time,event) ~ group, data = df.exPeron)
plot(KM.exPeron, col = c("blue","orange"), confint = FALSE)


e.exPeron <- BuyseTest(group ~ tte(time, event),
                       model.tte = KM.exPeron, method.inference = "none", 
                       data = data.frame(time = c(4.5,6,1.5), event = c(1,1,0), group = c("E","E","C")),
                       trace = FALSE, keep.pairScore = TRUE)
getPairScore(e.exPeron)
##    index.C index.E favorable unfavorable neutral uninf weight
##      <num>   <num>     <num>       <num>   <num> <num>  <num>
## 1:       3       1      0.75        0.25       0     0      1
## 2:       3       2      1.00        0.00       0     0      1


##----------------------------------------------------------------------
### descriptivePlot.R ends here
