# June 25th 2020
# Hector Bahamonde. "Still for Sale: The Micro-Dynamics of Vote Selling in the United States, Evidence From a List Experiment". Acta Politica.

###############################################
# Figure 3
###############################################

## ---- barplot:data:control:treatment ----
# Load Data
load("dat_list.RData") # Load data

# Histogram for Item Count

## Create a factor variable to use in the plot
dat$treatment.f = factor(dat$treatment, levels = c(0,1), labels=c("Control", "Treatment"))


treat.cont.bar.plot.d = data.frame(
        ycount = c(
                c(table(factor(dat$treatment100, levels = c(0,1,2,3,4), labels=c("0", "1", "2", "3", "4")))),
                c(table(factor(dat$treatment500, levels = c(0,1,2,3,4), labels=c("0", "1", "2", "3", "4")))),
                c(table(factor(dat$control, levels = c(0,1,2,3), labels=c("0", "1", "2", "3"))))
        ),
        Condition = c(
                rep("Low Treatment ($ 100)", 5),
                rep("High Treatment ($ 500)", 5),
                rep("Control", 4)
        ),
        xcount = c(
                0,1,2,3,4,#
                0,1,2,3,4,#
                0,1,2,3#
        ),
        n.items = c(
                rep("NonSens", 4), "Sens",
                rep("NonSens", 4), "Sens",
                rep("NonSens", 4)
        )
)

# reorder factor var to have low treatment BEFORE high treatment.
treat.cont.bar.plot.d$Condition = factor(treat.cont.bar.plot.d$Condition,levels(treat.cont.bar.plot.d$Condition)[c(1,3,2)])


# calculate percentages
## plyr is dicontinued/retired as it January 2020. Hence, I'll be installing from source
# install.packages("https://cran.r-project.org/src/contrib/plyr_1.8.5.tar.gz", repos=NULL, type="source")
library(plyr)

treat.cont.bar.plot.d = ddply(treat.cont.bar.plot.d, .(Condition), transform, percent = ycount/sum(ycount) * 100)
treat.cont.bar.plot.d = ddply(treat.cont.bar.plot.d, .(Condition), transform, pos = (cumsum(ycount) - 0.5 * ycount))
treat.cont.bar.plot.d$label = paste0(sprintf("%.0f", treat.cont.bar.plot.d$percent), "%")


# Plot
if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(ggplot2)

barplot.descriptive.plot = 
        ggplot(treat.cont.bar.plot.d, 
               aes(x=xcount, 
                   y = ycount,
                   fill=n.items)) + 
        geom_bar(stat = "identity") + 
        facet_grid(.~ Condition) + 
        xlab("Number of Items") + 
        ylab("Frequency") +
        theme_bw() +
        scale_fill_manual(values=c("forestgreen", "forestgreen")) +
        theme(axis.text.y = element_text(size=7), 
              axis.text.x = element_text(size=7), 
              axis.title.y = element_text(size=7), 
              axis.title.x = element_text(size=7), 
              legend.text=element_text(size=7), 
              legend.title=element_text(size=7),
              plot.title = element_text(size=7),
              strip.text.x = element_text(size = 7),
              legend.position="none") +
        geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2)
## ----



##############################
# difference in means test
##############################

## ---- ttests

# low
mean.t.low.d = c(
        # 0
        rep(treat.cont.bar.plot.d$xcount[treat.cont.bar.plot.d$xcount==0 & treat.cont.bar.plot.d$Condition=="Low Treatment ($ 100)"],treat.cont.bar.plot.d$ycount[treat.cont.bar.plot.d$xcount==0 & treat.cont.bar.plot.d$Condition=="Low Treatment ($ 100)"]), 
        # 1
        rep(treat.cont.bar.plot.d$xcount[treat.cont.bar.plot.d$xcount==1 & treat.cont.bar.plot.d$Condition=="Low Treatment ($ 100)"],treat.cont.bar.plot.d$ycount[treat.cont.bar.plot.d$xcount==1 & treat.cont.bar.plot.d$Condition=="Low Treatment ($ 100)"]),
        # 2
        rep(treat.cont.bar.plot.d$xcount[treat.cont.bar.plot.d$xcount==2 & treat.cont.bar.plot.d$Condition=="Low Treatment ($ 100)"],treat.cont.bar.plot.d$ycount[treat.cont.bar.plot.d$xcount==2 & treat.cont.bar.plot.d$Condition=="Low Treatment ($ 100)"]),
        # 3
        rep(treat.cont.bar.plot.d$xcount[treat.cont.bar.plot.d$xcount==3 & treat.cont.bar.plot.d$Condition=="Low Treatment ($ 100)"],treat.cont.bar.plot.d$ycount[treat.cont.bar.plot.d$xcount==3 & treat.cont.bar.plot.d$Condition=="Low Treatment ($ 100)"]),
        # 4
        rep(treat.cont.bar.plot.d$xcount[treat.cont.bar.plot.d$xcount==4 & treat.cont.bar.plot.d$Condition=="Low Treatment ($ 100)"],treat.cont.bar.plot.d$ycount[treat.cont.bar.plot.d$xcount==4 & treat.cont.bar.plot.d$Condition=="Low Treatment ($ 100)"])
)

mean.t.low = round(mean(mean.t.low.d),3)


## high 
mean.t.high.d = c(
        # 0
        rep(treat.cont.bar.plot.d$xcount[treat.cont.bar.plot.d$xcount==0 & treat.cont.bar.plot.d$Condition=="High Treatment ($ 500)"],treat.cont.bar.plot.d$ycount[treat.cont.bar.plot.d$xcount==0 & treat.cont.bar.plot.d$Condition=="High Treatment ($ 500)"]), 
        # 1
        rep(treat.cont.bar.plot.d$xcount[treat.cont.bar.plot.d$xcount==1 & treat.cont.bar.plot.d$Condition=="High Treatment ($ 500)"],treat.cont.bar.plot.d$ycount[treat.cont.bar.plot.d$xcount==1 & treat.cont.bar.plot.d$Condition=="High Treatment ($ 500)"]),
        # 2
        rep(treat.cont.bar.plot.d$xcount[treat.cont.bar.plot.d$xcount==2 & treat.cont.bar.plot.d$Condition=="High Treatment ($ 500)"],treat.cont.bar.plot.d$ycount[treat.cont.bar.plot.d$xcount==2 & treat.cont.bar.plot.d$Condition=="High Treatment ($ 500)"]),
        # 3
        rep(treat.cont.bar.plot.d$xcount[treat.cont.bar.plot.d$xcount==3 & treat.cont.bar.plot.d$Condition=="High Treatment ($ 500)"],treat.cont.bar.plot.d$ycount[treat.cont.bar.plot.d$xcount==3 & treat.cont.bar.plot.d$Condition=="High Treatment ($ 500)"]),
        # 4
        rep(treat.cont.bar.plot.d$xcount[treat.cont.bar.plot.d$xcount==4 & treat.cont.bar.plot.d$Condition=="High Treatment ($ 500)"],treat.cont.bar.plot.d$ycount[treat.cont.bar.plot.d$xcount==4 & treat.cont.bar.plot.d$Condition=="High Treatment ($ 500)"])
)

mean.t.high = round(mean(mean.t.high.d),3)

# control
mean.control.d = c(
        # 0
        rep(treat.cont.bar.plot.d$xcount[treat.cont.bar.plot.d$xcount==0 & treat.cont.bar.plot.d$Condition=="Control"],treat.cont.bar.plot.d$ycount[treat.cont.bar.plot.d$xcount==0 & treat.cont.bar.plot.d$Condition=="Control"]), 
        # 1
        rep(treat.cont.bar.plot.d$xcount[treat.cont.bar.plot.d$xcount==1 & treat.cont.bar.plot.d$Condition=="Control"],treat.cont.bar.plot.d$ycount[treat.cont.bar.plot.d$xcount==1 & treat.cont.bar.plot.d$Condition=="Control"]),
        # 2
        rep(treat.cont.bar.plot.d$xcount[treat.cont.bar.plot.d$xcount==2 & treat.cont.bar.plot.d$Condition=="Control"],treat.cont.bar.plot.d$ycount[treat.cont.bar.plot.d$xcount==2 & treat.cont.bar.plot.d$Condition=="Control"]),
        # 3
        rep(treat.cont.bar.plot.d$xcount[treat.cont.bar.plot.d$xcount==3 & treat.cont.bar.plot.d$Condition=="Control"],treat.cont.bar.plot.d$ycount[treat.cont.bar.plot.d$xcount==3 & treat.cont.bar.plot.d$Condition=="Control"])
)

mean.control = round(mean(mean.control.d),3)

# t-tests

## diff in means
diff.mean.high.control = round(as.numeric(mean(mean.t.high.d) - mean(mean.control.d)), 3)
diff.mean.low.control = round(as.numeric(mean(mean.t.low.d) - mean(mean.control.d)), 3)

## p-values
high.control.t.test.pvalue = round(as.numeric(t.test(mean.t.high.d,mean.control.d)$'p.value'),3)
low.control.t.test.pvalue = round(as.numeric(t.test(mean.t.low.d,mean.control.d)$'p.value'),3)

## t-statistic
high.control.t.test.t = round(as.numeric(t.test(mean.t.high.d,mean.control.d)$'statistic'),3)
low.control.t.test.t = round(as.numeric(t.test(mean.t.low.d,mean.control.d)$'statistic'),3)

## df
high.control.t.test.df = round(as.numeric(t.test(mean.t.high.d,mean.control.d)$'parameter'),0)
low.control.t.test.df = round(as.numeric(t.test(mean.t.low.d,mean.control.d)$'parameter'),0)

## ----




## ---- barplot:figure:control:treatment
# use this to explain plot in the paper
barplot.descriptive.plot
barplot.descriptive.plot.note <- paste(
        "{\\bf Frequency and Percentages of Subjects Declaring How Many (if any) Illegal Things They Would Do}.",
        "\\\\\\hspace{\\textwidth}", 
        "{\\bf Note}: Notice that the X-axis denotes the number of items, not which ones. Percentages show proportions per condition.",
        "\n")
## ----




######################################################
# Difference in Means
######################################################

## ---- diff:in:means ----

# customization of ictreg
options(scipen=999)
method = as.character("ml")
maxIter = as.numeric(200000)
options(digits=2)
ci.level  = 0.95



# Loading the Data
load("dat_list.RData") # Load data

# Constructing DFs

## subseting leaving only LOW/HIGH treatment condition
dat.low <- dat[ which(dat$treatlow==1 | dat$treatment==0), ] 
dat.high <- dat[ which(dat$treatlow==0 | dat$treatment==0), ] 

## Function to drop missing data
completeFun <- function(data, desiredCols) {
        completeVec <- complete.cases(data[, desiredCols])
        return(data[completeVec, ])
        }

## dropping rows with missing data
dat.low = completeFun(dat.low, c("socideo", "partyid","income.n","educ.n"))

## dropping rows with missing data
dat.high = completeFun(dat.high, c("socideo", "partyid","income.n","educ.n"))

# Total Sample Size
total.sample.size = as.character(formatC(c(nrow(dat.high) + nrow(dat.low)), format="d", big.mark=","))


# Difference in means 
if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(list)

# estimating diff in means low treatment
dif.means.low <- ictreg(ycount ~ 1, data = dat.low, treat = "treatment", J=3, method = "lm")
dif.means.low.predicted = predict(dif.means.low, avg = T, interval = "confidence")

# estimating diff in means high treatment
dif.means.high <- ictreg(ycount ~ 1, data = dat.high, treat = "treatment", J=3, method = "lm")
dif.means.high.predicted = predict(dif.means.high,avg = T, interval = "confidence")

# estimaging diff in means control
dif.means.control <- ictreg(ycount ~ 1, data = dat, treat = "treatment", J=3, method = "lm")
dif.means.control.high <- ictreg(ycount ~ 1, data = dat.high, treat = "treatment", J=3, method = "lm")
dif.means.control.low <- ictreg(ycount ~ 1, data = dat.low, treat = "treatment", J=3, method = "lm")

options(digits=4)
dif.means.control.mean = as.numeric(dif.means.control$par.control, length=5)
dif.means.control.se = as.numeric(dif.means.control$se.control, length=5)
dif.means.control.upr = dif.means.control.mean + 1.96*dif.means.control.se
dif.means.control.lrw = dif.means.control.mean - 1.96*dif.means.control.se

# difference in means low treatment
dif.in.means.low.mean = as.numeric(dif.means.low.predicted$fit[1], length=5) # estimated mean
dif.in.means.low.lwr = as.numeric(dif.means.low.predicted$fit[2], length=5) # lower bound
dif.in.means.low.upr = as.numeric(dif.means.low.predicted$fit[3], length=5) # upper bound

# difference in means high treatment
dif.in.means.high.mean = as.numeric(dif.means.high.predicted$fit[1], length=5) # estimated mean
dif.in.means.high.lwr = as.numeric(dif.means.high.predicted$fit[2], length=5) # lower bound
dif.in.means.high.upr = as.numeric(dif.means.high.predicted$fit[3], length=5) # upper bound

# construct DF
diff.means.plot.d = data.frame(
        c.1.3 = c("Difference\nin Means\n(without covariates)", "Difference\nin Means\n(without covariates)"),
        fit = c(dif.in.means.low.mean, dif.in.means.high.mean),
        upr = c(dif.in.means.low.upr,dif.in.means.high.upr),
        lwr = c(dif.in.means.low.lwr,dif.in.means.high.lwr),
        Condition = c("Low ($100)", "High ($500)")
)
## ---- 


###########################################################
# FIGURE 4
###########################################################

## ---- list:analysis:covariates ----
 
load("dat_list.RData") # Load data

##############
# List Low Condition
##############
if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(list)

list.low <- ictreg(ycount ~ socideo + partyid + income.n + educ.n,
                   data = dat.low, 
                   treat = "treatment", 
                   J=3, 
                   method = method)


## this is to construct a table later
n.draws = 200000
coeffs.treat.list.low = as.data.frame(summary(list.low, n.draws = n.draws)["par.treat"])[1:10,]
se.treat.list.low = as.data.frame(summary(list.low, n.draws = n.draws)["se.treat"])[1:10,]
coeffs.cont.list.low = as.data.frame(summary(list.low, n.draws = n.draws)["par.control"])[1:10,]
se.cont.list.low = as.data.frame(summary(list.low, n.draws = n.draws)["se.control"])[1:10,]



##############
# List High Condition
##############

if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(list)

list.high <- ictreg(ycount ~ socideo + partyid + income.n + educ.n, 
                    data = dat.high, 
                    treat = "treatment", 
                    J=3, 
                    method = method)

## this is to construct a table later
coeffs.treat.list.high = as.data.frame(summary(list.high, n.draws = n.draws)["par.treat"])[1:10,]
se.treat.list.high = as.data.frame(summary(list.high, n.draws = n.draws)["se.treat"])[1:10,]
coeffs.cont.list.high = as.data.frame(summary(list.high, n.draws = n.draws)["par.control"])[1:10,]
se.cont.list.high = as.data.frame(summary(list.high, n.draws = n.draws)["se.control"])[1:10,]
## ----

##############
# Direct: High Condition
##############

## ---- list:analysis:social:desirability:data
direct.q.high <- glm(directquestion ~ socideo + partyid + income.n + educ.n, 
                     data = dat.high, 
                     family = binomial("logit"))


##############
# Direct: Low Condition
##############
direct.q.low <- glm(directquestion ~ socideo + partyid + income.n + educ.n, 
                    data = dat.low, 
                    family = binomial("logit"))

avg.pred.social.desirability.high <- predict.ictreg(list.high, direct.glm = direct.q.high, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.social.desirability.low <- predict.ictreg(list.low, direct.glm = direct.q.low, se.fit = TRUE, interval = "confidence", level = ci.level)


### DF for individual prediction: High Condition
socdes.p.high = data.frame(avg.pred.social.desirability.high$fit, 
                           avg.pred.social.desirability.high$se.fit,
                           c(1:3),
                           Condition = rep("High ($500)"), 3)


socdes.p.high$c.1.3 = as.factor(socdes.p.high$c.1.3)
socdes.p.high$c.1.3 <- factor(socdes.p.high$c.1.3, labels = c("List\nExperiment\n(with covariates)", "Direct\nQuestion\n(with covariates)", "Social\nDesirability\n(with covariates)"))
socdes.p.high <- socdes.p.high[c("fit", "lwr", "upr", "Condition", "c.1.3")]

### DF for individual prediction: Low Condition
socdes.p.low = data.frame(
        avg.pred.social.desirability.low$fit, 
        avg.pred.social.desirability.low$se.fit,
        c(1:3),
        Condition = rep("Low ($100)"), 3
        )

socdes.p.low$c.1.3 = as.factor(socdes.p.low$c.1.3)
socdes.p.low$c.1.3 <- factor(socdes.p.low$c.1.3, labels = c("List\nExperiment\n(with covariates)", "Direct\nQuestion\n(with covariates)", "Social\nDesirability\n(with covariates)"))
socdes.p.low <- socdes.p.low[c("fit", "lwr", "upr", "Condition", "c.1.3")]

### Rbinding both DF's
socdes.p.high.low = rbind(socdes.p.high, socdes.p.low)
rownames(socdes.p.high.low) <- NULL


## Appending with diff in means dataframe
socdes.p.high.low.diff.in.means = rbind(diff.means.plot.d,socdes.p.high.low)


### Plot
if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(ggplot2)

soc.des.plot = ggplot(socdes.p.high.low.diff.in.means, 
                      aes(c.1.3, fit, colour = Condition)) + 
        theme_bw() +
        scale_colour_grey() +
        xlab("") + 
        ylab("Estimated Proportion") +
        geom_hline(yintercept=0, colour = "red", linetype = "dashed", size = 0.2) +
        geom_pointrange(aes(
                x = c.1.3,
                ymin = lwr, 
                ymax = upr), position = position_dodge(width = 0.25)) +
        theme(axis.text.y = element_text(size=7), 
              axis.text.x = element_text(size=7), 
              axis.title.y = element_text(size=7), 
              axis.title.x = element_text(size=7), 
              legend.text=element_text(size=7), 
              legend.title=element_text(size=7),
              plot.title = element_text(size=7),
              legend.position="bottom")

## ----


## ---- list:analysis:social:desirability:plot
### Plot
soc.des.plot
soc.des.plot.note <- paste(
        "{\\bf List Experiment Data: Declared and Predicted Vote-Sellers}.",
        "\\\\\\hspace{\\textwidth}", paste("{\\bf Note}: The figure summarizes \\autoref{tab:t:test} by showing the simple difference in means (without covariates). It also shows the proportion of declared (``Direct Question'') and predicted (``List Experiment'') hypothetical vote sellers, and the difference (``Social Desirability''). The three sets of main estimates were obtained via a multivariate procedure (including covariates). Combining both ``low'' and ``high'' treatments," , paste(round(mean(c(round(socdes.p.high.low.diff.in.means$fit[socdes.p.high.low.diff.in.means$c.1.3=="List\nExperiment\n(with covariates)" & socdes.p.high.low.diff.in.means$Condition=="High ($500)"]*100,0), round(socdes.p.high.low.diff.in.means$fit[socdes.p.high.low.diff.in.means$c.1.3=="List\nExperiment\n(with covariates)" & socdes.p.high.low.diff.in.means$Condition=="Low ($100)"]*100,0))),0), "\\%", sep=""), "would be willing to sell their votes. And of those who answered affirmatively when asked directly", paste("(", round(mean(c(socdes.p.high.low.diff.in.means$fit[socdes.p.high.low.diff.in.means$c.1.3=="Direct\nQuestion\n(with covariates)" & socdes.p.high.low.diff.in.means$Condition=="High ($500)"]*100, socdes.p.high.low.diff.in.means$fit[socdes.p.high.low.diff.in.means$c.1.3=="Direct\nQuestion\n(with covariates)" & socdes.p.high.low.diff.in.means$Condition=="Low ($100)"]*100)),0), "\\%)", sep = ""), " an estimated additional ", paste(round(mean(c(socdes.p.high.low.diff.in.means$fit[socdes.p.high.low.diff.in.means$c.1.3=="Social\nDesirability\n(with covariates)" & socdes.p.high.low.diff.in.means$Condition=="High ($500)"]*100, socdes.p.high.low.diff.in.means$fit[socdes.p.high.low.diff.in.means$c.1.3=="Social\nDesirability\n(with covariates)" & socdes.p.high.low.diff.in.means$Condition=="Low ($100)"]*100)),0), "\\%", sep = ""), "lied about it. ``Liars'' answer the direct question negatively, but they are likely sellers. The figure shows 95\\% confidence intervals. There are two arbitrarily ``low'' and ``high'' vote-selling prices. The reason for having both was to control for possible price elasticities. The figure suggests some small differences that are not statistically significant. Consequently, these arbitrary pricing decisions do not threaten the experimental design."
                                           )
        )
## ----


#####################################################
# Table A1
#####################################################

## ---- design:effects:data
# define parameters
j = 3
alpha = 0.05

# Design effects [High Treatment]
ycount.high <- dat.high$ycount
treatment.high <- dat.high$treatment
design.high <- ict.test(ycount.high, treatment.high, J=j, gms = TRUE, n.draws = n.draws, alpha =alpha, pi.table = TRUE)


# Design effects [Low Treatment]
ycount.low <- dat.low$ycount
treatment.low <- dat.low$treatment
design.low <- ict.test(ycount.low, treatment.low, J=j, gms = TRUE, n.draws = n.draws, alpha =alpha, pi.table = TRUE)

## ----







######################################################
# FIGURE 5
######################################################


## ---- predictions:independent:variables:data ----
# 2 socideo (ok)

load("dat_list.RData") # Load data


socideoVC <- socideoC <- socideoM <- socideoL <- socideoVL <- dat
socideoVL <- dat[which(dat$socideo == levels(dat$socideo)[1]),]
socideoL <- dat[which(dat$socideo == levels(dat$socideo)[2]),]
socideoM <- dat[which(dat$socideo == levels(dat$socideo)[3]),]
socideoC <- dat[which(dat$socideo == levels(dat$socideo)[4]),]
socideoVC <- dat[which(dat$socideo == levels(dat$socideo)[5]),]

if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(list)



## low
avg.pred.socideoVL.low <- predict.ictreg(list.low, newdata = socideoVL, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.socideoL.low <- predict.ictreg(list.low, newdata = socideoL, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.socideoM.low <- predict.ictreg(list.low, newdata = socideoM, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.socideoC.low <- predict.ictreg(list.low, newdata = socideoC, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.socideoVC.low <- predict.ictreg(list.low, newdata = socideoVC, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)

socideo.p.low = data.frame(
        t(avg.pred.socideoVL.low$fit), 
        t(avg.pred.socideoL.low$fit), 
        t(avg.pred.socideoM.low$fit), 
        t(avg.pred.socideoC.low$fit), 
        t(avg.pred.socideoVC.low$fit)
)

socideo.p.low = data.frame(t(socideo.p.low))
Significance = as.numeric(ifelse(sign(socideo.p.low$lwr) == sign(socideo.p.low$upr), 1,0))
socideo.p.low["Significance"] <- Significance
socideo.p.low$Significance[socideo.p.low$Significance==1] <- "Yes"
socideo.p.low$Significance[socideo.p.low$Significance==0] <- "No"
socideo.p.low$socioideo = as.factor(c(1:5))
socideo.p.low$socioideo <- factor(socideo.p.low$socioideo, levels = c(1:5), labels = c("Very \n Liberal", "Liberal", "Moderate", "Conservative", "Very \n Conservative"))
socideo.p.low$'Experimental Condition' <- "Low"



## high
avg.pred.socideoVL.high <- predict.ictreg(list.high, newdata = socideoVL, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.socideoL.high <- predict.ictreg(list.high, newdata = socideoL, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.socideoM.high <- predict.ictreg(list.high, newdata = socideoM, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.socideoC.high <- predict.ictreg(list.high, newdata = socideoC, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.socideoVC.high <- predict.ictreg(list.high, newdata = socideoVC, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)

socideo.p.high = data.frame(
        t(avg.pred.socideoVL.high$fit), 
        t(avg.pred.socideoL.high$fit), 
        t(avg.pred.socideoM.high$fit), 
        t(avg.pred.socideoC.high$fit), 
        t(avg.pred.socideoVC.high$fit)
)

socideo.p.high = data.frame(t(socideo.p.high))
Significance = as.numeric(ifelse(sign(socideo.p.high$lwr) == sign(socideo.p.high$upr), 1,0))
socideo.p.high["Significance"] <- Significance
socideo.p.high$Significance[socideo.p.high$Significance==1] <- "Yes"
socideo.p.high$Significance[socideo.p.high$Significance==0] <- "No"
socideo.p.high$socioideo = as.factor(c(1:5))
socideo.p.high$socioideo <- factor(socideo.p.high$socioideo, levels = c(1:5), labels = c("Very \n Liberal", "Liberal", "Moderate", "Conservative", "Very \n Conservative"))
socideo.p.high$'Experimental Condition' <- "High"


# Plot: High and Low
## merge both datasets
socideo.plot.d = data.frame(rbind(socideo.p.high, socideo.p.low))
names(socideo.plot.d)[names(socideo.plot.d) == "Experimental.Condition"] <- "Experimental Condition"

## plots
if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(ggplot2)

socio.plot = ggplot(socideo.plot.d, 
                    aes(sign, fit, colour = `Experimental Condition`)) + 
        theme_bw() +
        scale_colour_grey() +
        xlab("") + 
        ggtitle("Ideology") + 
        ylab("Estimated Proportion") +
        geom_hline(yintercept=0, colour = "red", linetype = "dashed", size = 0.9) +
        geom_pointrange(aes(
                x = socideo.plot.d$socioideo,
                ymin = socideo.plot.d$lwr, 
                ymax = socideo.plot.d$upr), 
                position = position_dodge(width = 0.25)) +
        theme(axis.text.y = element_text(size=13), 
              axis.text.x = element_text(size=13), 
              axis.title.y = element_text(size=13), 
              axis.title.x = element_text(size=13), 
              legend.text=element_text(size=13), 
              legend.title=element_text(size=13),
              plot.title = element_text(size=14, face = "bold", hjust = 0.5),
              legend.position="bottom")



######################################################
# 3 partyid (ok)
load("dat_list.RData") # Load data


partyidD <- partyidR <- partyidI <- partyidSE <- dat

partyidD <- dat[which(dat$partyid == levels(dat$partyid)[1]),]
partyidR <- dat[which(dat$partyid == levels(dat$partyid)[2]),]
partyidI <- dat[which(dat$partyid == levels(dat$partyid)[3]),]
partyidSE <- dat[which(dat$partyid == levels(dat$partyid)[4]),]

if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(list)


## low
avg.pred.partyidD.low  <- predict.ictreg(list.low, newdata = partyidD, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.partyidR.low  <- predict.ictreg(list.low, newdata = partyidR, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.partyidI.low  <- predict.ictreg(list.low, newdata = partyidI, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.partyidSE.low  <- predict.ictreg(list.low, newdata = partyidSE, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)

partyid.p.low = data.frame(
        t(avg.pred.partyidD.low$fit), 
        t(avg.pred.partyidR.low$fit), 
        t(avg.pred.partyidI.low$fit), 
        t(avg.pred.partyidSE.low$fit)
)

partyid.p.low = data.frame(t(partyid.p.low))
Significance = as.numeric(ifelse(sign(partyid.p.low$lwr) == sign(partyid.p.low$upr), 1,0))
partyid.p.low["Significance"] <- Significance
partyid.p.low$Significance[partyid.p.low$Significance==1] <- "Yes"
partyid.p.low$Significance[partyid.p.low$Significance==0] <- "No"
partyid.p.low$partyid = as.factor(c(1:4))
partyid.p.low$partyid <- factor(partyid.p.low$partyid, levels = c(1:4), labels = c("Democrat", "Republican", "Independent", "Something\nElse"))
partyid.p.low$'Experimental Condition' <- "Low"


## high
avg.pred.partyidD.high  <- predict.ictreg(list.high, newdata = partyidD, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.partyidR.high  <- predict.ictreg(list.high, newdata = partyidR, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.partyidI.high  <- predict.ictreg(list.high, newdata = partyidI, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.partyidSE.high  <- predict.ictreg(list.high, newdata = partyidSE, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)

partyid.p.high = data.frame(
        t(avg.pred.partyidD.high$fit), 
        t(avg.pred.partyidR.high$fit), 
        t(avg.pred.partyidI.high$fit), 
        t(avg.pred.partyidSE.high$fit)
)

partyid.p.high = data.frame(t(partyid.p.high))
Significance = as.numeric(ifelse(sign(partyid.p.high$lwr) == sign(partyid.p.high$upr), 1,0))
partyid.p.high["Significance"] <- Significance
partyid.p.high$Significance[partyid.p.high$Significance==1] <- "Yes"
partyid.p.high$Significance[partyid.p.high$Significance==0] <- "No"
partyid.p.high$partyid = as.factor(c(1:4))
partyid.p.high$partyid <- factor(partyid.p.high$partyid, levels = c(1:4), labels = c("Democrat", "Republican", "Independent", "Something\nElse"))
partyid.p.high$'Experimental Condition' <- "High"


# Plot: High and Low
## merge both datasets
partyid.plot.d = data.frame(rbind(partyid.p.high, partyid.p.low))
names(partyid.plot.d)[names(partyid.plot.d) == "Experimental.Condition"] <- "Experimental Condition"


## plot

if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(ggplot2)


partyid.plot = ggplot(partyid.plot.d, 
                      aes(sign, fit, colour = `Experimental Condition`)) + 
        theme_bw() +
        scale_colour_grey() +
        xlab("") + 
        ylab("Estimated Proportion") +
        geom_hline(yintercept=0, colour = "red", linetype = "dashed", size = 0.9) +
        ggtitle("Party Identification") +
        geom_pointrange(aes(
                x = partyid.plot.d$partyid,
                ymin = partyid.plot.d$lwr, 
                ymax = partyid.plot.d$upr), 
                position = position_dodge(width = 0.25)) +
        theme(axis.text.y = element_text(size=13), 
              axis.text.x = element_text(size=13), 
              axis.title.y = element_text(size=13), 
              axis.title.x = element_text(size=13), 
              legend.text=element_text(size=13), 
              legend.title=element_text(size=13),
              plot.title = element_text(size=14, face = "bold", hjust = 0.5),
              legend.position="bottom")



######################################################
# 6 educ (ok)
load("dat_list.RData") # Load data


educ.SHS <- educ.HS <- educ.T <- educ.SC <- educ.AD <- educ.BD <- educ.GS <- dat

educ.SHS = dat[which(dat$educ == levels(dat$educ)[1]),]
educ.HS = dat[which(dat$educ == levels(dat$educ)[2]),]
educ.T = dat[which(dat$educ == levels(dat$educ)[3]),]
educ.SC = dat[which(dat$educ == levels(dat$educ)[4]),]
educ.AD = dat[which(dat$educ == levels(dat$educ)[5]),]
educ.BD = dat[which(dat$educ == levels(dat$educ)[6]),]
educ.GS = dat[which(dat$educ == levels(dat$educ)[7]),]

if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(list)

## low

avg.pred.educ.SHS.low = predict.ictreg(list.low, newdata = educ.SHS, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.educ.HS.low = predict.ictreg(list.low, newdata = educ.HS, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.educ.T.low = predict.ictreg(list.low, newdata = educ.T, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.educ.SC.low = predict.ictreg(list.low, newdata = educ.SC, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.educ.AD.low = predict.ictreg(list.low, newdata = educ.AD, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.educ.BD.low = predict.ictreg(list.low, newdata = educ.BD, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.educ.GS.low = predict.ictreg(list.low, newdata = educ.GS, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)

educ.p.low = data.frame(
        t(avg.pred.educ.SHS.low$fit), 
        t(avg.pred.educ.HS.low$fit), 
        t(avg.pred.educ.T.low$fit), 
        t(avg.pred.educ.SC.low$fit), 
        t(avg.pred.educ.AD.low$fit), 
        t(avg.pred.educ.BD.low$fit), 
        t(avg.pred.educ.GS.low$fit)
)

educ.p.low = data.frame(t(educ.p.low))
Significance = as.numeric(ifelse(sign(educ.p.low$lwr) == sign(educ.p.low$upr), 1,0))
educ.p.low["Significance"] <- Significance
educ.p.low$Significance[educ.p.low$Significance==1] <- "Yes"
educ.p.low$Significance[educ.p.low$Significance==0] <- "No"
educ.p.low$education = as.factor(c(1:7))
educ.p.low$education <- factor(educ.p.low$education, levels = c(1:7), labels = c("Some \n High  \n School", "High \n  School \n  Graduate", "Technical \n  School", "Some  \n College", "Associate \n  Degree", "Bachelor's \n  Degree", "Graduate \n  School")
)
educ.p.low$'Experimental Condition' <- "Low"


## high
avg.pred.educ.SHS.high = predict.ictreg(list.high, newdata = educ.SHS, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.educ.HS.high = predict.ictreg(list.high, newdata = educ.HS, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.educ.T.high = predict.ictreg(list.high, newdata = educ.T, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.educ.SC.high = predict.ictreg(list.high, newdata = educ.SC, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.educ.AD.high = predict.ictreg(list.high, newdata = educ.AD, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.educ.BD.high = predict.ictreg(list.high, newdata = educ.BD, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
avg.pred.educ.GS.high = predict.ictreg(list.high, newdata = educ.GS, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)

educ.p.high = data.frame(
        t(avg.pred.educ.SHS.high$fit), 
        t(avg.pred.educ.HS.high$fit), 
        t(avg.pred.educ.T.high$fit), 
        t(avg.pred.educ.SC.high$fit), 
        t(avg.pred.educ.AD.high$fit), 
        t(avg.pred.educ.BD.high$fit), 
        t(avg.pred.educ.GS.high$fit)
)

educ.p.high = data.frame(t(educ.p.high))
Significance = as.numeric(ifelse(sign(educ.p.high$lwr) == sign(educ.p.high$upr), 1,0))
educ.p.high["Significance"] <- Significance
educ.p.high$Significance[educ.p.high$Significance==1] <- "Yes"
educ.p.high$Significance[educ.p.high$Significance==0] <- "No"
educ.p.high$education = as.factor(c(1:7))
educ.p.high$education <- factor(educ.p.high$education, levels = c(1:7), labels = c("Some \n High  \n School", "High \n  School \n  Graduate", "Technical \n  School", "Some  \n College", "Associate \n  Degree", "Bachelor's \n  Degree", "Graduate \n  School")
)
educ.p.high$'Experimental Condition' <- "High"

# Plot: High and Low
## merge both datasets
educ.plot.d = data.frame(rbind(educ.p.high, educ.p.low))
names(educ.plot.d)[names(educ.plot.d) == "Experimental.Condition"] <- "Experimental Condition"


## plots
if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(ggplot2)


educ.plot = ggplot(educ.plot.d, 
                   aes(sign, fit, colour = `Experimental Condition`)) + 
        theme_bw() +
        scale_colour_grey() +
        xlab("") + 
        ylab("Estimated Proportion") +
        geom_hline(yintercept=0, colour = "red", linetype = "dashed", size = 0.9) +
        ggtitle("Education") +
        geom_pointrange(aes(
                x = educ.plot.d$education,
                ymin = educ.plot.d$lwr, 
                ymax = educ.plot.d$upr), 
                position = position_dodge(width = 0.25)) +
        theme(axis.text.y = element_text(size=13), 
              axis.text.x = element_text(size=13), 
              axis.title.y = element_text(size=13), 
              axis.title.x = element_text(size=13), 
              legend.text=element_text(size=13), 
              legend.title=element_text(size=13),
              plot.title = element_text(size=14, face = "bold", hjust = 0.5),
              legend.position="bottom")


######################################################
# 6 income.n (ok)

load("dat_list.RData") # Load data

income.1 <- income.2 <- income.3 <- income.4 <- income.5 <- income.6 <- income.7 <- income.8 <- income.9 <- income.10 <- income.11 <- income.12 <- income.13 <- income.14 <- dat

income.1  = dat[which(dat$income.n == 1),]
income.2  = dat[which(dat$income.n == 2),]
income.3  = dat[which(dat$income.n == 3),]
income.4  = dat[which(dat$income.n == 4),]
income.5  = dat[which(dat$income.n == 5),]
income.6  = dat[which(dat$income.n == 6),]
income.7  = dat[which(dat$income.n == 7),]
income.8  = dat[which(dat$income.n == 8),]
income.9  = dat[which(dat$income.n == 9),]
income.10 = dat[which(dat$income.n == 10),] 
income.11 = dat[which(dat$income.n == 11),] 
income.12 = dat[which(dat$income.n == 12),] 
income.13 = dat[which(dat$income.n == 13),] 
income.14 = dat[which(dat$income.n == 14),] 


if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(list)


## low
income.low.1  = predict.ictreg(list.low, newdata = income.1, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.low.2  = predict.ictreg(list.low, newdata = income.2, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.low.3  = predict.ictreg(list.low, newdata = income.3, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.low.4  = predict.ictreg(list.low, newdata = income.4, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.low.5  = predict.ictreg(list.low, newdata = income.5, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.low.6  = predict.ictreg(list.low, newdata = income.6, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.low.7  = predict.ictreg(list.low, newdata = income.7, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.low.8  = predict.ictreg(list.low, newdata = income.8, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.low.9  = predict.ictreg(list.low, newdata = income.9, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.low.10 = predict.ictreg(list.low, newdata = income.10, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level) 
income.low.11 = predict.ictreg(list.low, newdata = income.11, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level) 
income.low.12 = predict.ictreg(list.low, newdata = income.12, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level) 
income.low.13 = predict.ictreg(list.low, newdata = income.13, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level) 
income.low.14 = predict.ictreg(list.low, newdata = income.14, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level) 


income.p.low = data.frame(
        t(income.low.1$fit), 
        t(income.low.2$fit), 
        t(income.low.3$fit), 
        t(income.low.4$fit), 
        t(income.low.5$fit), 
        t(income.low.6$fit), 
        t(income.low.7$fit), 
        t(income.low.8$fit), 
        t(income.low.9$fit), 
        t(income.low.10$fit), 
        t(income.low.11$fit), 
        t(income.low.12$fit), 
        t(income.low.13$fit), 
        t(income.low.14$fit)
)

income.p.low = data.frame(t(income.p.low))
Significance = as.numeric(ifelse(sign(income.p.low$lwr) == sign(income.p.low$upr), 1,0))
income.p.low["Significance"] <- Significance
income.p.low$Significance[income.p.low$Significance==1] <- "Yes"
income.p.low$Significance[income.p.low$Significance==0] <- "No"
income.p.low$income = as.factor(c(1:14))
income.p.low$income <- factor(income.p.low$income, levels = c(1:14), labels = c("Less\nthan\n$20,000", "$20,000\nto\n$24,999", "$25,000\nto\n$29,999", "$30,000\nto\n$34,999", "$35,000\nto\n$39,999", "$40,000\nto\n$49,999", "$50,000\nto\n$59,999", "$60,000\nto\n$74,999", "$75,000\nto\n$84,999", "$85,000\nto\n$99,999", "$100,00\n to\n $124,999", "$125,00\n to\n $149,999", "$150,00\n to\n $174,999", "$175,000\nor\nmore")
)
income.p.low$'Experimental Condition' <- "Low"


## high
income.high.1  = predict.ictreg(list.high, newdata = income.1, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.high.2  = predict.ictreg(list.high, newdata = income.2, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.high.3  = predict.ictreg(list.high, newdata = income.3, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.high.4  = predict.ictreg(list.high, newdata = income.4, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.high.5  = predict.ictreg(list.high, newdata = income.5, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.high.6  = predict.ictreg(list.high, newdata = income.6, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.high.7  = predict.ictreg(list.high, newdata = income.7, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.high.8  = predict.ictreg(list.high, newdata = income.8, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.high.9  = predict.ictreg(list.high, newdata = income.9, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level)
income.high.10 = predict.ictreg(list.high, newdata = income.10, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level) 
income.high.11 = predict.ictreg(list.high, newdata = income.11, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level) 
income.high.12 = predict.ictreg(list.high, newdata = income.12, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level) 
income.high.13 = predict.ictreg(list.high, newdata = income.13, avg = TRUE, se.fit = TRUE, interval = "confidence", level = ci.level) 
income.high.14 = predict.ictreg(list.high, newdata = income.14, avg = TRUE, se.fit = TRUE, interval = "confidence" ,level = ci.level) 


income.p.high = data.frame(
        t(income.high.1$fit), 
        t(income.high.2$fit), 
        t(income.high.3$fit), 
        t(income.high.4$fit), 
        t(income.high.5$fit), 
        t(income.high.6$fit), 
        t(income.high.7$fit), 
        t(income.high.8$fit), 
        t(income.high.9$fit), 
        t(income.high.10$fit), 
        t(income.high.11$fit), 
        t(income.high.12$fit), 
        t(income.high.13$fit), 
        t(income.high.14$fit)
)

income.p.high = data.frame(t(income.p.high))
Significance = as.numeric(ifelse(sign(income.p.high$lwr) == sign(income.p.high$upr), 1,0))
income.p.high["Significance"] <- Significance
income.p.high$Significance[income.p.high$Significance==1] <- "Yes"
income.p.high$Significance[income.p.high$Significance==0] <- "No"
income.p.high$income = as.factor(c(1:14))
income.p.high$income <- factor(income.p.high$income, levels = c(1:14), labels = c("Less\nthan\n$20,000", "$20,000\nto\n$24,999", "$25,000\nto\n$29,999", "$30,000\nto\n$34,999", "$35,000\nto\n$39,999", "$40,000\nto\n$49,999", "$50,000\nto\n$59,999", "$60,000\nto\n$74,999", "$75,000\nto\n$84,999", "$85,000\nto\n$99,999", "$100,00\n to\n $124,999", "$125,00\n to\n $149,999", "$150,00\n to\n $174,999", "$175,000\nor\nmore")
)
income.p.high$'Experimental Condition' <- "High"


# Plot: High and Low
## merge both datasets
income.plot.d = data.frame(rbind(income.p.high, income.p.low))
names(income.plot.d)[names(income.plot.d) == "Experimental.Condition"] <- "Experimental Condition"


## plots
if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(ggplot2)

income.plot = ggplot(income.plot.d, 
                     aes(sign, fit, colour = `Experimental Condition`)) + 
        theme_bw() +
        scale_colour_grey() +
        xlab("") + 
        ylab("Estimated Proportion") +
        geom_hline(yintercept=0, colour = "red", linetype = "dashed", size = 0.9) +
        geom_pointrange(aes(
                x = income.plot.d$income,
                ymin = income.plot.d$lwr, 
                ymax = income.plot.d$upr), 
                position = position_dodge(width = 0.25)) +
        ggtitle("Income") + 
        theme(axis.text.y = element_text(size=13), 
              axis.text.x = element_text(size=7, angle = 45, hjust=1), 
              axis.title.y = element_text(size=13), 
              axis.title.x = element_text(size=7), 
              legend.text=element_text(size=13), 
              legend.title=element_text(size=13),
              plot.title = element_text(size=14, face = "bold", hjust = 0.5),
              legend.position="bottom")


######################################################
#### MERGING ALL PLOTS

# load libraries
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(ggplot2,gridExtra)



# To force GGplots to share same legend.
grid_arrange_shared_legend <- function(...) {
        require(ggplot2)
        require(gridExtra)
        plots <- list(...)
        g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
        legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
        lheight <- sum(legend$height)
        grid.arrange(
                do.call(arrangeGrob, lapply(plots, function(x)
                        x + theme(legend.position="none"))),
                legend,
                ncol = 1,
                heights = grid::unit.c(unit(1, "npc") - lheight, lheight))
}

#### multiplot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
        library(grid)
        
        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)
        
        numPlots = length(plots)
        
        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
        }
        
        if (numPlots==1) {
                print(plots[[1]])
                
        } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                        # Get the i,j matrix positions of the regions that contain this subplot
                        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                        
                        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                        layout.pos.col = matchidx$col))
                }
        }
}
## ---- 

## ---- predictions:independent:variables:plot ----
grid_arrange_shared_legend(
        income.plot, 
        educ.plot,
        partyid.plot,
        socio.plot,
        ncol = 2, nrow = 2)

predictions.independent.variables.plot.note <- paste(
        "{\\bf List Experiment: Covariates Used to Estimate Likely Vote-Sellers}.",
        "\\\\\\hspace{\\textwidth}", 
        paste(paste("{\\bf Note}: These variables were used in the multivariate statistical model to estimate individual-level probabilities of vote-selling. The figure shows the predicted probabilities and their corresponding ", ci.level*100, "\\% confidence intervals for", sep = ""), paste("income, education, party identification, and ideology. Since the vote-selling prices were set arbitrarily, the reason for two experimental conditions (``low'' and ``high'') was to control for possible price elasticities. While there are some perceptible changes, they are not statistically significant. Consequently, these arbitrary decisions do not threaten the identification strategy.")),
        "\n")
## ---- 



###############################################
# Table 1
# Covariate Balance Statistical Tests
###############################################


## ---- cov:balance:multinomial:data ----
# data
load("dat_list.RData") # Load data
dat.multinomial <- dat


# dependent variable: treatment condition
dat.multinomial$condition <- as.factor(ifelse(dat.multinomial$treatlow == 1 & dat.multinomial$treatment == 1, "Low", 
                                    ifelse(dat.multinomial$treatlow == 0 & dat.multinomial$treatment == 1, "High",
                                           ifelse(dat.multinomial$treatment == 0, "Control", NA)
                                           )
                                    )
                                    )

# change reference cat
dat.multinomial$condition = relevel(dat.multinomial$condition, ref = "Control")

# Multinomial Model
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(nnet) # to fit model

multinomial.output = multinom(condition ~ as.numeric(socideo) + as.numeric(partyid) + income.n + educ.n, 
                              dat = dat.multinomial)
## ---- 



###############################################
# Figure 1
# Direct CLIENTELISM question plot from LAPOP
###############################################

## ---- lapop:bar:chart:data ----
# cat("\014")
# rm(list=ls())
load("/Users/hectorbahamonde/RU/Term5/Experiments_Redlawsk/Experiment/Data/LAPOP/datLAPOP.rdata")
clientelism = datLAPOP$clien1
clientelism <- factor(clientelism, labels = c("Often", "Sometimes", "Never"))
clientelism <- na.omit(clientelism)

if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(data.table)

clientelism = data.table(clientelism)

lapop.bar.chart.N = nrow(clientelism)

percentage.never = round(as.numeric(as.numeric(table(clientelism)["Never"]) * 100) / lapop.bar.chart.N, 1)
percentage.sometimes = round(as.numeric(as.numeric(table(clientelism)["Sometimes"]) * 100) / lapop.bar.chart.N, 1)
percentage.often = round(as.numeric(as.numeric(table(clientelism)["Often"]) * 100) / lapop.bar.chart.N, 1)



if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(ggplot2)
lapop.bar.chart.p = ggplot(clientelism, aes(clientelism)) + 
        geom_bar(stat="count", width = 0.5) + 
        theme_bw() +
        xlab("") + 
        ylab("Subjects (N)") +
        theme(axis.text.y = element_text(size=7), 
              axis.text.x = element_text(size=7), 
              axis.title.y = element_text(size=7), 
              axis.title.x = element_text(size=7), 
              legend.text=element_text(size=7), 
              legend.title=element_text(size=7),
              plot.title = element_text(size=7),
              legend.position="bottom") +
        scale_x_discrete(labels=c(
                paste("Often (", paste(percentage.often, "%", sep = ""), ")", sep = ""),
                paste("Sometimes (", paste(percentage.sometimes, "%", sep = ""), ")", sep = ""),
                paste("Never (", paste(percentage.never, "%", sep = ""), ")", sep = "")
        ))
## ----

## ---- lapop:bar:chart:plot ----
### calling plot
lapop.bar.chart.p
### defining legend, title and notes.
lapop.bar.chart.p.note <- paste(
        "{\\bf Frequency of Clientelism in the United States (2010)}.",
        "\\\\\\hspace{\\textwidth}", 
        paste("{\\bf Note}: Figure shows the frequency of survey respondents, N = ", paste(lapop.bar.chart.N, ".", sep = ""), sep = ""),
        "\\\\\\hspace{\\textwidth}", 
        paste("{\\bf Source}: \\href{https://www.vanderbilt.edu/lapop/usa/2010_United_States_Questionnaire.pdf}{LAPOP}, 2010 wave for the United States. Question is \\texttt{clien1}: ``In recent years and thinking about election campaigns, has a candidate or someone from a political party offered you something, like a favor, food, or any other benefit or object in return for your vote or support? Has this happened often, sometimes, or never?''"),
        "\n")
## ----


###############################################
# Vote-selling Pricing Test Analysis
###############################################

## ---- pricing:experiment:data ----
load("dat_list.RData") # Load data

# building auxiliary dataset
pricing.d = na.omit(data.frame(dat$pricecheap,dat$priceexpensive))
colnames(pricing.d)[1] <- "Cheap"
colnames(pricing.d)[2] <- "Expensive"


# reshaping dataset for density plots
if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(reshape2)

pricing.d<- melt(pricing.d)
# mean(pricing.d$value[pricing.d$variable=="Cheap"])

# Plot
if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(ggplot2)

price.plot = ggplot(pricing.d,aes(x=value)) + geom_density((aes(fill = variable)), alpha=0.25) + # colour = variable
        xlab("Price for your vote") + 
        ylab("Density") +
        theme_bw() +
        #scale_fill_grey() +
        theme(legend.position="bottom", legend.direction="horizontal")  +
        theme(axis.text.y = element_text(size=7), 
              axis.text.x = element_text(size=7), 
              axis.title.y = element_text(size=7), 
              axis.title.x = element_text(size=7), 
              legend.text=element_text(size=7), 
              legend.title=element_text(size=0),
              plot.title = element_text(size=7),
              legend.position="bottom",
              legend.key.size = unit(0.5,"cm"),
              legend.spacing.x = unit(0.3, 'cm'))



## getting intersecting price
pricing.d2 = data.frame(na.omit(data.frame(dat$pricecheap,dat$priceexpensive)))

# Getting the intersecting point.
## read data in the table
if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(gdata)
options(scipen=999)

price.plot.d = ggplot_build(price.plot)
price.plot.d = data.frame(price.plot.d$data[[1]])
price.plot.d <- price.plot.d[order(price.plot.d$group),] # sort by price.plot.d


price.plot.d.1 = subset(price.plot.d, group == 1)
price.plot.d.1 <- price.plot.d.1[order(price.plot.d.1$x),] # sort by x

price.plot.d.2 = subset(price.plot.d, group == 2)
price.plot.d.2 <- price.plot.d.2[order(price.plot.d.2$x),] # sort by x

### https://github.com/andrewheiss/reconPlots
### https://www.andrewheiss.com/blog/2017/09/15/create-supply-and-demand-economics-curves-with-ggplot2/
### this package was used on Dec 11 2017
### It looks very new, and the arguments might not work in the future
### Today it was compiling good. 
if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(devtools,Hmisc) # reconPlots needs, during installing process, Hmisc
install_github("andrewheiss/reconPlots")
library(reconPlots)


#### define x and y // package WONT work if these columns don't have same name
curve.1 = data.frame(price.plot.d.1$density, price.plot.d.1$x)
colnames(curve.1)[1] <- "y"
colnames(curve.1)[2] <- "x"

#### define x and y // package WONT work if these columns don't have same name
curve.2 = data.frame(price.plot.d.2$density, price.plot.d.2$x)
colnames(curve.2)[1] <- "y"
colnames(curve.2)[2] <- "x"

line_intersection <- curve_intersect(curve.1, curve.2)
## ---- 










################
#### ABSTRACT
################

## ---- abstract ----
fileConn <- file ("abstract.txt")
writeLines(c("Hello","World"), fileConn)
writeLines(paste("In nineteenth-century United States politics, vote buying was commonplace. Nowadays, vote buying seems to have declined.  The quantitative empirical literature emphasizes vote buying, ignoring the micro-dynamics of vote selling. We seem to know that vote buyers can no longer afford this strategy; however, we do not know what American voters would do if offered the chance to sell their vote. Would they sell, and at what price, or would they consistently opt out of vote selling? A novel experimental dataset representative at the national level comprises ", total.sample.size, "U.S. voters who participated in an online list experiment in 2016, and the results are striking: Approximately ", paste(round(((socdes.p.high.low.diff.in.means$fit[socdes.p.high.low.diff.in.means$c.1.3=="List\nExperiment\n(with covariates)"][1] + socdes.p.high.low.diff.in.means$fit[socdes.p.high.low.diff.in.means$c.1.3=="List\nExperiment\n(with covariates)"][2]) / 2) * 100, 0), "\\%", sep=""), "would sell their vote for a minimum payment of ", paste("\\$", round(mean(pricing.d$value[pricing.d$variable=="Cheap"]),0), ".", sep=""), "Democrats and Liberals are more likely to sell, while education or income levels do not seem to impact the likelihood of vote selling."), fileConn)
close(fileConn)
## ----






######################################################
# Figure A1
# Descriptive Maps
######################################################

## ---- us:map:plot ----
# Load Data
load("dat_list.RData") # Load data


## Map of Observations

# install package zipcode from source when I run this script during January 2020.
## zipcode was removed from R, and last version (bellow installed) is 1.0 (dated 2012).
install.packages("https://cran.r-project.org/src/contrib/Archive/zipcode/zipcode_1.0.tar.gz", repos=NULL, type="source")


# plyr is dicontinued/retired as it January 2020. Hence, I'll be installing from source
# install.packages("https://cran.r-project.org/src/contrib/plyr_1.8.6.tar.gz", repos=NULL, type="source")
library(plyr)

if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(zipcode,ggplot2,ggmap)
# install.packages("ggmap", type = "source")library(ggmap) # if it gives the following error ("Error: GeomRasterAnn was built with an incompatible version of ggproto"), install ggmap from source.



data(zipcode)
zipcode <- zipcode[c("zip", "latitude", "longitude")]
dat <- merge(dat,zipcode,by=c("zip"))
us <- c(left = -160, bottom = 15, right = -55, top = 50)
map <- get_stamenmap(us, zoom = 4, maptype = "toner-lite", scale = 2, format = "png")

levels(dat$partyid)[levels(dat$partyid)=="SomethingElse"] <- "Something Else"


ggmap(map) + geom_point(aes(
        x = longitude,
        shape=partyid,
        y = latitude), 
        alpha = .7,
        size = 0.8,
        #shape = 1,
        data = dat) +
        xlab("Longitude") + 
        ylab("Latitude") +
        theme_bw() +
        labs(color='') +
        #scale_colour_manual(values=c("blue", "red", "forestgreen", "cyan1")) +
        scale_color_grey() +
        theme_bw() +
        theme(axis.text.y = element_text(size=7), 
              axis.text.x = element_text(size=7), 
              axis.title.y = element_text(size=7), 
              axis.title.x = element_text(size=7), 
              legend.text=element_text(size=7), 
              legend.title=element_text(size=7),
              plot.title = element_text(size=7),
              legend.position="bottom")

## ----




##############
# Figure A2
# Plot Individual predictions (both High and Low conditions)
##############



## ---- list:analysis:individual:predictions:data  ----
### Individual posterior likelihoods of vote-selling (high)
list.high.predicted.2B <- predict.ictreg(list.high, se.fit = TRUE, interval= "confidence", avg = F, return.draws = T, level = ci.level)
list.high.predicted.2B$fit<-round(list.high.predicted.2B$fit, 2)
list.high.predicted.2B$se.fit<-round(list.high.predicted.2B$se.fit, 2)
indpred.p.high = data.frame(
        list.high.predicted.2B$fit, 
        list.high.predicted.2B$se.fit, 
        Significance = as.numeric(ifelse(sign(list.high.predicted.2B$fit$lwr) == sign(list.high.predicted.2B$fit$up), 1,0)))
indpred.p.high$Significance[indpred.p.high$Significance==1] <- "Yes"
indpred.p.high$Significance[indpred.p.high$Significance==0] <- "No"
names(indpred.p.high)[4] = "se.fit"
rownames(indpred.p.high) <- NULL
indpred.p.high.fit= indpred.p.high$fit

### Individual posterior likelihoods of vote-selling (low)
list.low.predicted.2B <- predict.ictreg(list.low, se.fit = TRUE, interval= "confidence", avg = F, return.draws = T, level = ci.level)
list.low.predicted.2B$fit<-round(list.low.predicted.2B$fit, 2)
list.low.predicted.2B$se.fit<-round(list.low.predicted.2B$se.fit, 2)
indpred.p.low = data.frame(
        list.low.predicted.2B$fit, 
        list.low.predicted.2B$se.fit, 
        Significance = as.numeric(ifelse(sign(list.low.predicted.2B$fit$lwr) == sign(list.low.predicted.2B$fit$upr), 1,0)))
indpred.p.low$Significance[indpred.p.low$Significance==1] <- "Yes"
indpred.p.low$Significance[indpred.p.low$Significance==0] <- "No"
names(indpred.p.low)[4] = "se.fit"
rownames(indpred.p.low) <- NULL
indpred.p.low.fit= indpred.p.low$fit


# load libraries
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(ggplot2,grid,gridExtra)




## Low
ind.pred.low.cond.plot = ggplot() + geom_pointrange(data=indpred.p.low, 
                                                    mapping =aes(
                                                            x=1:nrow(indpred.p.low), 
                                                            y=indpred.p.low$fit, 
                                                            ymin=indpred.p.low$lwr, 
                                                            ymax=indpred.p.low$upr, 
                                                            colour = Significance), 
                                                    size=0.25, 
                                                    alpha=.5) + 
        #theme(legend.position="none") + 
        geom_hline(yintercept=0, colour = "red", linetype = "dashed", size = 0.9) +
        xlab("Observations") + 
        ylab("Probability of Vote-Selling\n(Low Condition)") +
        #guides(colour=FALSE) + 
        theme_bw() +
        scale_colour_grey() +
        theme(axis.text.y = element_text(size=7), 
              axis.text.x = element_text(size=12), 
              axis.title.y = element_text(size=18), 
              axis.title.x = element_text(size=18), 
              legend.text=element_text(size=18), 
              legend.title=element_text(size=18),
              plot.title = element_text(size=7),
              legend.position="bottom")

## High
ind.pred.high.cond.plot = ggplot() + geom_pointrange(data=indpred.p.high, 
                                                     mapping =aes(
                                                             x=1:nrow(indpred.p.high), 
                                                             y=indpred.p.high$fit, 
                                                             ymin=indpred.p.high$lwr, 
                                                             ymax=indpred.p.high$upr, 
                                                             colour = indpred.p.high$Significance), 
                                                     size=0.25, 
                                                     alpha=.5) + 
        #theme(legend.position="none") + 
        geom_hline(yintercept=0, colour = "red", linetype = "dashed", size = 0.9) +
        xlab("Observations") + 
        ylab("Probability of Vote-Selling\n(High Condition)") +
        #guides(colour=FALSE) + 
        theme_bw() +
        scale_colour_grey() +
        theme(axis.text.y = element_text(size=7), 
              axis.text.x = element_text(size=12), 
              axis.title.y = element_text(size=18), 
              axis.title.x = element_text(size=18), 
              legend.text=element_text(size=18), 
              legend.title=element_text(size=18),
              plot.title = element_text(size=7),
              legend.position="bottom")

# computing the sample size of the list experiment (which also determines the sample size in the conjoint portion) for the paper
total.sample.size = as.character(formatC(c(nrow(indpred.p.high) + nrow(indpred.p.low)), format="d", big.mark=","))
# this here converts well into RNW if \Sexpr{} is not between $$ signs.

## merging the two plots
# load libraries
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(ggplot2,gridExtra)



# To force GGplots to share same legend.
grid_arrange_shared_legend <- function(...) {
        require(ggplot2)
        require(gridExtra)
        plots <- list(...)
        g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
        legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
        lheight <- sum(legend$height)
        grid.arrange(
                do.call(arrangeGrob, lapply(plots, function(x)
                        x + theme(legend.position="none"))),
                legend,
                ncol = 1,
                heights = grid::unit.c(unit(1, "npc") - lheight, lheight))
}

#### multiplot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
        library(grid)
        
        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)
        
        numPlots = length(plots)
        
        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
        }
        
        if (numPlots==1) {
                print(plots[[1]])
                
        } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                        # Get the i,j matrix positions of the regions that contain this subplot
                        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                        
                        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                        layout.pos.col = matchidx$col))
                }
        }
}


# USE THE VECTOR WITH INDIVUDUAL PREDICTIONS: High Condition
ind.pred.social.desirability.high <- predict.ictreg(list.high, se.fit = TRUE, interval= "confidence", avg = F, return.draws = T, level = ci.level)

ind.pred.social.desirability.high$fit<-round(ind.pred.social.desirability.high$fit, 10)
ind.pred.social.desirability.high$se.fit<-round(ind.pred.social.desirability.high$se.fit, 10)
ind.pred.social.desirability.high.d = data.frame(
        ind.pred.social.desirability.high$fit, 
        ind.pred.social.desirability.high$se.fit, 
        sign = ifelse(sign(ind.pred.social.desirability.high$fit$lwr) == sign(ind.pred.social.desirability.high$fit$upr), 1, 0))
names(ind.pred.social.desirability.high.d)[4] = "se.fit"
rownames(ind.pred.social.desirability.high.d) <- NULL


# USE THE VECTOR WITH INDIVUDUAL PREDICTIONS: Low Condition
ind.pred.social.desirability.low <- predict.ictreg(list.low, se.fit = TRUE, interval= "confidence", avg = F, return.draws = T, level = ci.level)

ind.pred.social.desirability.low$fit<-round(ind.pred.social.desirability.low$fit, 10)
ind.pred.social.desirability.low$se.fit<-round(ind.pred.social.desirability.low$se.fit, 10)
ind.pred.social.desirability.low.d = data.frame(
        ind.pred.social.desirability.low$fit, 
        ind.pred.social.desirability.low$se.fit, 
        sign = ifelse(sign(ind.pred.social.desirability.low$fit$lwr) == sign(ind.pred.social.desirability.low$fit$upr), 1, 0))
names(ind.pred.social.desirability.low.d)[4] = "se.fit"
rownames(ind.pred.social.desirability.low.d) <- NULL

# cbind regular DFs (with the two conditions) with the predictions
dat.low.with.predict = data.frame(cbind(dat.low, ind.pred.social.desirability.low.d))
dat.high.with.predict = data.frame(cbind(dat.high, ind.pred.social.desirability.high.d))
dat.with.predict = data.frame(rbind(dat.low.with.predict, dat.high.with.predict))
## ---- 



## ---- list:analysis:individual:predictions:plot  ----
grid_arrange_shared_legend(
        ind.pred.low.cond.plot, 
        ind.pred.high.cond.plot,
        ncol = 1, nrow = 2)

individual.predictions.plot.note <- paste(
        "{\\bf Individual Estimated Probabilities of Vote-Selling}.",
        "\\\\\\hspace{\\textwidth}", 
        paste(paste(paste(paste("{\\bf Note}: Figure shows the individual probabilities of vote-selling (N = ", total.sample.size, ")",  sep = ""), sep = ""), "under the ``low'' and ``high'' conditions. After fitting the model, and following the advice of \\textcite[]{Blair2012} and \\textcite[]{Imai2014a}, individual probabilities of vote-selling under the ``low'' and ``high'' conditions were estimated. ", paste("The figure also shows", paste(ci.level*100,"\\%", sep = ""), "confidence intervals.", sep = " "))),
        "\n")
## ----





