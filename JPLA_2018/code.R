## ---- loadings:data ----
load("dat.RData") # Load data


# Constructing Matched Set
set.seed(604)

if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(MatchIt,optmatch)

m.out <- matchit(large ~ wealth + munopp + polinv + pop.10,
                 discard = "both", 
                 method = "full",
                 data = dat,
                 verbose = F)


# print. <- print(m.out)
# sum.match = summary(m.out)


# Match Data
m.data <- match.data(m.out)


# Recode client1dummy after matching
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(car)

m.data$clien1dummy <- as.numeric(m.data$clien1dummy)
m.data$clien1dummy <- recode(m.data$clien1dummy, "1 = 0 ; 2 = 1")

# save matched dataset
save(m.data, file = "mdata.RData")


# Generating the Propensity Score 
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(CBPS)

fit <- CBPS(as.factor(wagehalf.4) ~  wealth + munopp + polinv + pop.10,
            data = dat, 
            iterations = 25000, 
            twostep = TRUE, 
            method = "over",
            ATT = 0,
            standardize = F)


## transform the weight var. // Attaching weights to DF // sorting for GEE models
dat$weights = as.numeric(fit$weights)
## ----


#####################################################################
###  M  O D E L S
#####################################################################

## ---- models ----
# load data
## Recode Before modeling
dat$clien1dummy <- as.numeric(dat$clien1dummy)

if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(car)

dat$clien1dummy <- recode(dat$clien1dummy, "1 = 0 ; 2 = 1")

# formulas 
model.m = formula(clien1dummy ~ wealth*munopp*large + pop.10 + urban + polinv + ing4 + vb3 + exc7 + ed)
model.gps = formula(clien1dummy ~ wealth*munopp*wagehalf.4 + pop.10 + urban + polinv + ing4 + vb3 + exc7 + ed + weights)


if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(Zelig)

model.m.s = zelig(model.m, 
                  model = "logit.gee",
                  id = "municipality", 
                  weights = "wt",
                  std.err = "san.se",
                  corstr = "exchangeable",
                  data = m.data, 
                  cite = F)




model.gps.s = zelig(model.gps, 
                    model = "logit.gee",
                    id = "municipality", 
                    weights = "wt",
                    std.err = "san.se",
                    corstr = "exchangeable",
                    data = dat,
                    cite = F)
## ----


#####################################################################
###  Table A3     
#####################################################################

## ---- tab:results:data ----

# Recode client1dummy after matching
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(car)

dat$clien1dummy <- as.numeric(dat$clien1dummy)
dat$clien1dummy <- recode(dat$clien1dummy, "1 = 0 ; 2 = 1")

# models

if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(texreg)

extract.geepack <- function(model) {
  s <- summary(model)
  names <- rownames(s$coef)
  co <- s$coef[, 1]
  se <- s$coef[, 2]
  pval <- s$coef[, 4]
  
  n <- nrow(model.frame(model))
  nclust <- length(s$geese$clusz)
  
  gof = c(n, nclust)
  gof.names = c("Num. obs.", "Num. clust.")
  
  tr <- createTexreg(
    coef.names = names,
    coef = co,
    se = se,
    pvalues = pval,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = rep(FALSE, length(gof))
  )
  return(tr)
}


if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(geepack)

# formulas 
model.m = formula(clien1dummy ~ wealth*munopp*large + pop.10 + urban + polinv + ing4 + vb3 + exc7 + ed)
model.gps = formula(clien1dummy ~ wealth*munopp*wagehalf.4 + pop.10 + urban + polinv + ing4 + vb3 + exc7 + ed + weights)


options(scipen=999)

model.m.t = extract.geepack(model.m.model <- geeglm(model.m,
                               family = binomial(link = "logit"), 
                               id = municipality, 
                               weights = wt,
                               std.err = "san.se",
                               corstr = "exchangeable",
                               data = m.data))

model.gps.t = extract.geepack(model.gps.model <- geeglm(model.gps,
                                 family = binomial(link = "logit"), 
                                 id = municipality, 
                                 weights = wt,
                                 std.err = "san.se",
                                 corstr = "exchangeable",
                                 data = dat))

custom.coef.names = c(
        "(Intercept)", 
        "Wealth Index", 
        "Municipal Opposition", 
        "High Poor Density", 
        "Municipal Population", 
        "Urban", 
        "Political Involvement", 
        "Support for Democracy", 
        "Party Id.", 
        "Perception of Corruption", 
        "Years of Education", 
        "Wealth Index * Municipal Opposition", 
        "Wealth Index * High Poor Density", 
        "Municipal Opposition * High Poor Density", 
        "Wealth Index * Municipal Opposition * High Poor Density", 
        "Density of the Poor", 
        "Wealth Index * Density of the Poor", 
        "Municipal Opposition * Density of the Poor", 
        "Wealth Index * Municipal Opposition * Density of the Poor")
## ----


## ---- tab:results:table ----
# table
texreg(
  c(model.m.t,model.gps.t), 
  caption = "Generalized Estimating Logistic Equations: Clientelism",
custom.model.names = c("Matched","Weighted"),
custom.coef.names = custom.coef.names,
omit.coef = "weights",
label = "tab:1",
custom.note = ("\\parbox{.65\\linewidth}{\\vspace{2pt}%stars. Clustered standard errors at the municipality level. First column shows the estimates using the matched dataset. Second column shows the estimates of the weighted model (the generalized propensity score was omitted in the table). Both models are logit GEE.}"),
fontsize = "scriptsize",
digits = 3,
center = TRUE,
no.margin = TRUE, 
float.pos = "h")
## ----




##########################################################################
### Figure 4
##########################################################################

## ---- plot:four:quadrants:d ----
set.seed(602); options(scipen=999)

N = 250

# simulation matched data
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(Zelig)


high.poor.lowcomp.m = data.frame(competition = rep("Low Competition", N), income = rep("Poor Individuals", N), x = sim(model.m.s, x = setx(model.m.s, cond = TRUE,large = max(m.data$large), wealth= quantile(m.data$wealth, .25), munopp = min(m.data$munopp)), num=N)$get_qi("ev"))
high.poor.highcomp.m = data.frame(competition = rep("High Competition", N),income = rep("Poor Individuals", N), x = sim(model.m.s, x = setx(model.m.s, cond = TRUE,large = max(m.data$large), wealth= quantile(m.data$wealth, .25), munopp = max(m.data$munopp)), num=N)$get_qi("ev"))
high.rich.lowcomp.m = data.frame(competition = rep("Low Competition", N),income = rep("Non-Poor Individuals", N), x= sim(model.m.s, x = setx(model.m.s, cond = TRUE,large = max(m.data$large), wealth= quantile(m.data$wealth, .75), munopp = min(m.data$munopp)), num=N)$get_qi("ev"))
high.rich.highcomp.m = data.frame(competition = rep("High Competition", N),income = rep("Non-Poor Individuals", N), x= sim(model.m.s, x = setx(model.m.s, cond = TRUE,large = max(m.data$large), wealth= quantile(m.data$wealth, .75), munopp = max(m.data$munopp)), num=N)$get_qi("ev"))
low.poor.lowcomp.m = data.frame(competition = rep("Low Competition", N),income = rep("Poor Individuals", N), x= sim(model.m.s, x = setx(model.m.s, cond = TRUE,large = min(m.data$large), wealth= quantile(m.data$wealth, .25), munopp = min(m.data$munopp)), num=N)$get_qi("ev"))
low.poor.highcomp.m = data.frame(competition = rep("High Competition", N),income = rep("Poor Individuals", N), x= sim(model.m.s, x = setx(model.m.s, cond = TRUE,large = min(m.data$large), wealth= quantile(m.data$wealth, .25), munopp = max(m.data$munopp)), num=N)$get_qi("ev"))
low.rich.lowcomp.m = data.frame(competition = rep("Low Competition", N),income = rep("Non-Poor Individuals", N), x= sim(model.m.s, x = setx(model.m.s, cond = TRUE,large = min(m.data$large), wealth= quantile(m.data$wealth, .75), munopp = min(m.data$munopp)), num=N)$get_qi("ev"))
low.rich.highcomp.m = data.frame(competition = rep("High Competition", N),income = rep("Non-Poor Individuals", N), x= sim(model.m.s, x = setx(model.m.s, cond = TRUE,large = min(m.data$large), wealth= quantile(m.data$wealth, .75), munopp = max(m.data$munopp)), num=N)$get_qi("ev"))



# simulation raw/GPS data
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(Zelig)

high.poor.lowcomp.gps = data.frame(competition = rep("Low Competition", N), income = rep("Poor Individuals", N), x = sim(model.gps.s, x = setx(model.gps.s, cond = TRUE, wagehalf.4 = quantile(dat$wagehalf.4, .75), wealth= quantile(dat$wealth, .25), munopp = min(dat$munopp)), num=N)$get_qi("ev"))
high.poor.highcomp.gps = data.frame(competition = rep("High Competition", N),income = rep("Poor Individuals", N), x = sim(model.gps.s, x = setx(model.gps.s, cond = TRUE, wagehalf.4 = quantile(dat$wagehalf.4, .75), wealth= quantile(dat$wealth, .25), munopp = max(dat$munopp)), num=N)$get_qi("ev"))
high.rich.lowcomp.gps = data.frame(competition = rep("Low Competition", N),income = rep("Non-Poor Individuals", N), x= sim(model.gps.s, x = setx(model.gps.s, cond = TRUE, wagehalf.4 = quantile(dat$wagehalf.4, .75), wealth= quantile(dat$wealth, .75), munopp = min(dat$munopp)), num=N)$get_qi("ev"))
high.rich.highcomp.gps = data.frame(competition = rep("High Competition", N),income = rep("Non-Poor Individuals", N), x= sim(model.gps.s, x = setx(model.gps.s, cond = TRUE, wagehalf.4 = quantile(dat$wagehalf.4, .75), wealth= quantile(dat$wealth, .75), munopp = max(dat$munopp)), num=N)$get_qi("ev"))
low.poor.lowcomp.gps = data.frame(competition = rep("Low Competition", N),income = rep("Poor Individuals", N), x= sim(model.gps.s, x = setx(model.gps.s, cond = TRUE, wagehalf.4 = quantile(dat$wagehalf.4, .25), wealth= quantile(dat$wealth, .25), munopp = min(dat$munopp)), num=N)$get_qi("ev"))
low.poor.highcomp.gps = data.frame(competition = rep("High Competition", N),income = rep("Poor Individuals", N), x= sim(model.gps.s, x = setx(model.gps.s, cond = TRUE, wagehalf.4 = quantile(dat$wagehalf.4, .25), wealth= quantile(dat$wealth, .25), munopp = max(dat$munopp)), num=N)$get_qi("ev"))
low.rich.lowcomp.gps = data.frame(competition = rep("Low Competition", N),income = rep("Non-Poor Individuals", N), x= sim(model.gps.s, x = setx(model.gps.s, cond = TRUE, wagehalf.4 = quantile(dat$wagehalf.4, .25), wealth= quantile(dat$wealth, .75), munopp = min(dat$munopp)), num=N)$get_qi("ev"))
low.rich.highcomp.gps = data.frame(competition = rep("High Competition", N),income = rep("Non-Poor Individuals", N), x= sim(model.gps.s, x = setx(model.gps.s, cond = TRUE, wagehalf.4 = quantile(dat$wagehalf.4, .25), wealth= quantile(dat$wealth, .75), munopp = max(dat$munopp)), num=N)$get_qi("ev"))





# data frame
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(Rmisc)

ci = .95

plot.d = data.frame(
        mean = c(
                as.numeric(CI(high.poor.lowcomp.gps$x, ci = ci)[2]),  # gps
                as.numeric(CI(high.poor.highcomp.gps$x, ci = ci)[2]),  # gps
                as.numeric(CI(high.rich.lowcomp.gps$x, ci = ci)[2]),  # gps
                as.numeric(CI(high.rich.highcomp.gps$x, ci = ci)[2]),  # gps
                as.numeric(CI(low.poor.lowcomp.gps$x, ci = ci)[2]),  # gps
                as.numeric(CI(low.poor.highcomp.gps$x, ci = ci)[2]),  # gps
                as.numeric(CI(low.rich.lowcomp.gps$x, ci = ci)[2]),  # gps
                as.numeric(CI(low.rich.highcomp.gps$x, ci = ci)[2]), # gps
                as.numeric(CI(high.poor.lowcomp.m$x, ci = ci)[2]),  # matched
                as.numeric(CI(high.poor.highcomp.m$x, ci = ci)[2]),  # matched
                as.numeric(CI(high.rich.lowcomp.m$x, ci = ci)[2]),  # matched
                as.numeric(CI(high.rich.highcomp.m$x, ci = ci)[2]),  # matched
                as.numeric(CI(low.poor.lowcomp.m$x, ci = ci)[2]),  # matched
                as.numeric(CI(low.poor.highcomp.m$x, ci = ci)[2]),  # matched
                as.numeric(CI(low.rich.lowcomp.m$x, ci = ci)[2]),  # matched
                as.numeric(CI(low.rich.highcomp.m$x, ci = ci)[2]) # matched
        ),
        upper = c(
                as.numeric(CI(high.poor.lowcomp.gps$x, ci = ci)[1]),  # gps
                as.numeric(CI(high.poor.highcomp.gps$x, ci = ci)[1]), # gps
                as.numeric(CI(high.rich.lowcomp.gps$x, ci = ci)[1]), # gps
                as.numeric(CI(high.rich.highcomp.gps$x, ci = ci)[1]), # gps
                as.numeric(CI(low.poor.lowcomp.gps$x, ci = ci)[1]), # gps
                as.numeric(CI(low.poor.highcomp.gps$x, ci = ci)[1]), # gps
                as.numeric(CI(low.rich.lowcomp.gps$x, ci = ci)[1]), # gps
                as.numeric(CI(low.rich.highcomp.gps$x, ci = ci)[1]), # gps
                as.numeric(CI(high.poor.lowcomp.m$x, ci = ci)[1]),  # matched
                as.numeric(CI(high.poor.highcomp.m$x, ci = ci)[1]), # matched
                as.numeric(CI(high.rich.lowcomp.m$x, ci = ci)[1]), # matched
                as.numeric(CI(high.rich.highcomp.m$x, ci = ci)[1]), # matched
                as.numeric(CI(low.poor.lowcomp.m$x, ci = ci)[1]), # matched
                as.numeric(CI(low.poor.highcomp.m$x, ci = ci)[1]), # matched
                as.numeric(CI(low.rich.lowcomp.m$x, ci = ci)[1]), # matched
                as.numeric(CI(low.rich.highcomp.m$x, ci = ci)[1]) # matched
        ),
        lower = c(
                as.numeric(CI(high.poor.lowcomp.gps$x, ci = ci)[3]),  # gps
                as.numeric(CI(high.poor.highcomp.gps$x, ci = ci)[3]), # gps
                as.numeric(CI(high.rich.lowcomp.gps$x, ci = ci)[3]), # gps
                as.numeric(CI(high.rich.highcomp.gps$x, ci = ci)[3]), # gps
                as.numeric(CI(low.poor.lowcomp.gps$x, ci = ci)[3]), # gps
                as.numeric(CI(low.poor.highcomp.gps$x, ci = ci)[3]), # gps
                as.numeric(CI(low.rich.lowcomp.gps$x, ci = ci)[3]), # gps
                as.numeric(CI(low.rich.highcomp.gps$x, ci = ci)[3]), # gps
                as.numeric(CI(high.poor.lowcomp.m$x, ci = ci)[3]),  # matched
                as.numeric(CI(high.poor.highcomp.m$x, ci = ci)[3]), # matched
                as.numeric(CI(high.rich.lowcomp.m$x, ci = ci)[3]), # matched
                as.numeric(CI(high.rich.highcomp.m$x, ci = ci)[3]), # matched
                as.numeric(CI(low.poor.lowcomp.m$x, ci = ci)[3]), # matched
                as.numeric(CI(low.poor.highcomp.m$x, ci = ci)[3]), # matched
                as.numeric(CI(low.rich.lowcomp.m$x, ci = ci)[3]), # matched
                as.numeric(CI(low.rich.highcomp.m$x, ci = ci)[3]) # matched
        ),
        Density = c(rep("High", 4), rep("Low", 4), rep("High", 4), rep("Low", 4)),
        Wealth = rep(c(rep("Poor Individual", 2), rep("Non-Poor Individual", 2)), 4),
        Competition = rep(c("Low Competition","High Competition"),8),
        Sample = c(rep("Weighted (GPS)", 8), rep("Matched", 8))
)

# plot
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(ggplot2)

plot.four.quadrants.plot = ggplot(plot.d, aes(Density, mean,
                   ymin = upper,
                   ymax=lower,
                   colour = Sample)) + 
  geom_errorbar(width=0.2) + 
  scale_colour_manual(values= c("#7fc97f","#ef3b2c")) + 
  facet_grid(Competition~Wealth) +
        ylab("Probability of being Targeted") + xlab("Density of the Poor") +
        theme_bw() + #theme(legend.position="none") +
        theme(strip.text.x = element_text(size = 8), 
              strip.text.y = element_text(size = 8), 
              axis.title=element_text(size=10), 
              legend.text = element_text(size = 8), 
              legend.title = element_text(size = 10),
              axis.text.y = element_text(size = 8),
              legend.position="top")  #+ 
       # scale_colour_discrete(name = "Sample")
## ----


## ---- plot:four:quadrants ----
plot.four.quadrants.plot
plot.four.quadrants.plot.legend = paste(
        "{\\bf Simulated Expected Values of Clientelism}.",
        "\\\\\\hspace{\\textwidth}", 
        paste("{\\bf Note}:", paste("After fitting the models shown in \\autoref{tab:1}, this figure shows the predicted probabilities of being targeted under different scenarios, with", paste(ci*100,"\\%", sep = ""), "confidence intervals.")),"Substantively, the figure emulates the theoretical predictions of \\autoref{tab:strategy:set}. Clientelism is higher when non-poor individuals are nested in poor groups (``high'' density of the poor) in highly contested municipalities (Q1), when non-poor individuals are nested in non-poor groups (``low'' density of the poor) in scarcely contested municipalities (Q3), when poor individuals are nested in poor areas in highly contested municipalities (Q2), and when poor individuals are nested in non-poor areas in scarcely contested municipalities (Q4). For every quadrant, estimates from both the matched and weighted datasets are shown. The idea is to show that the decision of dichotomizing the density of the poor variable at its median (\\autoref{fig:tgraph:plot}) gives substantively exact results than using the continuous version of that variable via the GPS analysis.",
        "\n")
## ----







######################################################
# Figure 1
######################################################

## ---- wealth:client:plot:d ----
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(ggplot2)

wealth.plot = ggplot() + geom_jitter(
  width = 4,
  height = .45, 
  alpha = 1/4,
  aes(
    y=as.factor(dat$clien1dummy), 
    x=as.numeric(dat$wealth), 
    colour=as.factor(dat$clien1dummy))) +
  xlab("Wealth Index") + 
  ylab("Offered Him/Her to Buy Vote") + 
  theme_bw()+
  scale_colour_manual(values= c("#7fc97f","#ef3b2c")) + 
  theme(strip.text.x = element_text(size = 8), 
        strip.text.y = element_text(size = 8), 
        axis.title=element_text(size=10), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        legend.position="none")
## ---- 



## ---- wealth:client:plot ----
wealth.plot
wealth.plot.legend <- paste(
        "{\\bf Individual Wealth and Vote-Buying in Brazil}.",
        "\\\\\\hspace{\\textwidth}", 
        "{\\bf Note}: Following the advice of \\textcite{Cordova2008} and \\textcite{Cordova2009,Cordova2010}, different socio-economic variables in \\textcite{LAPOP2010} dataset were used to construct a relative wealth index. With this information, in addition to the frequency of clientelism question (\\texttt{clien1}), the figure shows that clientelist brokers target individuals at all levels of income.",
        "\n")
## ---- 


##########################
# Figure A4
##########################


## ---- pol.inv:pop.size:plot:d ----
# simulation
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(Zelig)

set.seed(602); options(scipen=999)


# low 
model.m.s.low = data.frame(
        zelig_qi_to_df(sim(model.m.s,
    x = setx(model.m.s, cond = TRUE,
             large = min(m.data$large), 
             polinv:large,
             polinv = min(m.data$polinv):max(m.data$polinv)), 
    num=700)))


model.m.s.low <- model.m.s.low[, c("polinv", "expected_value")]

df.low = data.frame(
        mean = data.frame(aggregate(cbind(expected_value)~polinv, data=model.m.s.low, FUN=CI)[,2])[,2], # mean
        Upper = data.frame(aggregate(cbind(expected_value)~polinv, data=model.m.s.low, FUN=CI)[,2])[,1], # upper
        Lower = data.frame(aggregate(cbind(expected_value)~polinv, data=model.m.s.low, FUN=CI)[,2])[,3], # lower
        Type = rep("Low Density", max(m.data$polinv)+1),
        Opposition = min(m.data$polinv):max(m.data$polinv)
        )
        



# high
model.m.s.high = data.frame(
        zelig_qi_to_df(sim(model.m.s,
    x = setx(model.m.s, cond = TRUE,
             large = max(m.data$large), 
             polinv:large,
             polinv = min(m.data$polinv):max(m.data$polinv)), 
    num=700))) ; 

model.m.s.high <- model.m.s.high[, c("polinv", "expected_value")]


df.high = data.frame(
        mean = data.frame(aggregate(cbind(expected_value)~polinv, data=model.m.s.high, FUN=CI)[,2])[,2], # mean
        Upper = data.frame(aggregate(cbind(expected_value)~polinv, data=model.m.s.high, FUN=CI)[,2])[,1], # upper
        Lower = data.frame(aggregate(cbind(expected_value)~polinv, data=model.m.s.high, FUN=CI)[,2])[,3], # lower
        Type = rep("High Density", max(m.data$polinv)+1),
        Opposition = min(m.data$polinv):max(m.data$polinv)
        )


# polinv 
model.m.s.polinv = data.frame(
        zelig_qi_to_df(sim(model.m.s, x = setx(model.m.s, cond = TRUE,
               large:munopp,
               polinv = min(m.data$polinv):max(m.data$polinv)), 
      num=300)))



model.m.s.polinv <- model.m.s.polinv[, c("polinv", "expected_value")]


model.m.s.polinv = data.frame(
        mean = data.frame(aggregate(cbind(expected_value)~polinv, data=model.m.s.polinv, FUN=CI)[,2])[,2], # mean
        Upper = data.frame(aggregate(cbind(expected_value)~polinv, data=model.m.s.polinv, FUN=CI)[,2])[,1], # upper
        Lower = data.frame(aggregate(cbind(expected_value)~polinv, data=model.m.s.polinv, FUN=CI)[,2])[,3], # lower
        Type = rep("Political Involvement", max(m.data$polinv)+1),
        Opposition = min(m.data$polinv):max(m.data$polinv)
)




### combined 3 df's
polinv.d= rbind(df.high, df.low,model.m.s.polinv)


### plot
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(ggplot2)



p1= ggplot(polinv.d, aes(x=Opposition, y=mean, colour=Type)) + 
  stat_smooth() + 
  geom_ribbon(aes(ymin=Lower, ymax=Upper, linetype=NA), alpha=0.2) +
  scale_color_manual(values=c("#ef3b2c","#7fc97f","#ffff33")) +
  stat_smooth(aes(x=Opposition,y=mean)) +
  xlab("Political Involvement") + ylab("Expected Value of Clientelism") + 
  theme_bw() + 
  theme(legend.position="top", legend.title=element_blank(), legend.key = element_rect())

##########################
#  LARGE * POP
##########################

## low 
gee.dich.m.2.low = data.frame(
        zelig_qi_to_df(sim(model.m.s,x = setx(model.m.s, cond = TRUE,
                                              large = min(m.data$large), 
                                              pop.10 = min(m.data$pop.10):max(m.data$pop.10)), 
                           num=300)))

gee.dich.m.2.low <- gee.dich.m.2.low[, c("pop.10", "expected_value")]


gee.dich.m.2.low = data.frame(
        mean = data.frame(aggregate(cbind(expected_value)~pop.10, data=gee.dich.m.2.low, FUN=CI)[,2])[,2], # mean
        Upper = data.frame(aggregate(cbind(expected_value)~pop.10, data=gee.dich.m.2.low, FUN=CI)[,2])[,1], # upper
        Lower = data.frame(aggregate(cbind(expected_value)~pop.10, data=gee.dich.m.2.low, FUN=CI)[,2])[,3], # lower
        Type = rep("Low Density", max(m.data$pop.10)),
        Population = min(m.data$pop.10):max(m.data$pop.10)
        )




## high
gee.dich.m.2.high = data.frame(
        zelig_qi_to_df(sim(model.m.s, x = setx(model.m.s, cond = TRUE,
               large = max(m.data$large), 
               pop.10 = min(m.data$pop.10):max(m.data$pop.10)),num=300)))


gee.dich.m.2.high <- gee.dich.m.2.high[, c("pop.10", "expected_value")]


gee.dich.m.2.high = data.frame(
        mean = data.frame(aggregate(cbind(expected_value)~pop.10, data=gee.dich.m.2.high, FUN=CI)[,2])[,2], # mean
        Upper = data.frame(aggregate(cbind(expected_value)~pop.10, data=gee.dich.m.2.high, FUN=CI)[,2])[,1], # upper
        Lower = data.frame(aggregate(cbind(expected_value)~pop.10, data=gee.dich.m.2.high, FUN=CI)[,2])[,3], # lower
        Type = rep("High Density", max(m.data$pop.10)),
        Population = min(m.data$pop.10):max(m.data$pop.10)
        )


## pop.10
df.pop.alone = data.frame(
        zelig_qi_to_df(sim(model.m.s,x = setx(model.m.s, cond = TRUE,
                                              pop.10 = min(m.data$pop.10):max(m.data$pop.10)),num=300)))



df.pop.alone <- df.pop.alone[, c("pop.10", "expected_value")]

df.pop.alone = data.frame(
        mean = data.frame(aggregate(cbind(expected_value)~pop.10, data=df.pop.alone, FUN=CI)[,2])[,2], # mean
        Upper = data.frame(aggregate(cbind(expected_value)~pop.10, data=df.pop.alone, FUN=CI)[,2])[,1], # upper
        Lower = data.frame(aggregate(cbind(expected_value)~pop.10, data=df.pop.alone, FUN=CI)[,2])[,3], # lower
        Type = rep("Population Size", max(m.data$pop.10)),
        Population = min(m.data$pop.10):max(m.data$pop.10)
)


### combined 3 df's
pop.d= rbind(gee.dich.m.2.low, gee.dich.m.2.high,df.pop.alone)



### plot
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(ggplot2)


p2= ggplot(pop.d, aes(x=Population, y=mean, colour=Type)) + 
  stat_smooth() + 
  geom_ribbon(aes(ymin=Lower, ymax=Upper, linetype=NA), alpha=0.2) +
  scale_color_manual(values=c("#7fc97f","#ef3b2c","#ffff33"), breaks = c('High Density','Low Density','Population Size')) +
  stat_smooth(aes(x=Population,y=mean)) +
  xlab("Municipal Population Size") + ylab("Probability of being Targeted") + 
  theme_bw() + 
  theme(axis.title.y=element_text(colour="white"), legend.position="top", legend.title=element_blank(), legend.key = element_rect())



if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(cowplot)

pol.inv.pop.size.plot = plot_grid(p1,p2,  nrow = 1, labels = "auto")
## ----



## ---- pol.inv:pop.size:plot ----
pol.inv.pop.size.plot
pol.inv.pop.size.plot.legend <- paste(
        "{\\bf Simulated Expected Probability of being Targeted: Political Involvement and Population Size}.",
        "\\\\\\hspace{\\textwidth}", 
        "{\\bf Note}: Using the estimations in \\autoref{tab:1}, the figure shows the probability of being targeted at different values of political involvement (a) and population size at the municipal level (b). The figure suggests that being nested in high-poor density areas contributes substantially more to explaining clientelism.",
        "\n")
## ---- 


#########
# Table A1
#########


## ---- tab:sum:stats:r:table ----
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(stargazer)

#library(stargazer, quietly = T)
stargazer(dat.r, 
          summary=T, 
          title = "Summary Statistics: Raw Sample.",
          label = "sumtab:raw",
          type = "latex", # 'text' for word // 'latex' for latex.
          font.size = "scriptsize",
          style= "apsr",
          covariate.labels=labels.r,
          table.placement = "h",
          notes.align = "c"
          )
## ----



############
# Figure 3
############


## ---- municipality:sample:plot ----

## df of matched set
municipality.m = data.frame(
  Municipality = as.factor(m.data$municipality),
  Sample = c(rep("Matched", length(m.data$municipality))
             )
)

## df of raw set
municipality.r = data.frame(
  Municipality = as.factor(dat$municipality),
  Sample = c(rep("Raw", length(dat$municipality))
  )
)

## rbinding the two of them
municipality.d = data.frame(rbind(municipality.m, municipality.r))
municipality.d = data.frame(table(municipality.d))

if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(ggplot2)

#mun.p1 = 
  ggplot(municipality.d, aes(x = Municipality, y = Freq, fill = Sample)) + 
    geom_bar(stat = "identity", position=position_dodge()) + 
    coord_flip() +
    scale_fill_manual(values= c("#7fc97f","#ef3b2c")) + 
    xlab("") + 
    ylab("Frequency") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          legend.key = element_rect(colour = NA, fill = NA, size = 0.5))
  
## ----

 

## ---- municipality:income:large:plot:matched:data ----
## HIGH df
high.d = data.frame(
  Municipality = as.factor(m.data$municipality[m.data$large == 1]),
  Density = c(rep("High", length(m.data$municipality[as.numeric(m.data$large == 1)]))),
  Wealth = as.numeric(m.data$wealth)[m.data$large == 1]
)

## LOW df
low.d = data.frame(
  Municipality = as.factor(m.data$municipality[m.data$large == 0]),
  Density = c(rep("Low", length(m.data$municipality[as.numeric(m.data$large == 0)]))),
  Wealth = as.numeric(m.data$wealth)[m.data$large == 0]
)

## rbinding the two of them
density.d = data.frame(rbind(high.d, low.d))

if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(ggplot2)

municipality.income.large.plot.matched.plot = 
  ggplot(density.d, aes(factor(Municipality), fill = Density)) + 
    geom_bar() + 
    scale_fill_manual(values= c("#ef3b2c","#7fc97f")) + 
    geom_point(data=density.d, 
                   position = position_jitter(width = 0.22, height = 5), 
                   size = I(1),
                   aes(
                           x=as.factor(Municipality), 
                           y=Wealth*10,
                           alpha=Wealth)) + 
    coord_flip() +
    xlab("") + 
    ylab("Frequency") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.key = element_rect(colour = NA, fill = NA, size = 0.5))
## ----


## ---- municipality:income:large:plot:matched:plot ----
municipality.income.large.plot.matched.plot 
municipality.income.large.plot.matched.plot.legend  <- paste(paste("{\\bf Distribution of Observations by Municipality, Wealth Index and Density of the Poor}."),
                                                             "\\\\\\hspace{\\textwidth}", 
                                                             paste("{\\bf Note}: The figure shows the municipalities in the analyses (matched set). For every municipality, the figure shows (1) the number of inhabitants (Y-axis), (2) whether the municipality is considered having a high or low density of the poor. High-density municipalities have more than half of their inhabitants living on less than half of the minimum wage. The figure also shows (3) individual wealth indexes." ), 
                                                                                "\n")
## ----





############################################################
## Figure 2
############################################################


## ---- tgraph:plot:d ----
dens <- density(m.data$wagehalf)
df <- data.frame(x=dens$x, y=dens$y)
df$quant <- factor(findInterval(df$x, median(m.data$wagehalf)))


if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(ggplot2)
tgraph.plot = ggplot(df, aes(x,y)) +
        geom_line() + 
        geom_ribbon(aes(ymin=0, ymax=y, fill=quant, alpha = 0.5)) + 
        xlab("Percentage of People Living on\nLess than Half of the Minimum Wage") + 
        ylab("Density") + 
        # scale_fill_brewer(palette="Greens") +
  scale_fill_manual(values= c("#7fc97f","#ef3b2c")) + 
  theme_bw() +
        theme(
                legend.position="none", 
                axis.title.y = element_text(size = 10),
                axis.title.x = element_text(size = 10)) +
        geom_segment(aes(
                x = (wagehalf=median(m.data$wagehalf)), 
                y = 0, 
                xend = (wagehalf=median(m.data$wagehalf)), 
                yend = df$y[df$x == max(df$x[df$x<=median(m.data$wagehalf)])]), 
                linetype="dotted", 
                size=.5, 
                colour = "red")
## ----

## ---- tgraph:plot ----
tgraph.plot
tgraph.plot.legend <- paste(paste("{\\bf Distribution of the Density of the Poor}."),
                                  "\\\\\\hspace{\\textwidth}", 
                            paste("{\\bf Note}: Employing Brazilian census data from the \\href{http://www.ibge.gov.br}{IBGE} (2010), the figure shows  the percentage of individuals who live on less than half of the minimum wage in a given municipality. While individual income is measured using the relative wealth index (in \\autoref{fig:wealth:client:plot}), the variable plotted here is used to measure economic development at the group level. Due to statistical reasons explained in the paper, the variable had to be dichotomized at its median", 
                                  paste("(",round(as.numeric(median(m.data$wagehalf)), 2),"\\%",").", sep = ""), "However, in separate statistical analyses shown in \\autoref{tab:1} (weighted model), the variable is used without dichotomizing it, showing the same results.", sep = " "), 
                            "\n")
## ---- 


######################################################
#  FIGURE A3
######################################################


## ---- balance:plot ----
# [balance:plot]
plot(m.out, type="hist")
## ----


