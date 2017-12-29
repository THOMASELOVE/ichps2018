## Toy example Script
## Libraries used: Epi, survival, sm, Matching, lme4, reshape2, epibasix, Hmisc, lattice and survey
## Assumes toy.csv has been imported (200 observations, 11 variables)

## Preliminary Data Cleanup

summary(toy)

## Re-expressing Binary Variables
toy$treated.f <- factor(toy$treated, levels=c(1,0),labels=c("Treated","Control"))
toy$covB.f <- factor(toy$covB, levels=c(1,0), labels=c("Has B", "No B"))
toy$out2.f <- factor(toy$out2.event, levels=c("Yes","No"), labels=c("Event Occurred", "No Event"))
toy$out2 <- as.numeric(toy$out2.event)-1 # subtracting 1 at the end changes the default 1/2 code to 0/1

## Sanity Checks
table(toy$treated.f, toy$treated)
table(toy$covB.f, toy$covB)
table(toy$out2.f, toy$out2.event)
table(toy$out2, toy$out2.event)
table(toy$out2, toy$out2.f)

## Re-expressing the Multi-Categorical Variable
toy$covF.Low <- as.numeric(toy$covF=="1-Low")
toy$covF.Middle <- as.numeric(toy$covF=="2-Middle")
toy$covF.High <- as.numeric(toy$covF=="3-High")

## Sanity Checks
table(toy$covF, toy$covF.Low)
table(toy$covF, toy$covF.Middle)
table(toy$covF, toy$covF.High)

## We have three transformations to execute for the covariates
## A squared, plus B-C and B-D interactions
## Must use cov B in a numeric (0,1) form to build product terms with C and D
toy$Asqr <- toy$covA^2
toy$BC <- toy$covB*toy$covC
toy$BD <- toy$covB*toy$covD

## Now, toy should contain 200 observations, 21 variables
names(toy)
dim(toy)

## Question A. Ignoring covariate information, what is the estimated effect
## ... on Outcome 1 [a continuous outcome]
by(toy$out1.cost, toy$treated.f, summary)
unadj.out1 <- lm(out1.cost ~ treated, data=toy)
summary(unadj.out1); confint(unadj.out1) ## provides treated effect and confidence interval estimates

## ... on Outcome 2 [a binary outcome]
table(toy$treated.f, toy$out2.f)
library(Epi)
twoby2(table(toy$treated.f, toy$out2.f)) ## provides risk difference and odds ratio estimates

## To get a logistic regression model version of this, we use
unadj.out2 <- glm(out2 ~ treated, data=toy, family=binomial())
summary(unadj.out2)
exp(coef(unadj.out2)) # produces odds ratio estimate
exp(confint(unadj.out2)) # produced 95% CI for odds ratio

## ... on Outcome 3 [a time-to-event outcome]
## Patients with out2.event=No are right-censored, those with out2.event=Yes have their times to event observed
## Fit a simple unadjusted Cox proportional hazards model
## predicting time to event (with event=Yes indicating non-censored cases) based on treatment group (treated)
library(survival)
unadj.out3 <- coxph(Surv(out3.time, out2.event=="Yes") ~ treated, data=toy)
summary(unadj.out3) ## exp(coef) section indicates relative risk estimate and 95% CI
## Check proportional hazards assumption, at least a little bit - would like this p value to be non-significant
cox.zph(unadj.out3)
plot(cox.zph(unadj.out3), var="treated")

## Question B. Fitting the Propensity Score Model, then plotting it simply
psmodel <- glm(treated ~ covA + covB + covC + covD + covE + covF + Asqr + BC + BD, family=binomial(), data=toy)
toy$ps <- psmodel$fitted
toy$linps <- psmodel$linear.predictors

by(toy$ps, toy$treated.f, summary)

plot(toy$ps ~ toy$treated.f, ylab="Propensity for Treatment", xlab="")

## for a fancier plot to compare the distributions of PS across treatment groups, we could use
library(sm)
sm.density.compare(toy$ps, toy$treated.f, xlab="Propensity for Treatment", main="Propensity Score Comparison", xlim=c(0,1), col=c("red", "dark green"), lty=1)
legend("topleft", legend=levels(toy$treated.f), lty=c(1,2), lwd=1, col=c("red","dark green"), text.col=c("red","dark green"), bty="n")

## The approach above automatically places the legend on the top left of the plot. To place the legend
## yourself, comment out the legend line above, in favor of ...
## legend(locator(1), legend=levels(toy$treated.f), lty=c(1,2), lwd=1, col=c("red","dark green"), text.col=c("red","dark green"), bty="n")

## Checking Rubin's Rules Before Propensity Score Adjustment

## Rubin's Rule 1 - calculate the (absolute value of the) standardized difference of the linear propensity score
## We want this value to be close to 0, and certainly less than 50 in order to push forward to outcomes analysis
## Before Propensity Adjustment, Rubin's Rule 1 summary is...
rubin1.unadj <- with(toy, abs(100*(mean(linps[treated==1])-mean(linps[treated==0]))/sd(linps)))
rubin1.unadj

## Rubin's Rule 2 - calculate the ratio of variances of the linear propensity score comparing treated to control
## We want this value to be near 1, between 4/5 and 5/4, ideally, but certainly between 1/2 and 2.

## Prior to Propensity Score Adjustment, Rubin's Rule 2 summary is ...
rubin2.unadj <-with(toy, var(linps[treated==1])/var(linps[treated==0]))
rubin2.unadj

## Rubin's Rule 3 - ratio of variances of regression residuals for each covariate included in the propensity model
## comparing treated to control - again looking for near 1, between 4/5 and 5/4, ideally, definitely between 1/2 and 2.

## General function rubin3 to help calculate Rubin's Rule 3
rubin3 <- function(data, covlist, linps) {
  covlist2 <- as.matrix(covlist)
  res <- NA
  for(i in 1:ncol(covlist2)) {
    cov <- as.numeric(covlist2[,i])
    num <- var(resid(lm(cov ~ data$linps))[data$treated==1])
    den <- var(resid(lm(cov ~ data$linps))[data$treated==0])
    res[i] <- round(num/den, 3)
  }
  names(res) <- names(covlist)
  print(res)
}

## Prior to propensity adjustment, Rubin's Rule 3 summaries are:
rubin3.unadj <- rubin3(data=toy, covlist=toy[c("covA", "covB", "covC", "covD", "covE", "covF.Middle", "covF.High", "Asqr","BC", "BD")])

## Building a dotplot of the Rubin's Rule 3 Variance Ratios:
d0 <- rubin3(data=toy, covlist=toy[c("covA", "covB", "covC", "covD", "covE", "covF.Middle", "covF.High", "Asqr","BC", "BD")])
d <- sort(d0)
low <- min(min(d), 0.45)
high <- max(max(d), 2.05)
  
dotchart(d, pch=15, col="black", main="Rubin's Rule 3 Results (Unadjusted)", xlab="Residual Variance Ratio", xlim=c(low, high))
 abline(v=1, lty=1)
 abline(v=0.8, lty=2, lwd=2, col="blue")
 abline(v=1.25, lty=2, lwd=2, col="blue")
 abline(v=0.5, lty=2, lwd=2, col="red")
 abline(v=2, lty=2, lwd=2, col="red")

rm(d, d0, low, high)

## Question C. Use 1:1 greedy matching to match all 70 treated to 70 unique control patients
## on the linear propensity scores. We'll break ties at random, as well.
## Then check balance (and plot standardized differences and variance ratios) appropriately
## including the raw and linear propensity scores in plots (usually unnecessary to show both)
library(Matching)
X <- toy$linps ## matching on the linear propensity score
Tr <- as.logical(toy$treated)
match1 <- Match(Tr=Tr, X=X, M = 1, replace=FALSE, ties=FALSE)
summary(match1)

## Next, we need to assess balance imposed by match1 on covariates (and PS and the linear PS).
mb1 <- MatchBalance(treated ~ covA + covB + covC + covD + covE + covF + Asqr + BC + BD + ps + linps, data=toy, match.out = match1, nboots=500)
covnames <- c("covA", "covB", "covC", "covD", "covE", "covF - Middle", "covF - High", "A^2","B*C", "B*D", "raw PS", "linear PS")

## Extract standardized differences
pre.szd <- NULL; post.szd <- NULL
for(i in 1:length(covnames)) {
  pre.szd[i] <- mb1$BeforeMatching[[i]]$sdiff.pooled
  post.szd[i] <- mb1$AfterMatching[[i]]$sdiff.pooled
}

## Basic Standardized Difference Table
temp <- data.frame(pre.szd, post.szd, row.names=covnames)
print(temp, digits=3)

## Absolute Standardized Difference Plot
temp <- data.frame(pre.szd, post.szd, row.names=covnames)
tempsort <- temp[with(temp, order(abs(pre.szd))),]
high <- max(max(abs(pre.szd)), max(abs(post.szd)), 0.1)

dotchart(abs(tempsort$pre.szd), pch="", xlim=c(0, 1.05*high), labels=row.names(tempsort), main="Absolute Standardized Difference Plot", xlab="Absolute Standardized Difference (%)")
points(abs(tempsort$pre.szd), seq(1:length(tempsort$pre.szd)), pch=15, col="blue", cex=1.2)
points(abs(tempsort$post.szd), seq(1:length(tempsort$post.szd)), pch=19, col="red", cex=1.2)
abline(v=0, lty=1)
abline(v=10, lty=2, col="purple")
legend("bottomright", legend = c("Before Matching", "After Matching"), col=c("blue", "red"), text.col=c("blue", "red"), bty="o", pch = c(15, 19))

## Or, if you prefer, plotting the standardized differences (not absolute values)
temp <- data.frame(pre.szd, post.szd, row.names=covnames)
tempsort <- temp[with(temp, order(pre.szd)), ]
low <- min(min(pre.szd), min(post.szd), -0.1)
high <- max(max(pre.szd), max(post.szd), 0.1)

dotchart(tempsort$pre.szd, xlim=c(1.05*low, 1.05*high), pch="", labels=row.names(tempsort), main="Standardized Difference Plot", xlab="Standardized Difference (%)")
points(tempsort$pre.szd, seq(1:length(tempsort$pre.szd)), pch=15, col="blue", cex=1.2)
points(tempsort$post.szd, seq(1:length(tempsort$post.szd)), pch=19, col="red", cex=1.2)
abline(v=0, lty=1)
abline(v=10, lty=2, col="purple")
abline(v=-10, lty=2, col="purple")
legend("bottomright", legend = c("Before Matching", "After Matching"), col=c("blue", "red"), text.col=c("blue", "red"), bty="o", pch = c(15, 19))

## Extract variance ratios
pre.vratio <- NULL; post.vratio <- NULL
for(i in 1:length(covnames)) {
  pre.vratio[i] <- mb1$BeforeMatching[[i]]$var.ratio
  post.vratio[i] <- mb1$AfterMatching[[i]]$var.ratio
}

## Table of Variance Ratios
temp <- data.frame(pre.vratio, post.vratio, row.names=covnames)
print(temp, digits=2)

## Variance Ratio Plot
temp <- data.frame(pre.vratio, post.vratio, row.names=covnames)
tempsort <- temp[with(temp, order(pre.vratio)), ]
low <- min(min(pre.vratio), min(post.vratio))
high <- max(max(pre.vratio), max(post.vratio))

dotchart(tempsort$pre.vratio, xlim=c(0.95*low, 1.05*high), pch="", labels=row.names(tempsort), main="Plot of Variance Ratios", xlab="Treatment Variance / Control Variance")
points(tempsort$pre.vratio, seq(1:length(tempsort$pre.vratio)), pch=8, col="black", cex=1.2)
points(tempsort$post.vratio, seq(1:length(tempsort$post.vratio)), pch=7, col="magenta", cex=1.2)
abline(v=1, lty=1)
abline(v=3/4, lty=2, col="brown")
abline(v=4/3, lty=2, col="brown")
legend("topleft", legend = c("Before Matching", "After Matching"), col=c("black", "magenta"), text.col=c("black", "magenta"), bty="o", pch = c(8, 7))

## Finally, we'll create a new data frame, containing only the matched sample
matches <- factor(rep(match1$index.treated, 2))
toy.matchedsample <- cbind(matches, toy[c(match1$index.control, match1$index.treated),])

## Sanity Check
table(toy.matchedsample$treated.f) ## should be 70 treated and 70 control patients
head(toy.matchedsample)

## Checking Rubin's Rules
## Rubin's Rule 1 - calculate the (absolute value of the) standardized difference of the linear propensity score
## We want this value to be close to 0, and certainly less than 50 in order to push forward to outcomes analysis

## As a reminder, Prior to Matching, Rubin's Rule 1 summary is...
rubin1.unadj <- with(toy, abs(100*(mean(linps[treated==1])-mean(linps[treated==0]))/sd(linps)))
rubin1.unadj

## After Matching, same summary is ...
rubin1.match <- with(toy.matchedsample, abs(100*(mean(linps[treated==1])-mean(linps[treated==0]))/sd(linps)))
rubin1.match

## Rubin's Rule 2 - calculate the ratio of variances of the linear propensity score comparing treated to control
## We want this value to be near 1, between 4/5 and 5/4, ideally, but certainly between 1/2 and 2.

## Again, Prior to Matching, Rubin's Rule 2 summary is ...
rubin2.unadj <-with(toy, var(linps[treated==1])/var(linps[treated==0]))
rubin2.unadj

## After Matching, Rubin's Rule 2 summary is ...
rubin2.match <- with(toy.matchedsample, var(linps[treated==1])/var(linps[treated==0]))
rubin2.match

## Rubin's Rule 3 - ratio of variances of regression residuals for each covariate included in the propensity model
## comparing treated to control - again looking for near 1, between 4/5 and 5/4, ideally, definitely between 1/2 and 2.

## Prior to matching, Rubin's Rule 3 summaries are:
rubin3.unadj <- rubin3(data=toy, covlist=toy[c("covA", "covB", "covC", "covD", "covE", "covF.Middle", "covF.High", "Asqr","BC", "BD")])

## After matching, then, Rubin's Rule 3 summaries are:
rubin3.match <- rubin3(data=toy.matchedsample, covlist=toy.matchedsample[c("covA", "covB", "covC", "covD", "covE", "covF.Middle", "covF.High", "Asqr","BC", "BD")])

## Building a dotplot of the Rubin's Rule 3 Variance Ratios Pre- and Post-Match:
d.unadj <- rubin3(data=toy, covlist=toy[c("covA", "covB", "covC", "covD", "covE", "covF.Middle", "covF.High", "Asqr","BC", "BD")])
d.match <- rubin3(data=toy.matchedsample, covlist=toy.matchedsample[c("covA", "covB", "covC", "covD", "covE", "covF.Middle", "covF.High", "Asqr","BC", "BD")])

low <- min(d.unadj, d.match, 0.45)
high <- max(d.unadj, d.match, 2.05)

dotchart(d.unadj, pch=15, col="black", main="Rubin's Rule 3 Results for TOY Example", xlab="Residual Variance Ratio", xlim=c(low, high))
points(d.match, seq(1:length(d.match)), pch=19, col="purple", cex=1.2)
abline(v=1, lty=1)
abline(v=0.8, lty=2, lwd=2, col="blue")
abline(v=1.25, lty=2, lwd=2, col="blue")
abline(v=0.5, lty=2, lwd=2, col="red")
abline(v=2, lty=2, lwd=2, col="red")
legend("topright", legend=c("Before Matching", "After Matching"), col=c("black", "purple"), text.col=c("black", "purple"), bty="o", pch = c(15, 19))

## Question D. After matching, what is the estimated average causal effect of treatment?

## ... on Outcome 1 [a continuous outcome]
## First, we'll look at the automatic answer which can be obtained when matching
X <- toy$linps ## matching on the linear propensity score
Tr <- as.logical(toy$treated)
Y <- toy$out1.cost
match1 <- Match(Y=Y, Tr=Tr, X=X, M = 1, replace=FALSE, ties=FALSE)
summary(match1)

## As a brief interlude, we can also get an ATE estimate (as compared to the ATT above)
match1.ATE <- Match(Y=Y, Tr=Tr, X=X, M = 1, replace=FALSE, ties=FALSE, estimand="ATE")
summary(match1.ATE)
rm(match1.ATE)

## Next, we'll compare the continuous outcome between treated groups with regression on the matched sample
by(toy.matchedsample$out1.cost, toy.matchedsample$treated.f, summary)  ## description doesn't account for paired samples
t.test(toy.matchedsample$out1.cost~toy.matchedsample$treated.f, paired=TRUE) ## accounts for pairing

adj.m.out1 <- lm(out1.cost ~ treated + factor(matches), data=toy.matchedsample) 
## This takes the pairing into account, but treats pairing as a fixed, rather than random, factor
## and thus matches the paired t test but isn't totally satisfactory as a solution
summary(adj.m.out1); confint(adj.m.out1) 
## provides treated effect and overly wide confidence interval estimates, but too many of them
coef(adj.m.out1)["treated"] # point estimate for treated effect
confint(adj.m.out1)["treated",1] # lower limit of 95% CI
confint(adj.m.out1)["treated",2] # lower limit of 95% CI

## A more appropriate result comes from a mixed model where the matches are treated as a random factor
## but the treatment group is treated as a fixed factor
library(lme4)
toy.matchedsample$matches.f <- as.factor(toy.matchedsample$matches) ## Need to use matches as a factor in R here
mixedmodel.out1 <- lmer(out1.cost ~ treated + (1 | matches.f), data=toy.matchedsample)
summary(mixedmodel.out1); confint(mixedmodel.out1)

## ... on Outcome 2 [a binary outcome]
## First, we'll look at the automatic ATT answer which can be obtained when matching
X <- toy$linps ## matching on the linear propensity score
Tr <- as.logical(toy$treated)
Y <- toy$out2
match1 <- Match(Y=Y, Tr=Tr, X=X, M = 1, replace=FALSE, ties=FALSE)
summary(match1)

## As a brief interlude, we can also get an ATE estimate (as compared to the ATT above)
match1.ATE <- Match(Y=Y, Tr=Tr, X=X, M = 1, replace=FALSE, ties=FALSE, estimand="ATE")
summary(match1.ATE)
rm(match1.ATE)

## Since we have the matched sample, we can either do conditional logistic regression, accounting for the matched pairs
library(survival)
adj.m.out2 <- clogit(out2 ~ treated + strata(matches), data=toy.matchedsample)
summary(adj.m.out2)

## Or a McNemar test to obtain the odds ratio
library(reshape2)
wide.out2 <- dcast(toy.matchedsample, matches ~ treated.f, value.var="out2.f")
wide.out2$Treated <- factor(wide.out2$Treated, levels=c("Event Occurred","No Event"), labels=c("Treated, Event", "Treated, No Event"))
wide.out2$Control <- factor(wide.out2$Control, levels=c("Event Occurred","No Event"), labels=c("Ctrl, Event", "Ctrl, No Event"))
table(wide.out2$Treated, wide.out2$Control)
library(epibasix)
mcNemar(as.matrix(table(wide.out2$Treated, wide.out2$Control)))
## But this won't give an answer since we have < 30 discordant pairs.
## To force it to show us the answer anyway, we can use
mcNemar(as.matrix(table(wide.out2$Treated, wide.out2$Control)), force=TRUE)

## ... on Outcome 3 [a time to event]
## First, we'll again look at the automatic answer which can be obtained when matching
## note that this doesn't take into account the right censoring at all...
X <- toy$linps ## matching on the linear propensity score
Tr <- as.logical(toy$treated)
Y <- toy$out3.time
match1 <- Match(Y=Y, Tr=Tr, X=X, M = 1, replace=FALSE, ties=FALSE)
summary(match1)

## Since we have the matched sample, we should use a stratified Cox proportional hazards model to compare
## the treatment groups on our time-to-event outcome, while accounting for the matched pairs
## Results will be a relative hazard rate
library(survival)
adj.m.out3 <- coxph(Surv(out3.time, out2) ~ treated + strata(matches), data=toy.matchedsample)
summary(adj.m.out3)
cox.zph(adj.m.out3) # Quick check for proportional hazards assumption
plot(cox.zph(adj.m.out3), var="treated")

## OK - that's it for the matching, I'm going to clean up all of the loose variables I've generated so far here
## except for the propensity model psmodel, which I'll need again.
rm(temp, tempsort, toy.matchedsample, wide.out2)
rm(Tr, X, Y, adj.m.out1, adj.m.out2, adj.m.out3, covnames, high, i, low, match1)
rm(matches, mb1, mixedmodel.out1, post.szd, post.vratio, pre.szd, pre.vratio)
rm(unadj.out1, unadj.out2, unadj.out3)
rm(rubin1.match, rubin2.match, rubin3.match, d.match, d.unadj)

## Question E. Execute subclassification by quintile of the linear propensity score, then assess covariate balance

library(Hmisc)
toy$stratum <- cut2(toy$ps, g=5)
toy$quintile <- factor(toy$stratum, labels=1:5)
table(toy$stratum, toy$quintile) ## quick sanity check

## Quick numerical summaries of propensity scores by subclass...
tapply(toy$ps, toy$stratum, summary)

## Divide the sample into the five quintiles
quin1 <- subset(toy, quintile==1)
quin2 <- subset(toy, quintile==2)
quin3 <- subset(toy, quintile==3)
quin4 <- subset(toy, quintile==4)
quin5 <- subset(toy, quintile==5)

## We want to check the balance and propensity score overlap for each quintile.
## I'll start with a set of strip charts to look at overlap
par(mfrow=c(1,5))
stripchart(round(ps,2) ~ treated.f, data=quin1, main="Q1: Low PS", method="stack", cex=1, offset=1/2, col=c("blue", "red"), xlab="Propensity Score")
stripchart(round(ps,2) ~ treated.f, data=quin2, main="Q2: Moderate-Low PS", method="stack", cex=1, offset=1/2, col=c("blue", "red"), xlab="Propensity Score")
stripchart(round(ps,2) ~ treated.f, data=quin3, main="Q3: Moderate PS", method="stack", cex=1, offset=1/2, col=c("blue", "red"), xlab="Propensity Score")
stripchart(round(ps,2) ~ treated.f, data=quin4, main="Q4: Moderate-High PS", method="stack", cex=1, offset=1/2, col=c("blue", "red"), xlab="Propensity Score")
stripchart(round(ps,2) ~ treated.f, data=quin5, main="Q5: High PS", method="stack", cex=1, offset=1/2, col=c("blue", "red"), xlab="Propensity Score")
par(mfrow=c(1,1))

addmargins(table(toy$quintile, toy$treated.f))

## We'll need to be able to calculate standardized differences in this situation
## so I've created a szd function to do this - using the average denominator method 
szd <- function(covlist, g) {
  covlist2 <- as.matrix(covlist)
  g <- as.factor(g)
  res <- NA
  for(i in 1:ncol(covlist2)) {
    cov <- as.numeric(covlist2[,i])
    num <- 100*diff(tapply(cov, g, mean, na.rm=TRUE))
    den <- sqrt(mean(tapply(cov, g, var, na.rm=TRUE)))
    res[i] <- round(num/den,2)
  }
  names(res) <- names(covlist)   
  res
}

## Now get the standardized differences within each quintile, as well as overall
covs <- c("covA", "covB", "covC", "covD", "covE", "covF.Middle", "covF.High", "Asqr","BC", "BD", "ps", "linps")
d.q1 <- szd(quin1[covs], quin1$treated)
d.q2 <- szd(quin2[covs], quin2$treated)
d.q3 <- szd(quin3[covs], quin3$treated)
d.q4 <- szd(quin4[covs], quin4$treated)
d.q5 <- szd(quin5[covs], quin5$treated)
d.all <- szd(toy[covs], toy$treated)

d.q1
d.q2
d.q3
d.q4
d.q5
d.all

## One way to plot the data in a semi-useful way -> showing results of standardized difference plots by quintile and overall
library(lattice)
temp <- matrix(c(abs(d.all), abs(d.q1), abs(d.q2), abs(d.q3), abs(d.q4), abs(d.q5)), ncol=6)
dimnames(temp) <- list(names(d.q1), c("All Subjects", "Quintile 1", "Quintile 2", "Quintile 3", "Quintile 4", "Quintile 5"))
dotplot(temp, groups=F, xlab="Absolute Value of Standardized Differences within Quintiles", scales = list(alternating=3), as.table=TRUE)

## Another plotting approach that I like, but it gets large quickly when you have lots of covariates, is...
dotplot(t(temp), groups=F, xlab="Absolute Value of Standardized Differences within Quintiles", scales = list(alternating=3), as.table=TRUE)

## Checking Variance Ratios after Propensity Score Quintile Stratification
## This is a silly, toy example-specific function
vratio.toy <- function(data)
{
  covA.vr <- with(data, var(covA[treated==1])/var(covA[treated==0]))
  covB.vr <- with(data, var(covB[treated==1])/var(covB[treated==0]))
  covC.vr <- with(data, var(covC[treated==1])/var(covC[treated==0]))
  covD.vr <- with(data, var(covD[treated==1])/var(covD[treated==0]))
  covE.vr <- with(data, var(covE[treated==1])/var(covE[treated==0]))
  covF.Middle.vr <- with(data, var(covF.Middle[treated==1])/var(covF.Middle[treated==0]))
  covF.High.vr <- with(data, var(covF.High[treated==1])/var(covF.High[treated==0]))
  Asqr.vr <- with(data, var(Asqr[treated==1])/var(Asqr[treated==0]))
  BC.vr <- with(data, var(BC[treated==1])/var(BC[treated==0]))
  BD.vr <- with(data, var(BD[treated==1])/var(BD[treated==0]))
  ps.vr <- with(data, var(ps[treated==1])/var(ps[treated==0]))
  linps.vr <- with(data, var(linps[treated==1])/var(linps[treated==0]))
  res.vr <- round(c(covA.vr, covB.vr, covC.vr, covD.vr, covE.vr, covF.Middle.vr, covF.High.vr, Asqr.vr, BC.vr, BD.vr, ps.vr, linps.vr),3)
  names(res.vr) <- c("A", "B", "C", "D", "E", "F-middle", "F-high", "A^2", "BC", "BD", "ps", "linear ps")
  res.vr
  }

vratio.toy(toy)
vratio.toy(quin1)
vratio.toy(quin2)
vratio.toy(quin3)
vratio.toy(quin4)
vratio.toy(quin5)

rm(vratio.toy)

## Checking Rubin's Rules After Propensity Score Quintile Stratification
## As a reminder, Prior to Adjustment, Rubin's Rule 1 summary is...
rubin1.unadj <- with(toy, abs(100*(mean(linps[treated==1])-mean(linps[treated==0]))/sd(linps)))
rubin1.unadj

## After Subclassification, we can obtain the same summary within each of the five quintiles ...
rubin1.q1 <- with(quin1, abs(100*(mean(linps[treated==1])-mean(linps[treated==0]))/sd(linps)))
rubin1.q2 <- with(quin2, abs(100*(mean(linps[treated==1])-mean(linps[treated==0]))/sd(linps)))
rubin1.q3 <- with(quin3, abs(100*(mean(linps[treated==1])-mean(linps[treated==0]))/sd(linps)))
rubin1.q4 <- with(quin4, abs(100*(mean(linps[treated==1])-mean(linps[treated==0]))/sd(linps)))
rubin1.q5 <- with(quin5, abs(100*(mean(linps[treated==1])-mean(linps[treated==0]))/sd(linps)))
rubin1.sub <- c(rubin1.q1, rubin1.q2, rubin1.q3, rubin1.q4, rubin1.q5)
names(rubin1.sub)=c("Q1", "Q2", "Q3", "Q4", "Q5")

rubin1.sub

## As a reminder, Prior to Adjustment, Rubin's Rule 2 summary is...
rubin2.unadj <- with(toy, var(linps[treated==1])/var(linps[treated==0]))
rubin2.unadj

## After Subclassification, we can obtain the same summary within each of the five quintiles ...
rubin2.q1 <- with(quin1, var(linps[treated==1])/var(linps[treated==0]))
rubin2.q2 <- with(quin2, var(linps[treated==1])/var(linps[treated==0]))
rubin2.q3 <- with(quin3, var(linps[treated==1])/var(linps[treated==0]))
rubin2.q4 <- with(quin4, var(linps[treated==1])/var(linps[treated==0]))
rubin2.q5 <- with(quin5, var(linps[treated==1])/var(linps[treated==0]))
rubin2.sub <- c(rubin2.q1, rubin2.q2, rubin2.q3, rubin2.q4, rubin2.q5)
names(rubin2.sub)=c("Q1", "Q2", "Q3", "Q4", "Q5")

rubin2.sub

## Prior to propensity adjustment, recall that Rubin's Rule 3 summaries are:
covs <- c("covA", "covB", "covC", "covD", "covE", "covF.Middle", "covF.High", "Asqr","BC", "BD")
rubin3.unadj <- rubin3(data=toy, covlist=toy[covs])

## After subclassification, then, Rubin's Rule 3 summaries are:
rubin3.q1 <- rubin3(data=quin1, covlist=quin1[covs])
rubin3.q2 <- rubin3(data=quin2, covlist=quin2[covs])
rubin3.q3 <- rubin3(data=quin3, covlist=quin3[covs])
rubin3.q4 <- rubin3(data=quin4, covlist=quin4[covs])
rubin3.q5 <- rubin3(data=quin5, covlist=quin5[covs])

library(lattice)
temp <- matrix(c(rubin3.unadj, rubin3.q1, rubin3.q2, rubin3.q3, rubin3.q4, rubin3.q5), ncol=6)
names(rubin3.unadj) <- c("A", "B", "C", "D", "E", "F.Middle", "F.High", "Asqr","BC", "BD")
dimnames(temp) <- list(names(rubin3.unadj), c("All Subjects", "Quintile 1", "Quintile 2", "Quintile 3", "Quintile 4", "Quintile 5"))
dotplot(temp, groups=F, xlab="Residual Variance Ratios Within Quintiles", scales = list(alternating=3), as.table=TRUE)

## Question F. After subclassifying, what is the estimated average average causal effect of treatment?

## ... on Outcome 1 [a continuous outcome]
## First, find the estimated average causal effect (and standard error) within each quintile via linear regression
quin1.out1 <- lm(out1.cost ~ treated, data=quin1)
quin2.out1 <- lm(out1.cost ~ treated, data=quin2)
quin3.out1 <- lm(out1.cost ~ treated, data=quin3)
quin4.out1 <- lm(out1.cost ~ treated, data=quin4)
quin5.out1 <- lm(out1.cost ~ treated, data=quin5)

coef(summary(quin1.out1)); coef(summary(quin2.out1)); coef(summary(quin3.out1)); coef(summary(quin4.out1)); coef(summary(quin5.out1))

## Next, we find the mean of the five quintile-specific estimated regression coefficients
est.st <- (coef(quin1.out1)[2] + coef(quin2.out1)[2] + coef(quin3.out1)[2] + coef(quin4.out1)[2] + coef(quin5.out1)[2])/5
est.st

## To get the combined standard error estimate, we do the following
se.q1 <- summary(quin1.out1)$coefficients[2,2]
se.q2 <- summary(quin2.out1)$coefficients[2,2]
se.q3 <- summary(quin3.out1)$coefficients[2,2]
se.q4 <- summary(quin4.out1)$coefficients[2,2]
se.q5 <- summary(quin5.out1)$coefficients[2,2]

se.st <- sqrt((se.q1^2 + se.q2^2 + se.q3^2 + se.q4^2 + se.q5^2)*(1/25))

## Resulting 95% Confidence Interval for Average Causal Effect is ...
temp.result1 <- c(est.st, est.st - 1.96*se.st, est.st + 1.96*se.st)
names(temp.result1) <- c("Estimate", "Low 95% CI", "High 95% CI")
temp.result1

## ... on Outcome 2 [a binary outcome]
## First, find the estimated average causal effect (and standard error) within each quintile via logistic regression
quin1.out2 <- glm(out2 ~ treated, data=quin1, family=binomial())
quin2.out2 <- glm(out2 ~ treated, data=quin2, family=binomial())
quin3.out2 <- glm(out2 ~ treated, data=quin3, family=binomial())
quin4.out2 <- glm(out2 ~ treated, data=quin4, family=binomial())
quin5.out2 <- glm(out2 ~ treated, data=quin5, family=binomial())

coef(summary(quin1.out2)); coef(summary(quin2.out2)); coef(summary(quin3.out2)); coef(summary(quin4.out2)); coef(summary(quin5.out2))

## Next, we find the mean of the five quintile-specific estimated logistic regression coefficients
est.st <- (coef(quin1.out2)[2] + coef(quin2.out2)[2] + coef(quin3.out2)[2] + coef(quin4.out2)[2] + coef(quin5.out2)[2])/5
est.st ## this is the estimated log odds ratio
## And we exponentiate this to get the overall odds ratio estimate
exp(est.st)

## To get the combined standard error estimate, we do the following
se.q1 <- summary(quin1.out2)$coefficients[2,2]
se.q2 <- summary(quin2.out2)$coefficients[2,2]
se.q3 <- summary(quin3.out2)$coefficients[2,2]
se.q4 <- summary(quin4.out2)$coefficients[2,2]
se.q5 <- summary(quin5.out2)$coefficients[2,2]
se.st <- sqrt((se.q1^2 + se.q2^2 + se.q3^2 + se.q4^2 + se.q5^2)*(1/25))
## Of course, this standard error is also on the log odds ratio scale

## Resulting 95% Confidence Interval for Average Causal Effect (as an Odds Ratio) is ...
temp.result2 <- c(exp(est.st), exp(est.st - 1.96*se.st), exp(est.st + 1.96*se.st))
names(temp.result2) <- c("Estimate", "Low 95% CI", "High 95% CI")
temp.result2

## Alternatively, we can get the twoby2 results from the Epi library
library(Epi)
twoby2(table(quin1$treated.f, quin1$out2.f))
twoby2(table(quin2$treated.f, quin2$out2.f))
twoby2(table(quin3$treated.f, quin3$out2.f))
twoby2(table(quin4$treated.f, quin4$out2.f))
twoby2(table(quin5$treated.f, quin5$out2.f))

## ... on Outcome 3 [a time to event]
## Patients with out2.event = "Yes" are truly observed events, while those with 
## out2.event == "No" are censored before an event can happen to them
## The Cox model comparing treated to control, stratifying on quintile, is...
library(survival)
adj.s.out3 <- coxph(Surv(out3.time, out2) ~ treated + strata(quintile), data=toy)
summary(adj.s.out3) ## exp(coef) gives relative hazard associated with treatment
cox.zph(adj.s.out3) ## checking the proportional hazards assumption
plot(cox.zph(adj.s.out3), var="treated")

## Question G. Execute weighting by the inverse propensity score, then assess covariate balance
## ATT approach: Weight treated patients as 1; control patients as ps/(1-ps)
toy$wts1 <- ifelse(toy$treated==1, 1, toy$ps/(1-toy$ps))

## Two ways of building a better looking plot of the applied weights are 
library(lattice)
xyplot(wts1 ~ ps | treated.f, data=toy, ylab="ATT Weights", xlab="Propensity for Treatment", scales = list(alternating=3))
xyplot(wts1 ~ ps, groups=treated.f, data=toy, ylab="ATT Weights", xlab="Propensity for Treatment", auto.key=TRUE)

## A clever alternative way to do this ATT weighting is:
## toy$wts1a <- rep(1, nrow(toy))
## toy$wts1a[toy$treated==0] <- exp(predict(psmodel, subset(toy, treated==0)))

## ATE Approach: Weight treated patients by 1/ps; Control patients by 1/(1-PS)
toy$wts2 <- ifelse(toy$treated==1, 1/toy$ps, 1/(1-toy$ps))
plot(toy$wts2 ~ toy$ps)
library(lattice)
xyplot(wts2 ~ ps | treated.f, data=toy, ylab="ATE Weights", xlab="Propensity for Treatment", scales = list(alternating=3))
xyplot(wts2 ~ ps, groups=treated.f, data=toy, ylab="ATE Weights", xlab="Propensity for Treatment", auto.key=TRUE)

## Using the twang library to assess balance after weighting via a table
library(twang)
covlist <- c("covA", "covB", "covC", "covD", "covE", "covF", "Asqr","BC", "BD", "ps", "linps")

# for ATT weights
bal.wts1 <- dx.wts(x=toy$wts1, data=toy, vars=covlist, treat.var="treated", estimand="ATT")
bal.wts1
bal.table(bal.wts1)

# for ATE weights
bal.wts2 <- dx.wts(x=toy$wts2, data=toy, vars=covlist, treat.var="treated", estimand="ATE")
bal.wts2
bal.table(bal.wts2)

## Question H. After weighting, what is the estimated average average causal effect of treatment?
## ... on Outcome 1 [a continuous outcome]
library(survey)
toywt1.design <- svydesign(ids=~1, weights=~wts1, data=toy) # using ATT weights
toywt2.design <- svydesign(ids=~1, weights=~wts2, data=toy) # using ATE weights

adjout1.wt1 <- svyglm(out1.cost ~ treated, design=toywt1.design)
summary(adjout1.wt1); confint(adjout1.wt1)

adjout1.wt2 <- svyglm(out1.cost ~ treated, design=toywt2.design)
summary(adjout1.wt2); confint(adjout1.wt2)

## ... on Outcome 2 [a binary outcome]
library(survey)
toywt1.design <- svydesign(ids=~1, weights=~wts1, data=toy) # using ATT weights
toywt2.design <- svydesign(ids=~1, weights=~wts2, data=toy) # using ATE weights

## Build outcome model using quasibinomial, rather than binomial family
## first using ATT weights
adjout2.wt1 <- svyglm(out2 ~ treated, design=toywt1.design, family=quasibinomial())
summary(adjout2.wt1)
exp(summary(adjout2.wt1)$coef)
exp(confint(adjout2.wt1))

## then using ATE Weights
adjout2.wt2 <- svyglm(out2.event ~ treated, design=toywt2.design, family=quasibinomial())
summary(adjout2.wt2)
exp(summary(adjout2.wt2)$coef)
exp(confint(adjout2.wt2))

## ... on Outcome 3 [a time to event]
## Patients with out2.event = "Yes" are truly observed events, while those with 
## out2.event == "No" are censored before an event can happen to them
## The Cox model comparing treated to control, weighting by ATT weights (wt1), is...
library(survival)
adjout3.wt1 <- coxph(Surv(out3.time, out2) ~ treated, data=toy, weights=wts1)
summary(adjout3.wt1) ## exp(coef) output gives relative hazard of treated compared to control
cox.zph(adjout3.wt1); plot(cox.zph(adjout3.wt1), var="treated")

## And, weighting by ATE weights (wt2), we have...
library(survival)
adjout3.wt2 <- coxph(Surv(out3.time, out2) ~ treated, data=toy, weights=wts2)
summary(adjout3.wt2) ## exp(coef) output gives relative hazard of treated compared to control
cox.zph(adjout3.wt2); plot(cox.zph(adjout3.wt2), var="treated")

## Question I. After direct adjustment for the linear propensity score, what is the estimated average average causal effect of treatment?
## ... on Outcome 1 [a continuous outcome]
## Linear Regression model with linps added as a covariate
adj.reg.out1 <- lm(out1.cost ~ treated + linps, data=toy)
summary(adj.reg.out1); confint(adj.reg.out1) ## provides treated effect and confidence interval estimates

## ... on Outcome 2 [a binary outcome]
## Logistic Regression model with linps added as a covariate
adj.reg.out2 <- glm(out2 ~ treated + linps, data=toy, family=binomial())
summary(adj.reg.out2)
exp(coef(adj.reg.out2)) # produces odds ratio estimate
exp(confint(adj.reg.out2)) # produced 95% CI for odds ratio

## ... on Outcome 3 [a time-to-event outcome]
## Patients with out2.event=No are right-censored, those with out2.event=Yes have their times to event observed
## Fit a simple unadjusted Cox proportional hazards model
## predicting time to event (with event=Yes indicating non-censored cases) based on treatment group (treated)
## and now also the linear propensity score
library(survival)
adj.reg.out3 <- coxph(Surv(out3.time, out2) ~ treated + linps, data=toy)
summary(adj.reg.out3) ## exp(coef) section indicates relative risk estimate and 95% CI
## Check proportional hazards assumption, at least a little bit - would like this p value to be non-significant
cox.zph(adj.reg.out3)
plot(cox.zph(adj.reg.out3), var="treated")
plot(cox.zph(adj.reg.out3), var="linps")

## Question J. Double Robust Approach - Weighting + Adjustment, what is the estimated average average causal effect of treatment?
## ... on Outcome 1 [a continuous outcome]
library(survey)
toywt1.design <- svydesign(ids=~1, weights=~wts1, data=toy) # using ATT weights
toywt2.design <- svydesign(ids=~1, weights=~wts2, data=toy) # using ATE weights

adjout1.wt1 <- svyglm(out1.cost ~ treated + linps, design=toywt1.design)
summary(adjout1.wt1); confint(adjout1.wt1)

adjout1.wt2 <- svyglm(out1.cost ~ treated + linps, design=toywt2.design)
summary(adjout1.wt2); confint(adjout1.wt2)

## ... on Outcome 2 [a binary outcome]
library(survey)
toywt1.design <- svydesign(ids=~1, weights=~wts1, data=toy) # using ATT weights
toywt2.design <- svydesign(ids=~1, weights=~wts2, data=toy) # using ATE weights

## Build outcome model using quasibinomial, rather than binomial family
## first using ATT weights
adjout2.wt1 <- svyglm(out2 ~ treated + linps, design=toywt1.design, family=quasibinomial())
summary(adjout2.wt1)
exp(summary(adjout2.wt1)$coef)
exp(confint(adjout2.wt1))

## then using ATE Weights
adjout2.wt2 <- svyglm(out2.event ~ treated + linps, design=toywt2.design, family=quasibinomial())
summary(adjout2.wt2)
exp(summary(adjout2.wt2)$coef)
exp(confint(adjout2.wt2))

## ... on Outcome 3 [a time to event]
## Patients with out2.event = "Yes" are truly observed events, while those with 
## out2.event == "No" are censored before an event can happen to them
## The Cox model comparing treated to control, weighting by ATT weights (wt1), is...
library(survival)
adjout3.wt1 <- coxph(Surv(out3.time, out2) ~ treated + linps, data=toy, weights=wts1)
summary(adjout3.wt1) ## exp(coef) output gives relative hazard of treated compared to control
cox.zph(adjout3.wt1); plot(cox.zph(adjout3.wt1), var="treated")

## And, weighting by ATE weights (wt2), we have...
library(survival)
adjout3.wt2 <- coxph(Surv(out3.time, out2) ~ treated + linps, data=toy, weights=wts2)
summary(adjout3.wt2) ## exp(coef) output gives relative hazard of treated compared to control
cox.zph(adjout3.wt2); plot(cox.zph(adjout3.wt2), var="treated")

## Final version of toy in workspace has (200 observations, 27 variables)

## Final cleanup - leaving only the revised toy data set and the rubin3 and szd functions behind
rm(quin1, quin2, quin3, quin4, quin5, temp, covs, bal.wts1, bal.wts2, covlist)
rm(adj.reg.out1, adj.reg.out2, adj.reg.out3, adj.s.out3, adjout1.wt1, adjout1.wt2, adjout2.wt1, adjout2.wt2, adjout3.wt1, adjout3.wt2)
rm(d.all, d.q1, d.q2, d.q3, d.q4, d.q5, est.st, psmodel, quin1.out1, quin1.out2, quin2.out1, quin2.out2, quin3.out1, quin3.out2)
rm(quin4.out1, quin4.out2, quin5.out1, quin5.out2, se.q1, se.q2, se.q3, se.q4, se.q5, se.st)
rm(temp.result1, temp.result2, toywt1.design, toywt2.design)
rm(rubin1.q1, rubin1.q2, rubin1.q3, rubin1.q4, rubin1.q5, rubin1.sub, rubin1.unadj)
rm(rubin2.q1, rubin2.q2, rubin2.q3, rubin2.q4, rubin2.q5, rubin2.sub, rubin2.unadj)
rm(rubin3.q1, rubin3.q2, rubin3.q3, rubin3.q4, rubin3.q5, rubin3.unadj)
