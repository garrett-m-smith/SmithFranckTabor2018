# Analysis of Spring 2016 experiment "Decisions, decisions"
# Staub (2008, 2009)-style RSVP in which participants select either a singular or plural verb
# Predictors: subject type (Containments, Collections, Measures, Quantifiers), length (+/- PP)

library(lattice)
library(ggplot2)
library(lme4)
library(lsmeans)



d <- read.delim('./DecisionsData.txt')
d <- droplevels(subset(d, N1Type != 'training'))  # rm practice items
rownames(d) <- seq_len(nrow(d))
d$N1Type <- factor(d$N1Type, levels = c('containers', 'collections', 'measure', 'quant', 'filler'))
d$Subject <- as.factor(d$Subject)
d$Item <- as.factor(d$Item)
d$Slide1.RESP <- as.factor(d$Slide1.RESP)
levels(d$Slide1.RESP) <- c('Vsg', 'Vpl')

# basic plotting
d$LogRT <- log(d$Slide1.RT)
summary(d$Slide1.RT)
library(ggplot2)
ggplot(d[d$ItemType == 'test',], aes(x = Slide1.RT)) + geom_density() + xlim(0, 5000) + facet_grid(Length ~ N1Type)
ggplot(d[d$ItemType == 'test',], aes(N1Type, Slide1.RT, col = Length)) + stat_summary(fun.data = mean_cl_boot)
ggplot(d[d$ItemType == 'test',], aes(N1Type, LogRT, col = Length)) + stat_summary(fun.data = mean_cl_boot)

test <- droplevels(subset(d, N1Type != 'filler'))
library(lattice)
qqmath(~Slide1.RT | Subject, data = test)
qqmath(~LogRT | Subject, data = test)
nrow(test[test$Slide1.RT < 50 | test$Slide1.RT > 5000,]) / nrow(test)


# Rm'ing especially long RTs
library(plyr)
test.trimmed <- ddply(test, .(Subject), .fun = function(df) {
  return(droplevels(subset(df, LogRT < mean(LogRT) + 3 * sd(LogRT) &
                      LogRT > mean(LogRT) - 3 * sd(LogRT))))
}, .inform = T) # rm'ed 11 data points, about 0.5% of the data

# Recoding variables as numeric after Barr et al. (2013) and Levy (2014)
test.trimmed$LengthNum <- sapply(test.trimmed$Length, function(i) contr.sum(2)[i,])
N1TypeNum <- sapply(test.trimmed$N1Type, function(i) contr.sum(4)[i,])
test.trimmed$N1TypeNum1 <- N1TypeNum[1,]
test.trimmed$N1TypeNum2 <- N1TypeNum[2,]
test.trimmed$N1TypeNum3 <- N1TypeNum[3,]

# Looking at verb choice:
test.trimmed$ChoseSg <- ifelse(test.trimmed$Slide1.RESP == 'Vsg', 1, 0)
ggplot(test.trimmed, aes(x = N1Type, fill = Slide1.RESP)) + 
  geom_bar(position = 'dodge') + facet_grid(Length ~.)

my_summary <- function(x) {
  Mean <- mean(x)
  Upper <- Mean + 2 * sqrt((1 / length(x)) * Mean * (1 - Mean))
  Lower <- Mean - 2 * sqrt((1 / length(x)) * Mean * (1 - Mean))
  return(data.frame(y = Mean, ymax = Upper, ymin = Lower))
}

ggplot(test.trimmed, aes(x = N1Type, y = ChosePl, shape = Length)) + 
  stat_summary(fun.data = my_summary, size = 2, position = position_dodge(width = 0.2)) + 
  scale_y_continuous('Probability of choosing a plural verb', 
                     breaks = c(0.25, 0.5, 0.75, 1), limits = c(-0.005, 1.01)) +
  scale_x_discrete('Subject NP Type', 
                   labels = c('A box [of...]', 'A pile [of...]', 'A lot [of...]', 'Many [...]')) +
  scale_shape('N2 presence', labels = c('Present', 'Absent')) +
  theme_bw(base_size = 32) + theme(legend.position = c(0.8, 0.2))
ggsave('~/Dropbox/Structure and Flexibility (1)/Projects/Pseudopartitives/English/Productions/CUNY 2017/Talk/Figures/VerbChoice.pdf', height = 10, width = 12, units = 'in')


#### REPORTED: Re-done with vpl as reference for paper ####
test.trimmed$ChosePl <- ifelse(test.trimmed$ChoseSg == 1, 0, 1)
with(test.trimmed, table(N1Type, Length, ChosePl))

plchoice0 <- glmer(ChosePl ~ LengthNum + N1TypeNum1 + N1TypeNum2 + N1TypeNum3 +
                   LengthNum:N1TypeNum1 + LengthNum:N1TypeNum2 + LengthNum:N1TypeNum3 +
                   (1 + LengthNum + N1TypeNum1 + N1TypeNum2 + N1TypeNum3 +
                      LengthNum:N1TypeNum1 + LengthNum:N1TypeNum2 + LengthNum:N1TypeNum3 || Subject) +
                   (1 + LengthNum || Item), test.trimmed, family = binomial,
                   glmerControl(optCtrl = list(maxfun = 2e5)))
relgrad <- with(plchoice0@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
diag.vals <- getME(plchoice0, "theta")[getME(plchoice0, "lower") == 0]
any(diag.vals < 1e-6) # No singularity
# pars <- getME(plchoice0, c('theta', 'fixef'))
# plchoice0.restart <- update(plchoice0, start = pars)
# max(abs(with(plchoice0.restart@optinfo$derivs, solve(Hessian, gradient))))
# print(summary(plchoice0.restart), corr = F)

# Sigh...
# Recommended by ?convergence:
# source(system.file('utils', 'allFit.R', package = 'lme4'))
library(afex)
library(optimx)
plchoice.allFit <- allFit(plchoice0)
# summary(plchoice.allFit)
t(sapply(plchoice.allFit,fixef)) ## extract fixed effects
sapply(plchoice.allFit,logLik) ## log-likelihoods
sapply(plchoice.allFit,getME,"theta") ## theta parameters
!sapply(plchoice.allFit,inherits,"try-error") ## was fit OK?
lapply(plchoice.allFit, function(x) x@optinfo$conv$lme4$messages)
# bobyqa and nlminb gave no errors, so switching to nlminb for rewritten pub
plchoice0 <- update(plchoice0, control = glmerControl(optimizer = 'optimx', optCtrl = list(method = 'nlminb')))

# Doing LRTs
plnointeraction <- update(plchoice0, .~. - LengthNum:N1TypeNum1 - LengthNum:N1TypeNum2 - LengthNum:N1TypeNum3)
plnolength <- update(plchoice0, .~. - LengthNum)
plnonum <- update(plchoice0, .~. - N1TypeNum1 - N1TypeNum2 - N1TypeNum3)
lapply(list(plnointeraction, plnolength, plnonum, plchoice0), summary)
summary(plchoice0)

anova(plchoice0, plnolength)
anova(plchoice0, plnonum)
anova(plchoice0, plnointeraction)
summary(plchoice0)

plfactor <- glmer(ChosePl ~ Length * N1Type + (1 + Length * N1Type || Subject) +
                    (1 + Length || Item), test.trimmed, family = binomial,
                  glmerControl(optimizer = 'optimx', optCtrl = list(method = 'nlminb')))
summary(plfactor)

posthoc_consec <- lsmeans(plfactor, consec ~ N1Type | Length, adjust = 'none',
        lsm.options(disable.pbkrtest=TRUE, adjust = 'none'), type = 'response', round = 10)
posthoc_consec$contrasts
pvals <- c(0.0021, 0.0002, 0.0194, 0.0145, 0.0001, 0.0002)
p.adjust(pvals, method = 'bonf', n = 6)

posthoc_len <- lsmeans(plfactor, pairwise ~ Length | N1Type, adjust = 'none', type = 'response',
       lsm.options(disable.pbkrtest=TRUE, adjust = 'none'))#, type = 'response')
round(p.adjust(summary(posthoc_len$contrasts)$p.value, 'bonf', n = 4), 6)

lsmeans(plfactor, consec ~ N1Type, adjust = 'bonf',
        lsm.options(disable.pbkrtest = T), type = 'response')


#### Plots for publication ####
# Calculating mean proportions and SEs
props.subj <- aggregate(ChosePl ~ Length + N1Type + Subject, test.trimmed, mean)
props <- aggregate(ChosePl ~ Length + N1Type, props.subj, mean)
ses <- aggregate(ChosePl ~ Length + N1Type, props.subj, function(x) {
  sqrt((1 / length(x)) * mean(x) * (1 - mean(x)))
})
cil <- aggregate(ChosePl ~ Length + N1Type, props.subj, function(x) {
  m <- mean(x)
  se <- sqrt((1 / length(x)) * mean(x) * (1 - mean(x)))
  return(m - qnorm(0.975) * se)
})
ciu <- aggregate(ChosePl ~ Length + N1Type, props.subj, function(x) {
  m <- mean(x)
  se <- sqrt((1 / length(x)) * mean(x) * (1 - mean(x)))
  return(m + qnorm(0.975) * se)
})
props$SE <- ses$ChosePl
props$CIl <- cil$ChosePl
props$CIu <- ciu$ChosePl
props[order(props$Length),]
levels(props$Length)

ggplot(props, aes(x = N1Type, y = ChosePl)) +
  geom_pointrange(aes(ymax = CIu, ymin = CIl, shape = Length),
                  position = position_dodge(width = 0.25), size = 0.75) +
  scale_shape_manual('Modifier\npresence', values = c(16, 15)) + 
  scale_x_discrete(labels = c('Containments', 'Collections', 'Measures', 'Quantifiers'),
                   name = 'N1Type') + ylab('Probability of choosing plural verb') +
  # theme_bw(base_size = 12, base_family = 'Arial Bold') +
  theme_classic(base_size = 12, base_family = 'Arial') +
  theme(panel.grid.minor = element_blank())


#### Simulation data analysis ####
test.trimmed$DataType <- 'Human'
test.trimmed$DataType <- factor(test.trimmed$DataType, levels=c('Human','Model'))
sim_data = expand.grid(N1Type = unique(test.trimmed$N1Type), Length = unique(test.trimmed$Length), ModelData = 0, HumanData = 0)
sim_data$ModelData = c(0.627, 0.002, 1.0, 0.135, 0.842, 0.163, 1.0, 0.611)
sim_data$HumanData = tapply(test.trimmed$ChosePl, c(test.trimmed$N1Type, test.trimmed$Length), FUN=mean)
sim_data$DataType = 'Model'
sim_data$DataType <- factor(sim_data$DataType, levels=c('Human','Model'))

props$DataType='Human'
props$DataType=factor(props$DataType, levels=c('Human','Model'))

# With subject CIs, revised for response to reviewers
ggplot(props, aes(x = N1Type, y = ChosePl, shape = Length)) + 
  geom_pointrange(aes(ymax = CIu, ymin = CIl, shape = Length),
                  position = position_dodge(width = 0.25), size = 0.5) +
  scale_y_continuous('Proportion of plural verb choices', 
                     breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(-0.01, 1.02)) +
  scale_x_discrete('Subject NP Type', 
                   labels = c('Containers', 'Collections', 'Measure\nPhrases', 'Quantifiers')) +
  scale_shape('+/-N2', labels = c('+N2', '-N2')) +
  theme_classic(base_size = 12) + theme(legend.position = c(0.85, 0.3))
ggsave('~/Dropbox/Structure and Flexibility (1)/Projects/Pseudopartitives/English/Productions/Journal Submission/Resubmission/Figures/Exp2DataRevised.tiff', height = 4, width = 5.5, units = 'in', dpi=300)
