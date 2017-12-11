# Cleaned-up analysis script for Smith, Franck, & Tabor pseudopartitives paper
# Experiment 1: notional plurality norms

library(ggplot2)
library(lme4)

# Assumes data file is in the current working directory
# The variable 'MoreThings' is 1 if the participant indicated 
# that the phrase referred to more than one thing, and 0 for one thing.
exp1data <- read.csv('Exp1Data.csv')
# Reordering factor levels; 'FullCont' refers to Containers
exp1data$Type <- factor(exp1data$Type, levels = c('FullCont', 'Coll', 'Meas', 'Quant'))

# Setting up backwards difference contrasts
# These compare level n of a factor to level n-1
# See http://www.ats.ucla.edu/stat/r/library/contrast_coding.htm
contrasts(exp1data$Type) = matrix(c(-3/4, 1/4, 1/4, 1/4, -1/2, -1/2, 1/2, 1/2, 
                                    -1/4, -1/4, -1/4, 3/4), ncol = 3)

# The logistic regression model reported in the paper.
# This model ignores dependencies in the data, though...
exp1mod <- glm(MoreThings ~ Type, family = binomial, data = exp1data, na.action = na.exclude)
summary(exp1mod)
# Deviance Residuals: 
# Min       1Q   Median       3Q      Max  
# -2.2073  -0.8971   0.4279   0.7475   1.4865  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  1.26073    0.11564  10.902  < 2e-16 ***
#   Type1        1.83477    0.24924   7.362 1.82e-13 ***
#   Type2        1.13645    0.32782   3.467 0.000527 ***
#   Type3        0.07587    0.38967   0.195 0.845633    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 750.91  on 639  degrees of freedom
# Residual deviance: 575.45  on 636  degrees of freedom
# AIC: 583.45
# 
# Number of Fisher Scoring iterations: 5
drop1(exp1mod, test = 'Chisq')
# MoreThings ~ Type
# Df Deviance    AIC    LRT  Pr(>Chi)    
# <none>      575.45 583.45                     
# Type    3   750.91 752.91 175.46 < 2.2e-16 ***

# Now trying a hierarchical logistic model using the default settings:
# ResponseID = participant ID, NP = which noun phrase
exp1default <- glmer(MoreThings ~ Type + (1 + Type | ResponseID) + (1 | NP),
                     exp1data, family=binomial, na.action=na.exclude)
# Warning message:
#   In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                  Model failed to converge with max|grad| = 0.104347 (tol = 0.001, component 1)
summary(exp1default) 
# Similar results to the GLM, except that the difference between Measures and Collections is marginal
sigtest <- drop1(exp1default, test='Chisq')
# Same type of warning message for the reduced model: max|grad| = 0.0511559 (tol = 0.001, component 1)
# Model:
#   MoreThings ~ Type + (1 + Type | ResponseID) + (1 | NP)
# Df    AIC    LRT   Pr(Chi)    
# <none>    464.44                     
# Type    3 478.30 19.866 0.0001809 ***

# Trying to simplify the random effects structure:
exp1nocorr <- glmer(MoreThings ~ Type + (1 | ResponseID) + (0 + Type | ResponseID) + (1 | NP),
                    exp1data, family=binomial, na.action=na.exclude)
# Same gradient issue: max|grad| = 0.0298319
summary(exp1nocorr) # This fails to converge in 10000 iterations and has max|grad| = 0.041496
# Upping number of iterations:
exp1nocorr <- glmer(MoreThings ~ Type + (1 | ResponseID) + (0 + Type | ResponseID) + (1 | NP),
                    exp1data, family=binomial, na.action=na.exclude,
                    control=glmerControl(optCtrl = list(maxfun=50000)))
# Still fails to converge: max|grad| = 0.0421596

# Trying out all optimizers in a hierarchical logistic model:
# See ?convergence and https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
# The all_fit function is in the afex package
library(optimx) # So we have access to all the optimizers
library(afex)
all.mods <- all_fit(exp1default)
# Output during fitting:
# bobyqa. : [OK]
# Nelder_Mead. : [OK]
# optimx.nlminb : [OK]
# optimx.L-BFGS-B : [OK]
# nloptwrap.NLOPT_LN_NELDERMEAD : [OK]
# nloptwrap.NLOPT_LN_BOBYQA : [OK]
# nmkbw. : [OK]
# Warning messages:
# 1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
# Model failed to converge with max|grad| = 0.0716725 (tol = 0.001, component 1)
# 2: In optimx.check(par, optcfg$ufn, optcfg$ugr, optcfg$uhess, lower,  :
# Parameters or bounds appear to have different scalings.
# This can cause poor performance in optimization. 
# It is important for derivative free methods like BOBYQA, UOBYQA, NEWUOA.
# 3: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
# Model failed to converge with max|grad| = 0.0146211 (tol = 0.001, component 1)
# 4: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
# Model failed to converge with max|grad| = 0.0146211 (tol = 0.001, component 1)
# 5: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
# Model failed to converge with max|grad| = 0.00331505 (tol = 0.001, component 1)

# Seeing where the messages came from:
lapply(all.mods, function(x) x@optinfo$conv$lme4$messages)
# All but built-in bobyqa, optimx:nlminb, and optimx:L_BFGS-B give gradient convergence warnings
# Checking (following ?all_fit)
t(sapply(all.mods, fixef)) # extract fixed effects: all pretty close
sapply(all.mods, logLik) # log-likelihoods: all pretty close
sapply(all.mods, getME, "theta") # theta parameters
!sapply(all.mods, inherits, "try-error") # was fit OK?

# Looking at the warning-less models
lapply(list(all.mods$bobyqa., all.mods$optimx.nlminb, all.mods$`optimx.L-BFGS-B`), summary)
# All quite similar. Nlminb has extra warnings, so I'll avoid that one

