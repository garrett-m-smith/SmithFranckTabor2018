# Analyzing syn./sem. tests on pseudopartitive experiment

library(afex)
library(lattice)
library(ggplot2)
library(IBEX.to.R) # downloaded from Google forum; by Anton Malko
library(dplyr)

# Reading in data; unusable data already removed
d <- read.csv('./Exp3Data.csv')

# Acclimatization phase
acclim <- droplevels(subset(d, type %in% paste('Acclimatization-', seq(9, 1), sep = '')))
with(acclim, table(sentence, answer))
aggregate(answer ~ sentence, acclim, FUN = mean)
ggplot(acclim, aes(sentence, answer)) + stat_summary(fun.data = 'mean_cl_boot') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Checking fillers
fillers <- droplevels(d[grepl('Filler', d$type, fixed = T),])
table(fillers$answer) # looks about right
aggregate(answer ~ sentence, fillers, FUN = mean) # 3 Pbn. w/ means > 5, leave in for now
mean(fillers$answer); sd(fillers$answer) # seems reasonable


# Looking at abstract N2s: 
abstract <- droplevels(subset(d, Test == 'Abstract Noun'))
ggplot(abstract, aes(Category, answer)) + stat_summary(fun.data = 'mean_cl_boot') + ylim(1, 7)+
  ggtitle('Abstract N2 experiment')
abstract_mod <- mixed(answer ~ Category + (1 + Category | subject) + (1 | item.id), abstract, method = 'LRT', cl = cluster)
abstract_mod # significant
summary(abstract_mod); levels(abstract$Category)
tapply(abstract$answer, abstract$Category, function (x) {
  return(c(mean(x), sd(x)))
})
pairs(lsmeans(abstract_mod, 'Category'), adjust = 'mvt')
ggplot(abstract, aes(answer, col = Category)) + geom_density() + facet_grid(~Category)

# Overflowing
overflowing <- droplevels(subset(d, Test == 'Overflowing'))
levels(overflowing$Category)
ggplot(overflowing, aes(Category, answer)) + stat_summary(fun.data = 'mean_cl_boot') + ylim(1, 7) +
  ggtitle('Overflowing experiment')
overflowing_mod <- mixed(answer ~ Category + (1 + Category | subject) + (1 | item.id), overflowing, method = 'LRT', cl = cluster)
overflowing_mod # doch significant!
summary(overflowing_mod); levels(overflowing$Category)
tapply(overflowing$answer, overflowing$Category, function (x) {
  return(c(mean(x), sd(x)))
})
ggplot(overflowing, aes(answer, col = Category)) + geom_density() + facet_grid(~Category)

# Breakableness
breakable <- droplevels(subset(d, Test == 'Breakableness'))
levels(breakable$Category)
ggplot(breakable, aes(Category, answer)) + stat_summary(fun.data = 'mean_cl_boot') + ylim(1, 7) +
  ggtitle('Breakableness experiment')
breakable_mod <- mixed(answer ~ Category + (1 + Category | subject) + (1 | item.id), breakable, method = 'LRT', cl = cluster)
breakable_mod # significant
summary(breakable_mod); levels(breakable$Category)
tapply(breakable$answer, breakable$Category, function (x) {
  return(c(mean(x), sd(x)))
})
