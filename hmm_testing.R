##### Sera HMM #####


library(tidyverse)
library(lubridate)
library(ggplot2)
source('Sera_functions.R')

# Set environment to EAT
Sys.setenv(TZ="Africa/Nairobi") 

rawData <- read.csv('./movdata/Sera_datasets_10Mar21.csv')
meta <- read.csv('./movdata/sera_metadata.csv')
rawData <- merge(rawData, meta, by = 'MovDataID')

##### 1. Clean trajectories #####
#' *Prep data and check fixes
#' *Round fixes to nearest hour
#' *Remove individuals with fix rates <95%

## Prep and check data
sera <- data.prep(rawData, collar.data = c('Class', 'releaseDate', 'cohort'))

meta.table <- sera %>% group_by(MovDataID, CalcID, Class, releaseDate, cohort) %>%
  tally()
View(meta.table)

# viz of raw fix times
p <- ggplot(sera, aes(x = Fixtime, y = minute(Fixtime))) + geom_point() + facet_wrap(.~MovDataID) +
  ylab("minutes from 0")
p

# Downsample to clean up immobility alerts and recheck result. 
sera <- sera %>%
  filter(minute(Fixtime) > 55 | minute(Fixtime) < 5 | between(minute(Fixtime), 25, 35))

#Downsample to 1-hour because wild collars are on hourly schedule
sera <- sera %>%
  filter(minute(Fixtime) > 55 | minute(Fixtime) < 5)

# check downsample
p <- ggplot(sera, aes(x = Fixtime, y = minute(Fixtime))) + geom_point() + facet_wrap(.~MovDataID) +
  ylab("minutes from 0")
p

## Clean Trajectories 

# Set burst function -- no relocs within 4 hours creates a seperate burst
foo <- function(dt) {
  return(dt > (6*3600))
}

# regularize traj
system.time({
  ele.traj <- traj.func(df = sera, tol = 1, fix.units = "hour", regular = TRUE, 
                        infolocs = sera[,c(1,6:8)]) 
})

# mean fix rate looks good
1 - mean(ele.traj$summary.df$fixRate)

# check output
# head(ele.traj$tracking.df)

##### 2. HMM dataframe #####
#' HMMs are fit using bursts as ID. Filter out bursts with low fix rates and <500 points. We also take the log step length.

library(momentuHMM)
# create burst filter list that meets HMM specs (>500 points in a burst, <5% missing data)
t <- filter.traj(df = ele.traj$summary.df, n.reloc = 500, p.missing = .05)
hmm.filter <- t$df

# all individuals are still included after filtering out bursts 
levels(hmm.filter$id)

# check data loss after filtering out bursts - minimal
1-sum(hmm.filter$nb.reloc)/dim(sera)[1]

# create HMM dataframe
sera.hmm <- filter(ele.traj$tracking.df, burst %in% hmm.filter$burst) %>%
  rename(ID = burst) # treat bursts as individuals for fitting
sera.hmm <- subset(sera.hmm, !is.na(sera.hmm$x))

# prep data for HMM
sera.step <- prepData(sera.hmm)
sera.log.step <- sera.step
sera.log.step$step <- log(sera.log.step$step + 0.001) # add constant for zero steps

##### 3. Fit HMMs #####
# state label
stateNames2 <- c("encamped", "exploratory")
stateNames3 <- c("encamped", "meandering", "dirwalk")

# set distributions
# distGam = list(step = "gamma", angle = "vm") # zero-inflated gamma distributions
distNorm = list(step = "norm", angle = "vm") # normal step distribution since we use log step length

# Null 2-state
par2 <- list(step = c(2, 6, 2, 1),
             angle = c(.2,.2)) 
system.time(
  mod2 <- fitHMM(data = sera.log.step, nbStates = 2, dist = distNorm,
                 Par0 = par2,
                 retryFits = 5,
                 formula = ~ (1:MovDataID), # add a random effect for individual
                 stateName = stateNames2,
                 modelName = "twoStep")
)

mod2

#' Some autocorrelation in the step lengths
plotPR(mod2)

# Null 3-state
par3 <- list(step = c(1, 4, 6, 1.5, 1, 1),
             angle = c(.2,.2,.2)) 
             
system.time(
  mod3 <- fitHMM(data = sera.log.step, nbStates = 3, dist = distNorm,
                 Par0 = par3,
                 retryFits = 1,
                 formula = ~ (1:MovDataID), # add a random effect for individual
                 stateName = stateNames3,
                 modelName = "threeState null")
)

# save 3-state model
saveRDS(mod3, 'mod3_noCovs.RDS')
mod3 <- readRDS('mod3_noCovs.RDS')

# model summary
mod3

# still some autocorrelation but not as bad as 2-state. Seems to be a two-day pattern
plotPR(mod3)

##### Activity Budgets #####

# TEMP - assign orphan/non-orphan
#sera.hmm <- merge(sera.hmm, meta, by = "MovDataID")

# apply viterbi algorithm to assign states
sera.hmm$viterbi <- as.numeric(viterbi(mod3)) 
sera.hmm$state <- sera.hmm$viterbi
levels(sera.hmm$state) <- c('encamped', 'meandering', 'dirwalk')

# create true start and end dates
sera <- sera %>%
  group_by(CalcID) %>%
  mutate(StartDate = min(Fixtime), EndDate = max(Fixtime))

#' Overall time budgets show similar but not identical activity between cohorts. 4 = wild
# pct budgets
prop.table(table(sera.hmm$viterbi, sera.hmm$cohort), margin = 2)

# also check sample sizes
sera.hmm %>% group_by(cohort) %>% tally()

#' Activity patterns generated using kernel density estimates also show similarly timed peaks across the 
#' three states. There is a slightly higher peak in meandering for orphans around mid-day. Given the 
#' similarity, it could be interesting to split orphan time between the two release stages. 

# plot
plot_budget(t = sera.hmm, facet = cohort~viterbi, title = 'budget by cohort (4 = wild)')

##### Activty Time Slices #####
# ref dates for each release - wild will be assessed over the same time period

# split by cohort - 4 = wild
split <- split(sera.hmm, sera.hmm$cohort)
as.Date(unique(sera.hmm$releaseDate))

release.date <- as.Date(c("2019-05-02", "2019-11-16", "2020-05-28", "2019-05-30"))

# assign periods
date_slice <- function(data, fixtime, release.date, slice.length){
  # create a date variable
  data$Day <- as.Date(fixtime)
  # assign slices based on 14 day periods from release date
  data$slice <- as.numeric(data$Day - release.date) %/% slice.length
  # record start dates of each slice
  data <- data %>% group_by(slice) %>%
    mutate(slice.start = min(Day))
}


# Create cohort slices - all dates are relative to the first cohort's release date. 
# Can change this by updating the release date variable in the function to make it relative to their own release date
cohort.1 <- date_slice(data = split[[1]], fixtime = split[[1]]$date, release.date = release.date[1], slice.length = 14)
plot_budget(cohort.1, facet = slice~viterbi, title = 'cohort 1')

cohort.2 <- date_slice(data = split[[2]], fixtime = split[[2]]$date, release.date = release.date[1], slice.length = 14)
plot_budget(cohort.2, facet = slice~viterbi, title = 'cohort 2')

cohort.3 <- date_slice(data = split[[3]], fixtime = split[[3]]$date, release.date = release.date[1], slice.length = 14)
plot_budget(cohort.3, facet = slice~viterbi, title = 'cohort 3')

cohort.4 <- date_slice(data = split[[4]], fixtime = split[[4]]$date, release.date = release.date[1], slice.length = 14)
plot_budget(cohort.4, facet = slice~viterbi, title = 'wild')


##### Overlap Tests #####
library(overlapping)

# overlap test function - used in the for loop
overlap.test <- function(df, comp0, comp1, boot.it = 1000) {
  require(tidyverse)
  require(lubridate)
  require(overlapping)
  d0 <- dplyr::right_join(df, as.data.frame(comp0))
  d1 <- dplyr::right_join(df, as.data.frame(comp1))
  d <- list(as.numeric(hour(d0$date)), as.numeric(hour(d1$date)))
  
  # overlap and bootstrap test
  t <- overlap(d, plot = T)
  boot <- boot.overlap(d, B = boot.it)
  boot$OVboot_stats
  
  return(list(d = df, bootstrap = boot))
  
}

## Prep cohort dataframes 
#filter out pre-release data based on release date of each cohort, and split by relative time slice
cohort.1 <- subset(cohort.1, as.Date(cohort.1$date) >= release.date[1])
cohort.1.split <- split(cohort.1, cohort.1$slice)

cohort.2 <- subset(cohort.2, as.Date(cohort.2$date) >= release.date[2])
cohort.2.split <- split(cohort.2, cohort.2$slice)

cohort.3 <- subset(cohort.3, as.Date(cohort.3$date) >= release.date[2])
cohort.3.split <- split(cohort.3, cohort.3$slice)

# combine to single nested list of orphan cohorts 
cohort.all.list <- list(cohort.1.split, cohort.2.split, cohort.3.split) # all cohorts after filtering out pre-release data

# **for wild, backload slice to comp with cohort 1's earlier release date 
wild.split <- split(cohort.4, cohort.4$slice)
t <- length(wild.split)
wild.split[[t+1]] <- wild.split[[2]]
wild.split[[t+1]]$slice <- 0
wild.split[[1]] <- wild.split[[2]]
wild.split[[1]]$slice <- 1
# Bind back together and check result
cohort.4 <- do.call(rbind, wild.split)
plot_budget(cohort.4, facet = slice~viterbi, title = 'wild') # check
# **

# pre-filter out slices with low amounts of data - adjust based on slice length
n <- 100
cohort.all.list <- lapply(cohort.all.list, function(x) Filter(function(y) nrow(y) >= n, x))

## Tactic comp loop 
comp.all <- list()
OV_results <- list()
for(k in 1:3){ # orphan cohort loop - to be comped always with wild cohort 4
  for(i in 1:3){ # viterbi state
    for(j in 1:length(cohort.all.list[[k]])){ # slice loop must update for each cohort
      #create comps
      #create comps -- !! How to reference the correct slice as it iterates through. Will only work for individual cohorts
      comp0 <- list(viterbi = as.integer(i), cohort = k, slice = unique(cohort.all.list[[k]][[j]]$slice)) # cohort number
      comp1 <- list(viterbi = as.integer(i), cohort = 4, slice = unique(cohort.all.list[[k]][[j]]$slice)) # comp is always wild cohort 4
      
      # run and store overlap output
      t <- overlap.test(df = rbind(cohort.all.list[[k]][[j]], cohort.4), comp0 = comp0, comp1 = comp1)
      name <- paste(as.integer(i), j,k,4, sep='-') #pattern = viterbi-slice-cohort-wild
      q <- quantile(t$bootstrap$OVboot_dist, probs = c(0.05, 0.95))
      
      # store OV stats, CIs, and comp info for sorting
      comp.all[[name]] <- t$bootstrap$OVboot_stats
      comp.all[[name]]$lwr <- as.numeric(q[1])
      comp.all[[name]]$upr <- as.numeric(q[2])
      comp.all[[name]]$state <- comp0$viterbi
      comp.all[[name]]$cohort <- comp0$cohort
      comp.all[[name]]$slice <- comp0$slice
      comp.all[[name]]$slice.start <- unique(cohort.all.list[[k]][[j]]$slice.start)
      comp.all[[name]]$comp <- paste(k, 4, sep='-')
      
      # store OV bootstrap results
      OV_results[[j]] <- t$bootstrap$OVboot_dist
    }
  }
}

# join results into single dataframe
t <- do.call(rbind, comp.all)

ggplot(t, aes(slice, estOV)) + geom_pointrange(aes(ymin = lwr, ymax = upr)) + facet_wrap(state~comp) 

# explore specific comps
k = 1
i = 1
j = 44

#create comps -- !! How to reference the correct slice as it iterates through. Will only work for individual cohorts
comp0 <- list(viterbi = as.integer(i), cohort = k, slice = unique(cohort.all.list[[k]][[j]]$slice)) # cohort number
comp1 <- list(viterbi = as.integer(i), cohort = 4, slice = unique(cohort.all.list[[k]][[j]]$slice)) # comp is always wild
t <- overlap.test(df = rbind(cohort.all.list[[k]][[j]], cohort.4), comp0 = comp0, comp1 = comp1)


##### Activity Time Budgets #####
library(DescTools)
#NOT RUN
t <- rbind(cohort.1, cohort.2, cohort.3, cohort.4) %>%
  group_by(cohort, slice, viterbi) %>%
  tally()
(t)

ag.budget <- t %>%
  group_by(cohort, slice) %>%
  summarize(pct = n/sum(n))
ag.budget$state <- t$viterbi
ag.budget

# with CIs
ag.budget <- t %>%
  group_by(cohort, slice) %>%
  mutate(prop = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,1],
         lwr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,2],
         upr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,3],
  ) %>%
  mutate(cohort = as.factor(cohort), viterbi = as.factor(viterbi))

ag.budget

write.csv(ag.budget, 'ac.time.budget.csv')

# plot by each cohort-wild comp
c.1 <- subset(ag.budget, cohort %in% c('1', '4'))
c.2 <- subset(ag.budget, cohort %in% c('2', '4'))
c.3 <- subset(ag.budget, cohort %in% c('3', '4'))

# plot w/ CIs
ggplot(c.1, aes(x = as.factor(slice), y = prop, group = cohort)) + 
  geom_pointrange(aes(ymin = lwr.ci, ymax = upr.ci, color = cohort)) + 
  geom_line(aes(color = cohort), linetype = "dashed") + 
  facet_wrap(.~viterbi) + ylab("% time spent in state") + xlab("14-day period")

# get absolute differences between movement phases
split <- split(ag.budget, ag.budget$cohort)
wild <- split[[4]]
split[[4]] <- NULL

for(i in 1:3){
  split[[i]] <- merge(split[[i]], wild, by = c("slice", 'viterbi'))
  split[[i]]$diff <- abs(split[[i]]$prop.x - split[[i]]$prop.y) # absolute difference
  split[[i]]$viterbi.y <- NULL
}

i = 2

ggplot(split[[i]], aes(x = as.factor(slice), y = diff, group = viterbi)) + geom_point() + 
  geom_line(linetype = "dashed") + 
  facet_wrap(.~viterbi)

# budgets over time
temp <- do.call(rbind, split)
ggplot(temp, aes(as.factor(slice), y = prop.x, group = viterbi)) + geom_point(aes(color = viterbi)) + 
  geom_line(aes(color = viterbi), linetype = 'dashed') +
  facet_wrap(cohort.x~viterbi)


# Old For Loop that does the comps by cohort
# ## Tactic comp loop - cohort 1
# comp <- list()
# OV_results <- list()
# # viterbi loop
# for(i in 1:3){
#   # cohort loop
#   for(j in 1:length(cohort.1.split)){
#     #create comps
#     comp0 <- list(viterbi = as.integer(i), cohort = 1) # cohort number
#     comp1 <- list(viterbi = as.integer(i), cohort = 4) # comp is always wild
#     
#     # run and store overlap output
#     t <- overlap.test(df = rbind(cohort.1.split[[j]], cohort.4), comp0 = comp0, comp1 = comp1)
#     name <- paste(as.integer(i), j, 4, sep='-')
#     q <- quantile(t$bootstrap$OVboot_dist, probs = c(0.05, 0.95))
#     
#     # store state and comp name
#     comp[[name]] <- t$bootstrap$OVboot_stats
#     comp[[name]]$lwr <- as.numeric(q[1])
#     comp[[name]]$upr <- as.numeric(q[2])
#     comp[[name]]$state <- as.integer(i)
#     comp[[name]]$comp <- paste(j, 4, sep='-')
#     
#     # store OV bootstrap results
#     OV_results[[j]] <- t$bootstrap$OVboot_dist
#   }
# }
# 
# t <- do.call(rbind, comp)
#   



