##### Sera HMM & Activity Budgets #####
library(tidyverse, quietly = TRUE)
library(lubridate, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(momentuHMM, quietly = TRUE)
library(overlapping, quietly = TRUE)
library(lsmeans, quietly = TRUE)
source('Sera_functions.R')

# Set environment to EAT
Sys.setenv(TZ="Africa/Nairobi") 

# split by cohort - 4 = wild
split <- readRDS('budget.plot.df.RDS')

#as.Date(unique(sera.hmm$releaseDate))
release.date <- as.Date(c("2019-05-02", "2019-11-16", "2020-05-28", "2019-05-30"))

# function to assign periods
date_slice <- function(data, fixtime, release.date, slice.length){
  # create a date variable
  data$Day <- as.Date(fixtime)
  # assign slices based on 14 day periods from release date
  data$slice <- as.numeric(data$Day - release.date) %/% slice.length
  # record start dates of each slice
  data <- data %>% group_by(slice) %>%
    mutate(slice.start = min(Day))
}

#' Plot activity over 14-day periods - all dates are relative to the first cohort's release date. 
#' Activity time density plots are faceted by relative period, showing pre and post release. Immediate use of exploratory state after release

# set period length
per <- 14

#' ```{r fig.width=6, fig.height=15}
cohort.1 <- date_slice(data = split[[1]], fixtime = split[[1]]$date, release.date = release.date[1], slice.length = per)
plot_budget(cohort.1, facet = slice~viterbi, title = 'cohort 1')
#' ```

cohort.2 <- date_slice(data = split[[2]], fixtime = split[[2]]$date, release.date = release.date[1], slice.length = per)
plot_budget(cohort.2, facet = slice~viterbi, title = 'cohort 2')

cohort.3 <- date_slice(data = split[[3]], fixtime = split[[3]]$date, release.date = release.date[1], slice.length = per)
plot_budget(cohort.3, facet = slice~viterbi, title = 'cohort 3') 

cohort.4 <- date_slice(data = split[[4]], fixtime = split[[4]]$date, release.date = release.date[1], slice.length = per)
plot_budget(cohort.4, facet = slice~viterbi, title = 'wild') 

