##### GME_Movement Functions #####
# Nathan Hahn


## Functions

# Prep data for traj cleaning
data.prep <- function(df, cols =  c("MovDataID", "CalcID", "Fixtime", "X", "Y"), collar.data,
                      rm.names = FALSE, names) {
  # select, reorder, and rename relevant columns
  require(dplyr)
  x <- df %>%
    dplyr::select(cols) %>%
    setNames(nm = c("MovDataID", "CalcID", "Fixtime", "X", "Y")) %>%
    mutate(CalcID = as.factor(CalcID))
  
  meta <- df %>%
    dplyr::select(collar.data)
  
  x <- as.data.frame(cbind(x, meta))
  
  # remove names (if needed)
  if (rm.names == TRUE) {
    y <- x[-which(x$MovDataID%in%names),]
    y <- as.POSIXct(strptime(y$Fixtime, format = "%Y-%m-%d %H:%M:%S", tz="Africa/Nairobi" ))
    
    return(as.data.frame(y))
  }
  
  else {
    x$Fixtime <- as.POSIXct(strptime(x$Fixtime, format = "%Y-%m-%d %H:%M:%S", tz="Africa/Nairobi" ))
    
    return(as.data.frame(x))
  }
  
}


# converter between latlon and xy
xyConv <- function(df, xy, CRSin = '+proj=longlat +datum=WGS84', 
                   CRSout, replace = TRUE) {
  
  df <- df[complete.cases(df[, xy]), ]
  coord <- data.frame(df[, xy])
  colnames(coord) <- c('x', 'y') 
  coord[, 1] <- as.numeric(coord[, 1])
  coord[, 2] <- as.numeric(coord[, 2])
  conv <- SpatialPoints(coordinates(coord),
                        proj4string = CRS(CRSin))
  conv <- spTransform(conv, CRS(CRSout))
  conv <- data.frame(conv)
  colnames(conv) <- c('x', 'y')
  df <- cbind(df, conv)
  df.T <- df %>% dplyr::select(-xy) %>%
    dplyr::rename(X = x, Y = y)
  
  if (replace == TRUE) {
    return(df.T)
  }
  
  else {
    return(df)
  }
  
}

# Apply adeHabitat to a tracking dataset. 
# tol + fix.units sets the tolerance, eg 30/"min" = 30 minutes
# use regular = TRUE to regularize relocation times on the hour/half hour. 
traj.func <- function(df, refdate = "2011-09-01 00:00:00", 
                      tol = 1, fix.units = "hour", regular, ...) {
  require(adehabitatLT)
  
  # set a reference date
  refda <- strptime(refdate, "%Y-%m-%d %H:%M:%S", tz="Africa/Nairobi")   #add ref date
  
  #ltraj object
  x <- as.ltraj(xy = df[,c("X","Y")], date = df$Fixtime, id = df$CalcID, ...)
  #Seperate into multiple bursts
  x <- cutltraj(x, "foo(dt)", nextr = TRUE)
  # create NAs based on hourly relocations
  x <- setNA(x, refda, tol, units = fix.units)
  
  if (regular == TRUE) {
    # "round" the relocs to make regular
    x <- sett0(x, refda, tol, units = fix.units)
  }
  
  else {
    x <- x
  }
  
  df.sum <- summary.ltraj(x) #dataframe of ltraj summary
  df.sum$fixRate <- df.sum$NAs/df.sum$nb.reloc
  df <- ld(x) #dataframe of tracking data
  
  return(list(ltraj = x, summary.df = df.sum, tracking.df = df))
  
}


# Filter trajectories given missing fix threshold (p.missing) and minimum relocs for a burst (n.reloc)
filter.traj <- function(df, n.reloc = 500, p.missing = .1) {
  df$pct <- df$NAs/df$nb.reloc
  keep <- df %>%
    filter(nb.reloc > n.reloc) %>%
    filter(pct < p.missing) %>%
    droplevels()
  
  to.keep <- keep$burst
  
  return(list(df = keep, ref = to.keep))
}


##### plot_budget #####
# function for plotting activity budgets
plot_budget <- function(t=t, facet, title){
  require(tidyverse)
  require(lubridate)
  ggplot(t, aes(hour(date))) +
    facet_grid(facet) +
    geom_density(aes(x = hour(date),
                     fill = factor(viterbi), colour = factor(viterbi)), alpha = 0.3, adjust = 1.5) +
    # add day/night lines - code for shading below
    geom_vline(aes(xintercept=6),
               color="dark grey", linetype="dashed", size=1) +
    geom_vline(aes(xintercept=18),
               color="dark grey", linetype="dashed", size=1) +
    
    # add colors
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"), name = c("state")) +
    scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73"), name = c("state")) +
    xlab("hour (0-23)") + ggtitle(title)
}

##### velocity #####
velocity <- function(df, dateTime) {
  require(rlist)
  # calculate time differences
  tdiff <- diff(dateTime)
  units(tdiff) <- "hours"
  tdiff <- as.numeric(tdiff)
  tdiff <- list.append(tdiff, 0) # add speed zero for last fix
  
  # step lengths
  step <- df$step 
  
  # calculate velocity at each time step
  speed <- as.vector(step/tdiff)
  
  # attach to the data.frame
  df$step <- speed
  return(df)
}


##### log.velocity #####
log_velocity <- function(df) {
  require(rlist)
  # calculate time differences
  tdiff <- diff(df$date)
  units(tdiff) <- "hours"
  tdiff <- as.numeric(tdiff)
  tdiff <- list.append(tdiff, 0) # add speed zero for last fix
  
  # step lengths
  step <- df$step 
  
  # calculate velocity at each time step
  speed <- as.vector(step/tdiff)
  
  # attach to the data.frame
  df$step <- log(speed + 0.001) # add constant for zero steps
  return(df)
}


##### mask.poly.raster #####
mask_poly_raster <- function(poly, raster) {
  # crop raster to polygon extent
  cr <- crop(raster, extent(poly), snap = "out")
  # mask the raster by mcp 
  lr <- mask(x = cr, mask = poly)
  return(lr)
}


##### withold #####
# withold data
withold <- function(x, cut, type) {
  n <- nrow(x)*cut
  
  if (type == "train") {
    y <- head(x, (nrow(x)-n))
  }
  
  if (type == "test") {
    y <- tail(x, n)
  }
  
  return(y)
}


##### cluster cutpoints #####
#To extract the cutoff between any two clusters from an mclust object
cluster_cutpoint<-function(object, comp, grph=T) {
  
  mu1<-object$parameters$mean[comp[1]]
  mu2<-object$parameters$mean[comp[2]]
  sd1<-sqrt(object$parameters$variance$sigmasq[comp[1]])
  sd2<-sqrt(object$parameters$variance$sigmasq[comp[2]])
  
  if(length(object$parameters$variance$sigmasq)>1) {sd2<-sqrt(object$parameters$variance$sigmasq[comp[2]])}
  
  ss<-seq(0,.5, 0.001)
  cut<-ss[max(which(pnorm(ss, mu1, sd1, lower.tail=F)>pnorm(ss, mu2, sd2)))]
  
  if (grph) {plot(ss, pnorm(ss, mu1, sd1, lower.tail=F),col="red", type="l", xlab="Mean Ag Use", ylab="Probability", main="Prob of being in cluster")
    lines(ss,  pnorm(ss, mu2, sd2), col="blue")
    abline(v=cut, lwd=2, lty=2)}
  
  return(cut)
  
}


