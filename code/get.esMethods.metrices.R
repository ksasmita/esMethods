#######
# Written by: Karen Sasmita on 10/20/2019
# Modified: 04/27/2020
#
# Modified: 10/10/2020 
#   Added bw input to peak-to-peak distance and peakiness. 
#
# Modified: 10/20/2021
#   Added t.window parameter to surprise index
#   Corrected rarity calculation: n should be sequence from 0 to n.overlap -1. 
#
# Modified: 10/25/2021
#   Changed: calculation of button press density distribution over time
#            calculate density of button press for from 0s to movie duration 
#            added padding to density estimate, to account for low probabilities (approaching 0) around the start and end of movie
#
#   Changed: how to determine normative boundaries in peak-to-peak distance calculation to account for average number of button presses of both groups
#
#######
library(seewave)#required to calculate rugosity (for peakiness)

##Peak-to-peak distance

get.peaktopeak.distance = function(sample1, sample2, mov.dur, adj.size, ave.bp.sample1, ave.bp.sample2, bw = 'SJ', dens.hz = 1){
  #####
  #Returns: calculate the average nearest distance between normative boundaries in sample1 and sample 2 and vice versa
  #
  #Parameter sample1: array of button press times for every participant in the first group
  #Precondition sample1: numeric
  #
  #Parameter sample2: array of button press times for every participant in the second group
  #Precondition sample2: numeric
  #
  #Parameter mov.dur: movie duration
  #Precondition mov.dur: numeric
  #
  #Parameter adj.size: size of bandwidth adjustment to be applied in the density function
  #Precondition adj.size: numeric
  #
  #Parameter ave.bp,sample1: average number of button presses taken from sample1
  #Precondition ave.bp.sample1: numeric
  #
  #Parameter ave.bp,sample2: average number of button presses taken from sample2
  #Precondition ave.bp.sample2: numeric
  #
  #Parameter bw: method for bandwidth setting for the density function
  #Precondition bw: type
  #
  #Parameter dens.hz: resolution of density time series in hz (samples/ second)
  #Precondition dens.hz: numeric
  #####
  
  #create temporary density distribution to get bandwidth of sample1 density distribution  
  temp <- density(sample1, bw = bw, adjust = adj.size, kernel = 'g', n = ceiling(mov.dur*dens.hz), from = 1/dens.hz, to = mov.dur)
  
  #add padding to time series 
  pad.ts <- ceiling(temp$bw * 2)
  
  #create density distributions of button presses with padding 
  sample1.dens <- density(sample1, bw = bw, adjust = adj.size, kernel = 'g', n = ceiling((mov.dur+2*pad.ts)*dens.hz), from = 1/dens.hz - pad.ts, to = mov.dur + pad.ts)
  sample2.dens <- density(sample2, bw = bw, adjust = adj.size, kernel = 'g', n = ceiling((mov.dur+2*pad.ts)*dens.hz), from = 1/dens.hz - pad.ts, to = mov.dur + pad.ts)
  
  #put values into dataframe for easy access 
  sample1.dens.df <- data.frame("Time" = sample1.dens$x, "Density" = sample1.dens$y, "peak" = 0)
  sample2.dens.df <- data.frame("Time" = sample2.dens$x, "Density" = sample2.dens$y, "peak" = 0)
  
  #Determine peaks and troughs for the density distribution. 
  #Peak = when density at time t is higher than density at t-1 and t+1. 
  #Trough = when density at t is lower than density at t-1 and t+1.
  for (row in seq(from = 2, to = nrow(sample1.dens.df)-1, by = 1)){
    if (sample1.dens.df$Density[row] > sample1.dens.df$Density[row+1] && sample1.dens.df$Density[row] > sample1.dens.df$Density[row-1]){
      sample1.dens.df$peak[row] <- 1}
    else if (sample1.dens.df$Density[row] < sample1.dens.df$Density[row+1] && sample1.dens.df$Density[row] < sample1.dens.df$Density[row-1]){
      sample1.dens.df$peak[row] <- -1} 
  }

  for (row in seq(from = 2, to = nrow(sample2.dens.df)-1, by = 1)){
    if (sample2.dens.df$Density[row] > sample2.dens.df$Density[row+1] && sample2.dens.df$Density[row] > sample2.dens.df$Density[row-1]){
      sample2.dens.df$peak[row] <- 1}
    else if (sample2.dens.df$Density[row] < sample2.dens.df$Density[row+1] && sample2.dens.df$Density[row] < sample2.dens.df$Density[row-1]){
      sample2.dens.df$peak[row] <- -1} 
  }
  
  ##get all peaks
  sample1.peakTimes <- sample1.dens.df$Time[sample1.dens.df$peak == 1]
  sample2.peakTimes <- sample2.dens.df$Time[sample2.dens.df$peak == 1]
  
  ##calculate the number of normative boundaries (minimum of c(ave number of bp sample1, sample1 number of peaks, sample2 number of peaks))
  n.norm.bp <- min(c(round(mean(c(ave.bp.sample1, ave.bp.sample2))), length(sample1.peakTimes), length(sample2.peakTimes)))
  
  #Obtain times when peak occurs
  #Peak times then sorted in ascending order. 
  sample1.normativeTimes <- sort(head(sample1.dens.df[order(-sample1.dens.df$peak, -sample1.dens.df$Density),],n.norm.bp)$Time)
  sample2.normativeTimes <- sort(head(sample2.dens.df[order(-sample2.dens.df$peak, -sample2.dens.df$Density),],n.norm.bp)$Time)
  
  #Calculate the average distance between every peak in sample 1 to the nearest peak in sample 2 by getting the minimum absolute distance between each peak of sample 1 with every peak in sample 2 and vice versa. 
  nearest.dist <- data.frame()
  for(i in 1:length(sample1.normativeTimes)){
    nearest.dist[i,1] <- min(abs(sample1.normativeTimes[i] - sample2.normativeTimes))
    nearest.dist[i,2] <- min(abs(sample2.normativeTimes[i] - sample1.normativeTimes))
  }

  return(mean(c(nearest.dist[,1], nearest.dist[,2])))
}

##Peakiness

get.peakiness = function(y.vect, mov.dur, sub.list, adj.size, bw = 'SJ', dens.hz=1){
  #####
  # Adapted from script written: programmed 06.09.16 by Khena Swallow
  # Changed: take bandwidth input 
  #
  # Returns: minimum rugosity, actual rugosity and peakiness (actual rugosity/minimum rugosity)
  # Minimum peakiness calculated:
  # Evenly distributes button presses through the movie in order to minimize overlap and maximize smoothness in the group time series.
  # The algorithm will count the total number of button presses across participants, then generate a sequence of evenly spaced button presses across the movie.
  # These evenly spaced button presses will then be used to calculate button press density over time, and the minimum smoothness of the time series.
  #
  #Parameter y.vect: array of button press times for every participant in sample
  #Precondition sample1: numeric
  #
  #Parameter mov.dur: movie duration (seconds)
  #Precondition mov.dur: numeric
  #
  #Parameter sub.list: array of matching size to y.vect that indicates subject which each button press time belongs to
  #Precondition sub.list: numeric
  #
  #Parameter adj.size: size of bandwidth adjustment to be applied in the density function
  #Precondition adj.size: numeric
  #
  #Parameter bw: method of setting bandwidth for density estimate 
  #Precondition: type
  #
  #Parameter dens.hz: resolution of density time series in hz (samples/ second)
  #Precondition dens.hz: numeric
  #
  #####
  # first, get the distribution of button presses for the minimum
  total.bps = length(y.vect)
  even.bps = seq(from = 0, to = mov.dur, length.out = total.bps) 
  
  #get density distribution for minimum and actual data
  temp = density(y.vect, bw = bw, adjust = adj.size, kernel = 'g', n = mov.dur*dens.hz, from = 1/dens.hz, to = mov.dur)
  
  pad.ts <- ceiling(temp$bw * 2)
  y.dens = density(y.vect, bw = bw, adjust = adj.size, kernel = 'g', n = ceiling((mov.dur+2*pad.ts)*dens.hz), from = 1/dens.hz - pad.ts, to = mov.dur + pad.ts)
  min.y.dens = density(even.bps, bw = bw, adjust = adj.size, kernel = 'g', n = ceiling((mov.dur+2*pad.ts)*dens.hz), from = 1/dens.hz - pad.ts, to = mov.dur + pad.ts)
  
  #calculate peakiness using rugo function
  min.rugo = rugo(min.y.dens$y)
  act.rugo = rugo(y.dens$y)
  peakiness = act.rugo/min.rugo
  
  return(c(min.rugo, act.rugo, peakiness))
}

##Agreement Index

get.agreement.index = function(sample.data, group.data){
  ###
  #Returns the agreement index - how correlated are each subject's segmentation pattern over time with the rest of the sample, adjusted for the number of button press made by each individual. (Kurby & Zacks, 2011)
  #
  #Parameter sample.data: array of presence of button press for every n-second of the movie length for a single participant (x)
  #Precondition sample.data: numeric
  #
  #Parameter group.data: array of average button presses for every n-second of the movie length for a group data excluding participant x. 
  #Precondition group.data: numeric
  ####
  
  #calculate correlation between sample and group timeseries 
  corr <- cor(sample.data, group.data)
  
  #calculate maximum possible correlation between sample and group timeseries. 
  #this is done by calculating the correlation of sample and group timeseries ordered in the same direction (i.e. if both timeseries are perfectly correlated) 
  maxcorr <- cor(sample.data[order(-sample.data)], group.data[order(-group.data)]) 
  
  #calculate minimum possible correlation between sample and group timeseries. 
  #this is done by calculating the correlation of sample and group timeseries ordered in opposite directions (i.e. if they are perfectly anticorrelated)
  mincorr <- cor(sample.data[order(sample.data)], group.data[order(-group.data)])
  
  #calculate agreement index 
  agreement.index <- (corr - mincorr)/(maxcorr-mincorr)
  
  return(c(mincorr, maxcorr, corr, agreement.index))
}

##Surprise Index

get.surprise.index = function(sub.bp, gp.bp, mov.dur, ave.bp, bw = "SJ", adj.size, t.window = 1, dens.hz=1){
  ###
  #Returns surprise index (adapted from Katori et al., 2018) = Degree that expected overlap between individual button press and group normative boundary is exceeded.  
  #
  #Parameter sub.bp = array of button press times for a single participant (x). 
  #Precondition sub.bp = numeric
  #
  #Parameter gp.bp = array of button press times for every individual (concatenated) within a group data excluding participant x.
  #Precondition gp.bp = numeric
  #
  #Parameter mov.dur = total movie duration in seconds 
  #Precondition mov.dur = numeric
  #
  #Parameter ave.bp = average button press of the group excluding participant x. 
  #Precondition ave.bp = numeric
  #
  #Parameter bw: method of setting bandwidth for density estimate 
  #Precondition: type
  #
  #Parameter adj.size: size of bandwidth adjustment to be applied in the density function
  #Precondition adj.size: numeric
  #
  #Parameter t.window = time window of overlap between individual button press time and normative boundary
  #Precondition t.window = numeric
  #
  #Parameter dens.hz: resolution of density time series in hz (samples/ second)
  #Precondition dens.hz: numeric
  #
  #####
  N <- mov.dur/t.window #total number of possible normative peaks 
  sub.bprate <- length(sub.bp)/N
  
  #Transform group button press to density over time 
  if (length(gp.bp) > 1){
    temp = density(gp.bp, bw = bw, adjust = adj.size, kernel = 'g', n = mov.dur*dens.hz, from = 1/dens.hz, to = mov.dur)
    
    pad.ts <- ceiling(temp$bw * 2)
    gp.dens = density(gp.bp, bw = bw, adjust = adj.size, kernel = 'g', n = ceiling((mov.dur+2*pad.ts)*dens.hz), from = 1/dens.hz - pad.ts, to = mov.dur + pad.ts)
    
    if(2.35/gp.dens$bw < t.window){
      warning('Density bandwidth is larger than t.window')
    } 
    
    gp.dens.df <- data.frame("Time" = gp.dens$x, "Density" = gp.dens$y, "peak" = 0)
    for (row in seq(from = 2, to = nrow(gp.dens.df)-1, by = 1)){
      if (gp.dens.df$Density[row] > gp.dens.df$Density[row+1] && gp.dens.df$Density[row] > gp.dens.df$Density[row-1]){
        gp.dens.df$peak[row] <- 1}
      else if (gp.dens.df$Density[row] < gp.dens.df$Density[row+1] && gp.dens.df$Density[row] < gp.dens.df$Density[row-1]){
        gp.dens.df$peak[row] <- -1}
    }
  } else {
    gp.dens.df <- data.frame("Time" = gp.bp, "Density" = 1, "peak" = 1)
  }
  
  #get peaks from button press density estimates
  gp.peakTimes <- gp.dens.df$Time[gp.dens.df$peak == 1]
  
  #get n of normative boundary 
  if(length(gp.peakTimes) > ave.bp){
    n.normativePeaks <- ave.bp
  } else if (length(gp.peakTimes) <= ave.bp){ #if n peak times < ave.bp
    n.normativePeaks <- length(gp.peakTimes) 
  }
  
  #get normative boundary times
  gp.normativeTimes <- sort(head(gp.dens.df[order(-gp.dens.df$peak, -gp.dens.df$Density),],n.normativePeaks)$Time)
  
  gp.boundaryrate <- length(gp.normativeTimes)/N
  
  get.overlap <- function(bp.time){
    return(any(abs(bp.time-gp.normativeTimes) <= t.window/2))
  }
  
  n <- sum(do.call(rbind, lapply(sub.bp, get.overlap))) #number of times that individual button press fall within t.window (t.window/2 before and t.window/2 after) of normative boundaries
  n.seq <- seq(from = n, to = N, by = 1)
  
  #calculate surprise index 
  lambda <- sub.bprate*gp.boundaryrate*N
  rarity <- sum(dpois(n.seq, lambda))
  sup.index <- -log(rarity, base = 2)
  
  return(sup.index)
}


#To get values to needed for calculation of predictive value and detection accuracy 

get.normative.hitRate_bin <- function(sub.bp, gp.bp,ave.bp, timeSeries){
  ###
  #Returns:
  #  1. normativeHit = hit rate (proportion of indiv button press overlapping with normative boundary) 
  #  2. normativeFA = false alarm rate (proportion of indiv button press non-overlapping with non-normative boundary)
  #  3. normativeMiss = misses (proportion of indiv button press non-overlapping with normative boundary)
  #  4. normativeCR = correct rejection (proportion of indiv button press overlapping with non-normative boundary)
  #  5. normativePPV = positive predictive value (indiv button press overlap with normative boundary/total indiv bp times)
  #  6. normativeFalsePPV = indiv button press overlap with non normative boundary / total indiv bp times
  #
  #Parameter sub.bp = array of presence of button press for n-second bin size of the movie for a single participant (x). 
  #Precondition sub.bp = numeric 
  #
  #Parameter gp.bp = array of average button presses for every n-second bin size of the movie  for the group data excluding participant x.
  #Precondition gp.bp = numeric 
  #
  #Parameter timeSeries = array of timeseries to match sub.bp and gp.bp 
  #Precondition timeSeries = numeric 
  #
  ###
  #
  #Determine timepoints of individual button press 
  sub.bpTimes <- timeSeries[sub.bp == 1]
  
  gp.normativeTimes <- head(gp.bp[order(-gp.bp$group.bp),],ave.bp)
  gp.nonBoundaryTimes <- timeSeries[!timeSeries %in% gp.normativeTimes$timeSeries]
  #Determine normative hit rate (sub.bpTimes overall with gp.bpTimes/ total gp.bpTimes)
  normativeHitRate <- sum(sub.bpTimes %in% gp.normativeTimes$timeSeries)/length(gp.normativeTimes$timeSeries)
  normativeFalseAlarm <- sum(sub.bpTimes %in% gp.nonBoundaryTimes)/length(gp.nonBoundaryTimes)
  normativeMiss <- (length(gp.normativeTimes$timeSeries)-sum(sub.bpTimes %in% gp.normativeTimes$timeSeries))/length(gp.normativeTimes$timeSeries)
  normativeCorrectRejection <- (length(gp.nonBoundaryTimes)-sum(sub.bpTimes %in% gp.nonBoundaryTimes))/length(gp.nonBoundaryTimes)
  normativePPV <- sum(sub.bpTimes %in% gp.normativeTimes$timeSeries)/length(sub.bpTimes)
  normativeFalsePPV <- sum(sub.bpTimes %in% gp.nonBoundaryTimes)/length(sub.bpTimes)
  
  return(data.frame('normativeHit' = normativeHitRate, 'normativeFA' = normativeFalseAlarm, "normativeMiss" = normativeMiss, "normativeCR" = normativeCorrectRejection, "normativePPV" = normativePPV, "normativeFalsePPV" = normativeFalsePPV))
}


