#######
# Written by: Karen Sasmita on 10/20/2019
# Modified: 04/27/2020
#
# Modified: 10/10/2020 
#   Added bw input to peak-to-peak distance and peakiness. 
#
#
# 
#######
library(seewave)

##Peak-to-peak distance

get.peaktopeak.distance = function(sample1, sample2, mov.dur, adj.size, ave.bp, bw){
  #####
  #Returns: average nearest distance between every peak of sample 1 to every peak of sample 2 
  #
  #Parameter sample1: array of button press times for every participant in sample1
  #Precondition sample1: numeric
  #
  #Parameter sample2: array of button press times for every participant in sample2
  #Precondition sample2: numeric
  #
  #Parameter mov.dur: movie duration
  #Precondition mov.dur: numeric
  #
  #Parameter adj.size: size of bandwidth adjustment to be applied in the density function
  #Precondition adj.size: numeric
  #
  #Parameter ave.bp: average button press taken from sample1
  #Precondition ave.bp: numeric
  #
  #Parameter bw: bandwidth to be applied to density function
  #Precondition bw: numeric
  #####
  
  #create density distributions of both samples 
  sample1.dens <- density(sample1, bw = bw, adjust = adj.size, kernel = 'g', n = mov.dur)
  sample2.dens <- density(sample2, bw = bw, adjust = adj.size, kernel = 'g', n = mov.dur)
  
  #put values into dataframe for easy access 
  sample1.dens.df <- data.frame("Time" = sample1.dens$x, "Density" = sample1.dens$y, "peak" = 0)
  sample2.dens.df <- data.frame("Time" = sample2.dens$x, "Density" = sample2.dens$y, "peak" = 0)
  
  #Determine peaks and throughs for the density distribution. 
  #Peak = when density at time t is higher than density at t-1 and t+1. 
  #Through = when density at t is lower than density at t-1 and t+1.
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
  
  ##calculate the number of normative boundaries (ave.bp of sample1 or total number of sample1 peaks (if sample1 peaks < ave.bp))
  if(length(sample1.peakTimes) <= length(sample2.peakTimes)){
    if (length(sample1.peakTimes) > ave.bp){
      n.norm.bp <- ave.bp
    } else if (length(sample1.peakTimes <= ave.bp)){
      n.norm.bp <- length(sample1.peakTimes)
    }
  } else if (length(sample1.peakTimes)>length(sample2.peakTimes)){
    if(length(sample2.peakTimes)>ave.bp){
      n.norm.bp <- ave.bp
    } else if (length(sample2.peakTimes <= ave.bp)){
      n.norm.bp <- length(sample2.peakTimes)
    }
  }
  
  
  #Obtain times when peak occurs, rounded to the next second. Where number of peaks = average button press in sample data. 
  #Peak times then sorted in ascending order. 
  sample1.normativeTimes <- sort(ceiling(head(sample1.dens.df[order(-sample1.dens.df$peak, -sample1.dens.df$Density),],n.norm.bp)$Time))
  sample2.normativeTimes <- sort(ceiling(head(sample2.dens.df[order(-sample2.dens.df$peak, -sample2.dens.df$Density),],n.norm.bp)$Time))
  
  #Calculate the average distance between every peak in sample 1 to the nearest peak in sample 2 by getting the minimum absolute distance between each peak of sample 1 with every peak in sample 2 and vice versa. 
  nearest.dist <- data.frame()
  for(i in 1:length(sample1.normativeTimes)){
    nearest.dist[i,1] <- min(abs(sample1.normativeTimes[i] - sample2.normativeTimes))
    nearest.dist[i,2] <- min(abs(sample2.normativeTimes[i] - sample1.normativeTimes))
  }
  
  sample1.normativeLength <- mean(diff(sample1.normativeTimes))
  
  near.dist <- data.frame(nearestDistance = mean(c(nearest.dist[,1], nearest.dist[,2])), 
                          sample_nPeaks = length(sample1.normativeTimes), ave_samplePeakInterval = mean(diff(sample1.normativeTimes)), sd_samplePeakInterval = sd(diff(sample1.normativeTimes)),
                          test_nPeaks = length(sample2.normativeTimes), ave_testPeakInterval = mean(diff(sample2.normativeTimes)), sd_testPeakInterval = sd(diff(sample2.normativeTimes)),
                          sample_norm_boundary = n.norm.bp)
  
  return(near.dist)
}

##Peakiness

get.peakiness = function(y.vect, mov.dur, sub.list, adj.size, bw){
  #####
  # Adapted from script written: programmed 06.09.16 by Khena Swallow
  # Changed: take bandwidth input 
  #
  # Returns: minimum and actual rugosity values (Peakiness not yet calculated from this code chunk. Peakiness reported in the manuscript was calculated in esMethods_ManuscriptAnalysis.rmd)
  # Minimum rugosity calculated:
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
  #Parameter adj.size: size of bandwidth adjustmet to be applied in the density function
  #Precondition adj.size: numeric
  #
  #Parameter bw: size of bandwidth to be applied to the density function 
  #Precondition: numeric
  #
  #####
  # first, get the distribution of button presses for the minimum
  total.bps = length(y.vect)
  even.bps = seq(from = 0, to = mov.dur, length.out = total.bps) 
  
  #get density distribution for minimum and actual data
  min.y.dens = density(even.bps, bw = bw, adjust = adj.size, kernel = 'g', n = mov.dur)
  y.dens = density(y.vect, bw = bw, adjust = adj.size, kernel = 'g', n = mov.dur)
  
  #calculate peakiness using rugo function
  min.rugo = rugo(min.y.dens$y)
  act.rugo = rugo(y.dens$y)
  
  return(c(min.rugo, act.rugo))
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

get.surprise.index = function(sub.bp, gp.bp, mov.dur, ave.bp, bw, adj.size){
  ###
  #Returns surprise index (adapted from Katori et al., ) = frequency that a subject's button press coincides with the rest of the group's button press controlled for individual number of button presses. 
  #
  #Parameter sub.bp = array of presence of button press for n-second bin size of the movie for a single participant (x). 
  #Precondition sub.bp = numeric
  #
  #Parameter gp.bp = array of average button presses for every n-second bin size of the movie  for the group data excluding participant x.
  #Precondition gp.bp = numeric
  #
  #Parameter timeSeries = time series of movie length spaced out based on bin size. 
  #
  #Paramete mov.dur = total movie duration in bins. 
  #Precondition mov.dur = numeric
  #
  #Parameter ave.bp = average button press of the group excluding participant x. 
  #Precondition ave.bp = numeric
  #
  #calculate proportion of subject button presses
  sub.bp.times <- unique(floor(sub.bp))
  sub.eventrate <- length(sub.bp)/mov.dur
  
  if (length(gp.bp) > 1){
    gp.dens = density(gp.bp, bw = bw, adjust = adj.size, kernel = 'g', n = mov.dur)
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
  
  gp.peakTimes <- gp.dens.df$Time[gp.dens.df$peak == 1]
  
  if(length(gp.peakTimes) > ave.bp){
    n.normativePeaks <- ave.bp
  } else if (length(gp.peakTimes) <= ave.bp){
    n.normativePeaks <- length(gp.peakTimes)
  }
  
  gp.normativeTimes <- sort(floor(head(gp.dens.df[order(-gp.dens.df$peak, -gp.dens.df$Density),],n.normativePeaks)$Time))
  
  gp.eventrate <- length(gp.normativeTimes)/mov.dur
  
  #calculate number of overlapping button presses (between subject and rest of sample)
  n <- seq(from = 0, to = sum(sub.bp.times %in% gp.normativeTimes), by = 1)
  
  #calculate surprise index 
  lambda <- sub.eventrate*gp.eventrate*mov.dur
  rarity <- 1-sum(exp(-lambda)*((lambda^n)/factorial(n)))
  sup.index <- -log(rarity, base = 2)
  
  return(sup.index)
}


#To get values to needed for calculation of predictive value and detection accuracy 

get.normative.hitRate_bin <- function(sub.bp, gp.bp,ave.bp, timeSeries){
  ###
  #Returns:
  #  1. normativeHit = hit rate (proportion of indiv button press overlapping with normative boundary) 
  #  2. normativeFA = false alarm rate (proportion of indiv button press non-overlapping with non-normative boundary)
  #  3. normativeMiss = misses (proportion of indiv button press non-overlapping with normavie boundary)
  #  4. normativeCR = correct rejection (proportion of indiv bon button press overlapping with non-normative boundary)
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
  #Determine normative hit rate (sub.bpTimes overal with gp.bpTimes/ total gp.bpTimes)
  normativeHitRate <- sum(sub.bpTimes %in% gp.normativeTimes$timeSeries)/length(gp.normativeTimes$timeSeries)
  normativeFalseAlarm <- sum(sub.bpTimes %in% gp.nonBoundaryTimes)/length(gp.nonBoundaryTimes)
  normativeMiss <- (length(gp.normativeTimes$timeSeries)-sum(sub.bpTimes %in% gp.normativeTimes$timeSeries))/length(gp.normativeTimes$timeSeries)
  normativeCorrectRejection <- (length(gp.nonBoundaryTimes)-sum(sub.bpTimes %in% gp.nonBoundaryTimes))/length(gp.nonBoundaryTimes)
  normativePPV <- sum(sub.bpTimes %in% gp.normativeTimes$timeSeries)/length(sub.bpTimes)
  normativeFalsePPV <- sum(sub.bpTimes %in% gp.nonBoundaryTimes)/length(sub.bpTimes)
  
  return(data.frame('normativeHit' = normativeHitRate, 'normativeFA' = normativeFalseAlarm, "normativeMiss" = normativeMiss, "normativeCR" = normativeCorrectRejection, "normativePPV" = normativePPV, "normativeFalsePPV" = normativeFalsePPV))
}


