
##################
##################
# my plot functions

#' shortcut function for pdfs
#' @export
mypdf <- function(suffix, ...) pdf(paste0('plots/',Sys.Date(),'_',suffix,'.pdf'), ...)

#' shortcut function for pngs
#' @export
mypng <- function(suffix, ...) png(paste0('plots/',Sys.Date(),'_',suffix,'.png'), ...)

##################
##################
#' Fitting an HMM to an individual (shortcut function)
#' @param y activity; encouraged to be log(activity + 1)
#' @import depmixS4
#' @export
fit_depmix_individual <- function(y, nstates=3, family=gaussian(), ntimes=length(y), ...){
  depmixS4::depmix(y ~ 1, 
    nstates=nstates, 
    family=family,
    ntimes=ntimes, ...) %>% fit
}

#' @import depmixS4
#' @export
state_order <- function(mod, y){
  v <- depmixS4::viterbi(mod)
  v$y <- y
  (group_by(v, state) %>% 
    summarize(m = mean(y, na.rm=TRUE)))$m %>% order
}

#' Get object from HMM model to plot
#' @import depmixS4
#' @export
get_viterbi_with_columns_sorted_by_state <- function(mod, y, include_mean_in_level_name=FALSE){
  v <- depmixS4::viterbi(mod)
  if(length(unique(v$state)) < ncol(v)-1)
  stop('empty state')

  v$y <- y


  #now, change "S1, S2..." state names to "L1, L2,...," which are ordered by within-group mean.
  means <- group_by(v, state) %>% 
    summarize(m = signif(mean(y, na.rm=TRUE), 3))
  means <- means[order(means$m),]
  means$level <- paste0('L',1:nrow(means))
  #means

  # Replace S labels with L labels
  v$level <- NA
  for(j in 1:nrow(means)){
    Sj <- paste0('S',means$state[j])
    Lj <- means$level[j]
    if(include_mean_in_level_name) Lj <- paste0(Lj, ' (',means$m[j],')')
    v[[Lj]] <- v[[Sj]]
    v$level[v$state == means$state[j]] <- j

    #delete "state" labels
    v[[Sj]] <-NULL
  }
  v$state <- NULL 
  v$y <- NULL

  v
}

#' plot results from get_viterbi_with_columns_sorted_by_state
#' @import tidyr ggplot2
#' @export
plot_sorted_viterbi <- function(
  vit, #output of 
  y, 
  time_bin, 
  index_to_plot = 1:length(y),
  groups=rep(1,length(y)),
  ylab='probability',
  xlab='time bin',
  means_included_in_level_names = FALSE,
  filllab=paste0('states',
    ifelse(means_included_in_level_names,' (within-\nstate mean)',''))
  ){
  

  nstates <- length(unique(vit$level))
  vit$level <- NULL
  
  # To draw the top of a polygon, we go from left to right.
  top_left2right <- t(apply(vit,1,cumsum)) %>% as.data.frame #do not order states
  top_left2right$groups <- groups[index_to_plot]
  top_left2right$time_bin <- time_bin[index_to_plot]

  id_avgs <- data.frame(y=y, groups=groups)[index_to_plot,]%>%
    group_by(groups) %>%
    summarize(mean_y = mean(y, na.rm=TRUE))
  top_left2right$lab_id <- NA
  for(i in 1:nrow(top_left2right)){
    lab_ind <- which(id_avgs$groups == top_left2right$groups[i])
    top_left2right$lab_id[i] <- paste0(
      'id = ',
      id_avgs$groups[lab_ind],
      ';\n avg y = ',
      signif(id_avgs$mean_y[lab_ind], 3))
  }
  top_left2right$lab_id <- factor(top_left2right$lab_id)
  
  # For drawing the underside of polygons, we have to go backwards
  # along the plot, from right to left.
  bottom_right2left <- top_left2right[nrow(top_left2right):1,]
  bottom_right2left[,1] <- 0
  bottom_right2left[,2] <- top_left2right[nrow(top_left2right):1,1] #bottom edge of level 2 is the top edge of level 1
  bottom_right2left[,3] <- top_left2right[nrow(top_left2right):1,2]

  sfill <- scale_fill_brewer(palette ='RdYlBu')
  if(nstates==2) sfill <- scale_fill_manual(values=c('#fc8d59','#91bfdb'))
  rbind(top_left2right, bottom_right2left) %>%
    pivot_longer(c(1:nstates), names_to='level',values_to='prob') %>% 
    ggplot(aes(x = time_bin, y=prob, fill=level)) +
    geom_polygon(col='black') +
    facet_grid(lab_id~.) + theme_bw() +
    sfill +
    labs(y=ylab,x=xlab, fill=filllab)
  
  # r <- nrow(vit)
  # plot(c(),xlim=c(1,r),ylim=c(0,1),ylab='State probabilities')
  # polygon(c(1:r,r:1),c(vitc[,1],rep(0,r)),col='#fc8d59')
  # polygon(c(1:r,r:1),c(vitc[,2],vitc[r:1,1]),col='#ffffbf')
  # polygon(c(1:r,r:1),c(vitc[,3],vitc[r:1,2]),col='#99d594')
  # #http://colorbrewer2.org/#type=diverging&scheme=Spectral&n=3
  
}
####








### AF - this was initially Dmitri's code. AF added comments.
### return phase in radians from the amplitudes of cos = A and sin = B components
#### use A*cos(omega*t-phase) = A*cos(phase)*cos(omega*t) + A*sin(phase)*sin(omega*t)
phase_trans <- function(cos.coef, sin.coef){
  
  ### calculate amplitude A
  #amp <- sqrt(sin.coef^2 + cos.coef^2)

  ## calculate phase parameter, which depends on which quardrant it is in
  ## ultimately we want to constraint it to (0, 2*pi) range
  div <- sin.coef/cos.coef
  phaser <- NA
  ### if in quardrant 1
  if (cos.coef >= 0 & sin.coef > 0) {
    phaser <- atan(div)
  }
  ### if in quardrant 2
  if (cos.coef < 0 & sin.coef >= 0) {
    phaser <- atan(div) + pi
  }
  ### if in quardrant 3
  if (cos.coef <= 0 & sin.coef < 0) {
    phaser <- atan(div) + pi
  }
  ### if in quardrant 4
  if (cos.coef > 0 & sin.coef < 0) {
    phaser <- atan(div) 
  }
  if(is.na(phaser)) stop('phaser is NA')
  ### probably redundant safe check
  if (phaser < 0) {
    phaser <- phaser + (2 * pi)
  }
  if (phaser > (2 * pi)) {
    phaser <- phaser - (2 * pi)
  }
  
  names(phaser) <- NULL
  phaser
}


#' Fit a cosine wave
#' @param y vector representing a value that changes over time 
#' @param x time variable. must be in units of cycles (x ranging from 0 to 1 denotes 1 cycle). multiple cycles should be treated as "overlaid" on the same 0-1 interval.
#' @param plot_it If TRUE, a plot of the fit is shown
#' @export
fit_01_wave <- function(y, x, plot_it=FALSE,...){ 
  # Note, this function treats "tomorrow at 4pm" the same
  # as "today at 4pm"
  
  if(any(x<0 | x>1)) stop('x must be between 0 and 1')
  lmfit <- lm(y ~ sin(x*2*pi) + cos(x*2*pi))
  
  if(any(is.na(lmfit$coefficients)) |
    var(y) == 0){
    warning('NAs from lm fit procedure')
    noise_sd <- phase <- amplitude <- intercept <- NA
  }else{
    amplitude <- sqrt(sum(lmfit$coefficients[-1]^2))
    intercept <- lmfit$coefficients[1]
    phase <- phase_trans(
      sin.coef = lmfit$coefficients[2],
      cos.coef = lmfit$coefficients[3]) / (2*pi) #transform to 0,1 range 
    noise_sd <- sd(lmfit$residuals)
  }

  

  if(plot_it){
    #Note: The last few hours shown (evening) 
    #actually correspond to the *previous day*, although
    #this shift doesn't affect the amplitude.
    plot(x=x, y=y, ...)
    
    
    plotx <- seq(0,1, length=1000)
    ploty <- cbind(1,  
                   sin(plotx*2*pi), 
                   cos(plotx*2*pi)) %*% lmfit$coefficients
    lines(x=plotx, y=ploty, col='blue', lwd=2)
    
    xmax <- plotx[which(ploty==max(ploty))]
    xmin <- plotx[which(ploty==min(ploty))]
    segments(x0=c(xmax,xmin),
             x1=c(xmax,xmin), 
             y0=rep(lmfit$coefficients[1],2),
             y1=lmfit$coefficients[1] + c(1,-1)*amplitude, lty=3)
    abline(h=lmfit$coefficients[1], lty=2)
  }
  
  return(list(
    'intercept'=intercept, 
    'amplitude'=amplitude,
    'phase'=phase,
    'noise_sd'=noise_sd))
}


#' @export
min1 <- as.difftime(1, units='mins')

#' @export
hr1 <- as.difftime(1, units='hours')


# Function to find 10hr max
# x = vector to take max over
# posix = time variable
#' @export
fit_max10hr <- function(x, posix, min_per_epoch_num, plot_it=FALSE){ 
  if(length(x) != length(posix)) 'Each element of x must have a corresponding time value in posix'
  if(length(x) > 60*24 / min_per_epoch_num) stop('too many epochs')
  if(sum(!is.na(x)) < 0.8 *  60*24 / min_per_epoch_num) return(NA) #not enough non-NAs

  hr10 <- 10*60*min1
  valid_starts <- which(posix + hr10 - min1*min_per_epoch_num <= max(posix)) #don't count last minute, first min is inclusive
  L <- length(valid_starts)

  means <- rep(NA,L) #epochs to consider
  for(i in 1:L){
    start_time_i <- posix[valid_starts[i]]
    ind_i <- which(posix >= start_time_i &
                   posix <  start_time_i + hr10)
    means[i] <- mean(x[ind_i], na.rm=TRUE)
  }
  max_mean <- max(means, rm.NA=TRUE)
  if(plot_it){
    plot(x=posix[valid_starts],y=means, type='l',
         xlab='Start times', 
         ylab='')
    mtext('Activity count in 10hr window\n(sum of vector magnitude)',2,2)
    abline(h=max_mean, lty=2)
  }  
  return(max_mean)
} #Code here is unnecessarily slow, and not scalable, 
    # but not so bad for a few subjects. Memory requirements are not high.
    # could improve speed here, but probably not worth the writing
    # time as it runs in about 2.5 seconds / patient.
	# See next function.


#time_length = how long a time window. 1*hr1 for M10. hours(5) for L5.
#' @export
fit_minmax_window <- function(x, posix, minmax, time_length, plot_it=FALSE){ 
  
  valid_starts <- which(posix + time_length < max(posix)) #don't count last minute, first min is inclusive
  L <- length(valid_starts)

  means <- rep(NA,L) #epochs to consider
  for(i in 1:L){
    start_time_i <- posix[valid_starts[i]]
    ind_i <- which(posix >= start_time_i &
                   posix <  start_time_i + time_length)
    means[i] <- mean(x[ind_i], na.rm=TRUE)
  }
  minmax_mean <- minmax(means, na.rm=TRUE)
  if(plot_it){
    plot(x=posix[valid_starts],y=means, type='l',
         xlab='Start times', 
         ylab='')
    mtext('Activity in W window',2,2)
    abline(h=minmax_mean, lty=2)
  }  
  return(minmax_mean)
} #Code here is unnecessarily slow, and not scalable, 
    # but not so bad for a few subjects. Memory requirements are not high.
    # could improve speed here, but probably not worth the writing
    # time as it runs in about 2.5 seconds / patient.
  # See next function.

# #' @export
# Mode <- function(x) {
#   ux <- unique(x)
#   ux[which.max(tabulate(match(x, ux)))]
# }#https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
  #problem: in some cases this converts things (like difftime) to numeric
  #I would like to keep the formatting.
  
#' @export
Mode <- function(x){ 
  tt <- table(x)
  maxtt <- names(tt)[which.max(tt)]
  x[which(x==maxtt)[1]]
}


#' Get M10, L5 or other window summaries, in linear time
#' For L5 M10 - the smallest, or largest, window of size W, where W is in difftime units.
#' This faster version tracks a moving window, rather than summing over all elements in the window each time. Additional conserations are given to adjust for NAs.
# posix = time variable
#' @export
fit_fast_minmax_window <- function(
  x, posix,
  minmax, #max or min function
  time_length = 1*hr1, #in difftime
  return_type = 'value',
  nonzero_min_proportion = 0.1,
  max_gap = 2 * hr1 #largest time gap allowed.
  ){  

  # Consider a window x[a:(a+w)] of size w, and suppose we have previously calculated sum(x[a:(a+w)])
  # when we move the window 1 position to the right, from x[a:(a+w)] to x[(a+1):(a+w+1)], 
  # we add xa+w+1] and subtract x[a].


  #!! period seems to not respect daylight savings, and be easier to add to.
  # difftime does appear to respect this, and be hard to add to.

  ######################## Error checks part 1  
  if(any(is.na(posix))) stop('NA posix')
  if(length(x) != length(posix)) stop('Each element of x must have a corresponding time value in posix')
  if(length(x) <= 2){
    warning('<=2 values in x; returning NA')
    return(NA)
  }
  if(any(duplicated(posix))) stop('duplicated epochs')
  test_x <- x
  test_x[is.na(test_x)] <- 0
  if(mean(test_x > 0) < nonzero_min_proportion) return(NA) #not enough nonzeros
  ########################



  #replace posix with daylight savings adjusted version, in UTC, where there is no daylight savings.
  ftz <- force_tz(posix,'UTC')
  timelag <- ftz - posix
  if(any(timelag != timelag[1])){
    unadjusted_posix <- posix
    dposix <- as.period(diff(unadjusted_posix))
    cumsum_dposix <- dposix
    cumsum_dposix[2:length(dposix)] <- NA
    for(j in 2:length(dposix)) cumsum_dposix[j] <- cumsum_dposix[j-1] + dposix[j]
    daylight_savings_adjusted_posix <- c(ftz[1],ftz[1] + cumsum_dposix)
    posix <- daylight_savings_adjusted_posix
  }

  ######################## Error checks part 2
  dposix <- diff(posix)
  
  if(max(dposix) > max_gap){
    warning('gaps as large as ',max(dposix),' are present; returning NA')
    return(NA)
  }
  if(any(dposix<0)) stop('time must be ordered')
  if(any(dposix<=0)) stop('invalid times')
  if(diff(range(posix)) > hours(24) & all(timelag == timelag[1])) stop('too many epochs given no daylight savings')
  if(diff(range(posix)) > hours(25)) stop('too many epochs regardless of daylight savings')
  ########################
  
  lenx <- length(x)
      
  # window will go from posix[start_i] to posix[end_i],
  # averaging over non-na values in x[start_i:end_i].
  # note, the time difference posix[end_i] - posix[start_i] 
  # must be as high as possible while strictly less than time_length,
  # since end points do not have width = 0.
  # in other words, we sum over x for posix in [posix[start_i], posix[start_i] + time_length), 
  # where upper interval endpoint is closed.

  #Itialize a starting & ending point
  best_start <-
  start_i <- 1
  best_end <-
  end_i <- max(which(posix < posix[start_i] + time_length))
  if(length(end_i) == 0){
    warning('no valid intervals, returning NA'); return(NA)
  }
  sum_j <- sum(x[start_i:end_i], na.rm=TRUE)
  n_j <- sum(!is.na(x[start_i:end_i]))
  mean_j <- sum_j / n_j
  minmaxW <- mean_j
  if(end_i == lenx){
    warning('only 1 valid interval')
    return(mean_j)
  }

  j <- 1
  while(j < 3*lenx){
    #move up either window endpoint, depending on if there is room to grow.
    
    #print(c(start_i,end_i,posix[end_i] - posix[start_i]))
 
    if( (posix[end_i+1] - posix[start_i] < time_length)
        & end_i < lenx){ 
      #is there room to move up?
      #if so, extend
      end_i <- end_i + 1
      add_j <- x[end_i]
      if(is.na(add_j)){
        add_j <- 0
      }else{ n_j <- n_j + 1 }
      sum_j <- sum_j + add_j

    }else{
      #update minmaxW, and increase left index to move to next valid start.
      mean_j <- sum_j / n_j
      if(minmax(minmaxW, mean_j) == mean_j){
        best_start <- start_i
        best_end <- end_i
      }
      minmaxW <- minmax(minmaxW, mean_j)
      
      if(end_i==lenx) break()

      #else, move to next starting point
      sub_j <- x[start_i]
      start_i <- start_i+1
      if(is.na(sub_j)){
        sub_j <- 0
      }else{ n_j <- n_j - 1 }
      sum_j <- sum_j - sub_j

    }
    if(n_j <= 0) stop('zero count')

    j = j + 1
 
  }
  if(j == 3*lenx) warning('search concluded without an answer')

  if(return_type=='full_indeces') return(best_start:best_end)
  if(return_type=='endpoints') return(c('start'=best_start, 'end'=best_end))
  if(return_type=='value') return(minmaxW)
  if(return_type=='value+endpoints') return(c(
    'value'=minmaxW,'start'=best_start, 'end'=best_end))
  if(return_type=='full_values') return(x[best_start:best_end])

  stop('invalid return_type')
} 


#Test
if(FALSE){ #example / test
  #### even spaced
  n <- 60*24
  posix <- ymd_hms('2018/12/31 00:00:00') + 
      as.difftime(seq(1,n,length=n), units='mins')
  x <- rnorm(n) + sin(seq(-pi,pi, length=n))*2  
  x <- x + abs(min(x))
  plot(x=posix, y=x)

  #MAX
  system.time({
    fit_slow <- fit_minmax_window(x=x, posix=posix, minmax=max, time_length=8*hr1, plot_it=TRUE)
  })
  system.time({
    fit_fast <- fit_fast_minmax_window(x=x, posix=posix, minmax=max, time_length=8*hr1)
  })
  if(abs(fit_fast - fit_slow)  > 10^-10) stop('Error in fit_fast_minmax_window')
  
  #MIN
  system.time({
    fit_slow <- fit_minmax_window(x=x, posix=posix, minmax=min, time_length=8*hr1, plot_it=TRUE)
  })
  system.time({
    fit_fast <- fit_fast_minmax_window(x=x, posix=posix, minmax=min, time_length=8*hr1)
  })
  if(abs(fit_fast - fit_slow)  > 10^-10) stop('Error in fit_fast_minmax_window')


  #uneven spacing
  si <- sample(n, floor(n/4), replace=FALSE)
  si <- si[order(si)]
  hist(as.numeric(diff(posix[si])))
  #MAX
  system.time({
    fit_slow <- fit_minmax_window(x=x[si], posix=posix[si], minmax=max, time_length=8*hr1, plot_it=TRUE)
  })
  system.time({
    fit_fast <- fit_fast_minmax_window(x=x[si], posix=posix[si], minmax=max, time_length=8*hr1)
  })
  if(abs(fit_fast - fit_slow)  > 10^-10) stop('Error in fit_fast_minmax_window')
  
  #MIN
  system.time({
    fit_slow <- fit_minmax_window(x=x[si], posix=posix[si], minmax=base::min, time_length=8*hr1, plot_it=TRUE)
  })
  system.time({
    fit_fast <- fit_fast_minmax_window(x=x[si], posix=posix[si], minmax=base::min, time_length=8*hr1)
  })
  abs(fit_fast - fit_slow)
  (fit_slow)
  if(abs(fit_fast - fit_slow)  > 10^-10) stop('Error in fit_fast_minmax_window')

}












################################################################
################# ONLY FOR 1001, not 2001!!!! ##################

#' Shortcut function Aaron used briefly. 
#' 
#' Not maintained. (newer code should use Date rather than year + yday.
#' @param df is a data frame with columns: SubjectID, posix, and activity
#' @param hr_delay time since midnight when a "new day" starts. For example, for max10 hours, we may wish day-groups to start at 5am (hr_delay = 5) rather than midnight. Alternatively, for total sleep time, we may wish days to start at 6pm (hr_delay = 18) when computing total sleep time.
#' @export
#' @return a data frame with: SubjectID, group_year, group_yday, and [summary measure] 
group_by_shifted_day <- function(df, hr_delay = 0, drop_group_posix=TRUE,group_output=TRUE, uday=FALSE){
  
  if(any( c('group_yday','group_year','group_posix') %in% names(df) )) warning('group_yday, group_year, and/or group_posix present in df names, and will be overwritten. group_posix will be removed.')
  if(!is.numeric(hr_delay)) stop('invlaid hr_delay')
  
  if(hr_delay > 24 | hr_delay<0) stop('invlaid hr_delay')
  if(hr_delay > 12) hr_delay <- hr_delay - 24
  difftime_delay <- as.difftime(hr_delay, units='hours')

 if(group_output) df <- ungroup(df)

  out <- mutate(df,
      group_posix = posix - difftime_delay,
      group_year = year(group_posix),
      group_yday = yday(group_posix)
    )
      

  if(group_output) out <- group_by(out, SubjectID, group_year, group_yday)
  if(uday) out <- mutate(out, uday = paste0(group_year, group_yday))
  if(drop_group_posix | uday) out <- dplyr::select(out, -group_posix)
  
  out
}

if(FALSE){ #example / test
  test.time <- data.frame(SubjectID = NA, posix = ymd_hms('2018/12/31 12:00:00') + 
    as.difftime(0:24, units='hours')) %>%
    mutate(
      original_yday = yday(posix),
      original_year = year(posix))
  group_by_shifted_day(test.time, hr_delay = 0) %>% as.data.frame() #midnight switch
  group_by_shifted_day(test.time, hr_delay = 4) %>% as.data.frame() #4am switch
  group_by_shifted_day(test.time, hr_delay = 18) %>% as.data.frame()#6pm switch  
}
