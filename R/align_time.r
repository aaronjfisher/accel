# !diagnostics off


#' @export
varb <- function(x) var(c(x)) * (length(x)-1)/length(x) #biased variance

#' Shift time series forward, wrapping around (optionally with multiple channels)
#' @param x a matrix with column per signal channel
#' @export
shift_forward_by <- function(x, i, L = dim(x)[1]){
  #this corresponds to *subtracting* i units from a time series
  # (moving it forward in time)

  if(length(dim(x)) != 2){
    stop('x must be a 2d array, or matrix, with 1 column per signal channel')
  }

  if(i== 0) return(x)
  if(i > 0) return( x[ c((i+1):L, 1:i), ] ) 
  if(i < 0) return( x[c((L-(abs(i))+1):L, 1:(L-abs(i))), ] )
}


# returns amount to shift a
# a & b should be matrices, with one column per channel.
# These should correspond to *single* observations (with possibly many channels),
# not multiple observations
align_a2b <- function(a, b,
  shifts = 1:dim(a)[1] - round(dim(a)[1]/2)){

  if(length(dim(a)) != 2 | length(dim(b)) != 2){
    stop('a and b must be 2d arrays, or matrices, with 1 column per signal channel')
  }
  if(any(dim(a) != dim(b)))stop('dimensions of a and b must match')

  L <- length(shifts)
  K <- dim(a)[2] #number of channels

  shifted_as <-  array(NA, dim=c(dim(a)[1], L, K))
  dimnames(shifted_as) <- list(time=dimnames(a)[1], 'shift'=shifts, 'channels'=dimnames(a)[2])
  for(i in 1:L){
    shifted_as[,i,] <- shift_forward_by(a, shifts[i])
  }
  # image(shifted_as)
  # note, if we multiply by a "shift array," it is substantially slower.
  # the shifting is substantially slower than the actual alignment check!
  # it's not a bottleneck.
  
  agreement <- matrix(NA, L, K) #shifts by channels
  for(k in 1:K){ #for each channel
    agreement[,k] <- crossprod(shifted_as[,,k], b[,keep(k)]) #if index is length 1, you need to use "keep" to prevent it from dropping
  }
  shift_ind <- which.max(rowMeans(agreement)) #total up agreement over all channels (!!unweighted mean!!)
  shift_a_value <- shifts[shift_ind]
  
  names(shift_a_value) <- names(shift_ind)
  
  shift_a_value
}

#' @import RColorBrewer
set2 <- brewer.pal(8, 'Set2')

#' Align two time series (optionally with multiple channels)
#' @param m should be a n x T x K array, where n is the sample size, T is the period length, and K is the number of channels. A n x T matrix can also be used.
#' @export
#' @import keep
align_time <- function(m, 
  plot = FALSE,
  verbose = FALSE,
  max_iter = 20,
  get_lwr_bound = TRUE,
  max_combn = 50000,
  shifts = 1:(dim(m)[2]) - round((dim(m)[2])/2), #must be integers!
  warn_max_iter = TRUE){
  
  if(length(dim(m))==2){
    warning('m should be 3-dimensional; converting m to a 3d karray.')
    m <- karray(m, dim=c(dim(m),1))
  }
  m <- as.karray(m) #don't drop extra dimensions!!

  ##### Process m
  naT <- apply(m, 2, function(z) any(is.na(z)))
  if(any(naT)){
    stop('NAs in m')
  }

  n <- dim(m)[1]
  L = length(shifts)
  K <- dim(m)[3]
  
  overall_means <- apply(m, 3, mean)
  add_constants_to_channel <- function(arr, constants){
    for(channel in 1:dim(arr)[3]){
      arr[,,channel] <- arr[,,channel] + constants[channel]
    } 
    arr
  }
  m <- add_constants_to_channel(arr=m, constants= -overall_means)


  shift_m_values <- rep(0,n)
  shift_m_values_mat <- as.data.frame(t(shift_m_values))
  
  ##### Create objects for loop
  get_varb_mus <- function(arr){
    #!! unweighted mean
    apply(arr, 2:3, mean) %>% #get within time-channel means
    apply(. , 2, varb) %>% #get within channel variance
    mean #average accross channels
  } 
  get_varb_noise <- function(arr){
    #!! unweighted mean
    apply(arr, 2:3, varb) %>% #get within time-channel variance
    mean #average accross channels & time
  } 



  shiftDiffs <-
    muDiffs <-
    varbMus <- rep(NA, max_iter+1)
  muDiffs[1] <- 
    varbMus[1] <- get_varb_mus(m)
  shiftDiffs[1] <- n
  


  m_iter <- m
  ##### Iterate to find mu and time shifts
  if(verbose) pb <- timerProgressBar(0,max=n)
  for(iter in 2:max_iter){
    
    for(i in 1:n){
      means_minus_i <- apply(m_iter[setdiff(1:n,i),,], 2:3, mean)
      shift_m_values[i] <- align_a2b(m[i,,], means_minus_i, shifts=shifts)
      # Always shift from m, not m_iter, so shift has constant interpretation. 
         # Shift m *towards* m_iter.
      m_iter[i,,] <- shift_forward_by(m[i,,], shift_m_values[i])
      #don't update m! otherwise shift values will change.
    }
    shift_m_values_mat[iter,] <- shift_m_values
    
    varbMus[iter] <- get_varb_mus(m_iter)
    shiftDiffs[iter] <- sum(shift_m_values != shift_m_values_mat[iter-1,])

    if(verbose) setTimerProgressBar(pb, n - shiftDiffs[iter])
    if(shiftDiffs[iter] == 0) break
  }
  
  
  ### adjust for drift
  shift_m_values <- shift_m_values - round(mean(shift_m_values))
  for(i in 1:n){
    m_iter[i,,] <- shift_forward_by(m[i,,], shift_m_values[i])
  }
  
  mean_varb_resid <- list(
    random = varb(c(m)),
    init = get_varb_noise(m),
    proposed = get_varb_noise(m_iter),
    lwr_bound = NA
  )
  
  if(abs(get_varb_noise(m)       + get_varb_mus(m)      - varb(c(m))) > 10^-9|
      abs(get_varb_noise(m_iter) + get_varb_mus(m_iter) - varb(c(m))) > 10^-9 ){
    stop('variance does not add to one.') 
    #!! to do -- change varb to var? what did you miss here?
    #!! why doesn't var work? it should? n-1 correction changing for length of colmeans(m) versus m itself? Does it matter?
  }
  
  if(get_lwr_bound){
    ncb <- choose(n, 2)
    if(ncb <= max_combn) cb <- t(combn(n, 2))
    if(ncb >  max_combn){
      cb <- data.frame(
        V1 = sample(1:n, 5*max_combn, replace=TRUE),
        V2 = sample(1:n, 5*max_combn, replace=TRUE)
      ) %>%
        filter(V1 < V2) %>% 
        mutate(id = paste0(V1,'-',V2)) %>%  #must come after V1 < V2 sort!
        filter(!duplicated(id)) %>% 
        dplyr::select(-c('id'))
      if(nrow(cb) > max_combn) cb <- cb[1:max_combn,]

      warning('n-choose-2 is larger than max_combn (',(choose(n,2)),' > ',
        max_combn,'), so lwr_bound will be approximated.')
    }
    
    # Below, we shift only one of the 2 vectors (a -> a2), but this shift needs
    # to account for any combination, or difference, of allowed shifts.
    # Thus, we use the outer produce of differences
    lwr_shifts <- c(outer(shifts, shifts, function(a,b) a-b)) %>% unique
    lwr_shifts <- lwr_shifts[order(lwr_shifts)]

    rowDiffsProposed <-  
    rowDiffsBest <- rep(NA, nrow(cb))
    for(i in 1:nrow(cb)){
      a <-  m[cb[i,1],,]
      b <-  m[cb[i,2],,]
      shift_i <- align_a2b(a, b, shifts=lwr_shifts)
      a2 <- shift_forward_by(a, shift_i)
      rowDiffsBest[i] <- mean((a2-b)^2) #unweighted mean!! (for multi-chanel, this is an avg of variances)

      #Workcheck, tested below:
      miter1 <- m_iter[cb[i,1],,]
      miter2 <- m_iter[cb[i,2],,]
      rowDiffsProposed[i] <- mean((miter1-miter2)^2)
    }
    
    convert_rowdiffs_to_bound <- function(dd) {
      #convert to var, and then to biased var (varb).
      (1/2) * mean(dd) * (n-1)/n 
    }
    mean_varb_resid$lwr_bound <- convert_rowdiffs_to_bound(rowDiffsBest)

        
    if(ncb <= max_combn){
      if(abs(
        convert_rowdiffs_to_bound(rowDiffsProposed)  -
        get_varb_noise(m_iter))> 10^-10){
          stop('error in lower bound conversion')
      }
      if(any(unlist(mean_varb_resid) < mean_varb_resid$lwr_bound)) stop('lower bound error')
    }
    
    if( (ncb > max_combn) & mean_varb_resid$proposed < mean_varb_resid$lwr_bound) warning('Approximate lower bound is higher than returned solution.')
    
    
  }

  m_plus_means <- add_constants_to_channel(arr=m, constants= overall_means)
  m_iter_plus_means <- add_constants_to_channel(arr=m_iter, constants= overall_means)
  if(plot){
    par(mfrow=c(2,1))
    matplot(t(m_plus_means[,,1]), type='l', main=paste('unaligned\nvar explained =', round(varbMus[1]/varb(m[,,1]), digits=3)),
      ylab='channel 1', col=set2)
    lines(apply(m_plus_means[,,1],2,mean), lwd=3, col='gray35')
    matplot(t(m_iter_plus_means[,,1]), type='l',
      ylab='channel 1', col=set2,
      main=paste0('aligned\nvar explained = ', round(varbMus[iter]/varb(m), digits=3),'\n(',iter,' iterations)'))
    lines(apply(m_iter_plus_means[,,1], 2, mean), lwd=3, col='gray35')
  }
  
  if(any(diff(varbMus[1:iter]) < 0)) stop('variance explained by mean went down.')
  
  flag <- 'none'
  if(iter==max_iter){
    if(warn_max_iter) warning('max_iter reached')
    flag <- 'max_iter reached'
  }
  return(list(
    m_aligned = m_iter_plus_means,
    mu_aligned = colMeans(m_iter_plus_means),
    shift_m_values = shift_m_values,
    mean_varb_resid = mean_varb_resid,
    iter=iter,
    flag=flag))
}





#!! This batch approach does not seem to improve things! none of the extra steps I've added onto the baseline approach seem to help! I think there may be improvements to be made in the lwr_bound, but perhaps not in the search?
batch_align <- function(m, batch_iter=50, batch_size, get_lwr_bound=FALSE, ...){
  n <- dim(m)[1]
  if(batch_size >= n) stop('batch size too large.')
  
  shift_params <- matrix(0, batch_iter, n)
  varbs <- rep(NA, batch_iter)

  shift_params[1,] <- 0 #initialize
  varbs[1] <- mean(apply(m, 2:3, varb))
  m_update <- m
  
  standard_align <- align_time(m, get_lwr_bound = TRUE, ...)
  
  pbt <- timerProgressBar(0, batch_iter)
  for(i in 2:batch_iter){
    ind <- sample(n, batch_size)
    ai <- align_time(m_update[ind,,], get_lwr_bound = FALSE, ...)
    
    shift_params[i,ind] <- shift_params[i-1,ind] + ai$shift_m_values
    m_update[ind,,] <- ai$m_aligned
    varbs[i] <- mean(apply(m_update, 2:3, varb))
    setTimerProgressBar(pbt, i)
  }
  #plot(varbs, type='l')
  #matplot(shift_params, type='l')
  
  output <- 
  ai_final <- align_time(m_update, get_lwr_bound = get_lwr_bound, ...)
  
  shift_final <- (shift_params[batch_iter,] + ai_final$shift_m_values) %% ncol(m)
  shift_final <- (shift_final - round(mean(shift_final))) %% ncol(m)
  #ai_final$mean_varb_resid
  output$shift_m_values <- shift_final
  
  # output$mean_varb_resid$proposed / standard_align$mean_varb_resid$proposed
  output
}




