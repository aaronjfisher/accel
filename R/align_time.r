# !diagnostics off


#' @export
varb <- function(x) var(c(x)) * (length(x)-1)/length(x) #biased variance
shift_forward_by <- function(x, i, L = length(x)){
  #this corresponds to *subtracting* i units from a time series
  # (moving it forward in time)
  if(i== 0) return(x)
  if(i > 0) return(c(x[(i+1):L], x[1:i]))
  if(i < 0) return(c(x[(L-(abs(i))+1):L], x[1:(L-abs(i))]))
}


#returns amount to shift a
align_a2b <- function(a, b,
  shifts = 1:length(a) - round(length(a)/2)){

  L=length(shifts)
  shifted_as <-  matrix(NA, L, L)
  dimnames(shifted_as) <- list(time=names(a), 'shift'=shifts)
  for(i in 1:L){
    shifted_as[,i] <- shift_forward_by(a, shifts[i])
  }
  # image(shifted_as)
  # note, if we multiply by a "shift array," it is substantially slower.
  # the shifting is substantially slower than the actual alignment check!
  # it's not a bottleneck.
  
  agreement <-  t(shifted_as) %*% b
  shift_inds <- apply(agreement, 2, which.max)
  shift_a_values <- shifts[shift_inds]
  
  #hist(shift_a_values, breaks=20)
  # shift_a_values <- shift_a_values - round(mean(shift_a_values)) #to avoid drift; shifting all values doesn't change the fit. #!!! must be moved later
  
  names(shift_a_values) <- names(shift_inds)
  
  shift_a_values
}

#' @export
align_time <- function(m, 
  plot = FALSE,
  verbose = FALSE,
  max_iter = 20,
  get_lwr_bound = TRUE,
  max_combn = 50000,
  warn_max_iter = TRUE){
  
  ##### Process m
  naCol <- apply(m, 2, function(z) any(is.na(z)))
  m <- m[,!naCol]
  overall_mean <- mean(m)
  m <- m - overall_mean #de-mean
  
  L <- ncol(m)
  n <- nrow(m)
  shifts <- 1:L - round(L/2)
  shift_m_values <- rep(0,n)
  shift_m_values_mat <- as.data.frame(t(shift_m_values))
  
  ##### Create objects for loop
  
  shiftDiffs <-
    muDiffs <-
    varbMus <- rep(NA, max_iter+1)
  muDiffs[1] <- 
    varbMus[1] <- varb(colMeans(m)) 
  shiftDiffs[1] <- n
  
  m_iter <- m
  ##### Iterate to find mu and time shifts
  if(verbose) pb <- timerProgressBar(0,max=n)
  for(iter in 2:max_iter){
    
    for(i in 1:n){
      shift_m_values[i] <- align_a2b(m[i,], colMeans(m_iter[-i,]), shifts=shifts)
      # Always shift from m, not m_iter, so shift has constant interpretation. 
         # Shift m *towards* m_iter.
      m_iter[i,] <- shift_forward_by(m[i,], shift_m_values[i])
      #don't update m! otherwise shift values will change.
    }
    shift_m_values_mat[iter,] <- shift_m_values
    
    varbMus[iter] <- varb(colMeans(m_iter))
    shiftDiffs[iter] <- sum(shift_m_values != shift_m_values_mat[iter-1,])

    if(verbose) setTimerProgressBar(pb, n - shiftDiffs[iter])
    if(shiftDiffs[iter] == 0) break
  }
  
  
  ### adjust for drift
  shift_m_values <- shift_m_values - round(mean(shift_m_values))
  for(i in 1:n){
    m_iter[i,] <- shift_forward_by(m[i,], shift_m_values[i])
  }
  
  mean_varb_resid <- list(
    random = varb(c(m)),
    init = mean(apply(m,2,varb)),
    proposed = mean(apply(m_iter,2,varb)),
    lwr_bound = NA
  )
  if(abs(mean(apply(m,2,varb)) + varb(colMeans(m)) - varb(c(m))) > 10^-9|
      abs(mean(apply(m_iter,2,varb)) + varb(colMeans(m_iter)) - varb(c(m))) > 10^-9 ){
    stop('variance does not add to one.') 
    #!! to do -- change varb to var? what did you miss here?
    #!! why doesn't var work? it should? n-1 correction changing for length of colmeans(m) versus m itself? Does it matter?
  }
  
  if(get_lwr_bound){
    ncb <- choose(n, 2)
    if(ncb <= max_combn) cb <- t(combn(n, 2))
    if(ncb >  max_combn){
      cb <- data.frame(
        V1 = sample(1:n, choose(500, 2), replace=TRUE),
        V2 = sample(1:n, choose(500, 2), replace=TRUE)
      ) %>%
        filter(V1 < V2) %>% 
        mutate(id = paste0(V1,'-',V2)) %>%  #must come after V1 < V2 sort!
        filter(!duplicated(id)) %>% 
        select(-id)
    }
    
    rowDiffsBest <- rep(NA, nrow(cb))
    for(i in 1:nrow(cb)){
      a <-  m[cb[i,1],]
      b <-  m[cb[i,2],]
      shift_i <- align_a2b(a, b)
      a2 <- shift_forward_by(a, shift_i)
      rowDiffsBest[i] <- mean((a2-b)^2)
    }
    mean_varb_resid$lwr_bound <- (1/2) * mean(rowDiffsBest) * (n-1)/n #convert to var, and then to biased var (varb).
    
  }
  
  if(plot){
    par(mfrow=c(2,1))
    matplot(t(m + overall_mean), type='l', main=paste('unaligned\nvar explained =', round(varbMus[1]/varb(m), digits=3)),
      ylab='values', col=set2)
    lines(colMeans(m) + overall_mean, lwd=3, col='gray35')
    matplot(t(m_iter +overall_mean), type='l',
      ylab='values', col=set2,
      main=paste0('aligned\nvar explained = ', round(varbMus[iter]/varb(m), digits=3),'\n(',iter,' iterations)'))
    lines(colMeans(m_iter) +overall_mean, lwd=3, col='gray35')
  }
  
  if(any(diff(varbMus[1:iter]) < 0)) stop('variance explained by mean went down.')
  
  flag <- 'none'
  if(iter==max_iter){
    if(warn_max_iter) warning('max_iter reached')
    flag <- 'max_iter reached'
  }
  return(list(
    m_aligned = m_iter + overall_mean,
    mu_aligned = colMeans(m_iter) + overall_mean,
    shift_m_values = shift_m_values,
    mean_varb_resid = mean_varb_resid,
    iter=iter,
    flag=flag))
}





#!! This batch approach does not seem to improve things! none of the extra steps I've added onto the baseline approach seem to help! I think there may be improvements to be made in the lwr_bound, but perhaps not in the search?
batch_align <- function(m, batch_iter=50, batch_size, get_lwr_bound=FALSE, ...){
  browser()
  n <- nrow(m)
  if(batch_size >= n) stop('batch size too large.')
  
  shift_params <- matrix(0, batch_iter, n)
  varbs <- rep(NA, batch_iter)
  
  shift_params[1,] <- 0 #initialize
  varbs[1] <- mean(apply(m,2,varb))
  m_update <- m
  
  standard_align <- align_time(m, get_lwr_bound = TRUE, ...)
  
  pbt <- timerProgressBar(0, batch_iter)
  for(i in 2:batch_iter){
    ind <- sample(n, batch_size)
    ai <- align_time(m_update[ind,], get_lwr_bound = FALSE, ...)
    
    shift_params[i,ind] <- shift_params[i-1,ind] + ai$shift_m_values
    m_update[ind,] <- ai$m_aligned
    varbs[i] <- mean(apply(m_update,2,varb))
    setTimerProgressBar(pbt, i)
  }
  #plot(varbs, type='l')
  #matplot(shift_params, type='l')
  
  output <- 
  ai_final <- align_time(m_update, get_lwr_bound = get_lwr_bound, ...)
  
  shift_final <- (shift_params[batch_iter,] + ai_final$shift_m_values) %% ncol(m)
  shift_final <- (shift_final - round(mean(shift_final))) %% ncol(m)
  ai_final$mean_varb_resid
  output$shift_m_values <- shift_final
  
  # output$mean_varb_resid$proposed / standard_align$mean_varb_resid$proposed
  output
}




