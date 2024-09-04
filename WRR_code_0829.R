# Baseflow Separation
HYSEP <- function(Q, area_mi2, method=NULL){
  # R implementation of USGS HYSEP baseflow separation algorithms
  # as described in Pettyjohn & Henning (1979) and implemented 
  # in Sloto & Crouse (1996).
  #
  # Inputs:
  #   Q = discharge timeseries (no missing data) (any units are OK)
  #   area_mi2 = area in square miles
  #   method = HYSEP method; options are "fixed", "sliding", "local"
  #
  # Output:
  #   bf = baseflow timeseries, same length and units as Q
  
  ## package dependencies
  require(zoo)
  library(dplyr)
  require(magrittr)
  
  if (sum(is.na(Q))>0){
    stop(paste0(sum(is.na(Q)), " missing data points. You need to gap-fill or remove it."))
  }
  
  ## calculate interval width (2N*)
  N <- area_mi2^0.2
  int_width <- 2*floor((2*N)/2)+1  # nearest odd integer to 2N
  if (int_width<3)  int_width <- 3
  if (int_width>11) int_width <- 11
  
  ## calculation depends on method
  if (method=="fixed"){
    
    ## fixed interval of width 2Nstar
    n.ints <- ceiling(length(Q)/int_width)
    ints <- rep(seq(1,n.ints), each=int_width)[1:length(Q)]
    
    # build data frame
    df <- data.frame(int = ints,
                     Q = Q)
    
    # summarize by interval
    df <- df %>% 
      group_by(int) %>% 
      dplyr::summarize(bf = min(Q)) %>% 
      left_join(df, ., by="int")
    
    return(df$bf)
    
  } else if (method=="sliding"){
    
    ## sliding interval of width 2Nstar
    bf <- rollapply(Q, width=int_width, FUN=min, 
                    align="center", partial=T)
    
    return(bf)
    
  } else if (method=="local"){
    
    ## local minimum
    # local minima are points where Q is equal to the sliding interval minimum
    interval_min <- rollapply(Q, width=int_width, FUN=min,
                              align="center", partial=T)
    i_minima <- which(interval_min==Q)
    
    # interpolate between values using na.approx
    bf <- rep(NaN, length(Q))
    bf[i_minima] <- Q[i_minima]
    if (min(i_minima) != 1) bf[1] <- Q[1]*0.5
    bf <- as.numeric(zoo::na.approx(bf, na.rm=F))
    if (max(i_minima) != length(Q)) bf[(max(i_minima)+1):length(Q)] <- bf[max(i_minima)]
    
    # find any bf>Q and set to Q
    i_tooHigh <- which(bf>Q)
    bf[i_tooHigh] <- Q[i_tooHigh]
    
    return(bf)
    
  } else {
    
    # error
    stop("Wrong or missing method. Choose fixed, sliding, or local")
    
  }
  
}


UKIH <- function(Q, endrule="NA"){
  # R implementation of UKIH baseflow separation algorithm as
  # described in Piggott et al. (2005)
  #
  # Inputs:
  #   Q = discharge timeseries (no missing data) (any units are OK)
  #   endrule = what to do with endpoints, which will always have NAs?
  #     "NA" (default) = retain NAs
  #     "Q" = use Q of the first/last point
  #     "B" = use bf of the first/last point
  #
  # Output:
  #   bf = baseflow timeseries, same length and units as Q
  
  ## package dependencies
  require(zoo)
  require(dplyr)
  require(magrittr)
  
  ## fixed interval of width 5
  int_width <- 5
  n.ints <- ceiling(length(Q)/int_width)
  ints <- rep(seq(1,n.ints), each=int_width)[1:length(Q)]
  
  # build data frame
  df <- data.frame(int = ints,
                   day = seq(1,length(ints)),
                   Q = Q)
  
  # summarize by interval
  df <- 
    df %>% 
    group_by(int) %>% 
    summarize(Qmin = min(Q),
              n_int = sum(is.finite(int))) %>% 
    left_join(df, ., by="int") %>% 
    subset(n_int==int_width)
  
  # extract minimum Qmin for each interval; 
  # these are candidates to become turning points
  df.mins <- df[df$Q==df$Qmin, ]
  
  # if there are two minima for an interval  
  # (e.g. two days with same Q), choose the earlier one
  df.mins <- df.mins[!duplicated(df.mins$int), ]
  
  ## determine turning points, defined as:
  #    0.9*Qt < min(Qt-1, Qt+1)
  # do this using a weighted rolling min function
  df.mins$iQmin <- rollapply(df.mins$Qmin, width=3, align="center", 
                             fill=NA, FUN=function(z) which.min(z*c(1,0.9,1)))
  df.mins <- subset(df.mins, is.finite(iQmin))  # get rid of first/last point
  TP.day <- df.mins$day[df.mins$iQmin==2]
  TP.Qmin <- df.mins$Qmin[df.mins$iQmin==2]
  
  if (length(TP.day>1)){
    
    # linearly interpolate to length Q
    bf <- rep(NaN, length(Q))
    bf[TP.day] <- TP.Qmin
    bf <- as.numeric(zoo::na.approx(bf, na.rm=F))
    
    # need to fill in NAs?
    if (endrule=="Q"){
      # start
      bf[1:(TP.day[1]-1)] <- Q[1:(TP.day[1]-1)]
      
      # end
      bf[(TP.day[length(TP.day)]+1):length(Q)] <- 
        Q[(TP.day[length(TP.day)]+1):length(Q)]
      
    } else if (endrule=="B") {
      # start
      bf[1:(TP.day[1]-1)] <- bf[TP.day[1]]
      
      # end
      bf[(TP.day[length(TP.day)]+1):length(Q)] <- 
        bf[TP.day[length(TP.day)]]
      
    } else if (endrule != "NA") {
      
      stop("Invalid endrule")
      
    }
    
  } else {
    bf <- rep(0, length(Q))
  }
  
  # find any bf>Q and set to Q
  i_tooHigh <- which(bf>Q)
  bf[i_tooHigh] <- Q[i_tooHigh]
  return(bf)
}


Chapman = function(Q, k){
  # Inputs:
  #   Q = discharge timeseries (no missing data) (any units are OK)
  #   k = recession constant; this can be estimated with the function baseflow_RecessionConstant.
  #       
  # Output:
  #   bf = baseflow timeseries, same length and units as Q
  
  n <- length(Q)
  b <- Q
  bf[1] <- Q[1]
  for(i in 2:n){
    
    bf[i] <- (3*k-1) / (3-k) * bf[i-1] + (1-k) / (3-k) * (Q[i]+Q[i-1])
    
    if(bf[i] > Q[i]) bf[i] = Q[i]
  }
  return(bf)
}


BFlow <- function(Q, beta=0.925, passes=3){
  # R implementation of BFLOW baseflow separation algorithm as described in Arnold & Allen (1999).
  # This is the same as the original digital filter proposed by Lyne & Holick (1979) and tested in Nathan & McMahon (1990).
  #
  # It is called BFLOW because of this website: 
  #   http://www.envsys.co.kr/~swatbflow/USGS_GOOGLE/display_GoogleMap_for_SWAT_BFlow.cgi?state_name=indiana
  #
  # Inputs:
  #   Q = discharge timeseries (no missing data) (any units are OK)
  #   beta = filter parameter; recommended value 0.925 (Nathan & McMahon, 1990); 0.9-0.95 reasonable range
  #   passes = how many times to go through the data (3=default=forward/backward/forward)
  #       
  # Output:
  #   bf = baseflow timeseries, same length and units as Q
  
  ## package dependencies
  require(zoo)
  require(dplyr)
  require(magrittr)
  
  # Q for use in calculations
  bfP <- Q
  
  for (p in 1:passes){
    # figure out start and end
    if ((p %% 2)==0){
      # backward pass
      i.start <- length(Q)-1
      i.end   <- 1
      i.fill  <- length(Q)
      ts      <- -1
    } else {
      # forward pass
      i.start <- 2
      i.end   <- length(Q)
      i.fill  <- 1
      ts      <- 1
    }
    
    # make empty vector
    qf <- rep(NaN, length=length(Q))
    
    # fill in value for timestep that will be ignored by filter
    if (p==1){
      qf[i.fill] <- bfP[1]*0.5
    } else {
      qf[i.fill] <- max(c(0, (Q[i.fill]-bfP[i.fill])))
    }
    
    # go through rest of timeseries
    for (i in i.start:i.end){
      qf[i] <- (beta*qf[i-ts] + ((1+beta)/2)*(bfP[i]-bfP[i-ts]))
      
      # check to make sure not too high/low
      if (qf[i] > bfP[i]) qf[i] <- bfP[i]
      if (qf[i] < 0) qf[i] <- 0
    }
    
    # calculate bf for this pass
    bfP <- bfP-qf
    
    # when p==passes, return bfP
    if (p==passes){
      bf <- bfP
    }
    
  } # end of passes loop
  
  return(bf)
}


Lyne_Hollick = function(Q, k){
  # Inputs:
  #   Q = discharge timeseries (no missing data) (any units are OK)
  #   k = recession constant; this can be estimated with the function baseflow_RecessionConstant.
  #       
  # Output:
  #   bf = baseflow timeseries, same length and units as Q
  
  n <- length(Q)
  qs <- Q
  qs[1] <- Q[1]*0.2
  for(i in 2:n){
    qs[i] <- k * qs[i-1] + (1+k) / 2 * (Q[i] - Q[i-1])
    
    if(qs[i] > Q[i]) qs[i] <- Q[i]
    if(qs[i] < 0) qs[i] <- 0
  }
  bf <- Q - qs
  return(bf)
}


Chapman_Maxwell = function(Q, k){ 
    # Inputs:
    #   Q = discharge timeseries (no missing data) (any units are OK)
    #   k = recession constant; this can be estimated with the function baseflow_RecessionConstant.
    #       
    # Output:
    #   bf = baseflow timeseries, same length and units as Q
  
  n <- length(Q)
  bf <- Q
  bf[1] <- Q[1]
  for(i in 2:n){
    if(bf[i-1] < Q[i]) {
      bf[i] <- (k) / (2-k) * bf[i-1] + (1-k) / (2-k) * (Q[i])
    } else if( bf[i-1] > Q[i]){
      bf[i] = Q[i]
    }
  }
  return(bf)
}


Eckhardt <- function(Q, BFImax, k){
  # R implementation of Eckhardt (2005) baseflow separation algorithm.  
  #
  # Inputs:
  #   Q = discharge timeseries (no missing data) (any units are OK)
  #   BFImax = maximum allowed value of baseflow index; 
  #            the BFImax derived by the CMB method 
  #   k = recession constant; this can be estimated with the function baseflow_RecessionConstant.
  #       
  # Output:
  #   bf = baseflow timeseries, same length and units as Q
  
  # empty output vector
  bf <- rep(NaN, length(Q))
  
  # fill in initial value
  bf[1] <- Q[1]  # from Barlow 'Digital Filters' document
  
  # scroll through remaining values
  for (i in 2:length(Q)){
    # calculate bf using digital filter
    bf[i] <- (((1-BFImax)*k*bf[i-1]) + ((1-k)*BFImax*Q[i]))/(1-k*BFImax)
    
    # make sure 0 <= bf <= Q
    if (bf[i]<0)    bf[i] <- Q[i]*BFImax*0.9  # from Barlow 'Digital Filters' document
    if (bf[i]>Q[i]) bf[i] <- Q[i]
  }
  return(bf)
}


SC_baseflow <- function(df) {
  df[,1] <- as.Date(df[,1], tryFormats = "%m/%d/%Y" )
  
  df$Year <- lubridate::year(df$Date)
  df$SC <- as.numeric(df$C, na.rm = TRUE)
  
  df_year <- df %>%
    group_by(Year) %>%
    filter(n() >= 365) %>% 
    summarise(Percentile99 = quantile(SC, probs = 0.99, na.rm = TRUE))
  SCBF <- mean(df_year$Percentile99, na.rm = TRUE)
  return(SCBF)
}
CMB <- function(date, discharge, SC, SC_baseflow, 
                STAID="Unknown") {
  
  ## Start of code: initial processing
  STAID <- as.character(STAID[1L])
  discharge <- pmax(discharge, 0.000001 ) # Convert 0 to a small number
  if(any(is.na(discharge)))
    stop("Missing values in discharge vector.")
  
  n <- length( discharge )
  
  baseflow  <- numeric( n ) 
  quickflow <- numeric( n )
  
  SC_runoff <- quantile( SC, 0.01, na.rm = TRUE)  
  
  for(i in 1:n){
    
    baseflow[i] <- discharge[i] * ( SC[i] - SC_runoff ) / ( SC_baseflow - SC_runoff )
  }
  
  baseflow <- pmin( baseflow, discharge )  
  baseflow <- pmax( baseflow, 0.0 )
  baseflow = round( baseflow, 3)  
  
  quickflow <- discharge - baseflow
  
  retval <- data.frame( date=date, discharge=discharge, baseflow=round( baseflow, 3),
                        quickflow = round( quickflow, 3 ) )
  
  if(!is.null(STAID))
    attr(retval, "STAID") <- STAID  
  attr(retval, "type") <- "part"
  class(retval) <- c("CMB", "data.frame")  
  return(retval)
}


# Recession Constant
baseflow_RecessionConstant <- function(Q, UB_prc=0.95, method="Brutsaert", K=50){
  # Script to estimate baseflow recession constant.
  #
  # Inputs:
  #   Q = discharge timeseries (no missing data) (any units are OK)
  #   UB_prc = percentile to use for upper bound of regression
  #   method = method to use to calculate recession constant
  #     "Brutsaert" = Brutsaert (2008) WRR
  #   K = characteristic timescale of the catchment drainage process, 
  #       also commonly referred to as the storage coefficient;
  #       K = 45 ± 15 days is from Brutsaert (2008) WRR
  #       
  # Output:
  #   k = recession constant
  
  ## package dependencies
  require(quantreg)  # used for quantile regression

  # calculate lagged difference (dQ/dt) based on before/after point
  dQ_dt <- c(NaN, diff(Q, lag=2)/2, NaN)
  dQ_dt_left <- c(NaN, diff(Q))
  
  # screen data for which dQ_dt to calculate recession, based on rules in Brutsaert (2008) WRR Section 3.2
  which_negative <- which(dQ_dt < 0 & dQ_dt_left < 0 & Q > 0)
  which_positive <- which(dQ_dt >= 0)
  which_positive_with_buffer <- unique(c(which_positive-2, which_positive-1, which_positive,
                                         which_positive+1, which_positive+2, which_positive+3))  # 2 days before and 3 days after a positive or 0 value
  which_positive_with_buffer <- which_positive_with_buffer[which_positive_with_buffer > 0]  # get rid of negative indices; possible because of 2 days before
  which_keep <- which_negative[!(which_negative %in% which_positive_with_buffer)]  # get rid of points within buffer around flow increases
  which_keep <- which_keep[(which_keep-1) %in% which_keep]  # trim to dates with both the current and previous day retained
  
  # any data exist to fit?
  if (length(which_keep) >= K){
    
    # fit regression
    fit.qr <- rq(Q[which_keep] ~ 0+Q[which_keep-1], tau=UB_prc)  # force intercept to go through origin
    
    # extract constant
    k <- as.numeric(coef(fit.qr)[1])
    
  } else {
    k <- NaN
  }
  return(k)
    
}


























