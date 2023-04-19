################################################################################
# PGW Survival Function
################################################################################
#' PGW Survival Function
#' @param t: positive parameter
#' @param eta: scale parameter
#' @param  nu: shape parameter
#' @param delta: shape parameter
#' @param  log.p: return log survival (TRUE or FALSE)
#' @return PGW survival function
#' @export
spgw <- function(t, eta, nu, delta, log.p = FALSE){
  val <- 1 - ( 1 + (t/eta)^nu )^(1/delta)
  if(log.p) return(val) else return(exp(val))
}

################################################################################
# PGW Random Number Generation Function
################################################################################
#' Random Number Generation Function PGW distribution
#' @param n: sample size
#' @param eta: scale parameter
#' @param  nu: shape parameter
#' @param delta: shape parameter
#' @return vector of simulated values
#' @export
rpgw <- function(n, eta, nu, delta){
  p <- runif(n)
  out <- eta*(  ( 1 - log(1-p) )^delta - 1 )^(1/nu)
  return(as.vector(out))
}


################################################################################
#' Function to simulate survival times from a Kaplan-Meier estimator
#' It simulates survival times either from the original Kaplan-Meier estimator
#' or it fits a PGW distribution using least squares based on nlsp points
################################################################################
#' @param KM.obj : Kaplan-Meier object (a surv object)
#' @param times : times used to construct the Kaplan-Meier estimator
#' @param survprob : survival probabilities for the Kaplan-Meier estimator
#' @param n : number of observations
#' @param seed : seed used in the simulation
#' @param type.sim : "KM" for using the Kaplan-Meier estimator directly or "PGW" for using the PGW distribution fitted via least squares using nlsp points
#' @param nlsp : number of points to be used in the least squares fit (equally spaced between min and max)
#' @export 
KMSim <- function(KM.obj = NULL, times = NULL, survprob = NULL, n = 100, 
                  seed = 1234, type.sim = NULL, nlsp = NULL){
  
  # Fix the seed for simulation
  set.seed(seed)
  # Simulating n uniform observations
  u = runif(n)
  
  out <- vector()
  
  # Using the KM object  
  if(type.sim == "KM" & !is.null(KM.obj)){
    out <- as.vector(quantile(KM.obj, u) )$quantile
    names(out) <- 1:n
  }  
  
  # Using the KM based on times and survival probabilities  
  if(type.sim == "KM" & is.null(KM.obj)){
    minsurv <- min(survprob)
    maxsurv <- max(survprob)
    mintimes <- min(times)
    maxtimes <- max(times)
    
    # Kaplan-Meier estimator based on a step function
    stepf <- function(x) stepfun(x = times, y = c(1,survprob), right = FALSE)(x)
    for(i in 1:n){
      if(u[i] > maxsurv) out[i] <- mintimes
      if(u[i] < minsurv) out[i] <- Inf
      if(u[i] > minsurv & u[i] < maxsurv){
        # Inverting the step function
        tempf <- Vectorize(function(t) stepf(t) - u[i]  )
        
        out[i] <- uniroot(tempf, c(mintimes, maxtimes), maxiter = 10000)$root
      }
      
    }    
  }   
  
  # Using PGW distribution and KM object
  if(type.sim == "PGW" & !is.null(KM.obj)){
    # Required quantities
    survprob <- as.vector(KM.obj$surv)
    times = as.vector(KM.obj$time)
    index <- round(seq(from = 1, to = length(times), length.out = nlsp))
    vec <- times[index]
    svec <-   survprob[index]
    
    # Sum of squared errors
    lss <- function(par){
      pars = exp(par)
      sfit <- Vectorize(function(t) spgw(t, pars[1], pars[2], pars[3]))
      ssd <- sum( ( svec - sfit(vec))^2  )
      return(ssd)
    }
    
    # Optimisation step     
    OPT <- nlminb(c(0,0,0), lss, control = list(iter.max = 10000))
    LSE <- exp(OPT$par) 
    
    # Simulation step
    sim <- rpgw(n = n, eta = LSE[1], nu = LSE[2], delta = LSE[3])
    
    out <- list(LSE = LSE, sim = sim)
    
  }  
  
  # Using PGW distribution, times, and survival probabilities
  if(type.sim == "PGW" & is.null(KM.obj)){
    # Required quantities
    index <- round(seq(from = 1, to = length(times), length.out = nlsp))
    vec <- times[index]
    svec <- survprob[index]
    
    # Sum of squared errors  
    lss <- function(par){
      pars = exp(par)
      sfit <- Vectorize(function(t) spgw(t, pars[1], pars[2], pars[3]))
      ssd <- sum( ( svec - sfit(vec))^2  )
      return(ssd)
    }
    
    # Optimisation step
    OPT <- nlminb(c(0,0,0), lss, control = list(iter.max = 10000))
    LSE <- exp(OPT$par) 
    
    # Simulation step
    sim <- rpgw(n = n, eta = LSE[1], nu = LSE[2], delta = LSE[3])
    
    out <- list(LSE = LSE, sim = sim)
  }    
  
  return(out)
  
}