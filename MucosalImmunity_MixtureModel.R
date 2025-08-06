# Train mixture model to identify reinfection based on serial serology

# Load packages and source necessary files
# library(fitdistrplus)
library(tidyverse)

library(zoo)

theme_set(theme_classic())


# I: LOAD & FORMAT DATA--------------------------------------------------------

# Load data from training set
d0 <- readRDS("MucosalImmunity_PCR_training.rds")

# Create log (nodratio) variable
d0$logratio <- log(d0$nodratio)

# II: MIXTURE MODEL FUNCTIONS --------------------------

expit <- function(x){1/(1+exp(-x))}
logit <- function(x){log(x/(1-x))}

#### FUNCTIONS #####
# Function to return log likelihood of distribution and parameters given observed data
# param: vector of coefficients and distribution parameter names
# betas: vector of coefficient names for covariates
# gammas: vector of parameter names for generalized gamma distribution
# val: observed value
# reg_vars_p: value of coefficients for logistic regression for probability of being a responder
# reg_vars_g: value of coeffictients for gamma regression for probability of OD given responder status
# model_nonresp: distribution of values for nonresponders (no interval exposure)
# model_resp: distribution of values for responders (interval increase in titers reflective of exposure)
fitmixture_reg_gengamma<-function(param, val, reg_vars_p, reg_vars_g, model_nonresp, model_resp) {
  
  betas <- param[grepl("beta",names(param))] # grepl("beta", names(param)) gives logical vector with name matching "beta" in the vector "param"
  gammas <- param[grepl("gamma",names(param))] # same concept here, matching "gamma"
  
  # Mean ISR of seroresponders restricted to be greater than mean ISR of non-responders
  # So, set mean to mean_nonresp+exp(-gamma0)
  gammas[1] = param[["mean_nonresp"]]+exp(gammas[1]) # guessing from this that param is expressed in terms of log(E(x))
  
  if (length(betas)>1) { # If there are coefficients other than beta0
    p = expit(colSums(betas*reg_vars_p)) # constrain to be 0<.<1; p is prob of being a responder; multiplying beta x covar, summing up to get y predicted
  } else {
    p = expit(betas*reg_vars_p) # constrain to be 0<.<1
  }
  
  if (length(gammas)>1) {
    mu <- colSums(gammas*reg_vars_g) # mu is the mean; model for response
  } else {
    mu <- gammas*reg_vars_g #multiplying gamma x covar, summing up to get y predicted for mean response
  }
  
  if (model_resp=="gamma") {
    gamma_shape = exp(param[["sigma"]]) # in gamma dist E(x) = k*theta = alpha/beta; here defining alpha = exp(sigma)
    gamma_rate = exp(param[["sigma"]])/mu # this works out to beta = alpha/k*theta = exp(sigma)/mu
  } else if (model_resp=="lognormal") {
    logmean = log(mu) # mean of the log function is log(mu)
    logsd = exp(-param[["sigma"]]/2) # sd is 1/sqrt(alpha), such that variance = 1/alpha
  } else if (model_resp=="gengamma") {
    mu_gengamma = log(mu)
    sigma_gg=exp(-param[["sigma"]]/2)
    
    Q=param[["delta"]] # third parameter for generalized gamma
    if (abs(Q)<1e-4){Q=0}  # set to zero if sufficiently small
    
  } else if (model_resp=="normal") {
    sd_norm = exp(param[["sigma"]]) # for normal distribution sd = alpha; unstated is that mean = mu
  }
  
  # Seroresponse
  # Given distribution with the given parameters, probability of observed value
  if (model_resp=="gamma"){
    seroresponse_l <- dgamma(val,shape=gamma_shape,rate=gamma_rate)
  } else if (model_resp=="lognormal"){
    seroresponse_l <- dlnorm(val,meanlog=logmean,sdlog=logsd)
  } else if (model_resp=="gengamma"){
    seroresponse_l <- dgengamma(val,mu=mu_gengamma,sigma=sigma_gg,Q=Q)
  } else if (model_resp=="normal") {
    seroresponse_l <- dnorm(val,mean=mu,sd=sd_norm)
  }
  
  if (model_nonresp=="gamma") { # fewer options: just gamma or normal (no lognormal or gengamma)
    gamma_shape = exp(param[["b0"]])
    gamma_rate = exp(param[["b0"]])/param[["mean_nonresp"]] # same as above, sigma expressed as b0, mu expressed as mean_nonresp
    seroresponse_non <- dgamma(val,shape=gamma_shape,rate=gamma_rate)
  } else if (model_nonresp=="normal") {
    seroresponse_non <- dnorm(val, mean = param[["mean_nonresp"]], sd = exp(param[["b0"]]))
  }
  
  return(
    # neg log likelihood function # output is the likelihood of parameters given data and choice of distribution
    - sum(log(
      p * seroresponse_l +
        (1-p) * seroresponse_non
    ))
  )
}

# Function to optimize log likelihood
mixcens_optim <- function(data,y_name,cov_names_b,cov_names_g,mi,init_params,model_nonresp,model_resp) {
  
  # Data processing
  
  # Univariate and interaction variables
  cov_names_b_uni <- cov_names_b[!grepl(':',cov_names_b)] # subset of cov_names_b that does not match ":"
  cov_names_b_int <- setdiff(cov_names_b,cov_names_b_uni) # interaction terms
  
  cov_names_g_uni <- cov_names_g[!grepl(':',cov_names_g)]
  cov_names_g_int <- setdiff(cov_names_g,cov_names_g_uni)
  
  for (interaction_var in union(cov_names_b_int,cov_names_g_int)) { # values that have ":" for either cov_names_b or g (all interaction terms)
    
    var1 <- unlist(strsplit(interaction_var,':'))[1] # separate the interaction factors
    var2 <- unlist(strsplit(interaction_var,':'))[2]
    
    data[,interaction_var] <- data[,var1]*data[,var2] #multiply the interaction factors together -> single value
    
  }
  
  # Outcome variable
  y = data[,c(y_name)]
  
  # Regression covariates. Get into matrix and remove all with missing data
  numbetas = length(cov_names_b) # Number of covariates
  x_b = data[,c(cov_names_b)] # names of covariates
  
  if (numbetas>1) {x_b=t(x_b)} # transpose if in matrix form (i.e., if at least one covariate)
  
  numgammas = length(cov_names_g)
  x_g = data[,c(cov_names_g)]
  
  if (numgammas>1) {x_g=t(x_g)}
  
  # If no covariates, the outcome variable does not require transformation
  # value of covariates set to 1
  if (numbetas+numgammas==0) {
    
    y = as.numeric(y)
    y = y[!is.na(y)]
    reg_vars_b = reg_vars_g = rep(1,length(y))
    
    
  } else {
    
    covariates = data.matrix(rbind(y,x_b,x_g))
    covariates = covariates[,!is.na(apply(covariates,2,function(x) min(x)))]
    
    y = covariates[1,]
    covariates = covariates[-1,,drop=FALSE]
    
    if (numgammas==0) {
      reg_vars_b = rbind(rep(1,length(y)),
                         covariates)
      covariates <- NULL
    } else {
      reg_vars_b = rbind(rep(1,length(y)),
                         covariates[0:numbetas,])
      if (numbetas>0) {
        covariates = covariates[-(0:numbetas),,drop=FALSE] 
      } else {
        
      }
    }
    
    reg_vars_g = rbind(rep(1,length(y)),
                       covariates)
    
  }
  
  param <- init_params
  names = names(param) = c("b0","mean_nonresp","sigma","beta0","gamma0","delta")
  
  if (model_resp != "gengamma") {
    param = param[names(param) != "delta"]
    names = names[names != "delta"]
  }
  
  if (numbetas>0) {
    param <- c(param,rep(0,numbetas))
    names <- c(names,paste0("beta: ",cov_names_b))
    
  }
  
  if (numgammas>0) {
    param <- c(param,rep(0,numgammas))
    names <- c(names,paste0("gamma: ",cov_names_g))
  }
  
  names(param) <- names
  
  result=optim(param,fitmixture_reg_gengamma,val=y,reg_vars_p=reg_vars_b,reg_vars_g=reg_vars_g, # we are optimizing the likelihood function
               model_nonresp=model_nonresp,model_resp=model_resp,hessian=T,control=list(maxit=mi)
  )
  parfit=result$par #optimized parameters
  
  if(min(diag(result$hessian))>0) {
    ses <- sqrt(diag(solve(result$hessian)))
  } else {
    ses <- rep(NA,length(param))
  }
  cilowers = parfit-1.96*ses
  ciuppers = parfit+1.96*ses
  aic = 2*length(parfit)+2*result$value
  
  gammas = parfit[grepl('gamma',names(parfit))]
  gammas[1] = parfit[['mean_nonresp']]+exp(gammas[1])
  if (numgammas>0) {
    mu <- colSums(gammas*reg_vars_g)
  } else {
    mu <- gammas*reg_vars_g
  }
  
  if (model_resp=="gamma") {
    gamma_shape = exp(parfit[["sigma"]])
    gamma_rate = exp(parfit[["sigma"]])/mu
    y.signl1 = expit(parfit[["beta0"]])* # prob of being responder (y.signl is prob of seeing this value)
      dgamma(y,shape=gamma_shape,rate=gamma_rate) # Prob of seeing this value given distribution and parameters and responder
  } else if (model_resp=="lognormal") {
    logmean = log(mu)
    logsd = exp(-param[["sigma"]]/2)
    y.signl1 = expit(parfit[["beta0"]])*
      dlnorm(y,meanlog=logmean,sdlog=logsd)
  } else if (model_resp=="gengamma") {
    mu_gengamma = log(mu)
    sigma_gg=exp(-parfit[["sigma"]]/2)
    
    Q=param[["delta"]]
    if (abs(Q)<1e-4){Q=0}
    y.signl1 <- expit(parfit[["beta0"]])*dgengamma(y,mu=mu_gengamma,sigma=sigma_gg,Q=Q)
  } else if (model_resp=="normal") {
    sd_norm=exp(parfit[["sigma"]])
    y.signl1 <- expit(parfit[["beta0"]])*dnorm(y,mean=mu,sd=sd_norm)
  }
  
  if (model_nonresp=="gamma") {
    gamma_shape = exp(parfit[["b0"]])
    gamma_rate = exp(parfit[["b0"]])/parfit[["mean_nonresp"]]
    y.error = (1-expit(parfit[["beta0"]]))* # prob of being nonresponder
      dgamma(y,shape=gamma_shape,rate=gamma_rate)
  } else if (model_nonresp=="normal") {
    mean = parfit[["mean_nonresp"]]
    sd = exp(parfit[["b0"]])
    y.error = (1-expit(parfit[["beta0"]]))*dnorm(y, mean = parfit[["mean_nonresp"]], sd = exp(parfit[["b0"]]))
  }
  
  p.signal = y.signl1/(y.error+y.signl1)
  
  return(list(r=result,par=parfit,cilowers=cilowers,ciuppers=ciuppers,aic=aic,p.signal=as.numeric(p.signal))) # fit
  
}

# Function to choose next set of proposed parameters by disturbing current set according to tuning parameter
proposal <- function(cur_param,tuning_sd) { # Choose next parameters based on current ones by disturbing with tuning parameter
  # Try doubling the tuning parameter for beta parameters
  tuning_sds=rep(tuning_sd,length(cur_param))
  return(cur_param+rnorm(length(cur_param),0,tuning_sds)) # n, mean, sd
}
# log likelihood of prior parameters for normal distribution
norm_logprior <- function(param,prior_means,prior_sd) {
  sum(log(dnorm(param,mean=prior_means,sd=prior_sd)))
}

# log likelihood of prior parameters for uniform distribution
unif_logprior <- function(param,param.lower,param.upper){
  inlower = (param >= param.lower)
  inupper = (param <= param.upper)
  inlower = ifelse(is.na(inlower),1,inlower)
  inupper = ifelse(is.na(inupper),1,inupper)
  sum(log(inlower & inupper))
}

# Print to file
mycat = function(x,outcon){
  cat(
    paste0(paste(x, collapse=","),'\n')
    ,file = outcon
  )
}

# Function to run Metropolis-Hastings algorithm; parameters similar to optim function
run_mh_mixcens <- function(data,y_name,cov_names_b,cov_names_g,model_nonresp,model_resp,param_lowers,param_uppers,
                           tuning_sd.init,chain_len,mcmc_burnin,mcmc_subsample,
                           baseFile,negprior=F,prior_params=c()) {
  
  if (model_resp!="gengamma" & length(param_lowers) == 6) {stop("Length of uniform prior boundaries doesn't match assumed distribution")}
  
  # Data processing
  # Univariate and interaction variables
  cov_names_b_uni <- cov_names_b[!grepl(':',cov_names_b)]
  cov_names_b_int <- setdiff(cov_names_b,cov_names_b_uni)
  
  cov_names_g_uni <- cov_names_g[!grepl(':',cov_names_g)]
  cov_names_g_int <- setdiff(cov_names_g,cov_names_g_uni)
  
  for (interaction_var in union(cov_names_b_int,cov_names_g_int)) {
    
    var1 <- unlist(strsplit(interaction_var,':'))[1]
    var2 <- unlist(strsplit(interaction_var,':'))[2]
    
    data[,interaction_var] <- data[,var1]*data[,var2]
    
  }
  
  # Outcome variable
  y = data[,c(y_name)]
  
  # Regression covariates. Get into matrix and remove all with missing data
  numbetas = length(cov_names_b)
  x_b = data[,c(cov_names_b)]
  
  if (numbetas>1) {x_b=t(x_b)}
  
  numgammas = length(cov_names_g)
  x_g = data[,c(cov_names_g)]
  
  if (numgammas>1) {x_g=t(x_g)}
  
  if (numbetas+numgammas==0) {
    
    y = as.numeric(y)
    y = y[!is.na(y)]
    reg_vars_b = reg_vars_g = rep(1,length(y))
    
    
  } else {
    
    covariates = data.matrix(rbind(y,x_b,x_g))
    covariates = covariates[,!is.na(apply(covariates,2,function(x) min(x)))] # removing missing data
    
    y = covariates[1,]
    covariates = covariates[-1,,drop=FALSE]
    
    if (numgammas==0) {
      reg_vars_b = rbind(rep(1,length(y)),
                         covariates)
      covariates <- NULL
    } else {
      reg_vars_b = rbind(rep(1,length(y)),
                         covariates[0:numbetas,])
      if (numbetas>0) {
        covariates = covariates[-(0:numbetas),,drop=FALSE] 
      } else {
        
      }
    }
    
    reg_vars_g = rbind(rep(1,length(y)),
                       covariates)
    
  }
  
  currentLL = NaN
  # Initialize MCMC chain
  while (is.nan(currentLL) | is.infinite(currentLL)) { # while current log likelihood is infinity or 1/infinity
    
    if (negprior) {
      prior_means = c(prior_params[['b0']],prior_params[['mean_nonresp']]) # set priors for nonresponders
      param = rnorm(2,mean=prior_means,sd=prior_params[['prior_sd']]) # pick 2 random numbers w. mean and sd reflecting data
    } else {
      param = runif(2,min=param_lowers[1:2],max=param_uppers[1:2]) # pick from uniform dist
    }
    
    names = names(param) = c("b0","mean_nonresp")
    
    param <- c(param,runif(3,min=param_lowers[3:5],max=param_uppers[3:5])) # randomly pick sigma, beta, gamma from uniform distributions with diff bounds for each
    names = c(names,"sigma","beta0","gamma0")
    names(param)=names
    
    if (model_resp == "gengamma") {
      param = c(param,runif(1,min=param_lowers[6],max=param_uppers[6])) # randomly pick delta from uniform
      names = c(names,"delta")
      names = names(param)
    }
    
    
    if (numbetas>0) {
      param <- c(param,runif(numbetas,min=-1,max=1))
      names <- c(names,paste0("beta: ",cov_names_b))
      
    }
    
    if (numgammas>0) {
      param <- c(param,runif(numgammas,min=-1,max=1))
      names <- c(names,paste0("gamma: ",cov_names_g))
    }
    
    names(param) <- names
    
    cur_params <- param
    
    if (negprior) {
      priorLL = norm_logprior(cur_params[1:2],prior_means,prior_params[['prior_sd']]) # log likelihood for prior
    } else {
      priorLL = unif_logprior(cur_params[1:2],param_lowers[1:2],param_uppers[1:2])
    }
    # log of product of likelihoods for the set of parameters and the distribution of observed response values
    currentLL <- -fitmixture_reg_gengamma(cur_params, y, reg_vars_b, reg_vars_g, model_nonresp,model_resp) + 
      priorLL + #parameters b0 and mean_nonresp
      unif_logprior(cur_params[3:length(cur_params)], # parameters sigma, all beta and gamma coefficients
                    c(param_lowers[3:length(param_lowers)],rep(-1,numbetas),rep(-1,numgammas)),
                    c(param_uppers[3:length(param_uppers)],rep(1,numbetas),rep(1,numgammas)))
    
  }
  
  # Files to store results
  out = list()
  out$stepsize = paste0(baseFile,'_stepsize.txt')
  out$chain = paste0(baseFile,'_chain.csv')
  out$chainMeta = paste0(baseFile,'_chainMeta.csv')
  
  # open connections to output files
  out = lapply(out, function(x){ file(x,'w') })
  cat(
    paste("Initial step size:",tuning_sd.init)
    ,paste('Iteration','AcceptRate','NewTuningSD')
    ,file = out$stepsize
    ,sep = "\n"
  )
  
  i.stepsizeEval = 50000
  acceptRate.ideal = .2
  
  step_count<-0
  tuning_sd = tuning_sd.init
  sumAccept = 0
  
  for (i in 1:chain_len) {
    
    if((i %% i.stepsizeEval)==0){ # this will equal 0 when i = a multiple of 50k
      acceptRate = sumAccept/i.stepsizeEval # rate accepted in last 50k steps
      tuning_sd = tuning_sd * sqrt(1 + ((acceptRate-acceptRate.ideal)/acceptRate.ideal))
      mycat( c(i,acceptRate,tuning_sd) ,out$stepsize)
      sumAccept = 1 #reset summAccept to 1
    }
    
    # Propose new params
    if (i==1) {
      prop_params <- param #initial parms we set
    } else {
      prop_params <- proposal(cur_params,tuning_sd) # propose new params
      
    }
    
    # negll_mixture_censoring returns negative log likelihood
    if (negprior) {
      priorLL = norm_logprior(prop_params[1:2],prior_means,prior_params[['prior_sd']])
    } else {
      priorLL = unif_logprior(prop_params[1:2],param_lowers[1:2],param_uppers[1:2])
    }
    proposalLL <- -fitmixture_reg_gengamma(prop_params, y, reg_vars_b, reg_vars_g, model_nonresp,model_resp) + 
      priorLL +
      unif_logprior(prop_params[3:length(prop_params)],
                    c(param_lowers[3:length(param_lowers)],rep(-1,numbetas),rep(-1,numgammas)),
                    c(param_uppers[3:length(param_uppers)],rep(1,numbetas),rep(1,numgammas)))
    
    log_acc <-  proposalLL-currentLL # ratio of proposed vs prior
    
    if (!is.nan(log_acc)) {
      
      if (log(runif(1))<log_acc) {
        
        cur_params <- prop_params
        currentLL <- proposalLL
        sumAccept=sumAccept+1  # accept of meets certain ratio
        
      }
      
    }
    
    if (i>mcmc_burnin) {
      
      step_count <- step_count+1 # only count steps once we are out of the burn-in period
      
      if (step_count==mcmc_subsample) {
        
        mycat( c(i,currentLL) ,out$chainMeta)
        mycat( cur_params ,out$chain)
        
        step_count<-0
        
      }
      
    }
    
  } # End of for loop
  
} # end of mcmc function


# III: MCMC PARAMETERS-----------------------

tuning_sd = 0.1 
chain_len = 250000 # Chain length
mcmc_burnin = 125000 # Burn-in
mcmc_subsample=100 # Subsample

# Prior parameters for non-response
prior_sd = 0.1 # Variance of prior distribution - if low, seronegative distribution will stay close to prior

# Set parameters
model_resp = "lognormal" # Distribution of seroresponse
model_nonresp = "normal" # Distribution of non-response
y_name = "logratio" # Outcome variable to fit to 
suffix = '_noprior' # Prior: no prior, use negatives as prior, or a manual prior based on eyeballing negatives 
prior_lowers = c(-2,-2,-4,-4,-2) # b0, mean_nonresp, sigma, beta0, gamma0, delta
prior_uppers = c(2,2,4,4,4)


baseFile=paste0('mcmc','_',y_name,'_nr.',model_nonresp,'_r.',model_resp,suffix)

# Run MCMC
run_mh_mixcens(d0, y_name, cov_names_b = NULL, cov_names_g = NULL,
               model_nonresp, model_resp, prior_lowers, prior_uppers,
               tuning_sd, chain_len, mcmc_burnin, mcmc_subsample,
               paste0(baseFile), negprior=F, prior_params=c())
closeAllConnections()

