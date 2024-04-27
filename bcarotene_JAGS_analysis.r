################################################################
##                                                            ##
##  bcarotene_JAGS_analysis.r                                 ##
##                                                            ##
##    Demonstration file for beta-carotene analysis project   ##
##    Dr. Fletcher G.W. Christensen                           ##
##    3 May 2022                                              ##
##                                                            ##
################################################################

library(rjags)
library(coda)

bc_datafile <- "c:\\[Data]\\School\\UNM\\Teaching\\577\\bcarotene.csv"

setwd("d:\\[Data]\\School\\UNM\\Teaching\\577\\")


##  Set MCMC parameters.
chains <- 1           ##  One chain of MCMC values is obtained
burn_in <- 1000       ##  The first 1000 iterates from the MCMC chain are discarded
iterations <- 100000   ##  10,000 values are sampled from the posterior distribution
thin <- 1             ##  No thinning is performed


##  Import beta-carotene data an store it in new variable names for use with JAGS.
bc_data <- read.csv(bc_datafile)

ptid <- bc_data$ptid
month <- bc_data$month
bcarot <- bc_data$bcarot
vite <- bc_data$vite
dose <- bc_data$dose
age <- bc_data$age
male <- bc_data$male
bmi <- bc_data$bmi
chol <- bc_data$chol
cauc <- bc_data$cauc
vauc <- bc_data$vauc


##  Define new variables for JAGS to use.
n <- dim(bc_data)[1]
n_patients <- length( unique(ptid) )

intercept <- rep(1,n)
tx <- intercept - (month<4)


##  Define the design matrix for the model and specify the number of covariates, for use in JAGS.
x <- cbind(intercept,tx*dose,age)
n_covariates <- dim(x)[2]
tau_b <- 0.00001


##  Establish data and parameter lists for use with OpenBUGS; define function for generating initial parameter values.
data <- list( "bcarot"=bcarot, "n"=n, "n_patients"=n_patients, "n_covariates"=n_covariates,
              "ptid"=ptid, "x"=x, "tau_b"=tau_b )
inits <- function() {
  list( beta = rnorm( n_covariates, 0, 1 ),
        gamma = rnorm( n_patients, 0, 1 ),
        tau_bc = runif( 1, 0, 2 ),
        tau_g = runif( 1, 0, 2 ) )
}
parameters <- c( "beta", "gamma", "sigma_bc", "sigma_g" )



###############   CREATE MODEL   ###############

bc_modelstring<-"model
  {
  	for(i in 1:n){
	  	bcarot[i] ~ dnorm( mu[i], tau_bc )
		  mu[i] <- inprod( x[i,], beta[] ) + gamma[ ptid[i] ]
  	}
	  for(j in 1:n_patients){
  		gamma[j] ~ dnorm( 0, tau_g )
	  }
  	for(k in 1:n_covariates){
		  beta[k] ~ dnorm( 0, tau_b )
	  }
  	tau_bc ~ dgamma(0.001, 0.001)
	  tau_g ~ dgamma(0.001, 0.001)
  	sigma_bc <- pow(tau_bc,-0.5)
	  sigma_g <- pow(tau_g,-0.5)
  }"

beta_carotene.m <- jags.model( data=data, inits=inits, file=textConnection(bc_modelstring),
                               n.chains=chains )

beta_carotene.sim <- coda.samples(beta_carotene.m, parameters, n.iter=iterations,
                                  thin=thin, n.burn=burn_in)

## Assign variable names to posterior samples.
beta_carotene.iterates <- as.matrix(beta_carotene.sim)
beta_carotene.iterates <- as.data.frame(beta_carotene.iterates)

posterior_betas <- beta_carotene.iterates[,(1:n_covariates)]
posterior_gammas <- beta_carotene.iterates[,((n_covariates+1):(n_covariates+n_patients))]
posterior_sigmas <- beta_carotene.iterates[,((n_covariates+n_patients+1):(n_covariates+n_patients+2))]
