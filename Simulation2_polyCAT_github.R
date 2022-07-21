###############################################
### Simulation 2 script - Niek Frans; 2022  ###
### Empirical Priors in Polytomous CATs:    ###
### Risks and Rewards in Clinical Settings  ###
###############################################

# Set working directory to file path
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path)))

# Override run_local function to set prior_sd at a later stage
source("run_local2.R")
environment(run_local) <- asNamespace('mirtCAT')
assignInNamespace("run_local", run_local, ns = "mirtCAT",pos = "package:mirtCAT")

###############################################
### LOADING FUNCTIONS AND RELEVANT PACKAGES ###
###############################################

# Setting seed
set.seed(269513)

# Loading relevant packages
pkg <- c("mirt","mirtCAT")
sapply(pkg,library,character.only=TRUE,logical.return = TRUE)

# Function to set a prior in mirtCAT
# setting normal prior with person specific mean and sd = prior_sd

# prior_sd gets overwritten by latent_covariance or in this script by run_local2
# run_local2 is used to overwrite run_local function in order to set prior_sd after the response pattern has been generated
# Single_run_test script.R shows that prior_sd is indeed overwritten
custom_den <- function(Theta, prior_mean,...) dnorm(Theta, mean = prior_mean, sd = prior_sd) 

#mirtCAT,customUpdateThetas function, utilizing prior
customUpdate_Thetas <- function(design, person, test){
  mo          <- extract.mirtCAT(test, 'mo') # Extract CATirt model
  responses   <- extract.mirtCAT(person, 'responses') # extract response patterns
  
  # Person specific means
  pp          <- extract.mirtCAT(design, 'person_properties') # extract person_properties (prior mean values)
  ID          <- extract.mirtCAT(person, "ID") # row associated with person properties
  prior_mean  <- pp[ID,] # select prior mean associated with person row
  
  # Estimate theta using prior custom density
  tmp         <- fscores(mo, response.pattern = responses, custom_den = custom_den, prior_mean = prior_mean,  
                         method = 'MAP')
  
  # Update theta estimate and SE 
  person$Update_thetas(tmp[,'F1'],
                       tmp[,'SE_F1', drop=FALSE])
  
  invisible()
}


#############################
### SIMULATION CONDITIONS ###
#############################
stop_SE       <- .316         # stopping SE

# Load item bank from Hummelen et al. (2021) 
load("coef_AMPD.Rdata")
model    <- generate.mirt_object(coef_AMPD, itemtype='graded') # Generate mirt object
n_items  <- nrow(coef_AMPD)

# load distribution proportions P(Global|Theta) from data by Hummelen et al. (2021)
# (Table 3 in paper), tabobs2 in script
load("Table_3.Rdata")


  ##################
  ### TRUE THETA ###
  ##################
  n_persons     <- 5000                         # number of simulees
  theta_true    <- cbind(rnorm(n_persons,0,1))  # standard normal, based on data by Hummelen et al. (2021) 
  ID            <- 1:n_persons                  # simulee ID
  
  
  #########################
  ### RESPONSE PATTERNS ###
  #########################
  # Generate response pattern from item parameters and theta
              Y     <- generate_pattern(model, Theta = theta_true) 

            ######################     
            ### GENERATE PRIOR ###
            ######################
              # Prior location from observed frequency
              # Sample from global score (global) using prob of global score given true theta P(global|T)
              theta_global    <- ifelse(theta_true > 3.499, 3.499, # Cap off at 3.499 as there is no empirical data
                                ifelse(theta_true < -3.499, -3.499,theta_true)) # same for -3.499
  
              # Generate global score based on theta value
              global <- rep(NA, length(theta_global))
              for(i in 1:length(theta_global)){
                col <- match(round(theta_global[i]),colnames(tabobs2)) # match rounded theta to theta value in tabobs2
                global[i] <- sample(0:4,1,prob = prop.table(tabobs2[,col])) # sample global score 
              }
              table(global)


              # Calculate the mean theta value for each global score to use as prior mean
              prop1       <- prop.table(tabobs2,1)
              prior_means  <- rowSums(prop1 %*% as.numeric(colnames(prop1)))
              
              # Calculate the variance for each global score to use as prior precision
              vars <- c()
              for(i in 1:nrow(tabobs2)){
                vars[i] <- sum(tabobs2[i,] * (as.numeric(colnames(tabobs2)) - prior_means[i])^2)/sum(tabobs2[i,])
              }
              names(vars) <- rownames(tabobs2)
            
                  ###############
                  ### RUN CAT ###
                  ###############
                  
                  ### PRIOR CONDITION
                  # Run fixed precision CAT with person specific prior
                  # split response vector into different global values
                  # take as prior mean max(P(theta|global))
                  # take as prior standard deviation, entropy|global
                  # Run cat
                  # Patch output back together
                  
                  out <- list()
                  out_full <- list()
                  ID_x <- list()
                  
                  for( j in names(vars)){
                    # Split response vector by global score to run separate CATs with each prior
                    Y_x <- Y[global == j,]
                    
                    if(!all(is.na(Y_x))){ # if no respondents have a prior score skip
                    
                    # Take mean theta for a given global score as prior mean
                    person_prior_properties <- data.frame(prior_mean = rep(round(prior_means[j],1),sum(global == j)))
                    
                    # Take associated variance as prior variance
                    prior_sd <- sqrt(unname(vars[match(j,names(vars))]))
                    
                    # Set starting theta on prior location
                    theta_start <- matrix(person_prior_properties$prior_mean, ncol = 1, 
                                          nrow = nrow(person_prior_properties)) 
                    
                    # CAT design 
                    design      <- list(min_SEM= stop_SE,max_items=n_items, 
                                        thetas.start = theta_start,
                                        customUpdateThetas = customUpdate_Thetas, 
                                        person_properties=person_prior_properties) 
                    # run CAT with prior  
                    out[[j]]         <- mirtCAT(mo=model, local_pattern=Y_x, start_item= "Trule",
                                         method="MAP", criteria="MI", design=design)

                    # save ID's for merging later
                    ID_x[[j]] <- ID[global == j]
                  }
                  }
                  
                  ID_prior <- unname(unlist(ID_x))
                  

                  ### BASELINE CAT 
                  
                  # Run baseline CAT ~N(0,1) prior
                  prior_sd    <- 1    # Reset prior sd to 1

                  # Baseline ~N(0,1) prior CAT design
                  design_base <-  list(min_SEM= stop_SE,max_items=n_items) 
                  
                  # CAT output with standard normal prior
                  out_base    <- mirtCAT(mo=model, local_pattern=Y[ID_prior,], start_item= "Trule",
                                         method="MAP", criteria="MI", design=design_base)

  
                ############################
                ### SAVE RELEVANT OUTPUT ###
                ############################
                
                  # sort values according to CAT administration order
                  global      <- factor(global[ID_prior], levels = names(vars))
                  theta_true  <- theta_true[ID_prior]
                  
                  # Empty lists for CAT output
                  ni <- list(); theta_est <- list(); theta_se <- list()
                  prior_SD <- list(); prior_MEAN <- list()
                  item <- list()
                
                for(z in names(vars)){
                  if(table(global,exclude = NULL)[z] > 1){ # Use sapply if multiple respondents
                ni[[z]]        <- sapply(out[[z]], function(x) length(x$items_answered))  # test length
                theta_est[[z]] <- sapply(out[[z]], function(x) x$thetas)                  # theta estimates
                theta_se[[z]]  <- sapply(out[[z]], function(x) x$SE_thetas)               # theta SE
                
                prior_SD[[z]]   <- sapply(out[[z]],function(x) x$thetas_SE_history[1])  # prior sd
                prior_MEAN[[z]] <- sapply(out[[z]], function(x) x$thetas_history[1])    # prior location
                
                  }else if(table(global,exclude = NULL)[z] == 1){ # otherwise not
                    ni[[z]]        <- length(out[[z]]$items_answered)  # test length
                    theta_est[[z]] <- out[[z]]$thetas                  # theta estimates
                    theta_se[[z]]  <- out[[z]]$SE_thetas               # theta SE

                    prior_SD[[z]]   <- out[[z]]$thetas_SE_history[1]  # prior sd
                    prior_MEAN[[z]] <- out[[z]]$thetas_history[1]     # prior location
                    
                  }
                }

                ni <- unname(unlist(ni))
                theta_est <- unname(unlist(theta_est))
                theta_se <- unname(unlist(theta_se))

                ni_base        <- sapply(out_base, function(x) length(x$items_answered))  # test length N(0,1) prior
                theta_est_base <- sapply(out_base, function(x) x$thetas)                  # theta estimates
                theta_se_base  <- sapply(out_base, function(x) x$SE_thetas)               # theta SE


                prior_SD <- unname(unlist(prior_SD))        # prior sd empirical prior
                prior_MEAN <- unname(unlist(prior_MEAN))    # prior mean empirical prior

                ################################
                ### STORE OUTPUT IN CSV FILE ###
                ################################
                
                dat         <- data.frame(cbind(ID_prior,global, theta_true,
                                                prior_SD,prior_MEAN,
                                                theta_est, theta_se,ni,
                                                theta_est_base, theta_se_base, ni_base), 
                                          stringsAsFactors = F)
                
                names(dat)  <- c("ID","global_score","theta_true",
                                 "prior_sd","prior_mean",
                                 "theta_est","theta_se","test_length",
                                 "theta_est_base","theta_se_base","test_length_base")
  
                # Save as csv files
                write.table(dat,file = paste0("Sim_AMPD_single2.csv"), sep = ";", dec = ",", row.names = F)




