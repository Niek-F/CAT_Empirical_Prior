###############################################
### Simulation 1 script - Niek Frans; 2022  ###
### Empirical Priors in Polytomous CATs:    ###
### Risks and Rewards in Clinical Settings  ###
###############################################

# Override run_local function to set prior_sd at a later stage
source("run_local2.R")
environment(run_local) <- asNamespace('mirtCAT')
assignInNamespace("run_local", run_local, ns = "mirtCAT",pos = "package:mirtCAT")

###############################################
### LOADING FUNCTIONS AND RELEVANT PACKAGES ###
###############################################

# Setting seed
set.seed(269513)

pkg <- c("mirt","mirtCAT","truncnorm")
sapply(pkg,library,character.only=TRUE,logical.return = TRUE)

# Function to set a personalized prior in mirtCAT
# setting normal prior with person specific mean and sd = prior_sd
custom_den <- function(Theta, prior_mean,...)
  dnorm(Theta, mean = prior_mean, sd = prior_sd)

#mirtCAT,customUpdateThetas function, utilizing prior
customUpdate_Thetas <- function(design, person, test){
  mo          <- extract.mirtCAT(test, 'mo') # Extract CATirt model
  responses   <- extract.mirtCAT(person, 'responses') # extract response patterns
  
  # Person specific means
  pp          <- extract.mirtCAT(design, 'person_properties')# extract person_properties (prior mean values)
  ID          <- extract.mirtCAT(person, "ID") # row associated with person properties
  prior_mean  <- pp[ID,]# select prior mean associated with person row
  
  # Estimate theta using prior custom density  
  tmp         <- fscores(mo, response.pattern = responses,
                         custom_den = custom_den, prior_mean = prior_mean,
                         method = 'MAP')
  
  # Update theta estimate and SE 
  person$Update_thetas(tmp[,'F1'],
                       tmp[,'SE_F1', drop=FALSE])
  invisible()
}

##########################################################
### DEFINING TRUE THETA GRID AND SIMULATION CONDITIONS ###
##########################################################
Iter          <- 1            # Starting iteration number
Iter_max      <- 100           # Maximum number of iterations

n_grid        <- 100                                        # number of simulees in each theta grid point
theta_true    <- cbind(rep(seq(-1,3,.5),each = n_grid))     # theta values matched to item bank information
n_persons     <- length(theta_true)
ID            <- 1:n_persons                                # simulee ID
stop_SE       <- .316                                       # CAT stopping criterion SE(theta) < .316

repeat{ # Repeat the following for Iter_max iterations
  ######################################
  ### DEFINING SIMULATION PARAMETERS ###
  ######################################
  n_items_vec       <- c(30,60)             # Number of items in item bank
  prior_sd_vec      <- c(1, .707, .5)       # Prior precision
  bias_vec          <- c(-2, -1, 0, 1, 2)   # Degree of bias in prior location
  min_items_vec     <- c(1,2,3,4)           # minimum number of items constraint
  
  #######################################
  ### ITEM BANK AND RESPONSE PATTERNS ###
  #######################################
  # For each item bank
  for(i in 1:length(n_items_vec)){
    n_items <- n_items_vec[i]
    
    # Set item parameters
    a1 <- rtruncnorm(n_items, a = 1.5, b = 5, mean = 3.5, sd = 1.0)     # item discrimination parameters
    
    b  <-  rnorm(n_items, mean = 2.2, sd = .4)                          # fourth threshold parameters
    steps <- rlnorm(n_items * 3,meanlog = log(.75), sdlog = log(1.2))   # difference between threshold parameters
    
    d4 <- -a1*b                                         # fourth threshold
    d3 <- d4 - -a1*steps[1:n_items]                     # third threshold
    d2 <- d3 - -a1*steps[(n_items+1):(2*n_items)]       # second threshold
    d1 <- d2 - -a1*steps[(2*n_items+1):length(steps)]   # first threshold
    
    g  <- rep(NA, n_items)                    # Lower asymptote (guessing parameter), set to NA
    
    param          <- data.frame(cbind(a1,d1,d2,d3,d4,g), 
                                 row.names = paste0('item', 1:n_items))
    
    model <- generate.mirt_object(param, itemtype='graded') # Generate mirt object
    Y     <- generate_pattern(model, Theta = theta_true)    # Generate item response data Y
    
    ##########################################
    ### MINIMUM NUMBER OF ITEMS CONSTRAINT ###
    ##########################################
    for(z in 1:length(min_items_vec)){
      min_items <- min_items_vec[z] # select constraint on minimum number of items
      
      # CAT design generic default standard normal prior
      person_prior_properties_base0 <- data.frame(prior_mean = cbind(rep(0,n_persons)))
      theta_start_base0 <- matrix(person_prior_properties_base0$prior_mean, ncol = 1, nrow = n_persons)
      
      design_base0 <-  list(min_SEM= stop_SE,max_items=n_items,min_items = min_items, 
                            thetas.start = theta_start_base0,
                            customUpdateThetas = customUpdate_Thetas,
                            person_properties=person_prior_properties_base0) 
      
      # CAT design generic clinical prior
      person_prior_properties_base2 <- data.frame(prior_mean = cbind(rep(1,n_persons)))
      theta_start_base2 <- matrix(person_prior_properties_base2$prior_mean, ncol = 1, nrow = n_persons)
      
      design_base2 <-  list(min_SEM= stop_SE,max_items=n_items,min_items = min_items, 
                            thetas.start = theta_start_base2,
                            customUpdateThetas = customUpdate_Thetas,
                            person_properties=person_prior_properties_base2) 
      
      # Making sure prior sd is reset to 1 for each loop
      prior_sd <- 1
      
      # CAT output generic standard normal prior 
      out_base0    <- mirtCAT(mo=model, local_pattern=Y, start_item= "Trule",
                              method="MAP", criteria="MI", design=design_base0)
      
      # CAT output generic standard normal prior 
      out_base2    <- mirtCAT(mo=model, local_pattern=Y, start_item= "Trule",
                              method="MAP", criteria="MI", design=design_base2)
      
      ########################################
      ### DEGREE OF BIAS IN PRIOR LOCATION ###
      ########################################
      for(j in 1:length(bias_vec)){
        bias <- bias_vec[j]       # Select bias condition from the vector
        
        # Set person specific prior to theta_true + bias
        person_prior_properties <- data.frame(prior_mean = theta_true + bias)
        
        #######################    
        ### PRIOR PRECISION ###
        #######################
        for(q in 1:length(prior_sd_vec)){ 
          prior_sd <- prior_sd_vec[q]
          
          ###############
          ### RUN CAT ###
          ###############
          # setting starting values for theta as prior location
          theta_start <- matrix(person_prior_properties$prior_mean, ncol = 1, nrow = n_persons)
          
          # CAT design personalized prior
          design      <- list(min_SEM= stop_SE,max_items=n_items,min_items = min_items, 
                              thetas.start = theta_start,
                              customUpdateThetas = customUpdate_Thetas, 
                              person_properties=person_prior_properties) 

          # CAT output personalized prior          
          out         <- mirtCAT(mo=model, local_pattern=Y, start_item= "Trule",
                                 method="MAP", criteria="MI", design=design)
          
          ############################
          ### SAVE RELEVANT OUTPUT ###
          ############################
          # Personalized prior
          ni        <- sapply(out, function(x) length(x$items_answered))  # test length 
          theta_est <- sapply(out, function(x) x$thetas)                  # theta estimates
          theta_se  <- sapply(out, function(x) x$SE_thetas)               # theta standard error
          
          # Generic standard normal prior
          ni_base0        <- sapply(out_base0, function(x) length(x$items_answered))  # test length 
          theta_est_base0 <- sapply(out_base0, function(x) x$thetas)                  # theta estimates
          theta_se_base0  <- sapply(out_base0, function(x) x$SE_thetas)               # theta standard error
          
          # Generic clinical prior
          ni_base2        <- sapply(out_base2, function(x) length(x$items_answered))  # test length 
          theta_est_base2 <- sapply(out_base2, function(x) x$thetas)                  # theta estimates
          theta_se_base2  <- sapply(out_base2, function(x) x$SE_thetas)               # theta standard error
          
          ################################
          ### STORE OUTPUT IN CSV FILE ###
          ################################
          dat         <- data.frame(cbind(Iter,ID,n_items,bias, 
                                          prior_sd,min_items, theta_true, 
                                          theta_est, theta_se,ni,
                                          theta_est_base0, theta_se_base0,
                                          theta_est_base2, theta_se_base2,
                                          ni_base0,ni_base2), 
                                    stringsAsFactors = F)
          
          names(dat)  <- c("Iter","ID","bank_size","prior_bias","prior_sd",
                           "min_items","theta_true",
                           "theta_est","theta_se","test_length",
                           "theta_est_base0","theta_se_base0",
                           "theta_est_base2","theta_se_base2",
                           "test_length_base0","test_length_base2")
          
          # Save as csv files
          filename    <- paste("Sim_poly",Iter,"_I",n_items,"min",min_items,
                               "_bias",bias,"_sigma",prior_sd,".csv",sep = "")
          write.table(dat,file = paste0("Output APM comp/",filename), sep = ";",
                      dec = ",", row.names = F)
          
          
        }
      }
    }
  }
  
  # Stop simulating after 100 repetitions for each cell  
  Iter      <- Iter +1 # Increase iteration by 1 after each loop
  
  if(Iter > Iter_max){
    break
  }else{
    end_seed  <- .Random.seed
    save(end_seed, file = paste0("Output APM comp/","Random_seed_endIter",Iter-1,".R")) # Save seed 
    next
  }
}