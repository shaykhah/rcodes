# rcodes
#' gfr-rbf: Generate and fit the GFR model and its extensions GLM object.
#' @description Generalised functional responses (GFRs) describe the response of
#' organisms to different habitats as a function of the availability of all habitats
#'  within the organisms' reach.
#'
#' @param formula The basic glm formula involving just the main effects
#' @param data The data frame of response and explanatory variable data
#' @param family The likelihood function to be used for fitting
#' @param blockName The name of a blocking factor in the data frame indicating the
#'        membership of each row to a sampling instance
#' @param addexp Additional explanatory variables that represent values applying
#'         uniformly to an entire sampling instance
#' @param order The required number of basis functions
#' @param step Whether to apply stepwise model selection to reduce the terms in the
#'         model (default is FALSE)
#' @param blockName The name of the blocking factor in the dataframe
#' @param expFlag Sometimes (e.g. in use-availability designs) only the
#'        availability points are required to calculate the expectations
#' @return The fitted rbf-gfr model. Containing the following arguments:
#' @return \code{model} The glm model object fitted to the data
#' @return \code{I_{i,j}} The calculated I from the data
#' @return \code{order} The requested number of basis function
#' @return \code{formula} The basic formula (not containing interactions)
#' @export
#
#' @seealso {gfr} in HATOPO library
#######################################################################################
# Part 1: Training of model using the RBF-GFR function
#######################################################################################


rbf_gfr<-function(formula, data, family, blockName, expFlag=NULL, addexp,
                  order=1, step=FALSE, offset=NULL, Iij_order, temp_sqr_term, g)
{
  
  # Preparaing the data
  order_id <- paste0("data$", blockName)
  
  data<-data[order(eval(parse(text=order_id))), ]
  
  # Converting blockName as factor
  block <- as.factor(eval(parse(text=paste("data$",blockName))))
  
  formula <- as.formula(formula)
  
  #extracts the names of all the formula terms (including quadratics)
  vs<-c(attr(terms(formula),"term.labels"))
  
  # Extracts the response variable
  y<-all.vars(formula)[1]
  
  # Extracts the main effects
  mainVars<-all.vars(formula)[-1]
  
  # Creates dataframe that will be augmented with the I's
  newdata<-data
  
  # Getting the unique id of block
  cats<-as.integer(unique(block))
  
  # Calculating number of times to repeat the I_i_j
  repeat_times <- floor(nrow(data)/length(cats))
  
  # Generate a Gaussian mixture approximation of the habitat composition
  # of different id.
  dat<- data[mainVars]
  
  fa<-favail(dat, blocking= block, G=g)
  
  exps<-as.data.frame(matrix(vector(), nrow=length(cats), ncol=0,
                             dimnames=list(cats,c())))
  
  # I's
  if(length(expFlag)==0) expFlag=rep(T, length(blockName))
  expVars<-c()
  
  
  # making a dataframe to store quantile_vars for each main var across all Iij
  # row : main vars, col:
  quantile_vars_df <- data.frame(matrix(NA, nrow = 0, ncol=Iij_order+2))
  
  # Iterate over number of basis functions:
  
  
  # Iterate over each main variable for Iij calculation and dataframe updation
  for(i in 1:length(mainVars))
  {
    # Fetching the data of main variable from the dataset
    rawCol<-eval(parse(text=paste0("data$",mainVars[i],"[expFlag]")))
    
    # ordering the data
    rawCol <- rawCol[order(rawCol)]
    
    # Getting the quantiles of the variables
    quantile_vars <- quantile(rawCol, probs = seq(0,1,by=1/(Iij_order+1)))
    
    # Converting quantile into vector
    quantile_vars <- c(as.matrix(quantile_vars))
    
    # Assigning quantiles in dataframe to use in prediction function
    quantile_vars_df[nrow(quantile_vars_df)+1,] <- quantile_vars
    
    # Intitializing the vector to store the difference between quantiles
    quantile_diff <- vector(mode = 'numeric', length = length(quantile_vars)-1)
    
    # Iterating over each quantile except 1st one
    for(l in 2:length(quantile_vars)) {
      quantile_diff[l-1] <- quantile_vars[l]-quantile_vars[l-1]
    }
    
    # calculation of favail function, mean, variance and weight for the variable
    # Defining the matrix
    mu_mat <- matrix(0,nrow =length(cats),ncol = g)
    var_mat <- matrix(0,nrow =length(cats),ncol = 1 )
    weight_mat <- matrix(0,nrow =length(cats),ncol = g)
    
    # Iterate over each id
    for(bl in 1:length(cats)){
      # Variance calculation
      var_mat[bl,1]<- fa[[bl]]$var[i]
      
      # Iterater over each component of approximation
      for(k in 1:g){
        # Mean calcualtion
        mu_mat[bl,k]<- fa[[bl]]$mus[i,k]
        
        # Weight calcuation
        weight_mat[bl,k] <- fa[[bl]]$ws[k]
        
      }
    }
    
    # Intitializing the vector to store the max of two quantile differences
    quantile_max <- vector(mode = 'numeric', length = length(quantile_diff)-1)
    
    # Iterating over each differnce of quantile calculated above
    for(df in 1:(length(quantile_diff)-1)) {
      
      # Getting the max of two quantile_max
      quantile_max[df] <- max(quantile_diff[df],quantile_diff[df+1])
      
      
      i_j_m<- matrix(0, nrow =length(cats),ncol = g)
      
      # Iterate over each id
      for(bl in 1:length(cats)){
        # Iterate over each component
        for(k in 1:g){
          
          i_j_m[bl,k]<-weight_mat[bl,k]*(quantile_max[df]/sqrt(var_mat[bl]+quantile_max[df]^2))*exp(-(mu_mat[bl,k]-quantile_vars[df+1])^2/(2*(var_mat[bl]+quantile_max[df]^2)))
          
        }
      }
      
      Iij <- rowSums(i_j_m)
      
      # Making the name of main variable's Iij order
      exp_name <- paste0(mainVars[i],i,df, sep = "")
      expVars<-c(expVars, exp_name)
      
      exps[,paste(exp_name)]<-Iij
      
      mat_Ijm <- rep(Iij, each=repeat_times)
      
      # Appending the calculated in the dataset
      newdata[,paste(exp_name)] <- mat_Ijm
    }
    
    # Storing quantile_max in df for using with test data set

  }
  
  
  i<-0
  while(i<length(addexp)) {
    i<-i+1
    
    rawCol <- eval(parse(text=paste("data$",addexp[i])))
    
    # Getting the quantiles of the variables
    quantile_vars <- quantile(rawCol, probs = seq(0,1,by=1/(Iij_order+1)))
    
    # Converting quantile into vector
    quantile_vars <- c(as.matrix(quantile_vars))
    
    quantile_vars_df[nrow(quantile_vars_df)+1, ] <- quantile_vars
    
    # Intitializing the vector to store the difference between quantiles
    quantile_diff <- vector(mode = 'numeric', length = length(quantile_vars)-1)
    
    
    # Iterating over each quantile except 1st one
    for(l in 2:length(quantile_vars)) {
      quantile_diff[l-1] <- quantile_vars[l]-quantile_vars[l-1]
    }
    
    # calculation of favail function, mean, variance and weight for the variable
    mu_mat <- matrix(0,nrow =length(cats),ncol = g)
    var_mat <- matrix(0,nrow =length(cats),ncol = 1 )
    weight_mat <- matrix(0,nrow =length(cats),ncol = g)
    
    for(bl in 1:length(cats)){
      var_mat[bl,1]<- fa[[bl]]$var[i]
      
      for(k in 1:g){
        mu_mat[bl,k]<- fa[[bl]]$mus[i,k]
        weight_mat[bl,k] <- fa[[bl]]$ws[k]
        
      }
    }
    # Intitializing the vector to store the max of two quantile differences
    quantile_max <- vector(mode = 'numeric', length = length(quantile_diff)-1)
    
    # Iterating over each difference of quantile calculated above
    for(df in 1:(length(quantile_diff)-1)) {
      quantile_max[df] <- max(quantile_diff[df],quantile_diff[df+1])
      
      i_j_m<- matrix(0,nrow =length(cats),ncol = g)
      
      for(bl in 1:length(cats)){
        for(k in 1:g){
          
          i_j_m[bl,k]<-weight_mat[bl,k]*(quantile_max[df]/sqrt(var_mat[bl]+quantile_max[df]^2))*exp(-(mu_mat[bl,k]-quantile_vars[df+1])^2/(2*(var_mat[bl]+quantile_max[df]^2)))
          
        }
      }
      
      Iij <- rowSums(i_j_m)
      
      exp_name <- paste0(addexp[i],i,df, sep = "")
      
      exps[,paste(exp_name)]<-Iij
      
    }
    
    #    quantile_max_df[nrow(quantile_max_df)+1, ] <- quantile_max
    
    exp_name <- paste0('add_',addexp[i], sep = "")
    expVars<-c(expVars, exp_name)
    names(newdata)<-gsub(addexp[i], exp_name, names(newdata))
    
  }
  if(temp_sqr_term==TRUE) {
    # Adding squared term in the dataset
    newdata$temp2 <- data$temp^2
    
    # Adding the squared term in vector 'vs' that would go in main variable list
    vs <- c(vs, 'temp2')
    
  }
  # Generation of rbf-gfr formula
  eMain<-paste(vs, collapse="+")
  eExp<-paste(expVars, collapse="+")
  txt<-paste(y,"~(",eMain,")*(",eExp,")", sep="")
  
  
  newformula<-as.formula(txt)
  newformula<-update(newformula, newformula)
  
  # Model fitting and selection
  modFull<-glm(newformula, data=newdata, family=family)
  mod<-modFull
  
  if(step==TRUE)
  {
    modMains<-glm(formula(paste(y,"~",paste(mainVars, collapse="+"))), data=newdata, family=family)
    mod<-step(modFull, direction="both",scope=list(upper=modFull,lower=modMains))
  }
  return(list("model"=mod,"I's"=exps, "order"=order,
              "formula"=formula, "blockName"=blockName, "addexp"=addexp,
              "quantiles"= quantile_vars_df, "traindata"=data))
  
}

#######################################################################################
# Part 2:  Prediction using the trained RBF-GFR model
#######################################################################################

gfr.predict<-function(model_gfr, newdata, f_type,
                      se_fit, quantiles_df,
                      main_formula, order, addexp)
{
  # basis to be used for I's (Inherrited from model object)
  mod<-model_gfr$model
  
  # Ordering the newdata as per the id
  newdata <- newdata[order(newdata$id), ]
  
  # Extracts the main effects
  mainVars<-all.vars(main_formula)[-1]
  
  # Getting the main variables for favail function
  dat<- newdata[mainVars]
  
  # Getting the unique ids
  block <- unique(newdata$id)
  
  g <- par_G
  
  # Favail function
  fa<-favail(dat, blocking= newdata$id, G=g)
  
  # Defining a dataframe a store I's for the unique ids
  exps <- as.data.frame(matrix(vector(), nrow=length(block), ncol=0, dimnames=list(block,c())))
  
  # Calculating number of times to repeat the I_i_j
  repeat_times <- floor(nrow(newdata)/length(block))
  
  # Iterate over basis, it is kept as 1
  
  
  # Iterate over each main variable for Iij calculation and dataframe updation
  for(i in 1:length(mainVars))
  {
    # Converting quantile into vector
    quantile_vars <- as.numeric(quantiles_df[i, ])
    
    # Intitializing the vector to store the difference between quantiles
    quantile_diff <- vector(mode = 'numeric', length = length(quantile_vars)-1)
    
    # Iterating over each quantile except 1st one
    for(l in 2:length(quantile_vars)) {
      quantile_diff[l-1] <- quantile_vars[l]-quantile_vars[l-1]
    }
    
    # calculation of favail function, mean variance and weight for the variable
    # Defining the matrix
    mu_mat     <- matrix(0, nrow =length(block), ncol = g)
    var_mat    <- matrix(0, nrow =length(block), ncol = 1 )
    weight_mat <- matrix(0, nrow =length(block), ncol = g)
    
    # Iterate over each id
    for(bl in 1:length(block)){
      
      # Variance calculation
      var_mat[bl,1]<- fa[[bl]]$var[i]
      
      # Iterater over each component of approximation
      for(k in 1:g){
        # Mean calcualtion
        mu_mat[bl,k]<- fa[[bl]]$mus[i,k]
        
        # Weight calcuation
        weight_mat[bl,k] <- fa[[bl]]$ws[k]
        
      }
    }
    
    # Intitializing the vector to store the max of two quantile differences
    quantile_max <- vector(mode = 'numeric', length = length(quantile_diff)-1)
    
    # Iterating over each differnce of quantile calculated above
    for(df in 1:(length(quantile_diff)-1)) {
      
      # Getting the max of two quantile_max
      quantile_max[df] <- max(quantile_diff[df],quantile_diff[df+1])
      
      # Defining matrix to store the calculations
      i_j_m<- matrix(0, nrow =length(block),ncol = g)
      
      # Iterate over each id
      for(bl in 1:length(block)){
        # Iterate over each component
        for(k in 1:g){
          
          i_j_m[bl,k]<-weight_mat[bl,k]*(quantile_max[df]/sqrt(var_mat[bl]+quantile_max[df]^2))*exp(-(mu_mat[bl,k]-quantile_vars[df+1])^2/(2*(var_mat[bl]+quantile_max[df]^2)))
          
        }
      }
      
      Iij <- rowSums(i_j_m)
      
      # Making the name of main Iij 
      exp_name <- paste0(mainVars[i],i,df, sep = "")
      
      # Assigning Iij to basis df
      exps[,paste(exp_name)]<-Iij
      
      # Repeating Iij
      mat_Ijm <- rep(Iij, each=repeat_times)
      
      # Appending the calculated in the dataset
      newdata[,paste(exp_name)] <- mat_Ijm
    }
  }
  
  
  # Loop over the addexp: output of the loop go into only exps df
  i<-0
  while(i<length(addexp)) {
    i<-i+1
    
    # Converting quantile into vector
    quantile_vars <- as.numeric(quantiles_df[length(mainVars)+1,])
    
    # Intitializing the vector to store the difference between quantiles
    quantile_diff <- vector(mode = 'numeric', length = length(quantile_vars)-1)
    
    # Iterating over each quantile except 1st one
    for(l in 2:length(quantile_vars)) {
      quantile_diff[l-1] <- quantile_vars[l]-quantile_vars[l-1]
    }
    
    # calculation of favail function, mean variance and weight for the variable
    mu_mat <- matrix(0,nrow =length(block),ncol = g)
    var_mat <- matrix(0,nrow =length(block),ncol = 1 )
    weight_mat <- matrix(0,nrow =length(block),ncol = g)
    
    for(bl in 1:length(block)){
      var_mat[bl,1]<- fa[[bl]]$var[i]
      
      for(k in 1:g){
        mu_mat[bl,k]<- fa[[bl]]$mus[i,k]
        weight_mat[bl,k] <- fa[[bl]]$ws[k]
        
      }
    }
    # Intitializing the vector to store the max of two quantile differences
    quantile_max <- vector(mode = 'numeric', length = length(quantile_diff)-1)
    
    # Iterating over each differnce of quantile calculated above
    for(df in 1:(length(quantile_diff)-1)) {
      quantile_max[df] <- max(quantile_diff[df],quantile_diff[df+1])
      
      i_j_m<- matrix(0,nrow =length(block),ncol = g)
      
      for(bl in 1:length(block)) {
        for(k in 1:g) {
          
          i_j_m[bl,k]<-weight_mat[bl,k]*(quantile_max[df]/sqrt(var_mat[bl]+quantile_max[df]^2))*exp(-(mu_mat[bl,k]-quantile_vars[df+1])^2/(2*(var_mat[bl]+quantile_max[df]^2)))
          
        }
      }
      
      Iij <- rowSums(i_j_m)
      
      exp_name <- paste0(addexp[i],i,df, sep = "")
      
      exps[,paste(exp_name)]<-Iij
      
    }
    
    exp_name <- paste0('add_',addexp[i], sep = "")
    
    names(newdata)<-gsub(addexp[i], exp_name, names(newdata))
    
  }
  # Adding squared term in the dataset
  newdata$temp2 <- test$temp^2
  
  
  print(str(newdata))
  # Prediction using the trained model and test dataset
  predictions <-  predict(model_gfr, newdata = newdata,  type = c("response"))
  
  return(list("predictions"=predictions, "I's"=exps))
}


# Appling the first simulated data as an example: 
library(HATOPO)

# Assign formula, use variable names (food & temperature)
par_formula <- use ~food + temp

# Define the likelihood function to be used for fitting
par_family <- poisson

# name of the variable for blocking id, just name
par_block_id <- 'id'

# Additional variable to be considered in
par_addexp <- 'N'

# It is the order of model
par_basis<- 1

# to apply stepwise model selection
par_step <- FALSE

par_offset <- NULL

par_expFlag <- NULL

# number of basis functions
par_Iij_order <- 10

# Set the parameter to involve temp_sqaure term in calculation, set FALSE to exclude
temp_sqr_term <- TRUE
# The number of Gaussian components to be used for the approximating mixture
par_G <- 9

# Order data by id
habitatDat<-habitatDat[order(habitatDat$id),]  

# Create a training data set
set.seed(123)
ids<-sample(as.numeric(unique(habitatDat$id)),360)
train<-na.omit(habitatDat[habitatDat$id%in%ids,])
# Create a testing data set
test<-na.omit(habitatDat[!habitatDat$id%in%ids,])

# assign dataset,
par_data <- train
list_of_outputs <- rbf_gfr(par_formula,   par_data,      par_family,
                       par_block_id,  par_expFlag,   par_addexp,
                       par_basis,     par_step,      par_offset,
                       par_Iij_order, temp_sqr_term, par_G)


list_of_outputs$model


############################### Prediction Function Call ##############################

gfr_predict <- gfr.predict(list_of_outputs$model, test,
                           "response", se.fit,
                           list_of_outputs$quantiles, list_of_outputs$formula,
                           list_of_outputs$order, list_of_outputs$addexp)

# Processing the output of the prediction function to generate output file
predictions<- as.data.frame(gfr_predict[[1]])

# calculating the out-of-sample score
R2<- 1 - (sum((test$use - predictions)^2)/sum((test$use - mean(test$use))^2))


################################################################################################################################
############################### Using another way (step-by-step code) to apply the RBF-GFR model ##############################
################################################################################################################################



