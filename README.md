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


# Initialisation & libraries
require(HATOPO)

# Order data by id
habitatDat<-habitatDat[order(habitatDat$id),]  

# Create a training data set
set.seed(123)
ids<-sample(as.numeric(unique(habitatDat$id)),360)
train<-na.omit(habitatDat[habitatDat$id%in%ids,])
# Create a testing data set
test<-na.omit(habitatDat[!habitatDat$id%in%ids,])

  #######################
  data<- as.data.frame(train)
  #extract the main variables
  food.data<- subset.data.frame(data,select=c(food))
  temp.data<- subset.data.frame(data,select=c(temp))
  ### order each variable
  food.data<-food.data[order(food.data),]
  temp.data<-temp.data[order(temp.data),]
  
  # Getting the quantiles of the variables & Converting quantile into vector
  x<-  c(as.matrix(quantile(food.data,probs = seq(0,1,by=1/11))))
  y<- c(as.matrix(quantile(temp.data,probs = seq(0,1,by=1/11))))

  
  # calculation the variances of the basis functions using quantiles
  f.sd1<- max(x[3]-x[2],x[2]-x[1])
  f.sd2<- max(x[3]-x[2],x[4]-x[3])
  f.sd3<- max(x[4]-x[3],x[5]-x[4])
  f.sd4<- max(x[5]-x[4], x[6]-x[5])
  f.sd5<- max( x[6]-x[5],x[7]-x[6])
  f.sd6<- max(x[7]-x[6], x[8]-x[7])
  f.sd7<- max( x[8]-x[7], x[9]-x[8])
  f.sd8<- max(x[9]-x[8] ,x[10]-x[9])
  f.sd9<- max(x[10]-x[9],x[11]-x[10])
  f.sd10<- max(x[11]-x[10],x[12]-x[11])
  
  
  t.sd1<- max(y[3]-y[2],y[2]-y[1])
  t.sd2<- max(y[3]-y[2],y[4]-y[3])
  t.sd3<- max(y[4]-y[3],y[5]-y[4])
  t.sd4<- max(y[5]-y[4], y[6]-y[5])
  t.sd5<- max(y[6]-y[5], y[7]-y[6])
  t.sd6<- max(y[7]-y[6],y[8]-y[7])
  t.sd7<- max(y[8]-y[7], y[9]-y[8])
  t.sd8<- max( y[9]-y[8],y[10],y[9])
  t.sd9<- max( y[10],y[9], y[11]-y[10])
  t.sd10<- max( y[11]-y[10],y[12]-y[11])
  
  x<- x[2:11]
  y<- y[2:11]

  # Generate a Gaussian mixture approximation of the habitat composition of different id.
  
  dat<- subset.data.frame( data,select=c(food,temp))
  g<-9
  set.seed(123)
  fa<-favail(dat, blocking=data$id, G=g)
  
  ##### Extract the mean in Mixture Gaussian approximation  
    # Extracts unique values of blocking index
  jv<-unique(data$id) 
  jvmax<-length(jv)
    #extracts the names of all the formula terms
  vs<-colnames(dat)
  comp<- 1:g
  # Extracts the main effects
  mainVars<-length(vs)
  jvmax<- 1:length(jv)
  mainVars<- 1:length(vs)
  
  ### mean of the mixture Gaussian distribution of the first variable for all block and each component.
  m1<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      m1[i,k]<- fa[[i]]$mus[1,k]
      
    }
  }
  ### mean of the second variable for all block and each component 
  m2<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      m2[i,k]<- fa[[i]]$mus[2,k]
      
    }
  }

  #### Variance of the mixture Gaussian approximation
  ### variance of the first variable "food" of all block
  c1<- matrix(0,nrow =length(jvmax),ncol = 1 )
  
  for(i in 1:length(jvmax)){
    
    c1[i]<- fa[[i]]$var[1]
    c1
  }
  ### variance of all block for second variable 
  c2<- matrix(0,nrow =length(jvmax),ncol = 1 )
  
  for(i in 1:length(jvmax)){
    
    c2[i]<- fa[[i]]$var[2]
    c2
  }

  #### Matrix of weight of the mixture Gaussian approximation  
  ### weight for all component for each block
 
  w<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      w[i,k]<- fa[[i]]$ws[k]
    }
  }

  
  ### I for 'food' the first variable,j=1, and basis function m=1

  i11<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i11[i,k]<-w[i,k]*(f.sd1/sqrt(c1[i]+f.sd1^2))*exp(-(m1[i,k]-x[1])^2/(2*(c1[i]+f.sd1^2)))
      
      i11[i,k]
    }
  }
  s11<- as.matrix(rowSums(i11)) 
  ### I for 'food' the first variable,j=1, and basis function m=2
  i12<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i12[i,k]<-w[i,k]*(f.sd2/sqrt(c1[i]+f.sd2^2))*exp(-(m1[i,k]-x[2])^2/(2*(c1[i]+f.sd2^2)))
      i12[i,k]
    }
  }
  s12<- as.matrix(rowSums(i12)) 
  ### I for 'food' the first variable,j=1, and basis function m=3
  i13<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i13[i,k]<-w[i,k]*(f.sd3/sqrt(c1[i]+f.sd3^2))*exp(-(m1[i,k]-x[3])^2/(2*(c1[i]+f.sd3^2)))
      i13[i,k]
    }
  }
  s13<- as.matrix(rowSums(i13)) 
  ### I for 'food' the first variable,j=1, and basis function m=4
  i14<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i14[i,k]<-w[i,k]*(f.sd4/sqrt(c1[i]+f.sd4^2))*exp(-(m1[i,k]-x[4])^2/(2*(c1[i]+f.sd4^2)))
      i14[i,k]
    }
  }
  s14<- as.matrix(rowSums(i14)) 
  ### I for 'food' the first variable,j=1, and basis function m=5
  i15<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i15[i,k]<-w[i,k]*(f.sd5/sqrt(c1[i]+f.sd5^2))*exp(-(m1[i,k]-x[5])^2/(2*(c1[i]+f.sd5^2)))
      i15[i,k]
    }
  }
  s15<- as.matrix(rowSums(i15)) 
  ### I for 'food' the first variable,j=1, and basis function m=6
  i16<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i16[i,k]<-w[i,k]*(f.sd6/sqrt(c1[i]+f.sd6^2))*exp(-(m1[i,k]-x[6])^2/(2*(c1[i]+f.sd6^2)))
      i16[i,k]
    }
  }
  s16<- as.matrix(rowSums(i16)) 

  
  ### I for 'food' the first variable,j=1, and basis function m=7
  i17<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i17[i,k]<-w[i,k]*(f.sd7/sqrt(c1[i]+f.sd7^2))*exp(-(m1[i,k]-x[7])^2/(2*(c1[i]+f.sd7^2)))
      i17[i,k]
    }
  }
  s17<- as.matrix(rowSums(i17)) 
  ### I for 'food' the first variable,j=1, and basis function m=8
  i18<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i18[i,k]<-w[i,k]*(f.sd8/sqrt(c1[i]+f.sd8^2))*exp(-(m1[i,k]-x[8])^2/(2*(c1[i]+f.sd8^2)))
      i18[i,k]
    }
  }
  s18<- as.matrix(rowSums(i18)) 
  ### I for 'food' the first variable,j=1, and basis function m=9
  i19<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i19[i,k]<-w[i,k]*(f.sd9/sqrt(c1[i]+f.sd9^2))*exp(-(m1[i,k]-x[9])^2/(2*(c1[i]+f.sd9^2)))
      i19[i,k]
    }
  }
  s19<- as.matrix(rowSums(i19)) 
  ### I for 'food' the first variable,j=1, and basis function m=10
  i110<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i110[i,k]<-w[i,k]*(f.sd10/sqrt(c1[i]+f.sd10^2))*exp(-(m1[i,k]-x[10])^2/(2*(c1[i]+f.sd10^2)))
      i110[i,k]
    }
  }
  s110<- as.matrix(rowSums(i110)) 

  
  ### I for 'temp' the second variable,j=1, and basis function m=1
  i21<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i21[i,k]<-w[i,k]*t.sd1/sqrt(c2[i]+t.sd1^2)*exp(-(m2[i,k]-y[1])^2/(2*(c2[i]+t.sd1^2)))
      i21[i,k]
    }
  }
  s21<- as.matrix(rowSums(i21)) 
  ### I for 'temp' the second variable,j=2, and basis function m=2
  i22<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i22[i,k]<-w[i,k]*t.sd2/sqrt(c2[i]+t.sd2^2)*exp(-(m2[i,k]-y[2])^2/(2*(c2[i]+t.sd2^2)))
      i22[i,k]
    }
  }
  s22<- as.matrix(rowSums(i22)) 
  ### I for 'temp' the second variable,j=2, and basis function m=3
  i23<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i23[i,k]<-w[i,k]*t.sd3/sqrt(c2[i]+t.sd3^2)*exp(-(m2[i,k]-y[3])^2/(2*(c2[i]+t.sd3^2)))
      i23[i,k]
    }
  }
  s23<- as.matrix(rowSums(i23)) 
  ### I for 'temp' the second variable,j=2, and basis function m=4
  i24<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i24[i,k]<-w[i,k]*t.sd4/sqrt(c2[i]+t.sd4^2)*exp(-(m2[i,k]-y[4])^2/(2*(c2[i]+t.sd4^2)))
      i24[i,k]
    }
  }
  s24<- as.matrix(rowSums(i24)) 
  ### I for 'temp' the second variable,j=2, and basis function m=5
  i25<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i25[i,k]<-w[i,k]*t.sd5/sqrt(c2[i]+t.sd5^2)*exp(-(m2[i,k]-y[5])^2/(2*(c2[i]+t.sd5^2)))
      i25[i,k]
    }
  }
  s25<- as.matrix(rowSums(i25)) 
  ### I for 'temp' the second variable,j=2, and basis function m=6
  i26<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i26[i,k]<-w[i,k]*t.sd6/sqrt(c2[i]+t.sd6^2)*exp(-(m2[i,k]-y[6])^2/(2*(c2[i]+t.sd6^2)))
      i26[i,k]
    }
  }
  s26<- as.matrix(rowSums(i26)) 
  
  ### I for 'temp' the second variable,j=2, and basis function m=7
  i27<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i27[i,k]<-w[i,k]*t.sd7/sqrt(c2[i]+t.sd7^2)*exp(-(m2[i,k]-y[7])^2/(2*(c2[i]+t.sd7^2)))
      i27[i,k]
    }
  }
  s27<- as.matrix(rowSums(i27)) 
  ### I for 'temp' the second variable,j=2, and basis function m=8
  i28<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i28[i,k]<-w[i,k]*t.sd8/sqrt(c2[i]+t.sd8^2)*exp(-(m2[i,k]-y[8])^2/(2*(c2[i]+t.sd8^2)))
      i28[i,k]
    }
  }
  s28<- as.matrix(rowSums(i28)) 
  
  ### I for 'temp' the second variable,j=2, and basis function m=9
  i29<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i29[i,k]<-w[i,k]*t.sd9/sqrt(c2[i]+t.sd9^2)*exp(-(m2[i,k]-y[9])^2/(2*(c2[i]+t.sd9^2)))
      i29[i,k]
    }
  }
  s29<- as.matrix(rowSums(i29))  
  ### I for 'temp' the second variable,j=2, and basis function m=10
  i210<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i210[i,k]<-w[i,k]*t.sd10/sqrt(c2[i]+t.sd10^2)*exp(-(m2[i,k]-y[10])^2/(2*(c2[i]+t.sd10^2)))
      i210[i,k]
    }
  }
  s210<- as.matrix(rowSums(i210)) 
  
  ### Creat columns I11,I12, ....
  data$I11<- c(rep(c(s11[1:length(s11)]), each=500))
  data$I12<- c(rep(c(s12[1:length(s11)]), each=500))
  data$I13<- c(rep(c(s13[1:length(s11)]), each=500))
  data$I14<- c(rep(c(s14[1:length(s11)]), each=500))
  data$I15<- c(rep(c(s15[1:length(s11)]), each=500))
  data$I16<- c(rep(c(s16[1:length(s11)]), each=500))
  data$I17<- c(rep(c(s17[1:length(s11)]), each=500))
  data$I18<- c(rep(c(s18[1:length(s11)]), each=500))
  data$I19<- c(rep(c(s19[1:length(s11)]), each=500))
  data$I110<- c(rep(c(s110[1:length(s11)]), each=500))
  
  
  data$I21<- c(rep(c(s21[1:length(s11)]), each=500))
  data$I22<- c(rep(c(s22[1:length(s11)]), each=500))
  data$I23<- c(rep(c(s23[1:length(s11)]), each=500))
  data$I24<- c(rep(c(s24[1:length(s11)]), each=500))
  data$I25<- c(rep(c(s25[1:length(s11)]), each=500))
  data$I26<- c(rep(c(s26[1:length(s11)]), each=500))
  data$I27<- c(rep(c(s27[1:length(s11)]), each=500))
  data$I28<- c(rep(c(s28[1:length(s11)]), each=500))
  data$I29<- c(rep(c(s29[1:length(s11)]), each=500))
  data$I210<- c(rep(c(s210[1:length(s11)]), each=500))
  data$temp2<- data$temp^2

  #Interaction Terms
  data$x11<- data$food*data$I11
  data$x12<- data$food*data$I12
  data$x13<- data$food*data$I13
  data$x14<- data$food*data$I14
  data$x15<- data$food*data$I15
  data$x16<- data$food*data$I16
  data$x17<- data$food*data$I17
  data$x18<- data$food*data$I18
  data$x19<- data$food*data$I19
  data$x110<- data$food*data$I110

  data$x21<- data$temp*data$I21
  data$x22<- data$temp*data$I22
  data$x23<- data$temp*data$I23
  data$x24<- data$temp*data$I24
  data$x25<- data$temp*data$I25
  data$x26<- data$temp*data$I26
  data$x27<- data$temp*data$I27
  data$x28<- data$temp*data$I28
  data$x29<- data$temp*data$I29
  data$x210<- data$temp*data$I210
  
  data$y11<- data$temp*data$I11
  data$y12<- data$temp*data$I12
  data$y13<- data$temp*data$I13
  data$y14<- data$temp*data$I14
  data$y15<- data$temp*data$I15
  data$y16<- data$temp*data$I16
  data$y17<- data$temp*data$I17
  data$y18<- data$temp*data$I18
  data$y19<- data$temp*data$I19
  data$y110<- data$temp*data$I110
  
  data$y21<- data$food*data$I21
  data$y22<- data$food*data$I22
  data$y23<- data$food*data$I23
  data$y24<- data$food*data$I24
  data$y25<- data$food*data$I25
  data$y26<- data$food*data$I26
  data$y27<- data$food*data$I27
  data$y28<- data$food*data$I28
  data$y29<- data$food*data$I29
  data$y210<- data$food*data$I210
  
  data$p21<- data$temp2*data$I21
  data$p22<- data$temp2*data$I22
  data$p23<- data$temp2*data$I23
  data$p24<- data$temp2*data$I24
  data$p25<- data$temp2*data$I25
  data$p26<- data$temp2*data$I26
  data$p27<- data$temp2*data$I27
  data$p28<- data$temp2*data$I28
  data$p29<- data$temp2*data$I29
  data$p210<- data$temp2*data$I210
  
  data$l11<- data$temp2*data$I11
  data$l12<- data$temp2*data$I12
  data$l13<- data$temp2*data$I13
  data$l14<- data$temp2*data$I14
  data$l15<- data$temp2*data$I15
  data$l16<- data$temp2*data$I16
  data$l17<- data$temp2*data$I17
  data$l18<- data$temp2*data$I18
  data$l19<- data$temp2*data$I19
  data$l110<- data$temp2*data$I110
  
  data$z1<- data$food*data$N
  data$z2<- data$temp*data$N
  data$z3<- data$temp2*data$N
  
  
  #### Model fitting and selection
  rbf.gfr<- glm(data$use~data$food+data$temp+data$I11+data$I12+data$I13+data$I14+data$I15+data$I16+data$I17+data$I18+data$I19+data$I110+
                            data$I21+data$I22+data$I23+data$I24+data$I25+data$I26+data$I27+data$I28+data$I29+data$I210+
                            data$x11+data$x12+data$x13+data$x14+data$x15+data$x16+data$x17+data$x18+data$x19+data$x110+
                            data$x21+data$x22+data$x23+data$x24+data$x25+data$x26+data$x27+data$x28+data$x29+data$x210+
                            data$y11+data$y12+data$y13+data$y14+data$y15+data$y16+data$y17+data$y18+data$y19+data$y110+
                            data$temp2+data$N+data$z1+data$z2+data$z3+
                            data$p21+data$p22+data$p23+data$p24+data$p25+data$p26+data$p27+data$p28+data$p29+data$p210+
                            data$l11+data$l12+data$l13+data$l14+data$l15+data$l16+data$l17+data$l18+data$l19+data$l110+
                            data$y21+data$y22+data$y23+data$y24+data$y25+data$y26+data$y27+data$y28+data$y29+data$y210,
                          family = poisson)
  

  #######################################################################################
  # Part 2:  Prediction using the trained RBF-GFR model
  #######################################################################################
  data<- as.data.frame(test)
  data<-data[order(data$id),]
  
  
  
  # Generate a Gaussian mixture approximation of the habitat composition of different id.
  
  dat<- subset.data.frame( data,select=c(food,temp))
  g<-9
  set.seed(123)
  fa<-favail(dat, blocking=data$id, G=g)
  
  ##### Extract the mean in Mixture Gaussian approximation  
  # Extracts unique values of blocking index
  jv<-unique(data$id) 
  jvmax<-length(jv)
  #extracts the names of all the formula terms
  vs<-colnames(dat)
  comp<- 1:g
  # Extracts the main effects
  mainVars<-length(vs)
  jvmax<- 1:length(jv)
  mainVars<- 1:length(vs)
  
  ### mean of the mixture Gaussian distribution of the first variable for all block and each component.
  m1<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      m1[i,k]<- fa[[i]]$mus[1,k]
      
    }
  }
  ### mean of the second variable for all block and each component 
  m2<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      m2[i,k]<- fa[[i]]$mus[2,k]
      
    }
  }
  
  #### Variance of the mixture Gaussian approximation
  ### variance of the first variable "food" of all block
  c1<- matrix(0,nrow =length(jvmax),ncol = 1 )
  
  for(i in 1:length(jvmax)){
    
    c1[i]<- fa[[i]]$var[1]
    c1
  }
  ### variance of all block for second variable 
  c2<- matrix(0,nrow =length(jvmax),ncol = 1 )
  
  for(i in 1:length(jvmax)){
    
    c2[i]<- fa[[i]]$var[2]
    c2
  }
  
  #### Matrix of weight of the mixture Gaussian approximation  
  ### weight for all component for each block
  
  w<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      w[i,k]<- fa[[i]]$ws[k]
    }
  }
  
  
  ### I for 'food' the first variable,j=1, and basis function m=1
  
  i11<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i11[i,k]<-w[i,k]*(f.sd1/sqrt(c1[i]+f.sd1^2))*exp(-(m1[i,k]-x[1])^2/(2*(c1[i]+f.sd1^2)))
      
      i11[i,k]
    }
  }
  s11<- as.matrix(rowSums(i11)) 
  ### I for 'food' the first variable,j=1, and basis function m=2
  i12<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i12[i,k]<-w[i,k]*(f.sd2/sqrt(c1[i]+f.sd2^2))*exp(-(m1[i,k]-x[2])^2/(2*(c1[i]+f.sd2^2)))
      i12[i,k]
    }
  }
  s12<- as.matrix(rowSums(i12)) 
  ### I for 'food' the first variable,j=1, and basis function m=3
  i13<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i13[i,k]<-w[i,k]*(f.sd3/sqrt(c1[i]+f.sd3^2))*exp(-(m1[i,k]-x[3])^2/(2*(c1[i]+f.sd3^2)))
      i13[i,k]
    }
  }
  s13<- as.matrix(rowSums(i13)) 
  ### I for 'food' the first variable,j=1, and basis function m=4
  i14<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i14[i,k]<-w[i,k]*(f.sd4/sqrt(c1[i]+f.sd4^2))*exp(-(m1[i,k]-x[4])^2/(2*(c1[i]+f.sd4^2)))
      i14[i,k]
    }
  }
  s14<- as.matrix(rowSums(i14)) 
  ### I for 'food' the first variable,j=1, and basis function m=5
  i15<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i15[i,k]<-w[i,k]*(f.sd5/sqrt(c1[i]+f.sd5^2))*exp(-(m1[i,k]-x[5])^2/(2*(c1[i]+f.sd5^2)))
      i15[i,k]
    }
  }
  s15<- as.matrix(rowSums(i15)) 
  ### I for 'food' the first variable,j=1, and basis function m=6
  i16<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i16[i,k]<-w[i,k]*(f.sd6/sqrt(c1[i]+f.sd6^2))*exp(-(m1[i,k]-x[6])^2/(2*(c1[i]+f.sd6^2)))
      i16[i,k]
    }
  }
  s16<- as.matrix(rowSums(i16)) 
  
  
  ### I for 'food' the first variable,j=1, and basis function m=7
  i17<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i17[i,k]<-w[i,k]*(f.sd7/sqrt(c1[i]+f.sd7^2))*exp(-(m1[i,k]-x[7])^2/(2*(c1[i]+f.sd7^2)))
      i17[i,k]
    }
  }
  s17<- as.matrix(rowSums(i17)) 
  ### I for 'food' the first variable,j=1, and basis function m=8
  i18<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i18[i,k]<-w[i,k]*(f.sd8/sqrt(c1[i]+f.sd8^2))*exp(-(m1[i,k]-x[8])^2/(2*(c1[i]+f.sd8^2)))
      i18[i,k]
    }
  }
  s18<- as.matrix(rowSums(i18)) 
  ### I for 'food' the first variable,j=1, and basis function m=9
  i19<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i19[i,k]<-w[i,k]*(f.sd9/sqrt(c1[i]+f.sd9^2))*exp(-(m1[i,k]-x[9])^2/(2*(c1[i]+f.sd9^2)))
      i19[i,k]
    }
  }
  s19<- as.matrix(rowSums(i19)) 
  ### I for 'food' the first variable,j=1, and basis function m=10
  i110<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i110[i,k]<-w[i,k]*(f.sd10/sqrt(c1[i]+f.sd10^2))*exp(-(m1[i,k]-x[10])^2/(2*(c1[i]+f.sd10^2)))
      i110[i,k]
    }
  }
  s110<- as.matrix(rowSums(i110)) 
  
  
  ### I for 'temp' the second variable,j=1, and basis function m=1
  i21<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i21[i,k]<-w[i,k]*t.sd1/sqrt(c2[i]+t.sd1^2)*exp(-(m2[i,k]-y[1])^2/(2*(c2[i]+t.sd1^2)))
      i21[i,k]
    }
  }
  s21<- as.matrix(rowSums(i21)) 
  ### I for 'temp' the second variable,j=2, and basis function m=2
  i22<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i22[i,k]<-w[i,k]*t.sd2/sqrt(c2[i]+t.sd2^2)*exp(-(m2[i,k]-y[2])^2/(2*(c2[i]+t.sd2^2)))
      i22[i,k]
    }
  }
  s22<- as.matrix(rowSums(i22)) 
  ### I for 'temp' the second variable,j=2, and basis function m=3
  i23<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i23[i,k]<-w[i,k]*t.sd3/sqrt(c2[i]+t.sd3^2)*exp(-(m2[i,k]-y[3])^2/(2*(c2[i]+t.sd3^2)))
      i23[i,k]
    }
  }
  s23<- as.matrix(rowSums(i23)) 
  ### I for 'temp' the second variable,j=2, and basis function m=4
  i24<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i24[i,k]<-w[i,k]*t.sd4/sqrt(c2[i]+t.sd4^2)*exp(-(m2[i,k]-y[4])^2/(2*(c2[i]+t.sd4^2)))
      i24[i,k]
    }
  }
  s24<- as.matrix(rowSums(i24)) 
  ### I for 'temp' the second variable,j=2, and basis function m=5
  i25<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i25[i,k]<-w[i,k]*t.sd5/sqrt(c2[i]+t.sd5^2)*exp(-(m2[i,k]-y[5])^2/(2*(c2[i]+t.sd5^2)))
      i25[i,k]
    }
  }
  s25<- as.matrix(rowSums(i25)) 
  ### I for 'temp' the second variable,j=2, and basis function m=6
  i26<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i26[i,k]<-w[i,k]*t.sd6/sqrt(c2[i]+t.sd6^2)*exp(-(m2[i,k]-y[6])^2/(2*(c2[i]+t.sd6^2)))
      i26[i,k]
    }
  }
  s26<- as.matrix(rowSums(i26)) 
  
  ### I for 'temp' the second variable,j=2, and basis function m=7
  i27<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i27[i,k]<-w[i,k]*t.sd7/sqrt(c2[i]+t.sd7^2)*exp(-(m2[i,k]-y[7])^2/(2*(c2[i]+t.sd7^2)))
      i27[i,k]
    }
  }
  s27<- as.matrix(rowSums(i27)) 
  ### I for 'temp' the second variable,j=2, and basis function m=8
  i28<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i28[i,k]<-w[i,k]*t.sd8/sqrt(c2[i]+t.sd8^2)*exp(-(m2[i,k]-y[8])^2/(2*(c2[i]+t.sd8^2)))
      i28[i,k]
    }
  }
  s28<- as.matrix(rowSums(i28)) 
  
  ### I for 'temp' the second variable,j=2, and basis function m=9
  i29<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i29[i,k]<-w[i,k]*t.sd9/sqrt(c2[i]+t.sd9^2)*exp(-(m2[i,k]-y[9])^2/(2*(c2[i]+t.sd9^2)))
      i29[i,k]
    }
  }
  s29<- as.matrix(rowSums(i29))  
  ### I for 'temp' the second variable,j=2, and basis function m=10
  i210<- matrix(0,nrow =length(jvmax),ncol = length(comp) )
  for(i in 1:length(jvmax)){
    for(k in 1:length(comp)){
      i210[i,k]<-w[i,k]*t.sd10/sqrt(c2[i]+t.sd10^2)*exp(-(m2[i,k]-y[10])^2/(2*(c2[i]+t.sd10^2)))
      i210[i,k]
    }
  }
  s210<- as.matrix(rowSums(i210)) 
  
  ### Creat columns I11,I12, ....
  data$I11<- c(rep(c(s11[1:length(s11)]), each=500))
  data$I12<- c(rep(c(s12[1:length(s11)]), each=500))
  data$I13<- c(rep(c(s13[1:length(s11)]), each=500))
  data$I14<- c(rep(c(s14[1:length(s11)]), each=500))
  data$I15<- c(rep(c(s15[1:length(s11)]), each=500))
  data$I16<- c(rep(c(s16[1:length(s11)]), each=500))
  data$I17<- c(rep(c(s17[1:length(s11)]), each=500))
  data$I18<- c(rep(c(s18[1:length(s11)]), each=500))
  data$I19<- c(rep(c(s19[1:length(s11)]), each=500))
  data$I110<- c(rep(c(s110[1:length(s11)]), each=500))
  
  
  data$I21<- c(rep(c(s21[1:length(s11)]), each=500))
  data$I22<- c(rep(c(s22[1:length(s11)]), each=500))
  data$I23<- c(rep(c(s23[1:length(s11)]), each=500))
  data$I24<- c(rep(c(s24[1:length(s11)]), each=500))
  data$I25<- c(rep(c(s25[1:length(s11)]), each=500))
  data$I26<- c(rep(c(s26[1:length(s11)]), each=500))
  data$I27<- c(rep(c(s27[1:length(s11)]), each=500))
  data$I28<- c(rep(c(s28[1:length(s11)]), each=500))
  data$I29<- c(rep(c(s29[1:length(s11)]), each=500))
  data$I210<- c(rep(c(s210[1:length(s11)]), each=500))
  data$temp2<- data$temp^2
  
  #Interaction Terms
  data$x11<- data$food*data$I11
  data$x12<- data$food*data$I12
  data$x13<- data$food*data$I13
  data$x14<- data$food*data$I14
  data$x15<- data$food*data$I15
  data$x16<- data$food*data$I16
  data$x17<- data$food*data$I17
  data$x18<- data$food*data$I18
  data$x19<- data$food*data$I19
  data$x110<- data$food*data$I110
  
  data$x21<- data$temp*data$I21
  data$x22<- data$temp*data$I22
  data$x23<- data$temp*data$I23
  data$x24<- data$temp*data$I24
  data$x25<- data$temp*data$I25
  data$x26<- data$temp*data$I26
  data$x27<- data$temp*data$I27
  data$x28<- data$temp*data$I28
  data$x29<- data$temp*data$I29
  data$x210<- data$temp*data$I210
  
  data$y11<- data$temp*data$I11
  data$y12<- data$temp*data$I12
  data$y13<- data$temp*data$I13
  data$y14<- data$temp*data$I14
  data$y15<- data$temp*data$I15
  data$y16<- data$temp*data$I16
  data$y17<- data$temp*data$I17
  data$y18<- data$temp*data$I18
  data$y19<- data$temp*data$I19
  data$y110<- data$temp*data$I110
  
  data$y21<- data$food*data$I21
  data$y22<- data$food*data$I22
  data$y23<- data$food*data$I23
  data$y24<- data$food*data$I24
  data$y25<- data$food*data$I25
  data$y26<- data$food*data$I26
  data$y27<- data$food*data$I27
  data$y28<- data$food*data$I28
  data$y29<- data$food*data$I29
  data$y210<- data$food*data$I210
  
  data$p21<- data$temp2*data$I21
  data$p22<- data$temp2*data$I22
  data$p23<- data$temp2*data$I23
  data$p24<- data$temp2*data$I24
  data$p25<- data$temp2*data$I25
  data$p26<- data$temp2*data$I26
  data$p27<- data$temp2*data$I27
  data$p28<- data$temp2*data$I28
  data$p29<- data$temp2*data$I29
  data$p210<- data$temp2*data$I210
  
  data$l11<- data$temp2*data$I11
  data$l12<- data$temp2*data$I12
  data$l13<- data$temp2*data$I13
  data$l14<- data$temp2*data$I14
  data$l15<- data$temp2*data$I15
  data$l16<- data$temp2*data$I16
  data$l17<- data$temp2*data$I17
  data$l18<- data$temp2*data$I18
  data$l19<- data$temp2*data$I19
  data$l110<- data$temp2*data$I110
  
  data$z1<- data$food*data$N
  data$z2<- data$temp*data$N
  data$z3<- data$temp2*data$N
  
  
  # Processing the output of the predictions 
  predictions<-predict.glm(rbf.gfr,type = 'response',test)
  predictions<- pmax(predictions, 0.001)
  
  # calculating the out-of-sample score
  1 - (sum((test$use - predictions)^2)/sum((test$use - mean(test$use))^2))
  
