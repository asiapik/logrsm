# based on regRSM package (http://cran.r-project.org/web/packages/regRSM/)

#' Checks if the given number is whole number.
#' 
#' @param x number to be checked
#' @param tol precision
is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}

#' Compute RSM scores.
compute_scores = function(y, x, m, B, initial_weights = NULL){

  p = ncol(x)
  scores = numeric(p)
  ns = numeric(p)

  for(k in 1:B){
    submodel = sample(1:p,size=m,replace=FALSE,prob=initial_weights)
    lm1 = glm(y~x[,submodel], family = binomial)
    weights = as.numeric((summary(lm1)$coef[-1,3])^2)
    scores[submodel] =  scores[submodel] + weights
    ns[submodel] = ns[submodel] + 1
  }
  ns = ifelse(ns!=0,ns,1)
  scores = scores/ns

  return(scores)
}

new.logRSM <- function()
{

  logRSM=list(scores=NULL,model=NULL,time=list(user=0,system=0,elapsed=0),
      data_transfer=list(user=0,system=0,elapsed=0),
      coefficients=NULL, predError=NULL,input_data=list(x=NULL,y=NULL),
      control=list(selval=NULL,screening=NULL,m=NULL,B=NULL))

  attr(logRSM,"class")="logRSM"
  return(logRSM)
}

#' Main function - calculate logRSM
#'
#' @param x
#' @param y
#' @param m
#' @param B
#' @param store_data
#' @param initial_weights
#' @return 
#' 
#' \code{initial_weights} parameter could be an array with weights or
#' a function, that accepts two parameters: \code{class} and \code{data}
#' and returns weigths.
#' Two functions \link{calculate_weigths_with_cor} and \link{calculate_weigths_with_t}
#' are provided.
logRSM = function(y, x, yval = NULL, xval = NULL, m = NULL, B = NULL,
    store_data = FALSE, screening = NULL, initial_weights = NULL)
{
  
  # checsk, if initial weights is a function
  if (!is.null(initial_weights)) {
    # if it's a functtion - call it
    if (is.function(initial_weights)) {
      initial_weights = initial_weights(y, x)
    }
  }
  
  data_x = x;
  x = as.matrix(x)
  y = as.numeric(y)
  n = length(y)
  p = ncol(x)
  scores = NULL

  startTime <- proc.time()

  # Set default values of m and B
  if(is.null(m)){
    m = floor(min(n-1,p)/2)
  }else{
    if(m>(n-2)) stop("Parameter m cannot be larger than the number of observations minus two!")
    if(m<=0) stop("Parameter m must be a positive number!")
    if(!is.wholenumber(m)) stop("Parameter m must be a whole number!")
  }
  if(is.null(B)){
    B = 1000
  }else{
    if(B<=0) stop("Parameter B must be a positive number!")
    if(!is.wholenumber(B)) stop("Parameter B must be a whole number!")
  }

  #Check for screeneing
  if(!is.null(screening))
  {
    if((screening>=1)||(screening<=0)) stop("screening must be in (0,1)")

    iw =  compute_initial_weights(y,x)
    sel = which(iw>=quantile(iw,screening))
    if(m>length(sel)) stop('Parameter m cannot be larger than the number of attributes remaining after screening procedure!')
    x = x[,sel]
  }

  #RSM method esence
  d1=d2=proc.time()
  scores = compute_scores(y,x,m,B,initial_weights)

  #Set score 0, when variable is not selected by screeneing
  if(!is.null(screening)){
    scores1 = numeric(ncol(data_x))
    scores1[sel] = scores
    scores = scores1
  }

  selval = ifelse(!is.null(yval) && !is.null(xval),TRUE,FALSE)
  if(selval==TRUE){
    order1 = sort(scores,decreasing=TRUE,index.return=TRUE)$ix
    selected_model = select_finalmodel_qr(y,data_x,yval,xval,order1)
    model = selected_model$model
    coefficients =  as.numeric(selected_model$coefficients)
    predError = selected_model$predError
    informationCriterion = NULL
  }else{
    model = NULL
    coefficients =  NULL
    predError = NULL
    informationCriterion = NULL
  }

  stopTime <- proc.time()

  logRSM = new.logRSM()
  logRSM$scores = scores
  logRSM$model = model
  logRSM$time = stopTime-startTime
  logRSM$coefficients = coefficients
  logRSM$predError = predError
  logRSM$informationCriterion = informationCriterion
  logRSM$data_transfer = d2-d1
  if (store_data) {
    logRSM$input_data$x = data_x;
    logRSM$input_data$y = y 
  }

  logRSM$control$selval = selval
  logRSM$control$screening = screening
  logRSM$control$initial_weights =  initial_weights
  logRSM$control$m = m
  logRSM$control$B = B

  return(logRSM)
}

#' calculate weights using T test.
calculate_weigths_with_t = function(class, data) {
  result = numeric(ncol(data)); 
  for (i in 1:ncol(data)) {
    x = numeric();
    y = numeric();
    for (j in 1:length(class)) {
      if (class[j] == 0) {
        x[length(x) + 1] <- data[j, i];
      } else {
        y[length(y) + 1] <- data[j, i];
      }
    }
    tResult = t.test(x, y, paired = FALSE, alternative = 'two.sided', var.equal = TRUE)$statistic;
    result[i] = abs(tResult);
  }
  result = result / sum(result)
  return (result)
}

#' calculate weights using correlations.
calculate_weigths_with_cor = function(class, data) {
  initial_weights = as.numeric(cor(class, data))^2
  initial_weights = initial_weights/(sum(initial_weights))
  return(initial_weights)
}
