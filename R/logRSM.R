# based on regRSM package (http://cran.r-project.org/web/packages/regRSM/)

compute_initial_weights = function(y,x){
# The function returns initial weights.

  initial_weights = as.numeric(cor(y,x))^2
  initial_weights = initial_weights/(sum(initial_weights))
  return(initial_weights)
}

is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
# The function checks if the given number is whole number.

  abs(x - round(x)) < tol
}

compute_scores = function(y,x,m,B,initial_weights=NULL){
# The function returns RSM scores.

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
      control=list(selval=NULL,screening=NULL,init_weights=FALSE,m=NULL,B=NULL))

  attr(logRSM,"class")="logRSM"
  return(logRSM)
}

logRSM = function(y,x,yval=NULL,xval=NULL,m=NULL,B=NULL,
    store_data=FALSE,screening=NULL,init_weights=FALSE,thrs=NULL,penalty=NULL,initial_weights=NULL)
{
  if (init_weights) {
    if (!is.null(initial_weights)) {
      stop('init_weights cannot be TRUE if initial_weigths are provided')
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
  #Check for initial_weights
  if(init_weights){
    initial_weights = compute_initial_weights(y,x)
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
  if(store_data) { logRSM$input_data$x=data_x; logRSM$input_data$y=y }

  logRSM$control$selval = selval
  logRSM$control$screening = screening
  logRSM$control$init_weights =  init_weights
  logRSM$control$m = m
  logRSM$control$B = B

  return(logRSM)
}

calculate_initial_weigths_with_t = function(class, data) {
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

