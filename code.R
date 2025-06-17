## Gauss Lasso Function:
## lasso.fit() is called by fitlasso().  X.model is the full model matrix without an intercept column.
## This function uses a large search grid for lambda. 

lasso.fit<-function(X, Y){
  require(broom)
  require(glmnet)
  n = dim(X)[1]
  Yc=Y-mean(Y)
  Xc<-apply(X, 2, function(X) X - mean(X))
  vlength<-sqrt((1/n)*(diag(t(Xc)%*%Xc)))
  Xs<-t(t(Xc)*(1/(vlength)))
  log_lam = seq(from = -8, to=round(log(max(abs((t(Xs)%*%Yc)))),digits=1), by = 0.1)  
  lambda = exp(log_lam)
  fit<-glmnet(Xs, Yc, standardize = FALSE, intercept = FALSE, 
              lambda=lambda)
  return(fit)
}

fitlasso<-function(X.model, Y){
  fit <- lasso.fit(X.model,Y)
  gamma<-0.1*max(abs(fit$beta)) # data driven threshold
  betat <- fit$beta
  betat[abs(betat) < gamma] <- 0
  xm<-as.data.frame(X)
  xm<-cbind(Y,xm)
  SSE<-c()
  p<-c()
  
  for (j in 1:length(fit$lambda)) {
    coef<-betat[,j]
    sub<-which(abs(coef)>0)
    coef<-coef[sub]
    
    if (length(coef)==0){ 
      fit2<-lm(Y~1, xm)
      SSE[j]  <- deviance(fit2)
      p[j]    <- 1 
      
    } else {
      
      Eqn<-gsub(":","*",names(coef))
      Eqn <- paste("I(",Eqn,")")
      Eqn <- paste(Eqn,collapse="+")
      Eqn <- as.formula(paste("Y~",Eqn))
      fit2<-lm(Eqn,data=xm)
      SSE[j] <-deviance(fit2)
      p[j]   <-1+length(coef)
      
    }
  }
  
  info=n*log(SSE/n)+p*log(n) #BIC
  lambda.min <- fit$lambda[which.min(info)]
  minIndex <- which.min(info)
  lasso_nz <- betat[,minIndex] # non-zero lasso estimates
  activeL<-which(abs(lasso_nz)>0)
  return(activeL)
}

## Regression Best Subsets adapted from Chapter 6 of James et al. (2016)
## predict.regsubets() provides predted values for a regsubsets() object.

predict.regsubsets =function (object ,newdata ,id){
  form=as.formula (object$call [[2]])
  mat=model.matrix (form ,newdata )
  coefi =coef(object ,id=id)
  xvars =names (coefi )
  mat[,xvars ]%*% coefi
}

## Best subsets regression with CV or LOOCV
## Change k to n for LOOCV

library(leaps)
k=5
cv.errors =matrix (NA ,k,  params, dimnames =list(NULL, paste(1:params)))
for(j in 1:k){
  best.fit =regsubsets (Y~.,data=data[cvindex[[j]],],nvmax = params)
  
  for(i in 1:(dim(summary(best.fit)$which)[1])) { #dim(summary...) extracts the nvmax from best.fit when it is reduced.
    pred=predict(best.fit, data[-cvindex[[j]],], id=i)
    cv.errors[j,i]=sqrt(mean(((data$Y[-cvindex[[j]]])-pred)^2))
  }
}
p<-which.min(apply(cv.errors ,2, mean, na.rm=TRUE))[[1]]
reg.best<-regsubsets(Y~.,data=data,nvmax =p)#params)
best.fit_summary <- summary(reg.best)
model_id <- p
model_vars <- best.fit_summary$outmat[model_id, ]
variable_vector <- names(model_vars)[model_vars == "*"]
variable_vector <- gsub("\\.", ":", variable_vector)
if (variable_vector[1]=="(Intercept)"){
  variable_vector <- variable_vector[-1]
} else {
  variable_vector=variable_vector
}


## Best Subsets Regression with Little Bootstrap
## Adapted from Breiman (1996)
## X is the model matrix



library(leaps)

little_bootstrap_best_subsets_96 <- function(X, y, t, n_bootstrap) {
  
  n <- length(y)
  p <- ncol(X)
  
  
  full_model <- lm(y ~ ., data = as.data.frame(X))
  full_model_summary<-summary(full_model)
  sigma<-full_model_summary$sigma
  subsets_model <- regsubsets(y ~ ., data = as.data.frame(X), nvmax = p)
  subsets_summary <- summary(subsets_model)
  RSS_s<-numeric()
  
  for (i in 1:dim(subsets_summary$which)[1]){
    selected<-names(which(subsets_summary$which[i, -1]))
    selected<-gsub("\\`", "", selected)
    Xsub<-data.frame(X[,selected])
    Xsub$y<-y
    fitted_values<-lm(y~., data=data.frame(Xsub))$fitted.values
    RSS_s[i]<-sum((y-fitted_values)^2)
  }
  
  
  B<- matrix(0, nrow = n_bootstrap, ncol = dim(subsets_summary$which)[1])
  
  for (b in 1:n_bootstrap) {
    
    
    noise <- rnorm(n, mean = 0, sd = (t * sigma))
    y_perturbed <- y + noise
    
    
    for (k in 1:(dim(subsets_summary$which)[1])) {
      selected<-names(which(subsets_summary$which[k, -1]))
      selected<-gsub("\\`", "", selected)
      Xsub<-data.frame(X[,selected])
      Xsub$y<-y_perturbed
      fitted_values<-lm(y~., data=data.frame(Xsub))$fitted.values
      mu_hat<-fitted_values
      bias <-sum(noise*(mu_hat))*(1/t^2)
      B[b,k]<-bias
      
      
    }
    
  }
  
  

  Avg_B<-colMeans(B)
  ME=RSS_s+(2*Avg_B) 
  optimal_model_size <- which.min(ME)
  final_model <- regsubsets(y ~ ., data = as.data.frame(X), nvmax = optimal_model_size)
  
  return(list(optimal_model_size = optimal_model_size, final_model = final_model)) #optimal model size includes the
}



## Best Subsets Regression with Little Bootstrap with estimate of error from ridge regression
## X is the model matrix

library(devtools)
install_github("xliusufe/RidgeVar")
library(RidgeVar)

little_bootstrap_best_subset_screen <- function(X, y, t = 0.6, n_bootstrap = 25) {
  
  n <- length(y)
  p <- ncol(X)
  X_df<-as.data.frame(X)
  
  
  
  fit<-VAR_RR(y, X_df, eta=NULL, alpha=0.1)
  sigma<-sqrt(fit$sigma2)
  subsets_model <- regsubsets(y ~ ., data = as.data.frame(X), nvmax = n-1)
  subsets_summary <- summary(subsets_model)
  #RSS_s<-subsets_model$rss[1:dim(subsets_summary$which)[1]] #this is sus
  
  RSS_s<-numeric() 
  for (i in 1:dim(subsets_summary$which)[1]){
    selected<-names(which(subsets_summary$which[i, -1]))
    selected<-gsub("\\`", "", selected)
    Xsub<-data.frame(X[,selected])
    Xsub$y<-y
    fitted_values<-lm(y~., data=data.frame(Xsub))$fitted.values
    RSS_s[i]<-sum((y-fitted_values)^2)
  }
  
  B<- matrix(0, nrow = n_bootstrap, ncol =dim(subsets_summary$which)[1])
  
  
  for (b in 1:n_bootstrap) {
    
    # Perturb response
    noise <- rnorm(n, mean = 0, sd = (t * sigma))
    y_perturbed <- y + noise
    
    
    
    
    for (k in 1:(dim(subsets_summary$which)[1])) {
      selected<-names(which(subsets_summary$which[k, -1]))
      selected<-gsub("\\`", "", selected)
      Xsub<-data.frame(X[,selected])
      Xsub$y<-y_perturbed
      fitted_values<-lm(y~., data=data.frame(Xsub))$fitted.values
      mu_hat<-fitted_values
      bias <-sum(noise*(mu_hat))*(1/t^2)
      B[b,k]<-bias
      
      
    }
    
  }
  
  
  
  Avg_B<-colMeans(B)
  ME=RSS_s+(2*Avg_B) #1996 formula
  optimal_model_size <- which.min(ME)
  final_model <- regsubsets(y ~ ., data = as.data.frame(X), nvmax = optimal_model_size)
  
  return(list(optimal_model_size = optimal_model_size, final_model = final_model))
}
