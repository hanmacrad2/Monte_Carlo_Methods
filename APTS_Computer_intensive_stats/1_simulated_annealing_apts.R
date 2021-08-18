#APTS
library(MASS)


#Question 1: Density
get_fx = function(X){
  
  x_1 = X[1]
  x_2 = X[2]
  fx = (4 -2.1*x_1^2 + (x_1^4)/3)*x_1^2 + x_1*x_2 + 4*(x_2^2 -1)*x_2^2
  
  fx
}

#Apply to get fx
X <- matrix(0, nrow=2, ncol= 1)#
X[1,] = -3
X[2,] = -2
fx = get_fx(X)
fx

#fx no2
simulate_fx = function(len_x = 100, X){
  
  #Matrix
  count = 1 
  #fx = vector('numeric', length = len_x^2)
  fx = 0 
  while (count <= len_x) {
    for (x_1 in X[1,]) {
      for (x_2 in X[2,]){
        #fx[1, count] = x_1
        #fx[2, count] = x_2
        fx[count] = (4 -2.1*x_1^2 + (x_1^4)/3)*x_1^2 + x_1*x_2 + 4*(x_2^2 -1)*x_2^2
        
        count = count + 1
        
      }
    } 
  }
  
  
  fx
}

#apply to get fx
len_x = 100
X1 = seq(-3,3,length=len_x) 
X2 = seq(-2,2,length=len_x)
X <- matrix(0, nrow=2, ncol= len_x)#
X[1,] = X1
X[2,] = X2
fx = get_fx(X)

#Function v2
get_fx_matrix = function(len_x = 100){
  
  #Matrix
  X1 = seq(-3,3,length=len_x) 
  X2 = seq(-2,2,length=len_x) 
  fx = matrix(0, nrow = 3, ncol = len_x^2)
  count = 1 
  
  while (count <= len_x) {
    for (x_1 in X1) {
      for (x_2 in X2){
        fx[1, count] = x_1
        fx[2, count] = x_2
        fx[3, count] = (4 -2.1*x_1^2 + (x_1^4)/3)*x_1^2 + x_1*x_2 + 4*(x_2^2 -1)*x_2^2
        
        count = count + 1
        
      }
    } 
  }

  
  fx
}

#apply to get fx
fx_matrix = get_fx_matrix()
#Plot
scatterplot3js(fx[1,], fx[2,], fx[3,], phi = 40, theta = 20,
               color=rainbow(length(z)),
               colkey = FALSE,
               cex = .3,
               main = "f(x1, x2)")

#fx - Mode point
fx_max = max(fx[3,])
fx_max
x1_max = X1[which.max(fx[3,])]
x1_max
x2_max = X2[which.max(fx[3,])]
x2_max

#fx - Min point
fx_min = min(fx_matrix[3,])
fx_min
x1_min = fx_matrix[1,][which.min(fx_matrix[3,])]
x1_min
x2_min = fx_matrix[2,][which.min(fx_matrix[3,])]
x2_min



#************************************************************************************
#Simulated Annelaing Bivariate function
simulated_annealing_multi = function(N = 2000, x0=as.vector(c(0,0)), Sigma = diag(1,nrow =2), beta1 = 1, alpha = 1.001){
  
  X <- matrix(0, nrow=2, ncol= N)
  X[,1] <- x0
  x_current = X[,1]
  fx_current = get_fx(X)
  beta_current = beta1
  fn = 0
  count_na = 0 
  
  #Implement Markov Chain in loop 
  for (i in 1: N){
    
    x_new = x_current + mvrnorm(mu=as.vector(c(0,0)), Sigma=Sigma)
    print(x_new)
    fx_new = get_fx(x_new)
    #Alpha
    alpha_accept = (fx_new/fx_current)^beta_current #alpha -acceptance probability
    
    #Check
    if (is.na(alpha_accept)){
      count_na = count_na + 1
    }

    else if ((alpha_accept > runif(1)) & abs(x_new[1]) <= 3 & abs(x_new[2]) <= 2){ #Criterion for acceptance -
      x_current <- x_new #If holds set next sample in the chain as the proposed sample
      fx_current <- fx_new
    }
    X[i] = x_current #If not the next sample in the chain is again set to as the sample value at previous step
    fn[i] = fx_current
    beta_current = beta_current*alpha #Update temperature/value of beta. 
    
    
  }
  print('count_na')
  print(count_na)
  
  list(X, fn)
}

#Simulated annealing to find the minimum of a multi paramater function
simulated_annealing_minimum_multi_param = function(N = 10000, x0=as.vector(c(0,0)), Sigma = diag(1,nrow =2), beta1 = 1, alpha = 1.001){
  
  X <- matrix(0, nrow=2, ncol= N)
  X[,1] <- x0
  x_current = X[,1]
  fx_current = get_fx(X)
  beta_current = beta1
  fn = 0
  count_accept = 0
  count_na = 0 
  
  #Implement Markov Chain in loop 
  for (i in 1: N){
    
    x_new = x_current + mvrnorm(mu=as.vector(c(0,0)), Sigma=Sigma)
    #print(x_new)
    fx_new = get_fx(x_new)
    #Alpha
    #alpha_accept = (fx_new/fx_current)^beta_current #alpha -acceptance probability
    alpha_accept = exp(-beta_current*(fx_new - fx_current))#alpha -acceptance probability
    
    #Check
    if (is.na(alpha_accept)){
      count_na = count_na + 1
    }
    
    else if ((alpha_accept > runif(1)) & abs(x_new[1]) <= 3 & abs(x_new[2]) <= 2){ #Criterion for acceptance -
      x_current <- x_new #If holds set next sample in the chain as the proposed sample
      fx_current <- fx_new
      count_accept = count_accept + 1
    }
    X[i] = x_current #If not the next sample in the chain is again set to as the sample value at previous step
    fn[i] = fx_current
    beta_current = beta_current*alpha #Update temperature/value of beta. 
    
    
  }
  #Print rates 
  rate_accept = count_accept/N
  print('count accept')
  print(count_accept)
  print('Rate accept')
  print(rate_accept)
  
  print('count_na')
  print(count_na)
  
  list(X, fn, rate_accept)
}


#***************
#Implement
N = 10000
vec_sm = simulated_annealing_minimum_multi_param(N)
head(vec_sm[1])
#Plot trace of X
plot.ts(vec_sm[[1]][1,], 
        xlab = 't', ylab = 'X1', col = 'green',main = 'Trace Plot of X1', cex.lab=1.5)
plot.ts(vec_sm[[1]][2,], 
        xlab = 't', ylab = 'X2', col = 'red', main = 'Trace Plot of X2', cex.lab=1.5)

#Plot Trace of fx
plot.ts(vec_sm[[2]], 
        xlab = 't', ylab = 'f(Xt)', main = 'Trace Plot of f(Xt)', cex.lab=1.5, col = 'blue')

#Plot cumulative mean of x1, x2, f
par(mfrow=c(2,1))
x1_mean = cumsum(vec_sm[[1]][1,])/seq_along(vec_sm[[1]][1,])
plot(seq_along(vec_sm[[1]][1,]), x1_mean, xlab = 'Time', ylab = 'Mean of X1', col = 'green', main = paste("Mean of X1, sd of proposal = ", Sigma))

x2_mean = cumsum(vec_sm[[1]][2,])/seq_along(vec_sm[[1]][2,])
plot(seq_along(vec_sm[[1]][2,]), x2_mean, xlab = 'Time', ylab = 'Mean of X2', col = 'red', main = paste("Mean of X2, sd of proposal = ", Sigma))


Sigma = 1
fx_mean = cumsum(vec_sm[[2]])/seq_along(vec_sm[[2]])
plot(seq_along(vec_sm[[2]]), fx_mean, xlab = 'Time', ylab = 'R0', col = 'red', main = paste("Mean of MCMC chain, True R0 = ",r0_true, ", sd of proposal = ", Sigma))

#*******************************************************************
#*Investigate Tuning parameters
list_alpha = c(1.001, 1.01, 1.1, 2, 10)
list_sigma = c(0.001,0.01, 0.1,1,2, 5, 10,100)

#Function to inspect
results_plot_ts_rate_acc = function(list_sigma, list_alpha, N){
  #Params
  vec_mean_rate_acc = vector('numeric', length(list_sigma))
  par(mfrow=c(3,5))
  for (sigma_j in list_sigma){
    print('Sigma')
    print(sigma_j)
    rate_accept_sigma = 0
    for (alpha_i in list_alpha){
      print('Alpha')
      print(alpha_i)
      vec_sm = simulated_annealing_minimum_multi_param(N = N, x0=as.vector(c(0,0)), Sigma = diag(sigma_j, nrow =2), beta1 = 1, alpha = alpha_i)
      plot.ts(vec_sm[[2]],  xlab = 't', ylab = 'f(Xt)',
              cex.lab=1.5, col = 'orange')
      rate_accept_sigma = rate_accept_sigma + vec_sm[[3]]
    }
    vec_mean_rate_acc[sigma_j] = rate_accept_sigma/length(list_alpha)
  }
  vec_mean_rate_acc
}

#Apply
N = 100000
list_alpha = c(1.001, 1.01, 1.1, 2, 10)
list_sigma1 = c(1, 2, 5) #Probably be enough c(0.001, 0.01, 0.1) #c(10, 100) #c(1, 2, 5) #c(0.001, 0.01, 0.1) #c(1, 2, 5) # #c(1, 10, 100) 
rate_acc = results_plot_ts_rate_acc(list_sigma1, list_alpha, N)


#****************************************************************************************************************
#*Part II 
#*

get_fxII = function(X){
  
  x_1 = X[1]
  x_2 = X[2]
  fx = (4 -2.1*x_1^2 + (x_1^4)/3)*x_1^2 + x_1*x_2 + 4*(x_2^2 -1)*x_2^2
  
  fx
}