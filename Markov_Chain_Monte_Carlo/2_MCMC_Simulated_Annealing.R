#MC Assign II Part b: Simulated Annealing
par(mar=c(1,1,1,1))
par(oma=c(2,2,2,2)) #:D 

#Question 1: Density
get_fx = function(x){
  
  #Student number
  sn = c(1,9,8,3,1,6,2) #Other check sn = c(1,9,9,0,6,2,7)
  #Density
  a = sn[7]*exp(-sin((sn[1]*x^2)/(15-sn[1])) - ((x-3-(sn[2]*pi))^2)/(2*(5+sn[3])^2) )
  b = 2*(1+sn[7])*exp(-((x^2)/32))
  c = (10-sn[7])*exp( -cos((sn[4]*x^2)/(15+sn[4])) - ((x+3+(sn[5]*pi))^2)/(2*(5+sn[6])^2))
  fx = a+b+c
  fx
}

#fx - apply
x = seq(-50,80,length=1000) 
fx = get_fx(x)
length(fx)
norm_fx = normalize(fx, range = c(0,1))

#Global params#
N = 100000 #0
sd_sim_an = c(0.001, 0.01, 0.1, 1, 10, 100)
list_alpha = c(1.001, 1.01, 1.1, 2, 10)

#Function
sim_anneal_1 = function(N, x0 = 0, beta1 = 1, alpha = 1.001, sigma = 1){
  
  #Containers
  x = vector('numeric', N)
  fn = numeric(N)
  count_na = 0 
  
  #Variables current
  x_current = x0
  fx_current = get_fx(x_current)
  beta_current = beta1
  
  for (i in 1: N){
    x_new = x_current + rnorm(2)*sigma
    fx_new = get_fx(x_new)
    alpha_accept = -beta_current*(fx_new - fx_current)
    if (is.na(alpha_accept)){
      count_na = count_na + 1
    }
    else if (alpha_accept > log(runif(1))){
      x_current <- x_new
      fx_current <- fx_new
    }
    x[i] = x_current
    fn[i] = fx_current
    beta_current = beta_current*alpha
    
    
  }
  
  list(x, fn)
}

#Sim_anneal
sim_anneal = function(N, x0 = 0, beta1 = 1, alpha = 1.001, sigma = 1){
  
  #Vectors for storage
  x = vector('numeric', N)
  fn = numeric(N)
  count_na = 0 
  
  #Variables - current values
  x_current = x0
  fx_current = get_fx(x_current)
  beta_current = beta1
  
  #Implement Markov Chain in loop 
  for (i in 1: N){
    x_new = x_current + rnorm(1)*sigma #Normal random walk proposal  
    fx_new = get_fx(x_new)
    alpha_accept = (fx_new/fx_current)^beta_current #alpha -acceptance probability
    if (is.na(alpha_accept)){
      count_na = count_na + 1
    }
    else if (alpha_accept > runif(1)){ #Criterion for acceptance -
      x_current <- x_new #If holds set next sample in the chain as the proposed sample
      fx_current <- fx_new
    }
    x[i] = x_current #If not the next sample in the chain is again set to as the sample value at previous step
    fn[i] = fx_current
    beta_current = beta_current*alpha #Update temperature/value of beta. 
    
    
  }
  print('count_na')
  print(count_na)
  
  list(x, fn)
}

#***************
#Apply
vec_sm = sim_anneal(N2, sigma = 10)
head(vec_sm[1])
#Plot trace of X
plot.ts(vec_sm[[1]], 
        xlab = 't', ylab = 'Xt', main = 'Trace Plot of Xt', cex.lab=1.5)
#Plot Trace of fx
plot.ts(vec_sm[[2]], 
        xlab = 't', ylab = 'f(Xt)', main = 'Trace Plot of f(Xt)', cex.lab=1.5)

#************************************************************************************
#Tuning params
N = 5000
#1. Sigma - 0.1. All alpha
plot_fx_trace = function(list_sigma, list_alpha, N){
  par(mfrow=c(3,5))
  for (sigma_j in list_sigma){
    print('Sigma')
    print(sigma_j)
    for (alpha_i in list_alpha){
      print(alpha_i)
      vec_sm = sim_anneal(N = N, alpha = alpha_i, sigma = sigma_j)
      plot.ts(vec_sm[[2]],  xlab = 't', ylab = 'f(Xt)',
              cex.lab=1.5, col = 'purple') # main = 'Trace Plot',, main = substitute(paste('Trace Plot, sigma = ', sd_i), sd_i))
    }
  }
}

#Apply 
list_sigma1 =  c(1, 10, 100) #c(0.001, 0.01, 0.1) #c(1, 10, 100) #c(0.001, 0.01, 0.1)
plot_fx_trace(list_sigma1, list_alpha, N)




#***********************************************************************************
#Log form :D 
sim_anneal_log = function(N, x0 = 0, beta1 = 1, alpha = 1.001, sigma = 1){
  #Containers
  x = vector('numeric', N)
  fn = numeric(N)
  count_na = 0 
  count_accept = 0 
  
  #Variables current
  x_current = x0
  fx_current = get_fx(x_current)
  beta_current = beta1
  
  for (i in 1: N){
    x_new = x_current + rnorm(1)*sigma #Before: wrong: x_new = x_current + rnorm(2)*sigma
    fx_new = get_fx(x_new)
    alpha_accept = beta_current*(log(fx_new)-log(fx_current)) 
    if (is.na(alpha_accept)){
      count_na = count_na + 1
    }
    else if (alpha_accept > log(runif(1))){
      x_current <- x_new
      fx_current <- fx_new
      count_accept = count_accept + 1
    }
    x[i] = x_current
    fn[i] = fx_current
    beta_current = beta_current*alpha
    
    
  }
  #print('count accept')
  #print(count_accept)
  rate_accept = count_accept/N
  print('Rate accept')
  print(rate_accept)
  
  
  list(x, fn, rate_accept)
}


#Apply
vec_sm = sim_anneal_log(N2, sigma = 10)
#Plot trace of X
plot.ts(vec_sm[[1]], 
        xlab = 't', ylab = 'Xt', main = 'Trace Plot of Xt', cex.lab=1.5)
#Plot Trace of fx
plot.ts(vec_sm[[2]], 
        xlab = 't', ylab = 'f(Xt)', main = 'Trace Plot of f(Xt)', cex.lab=1.5)

#Diagnositcs
#sim_anneal_log = function(N, x0 = 0, beta1 = 1, alpha = 1.001, sigma = 1)

for (alpha_i in list_alpha){
  print(alpha_i)
  vec_sm = sim_anneal_log(N, x0 = 0, beta1 = 1, alpha = 1.001, sigma = 1)
  #print('Autocorr')
  #print(acf( vec_sm[[2]], plot = FALSE)[1][1])
  print('ASJD')
  print(asjd(vec_sm[[2]])) #, plot = FALSE)[1][1])
}

#Diagnostics SA a new
sim_anneal_an

for (alpha_i in list_alpha){
  print(alpha_i)
  vec_sm = sim_anneal_an(N, x0 = 0, beta1 = 1, alpha = 1.001, sigma = 1)
  print('ASJD')
  print(asjd(vec_sm[[2]]))
  print('Autocorr')
  print(acf( vec_sm[[2]], plot = FALSE)[1][1])

}



#************************************************************************************
#Tuning params

#1. Sigma - 0.1. All alpha
results_plot_ts_rate_acc = function(list_sigma, list_alpha, N){
  #Params
  vec_mean_rate_acc = vector('numeric', length(list_sigma))
  par(mfrow=c(3,5))
  for (sigma_j in list_sigma){
    print('Sigma')
    print(sigma_j)
    rate_accept_sigma = 0
    for (alpha_i in list_alpha){
      print(alpha_i)
      vec_sm = sim_anneal_log(N = N, alpha = alpha_i, sigma = sigma_j)
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
list_sigma1 = c(10, 100) #c(1, 2, 5) #c(0.001, 0.01, 0.1) #c(1, 2, 5) # #c(1, 10, 100) 
rate_acc = results_plot_ts_rate_acc(list_sigma1, list_alpha, N)














#************************************************************************************
#Modification

sim_anneal_an = function(N, x0 = 0, beta1 = 1, alpha = 1.001, sigma = 1){
  #Containers
  x = vector('numeric', N)
  fn = numeric(N)
  count_na = 0
  count_accept = 0 
  
  #Variables current
  x_current = x0
  fx_current = get_fx(x_current)
  beta_current = beta1
  
  for (i in 1: N){
    x_new = x_current + rnorm(1)*sigma
    fx_new = get_fx(x_new)
    
    #****
    alpha_accept_new = (fx_new^beta_current)/(fx_new^beta_current + fx_current^beta_current)
    
    #alpha_accept_new = beta_current*log(fx_new^beta_current)/(beta_current*log(fx_new) + beta_current*log(fx_current))
    
    #Try b*log(fx)
    
    if (is.na(alpha_accept_new)){ #is.na
      count_na = count_na + 1
    }
    else if (alpha_accept_new > runif(1)){ #log(runif)
      x_current <- x_new
      fx_current <- fx_new
      count_accept = count_accept + 1
    }
    x[i] = x_current
    fn[i] = fx_current
    beta_current = beta_current*alpha
    
    
  }
  print('count_acccept') #print('count_na')
  print(count_accept/N) #print(count_na)
  
  list(x, fn)
}

#***************
#Apply
vec_sm_an = sim_anneal_an(N2, x0 = 1, sigma = 20)
#Check

#Plot trace of X
plot.ts(vec_sm_an[[1]], 
        xlab = 't', ylab = 'Xt', main = 'Trace Plot of Xt', cex.lab=1.5)
#Plot Trace of fx
plot.ts(vec_sm_an[[2]], 
        xlab = 't', ylab = 'f(Xt)', main = 'Trace Plot of f(Xt)', cex.lab=1.5)

