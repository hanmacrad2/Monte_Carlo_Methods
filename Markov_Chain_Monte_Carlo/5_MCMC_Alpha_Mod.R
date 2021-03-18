#5_MCMC Modification

#Alpha

logR = sum(log(F_sample(y.prop))) -
  sum(log(F_sample(y[t-1])))  

x_new = random_walk

alpha_new = fx(x_new)/( fx(x_new) + fx(x[i-1]))

#Function
#Part 1: Random Walk

rwm_alpha_new = function(N, sigma, x0 = 2){
  
  #Params
  count_accept = 0
  X = vector('numeric', N)
  X[1] = x0
  
  #Markov chain
  for (t in 2:N){
    X[t] = X[t-1] + rnorm(1, mean = 0, sd= sigma) 
    u = runif(1)
    alpha_new = get_fx(X[t])/(get_fx(X[t]) +get_fx(X[t-1]))
    if(u > alpha_new) {
      X[t] = X[t-1]
      count_accept = count_accept + 1
    }
  }
  #print('count_accept:')
  #print(count_accept)
  X
}

#Apply
N2 = 50000
X_rw_an = rwm_alpha_new(N2, sigma = 6.4, x0 = 2)
#i. Trace Plots Plot trace
plot.ts(X_rw_an)
#Hist samples
hist(X_rw_an)  



#**********************************************
#Tuning Parameters
list_sd = list(1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4, 153,6)
list_sd = list(1.6, 3.2, 6.4, 12.8, 25.6, 51.2)


#i. Trace plot
plot_ts_samp_rwan = function(list_sd, N_samp){
  par(mfrow=c(2,3))
  for (sd_i in list_sd){
    print(sd_i)
    X = rwm_alpha_new(N = N_samp, sigma = sd_i, x0 = 2)
    plot.ts(X,  xlab = 't', ylab = 'Xt',
            cex.lab=2.5, col = 'blue') # main = 'Trace Plot',, main = substitute(paste('Trace Plot, sigma = ', sd_i), sd_i))
  }
  
}

#Apply: sd
N = 100000
N = 5000
list_sd = list(1.6, 3.2, 6.4, 12.8, 25.6, 51.2) #, 102.4, 153.6) #list_sd = list(1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4, 153.6, 204.8)
plot_ts_samp_rwan(list_sd, N)
plot_ts_samp_ss(c(102.4, 153.6), N)

#***********************************************
#Diagnostics
#list_sd = list(25.6, 51.2, 102.4)
N = 100000
list_sd = list(1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4, 153.6)

rwan_get_acf1 = function(list_sd, N_samp){
  list_acf1 <- c() #vector("list", length(list_sd))
  count = 1
  for (sd_i in list_sd){
    X = rwm_alpha_new(N = N_samp, sigma = sd_i, x0 = 2)
    list_acf1[count] = acf(X, plot = FALSE)[1][1]
    count = count + 1
    
  }
  list_acf1
} 


list_acf = rwan_get_acf1(list_sd, N)  
list_acf

#***************************
#Acceptance Rate
ar_rwm_alpha_new = function(N, sigma, x0 = 2){
  
  #Params
  count_accept = 0
  X = vector('numeric', N)
  X[1] = x0
  count_accept = N 
  
  #Loop
  for (t in 2:N){
    X[t] = X[t-1] + rnorm(1, mean = 0, sd= sigma)
    u = runif(1)
    alpha_new = get_fx(X[t])/(get_fx(X[t]) +get_fx(X[t-1]))
    if(u > alpha_new) {
      X[t] = X[t-1]
      count_accept = count_accept - 1
    }
  }
  count_accept/N
}


#Apply 
rwan_apply_ar = function(list_sd, N_samp) {
  list_ar <- c() #vector("list", length(list_sd))
  count = 1
  for (sd_i in list_sd){
    X 
    list_ar[count] = ar_rwm_alpha_new(N_samp, sigma = sd_i, x0 = 2)
    count = count + 1
  }
  list_ar
}

#Apply
rwan_apply_ar(list_sd, N)  



#ASJD
rn_apply_asjd = function(list_sd, N_samp) {
  list_asjd <- c() #vector("list", length(list_sd))
  count = 1
  for (sd_i in list_sd){
    X = rwm_alpha_new(N = N_samp, sigma = sd_i, x0 = 2) 
    list_asjd[count] = asjd(X)
    count = count + 1
  }
  list_asjd
}

list_asjd = rn_apply_asjd(list_sd, N)  
list_asjd



#***********************************************
#Assesment
#Ergodic Mean

rwn_erg_mean_samp = function(list_sd, N_samp){
  #Plot
  par(mfrow=c(2,3))
  for (sd_i in list_sd){
    print(sd_i)
    X = rwm_alpha_new(N = N_samp, sigma = sd_i, x0 = 2)
    x_mean = cumsum(X)/seq_along(X)
    print('Mean:')
    print(mean(x_mean[(length(x_mean)-1000):length(x_mean)]))
    print(mean(x_mean[(length(x_mean)-10):length(x_mean)]))
    print('...')
    plot(seq_along(x_mean), x_mean,
         xlab = 't', ylab = 'Xt', #main = 'Ergodic Mean',
         type = 'l', col = 'black',lwd = '2', cex = 1.5)
    #mtext('Ergodic Mean', outer = TRUE, )
    
  }
  
}

N = 5000
#Original
erg_mean_samp(c(25.6, 51.2, 102.4), N)

#New 
rwn_erg_mean_samp(c(25.6, 51.2, 102.4), N)



#**************
#Hist of samples 
rn_hist_samps = function(x, scale_fx, list_sd, N_samp){
  #Plot
  par(mfrow=c(2,3))
  for (sd_i in list_sd){
    print(sd_i)
    Xs = rwm_alpha_new(N = N_samp, sigma = sd_i, x0 = 2)
    hist(Xs, 50, main = '', prob = TRUE, col = 'lightcyan2') #400 --> 100
    lines(x, scale_fx, type = 'l', lwd=1.5, col = 'red')
    
    #lines(x, fx, xlab='x', ylab='f(x)', 
    #cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='red') #xlim = range(0,1)
    
  }
  
}

#Compare
N = 5000

#Original
hist_samps(c(25.6, 51.2, 102.4), N)
lines(x, range02(fx), type = 'l')

#New
rn_hist_samps(x, scale_fx, c(25.6, 51.2, 102.4), N)



#commented Code
#*******************************************************
#Part 1: Random Walk Metropolis

rwm_alpha_new = function(N, sigma, x0 = 2){
  
  #Vectors to store samples
  X = vector('numeric', N)
  X[1] = x0 #Initialise 1st sample in chain in order to start the chain running 
  count_accept = N #Count of number of accepted values 
  #Loop
  for (t in 2:N){
    X[t] = X[t-1] + rnorm(1, mean = 0, sd= sigma) #Symmetric random walk proposal
    u = runif(1)
    #Modified alpha
    alpha_new = get_fx(X[t])/(get_fx(X[t]) +get_fx(X[t-1]))
    
    if(u > alpha_new = get_fx(X[t])/(get_fx(X[t]) +get_fx(X[t-1]))) { #Criterion for rejection 
      X[t] = X[t-1] #If criterion for acceptance is not met the next sample in chain is set to current sample
      count_accept = count_accept - 1
    }
  }
  X #list(X, count_accept) #Return samples and count of accepted values 
}