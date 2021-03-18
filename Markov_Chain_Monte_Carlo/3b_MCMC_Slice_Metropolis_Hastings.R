#6. Metropolis Slice Sampler :D boom boom boom!

par(mfrow = c(2,2))
slice_samp_metropolis = function(N, sigma, x0 = 2){
  
  #x- markov chain
  x = vector('numeric', N)
  y = vector('numeric', N)
  x[1] = x0
  count_accept = N
  
  #Markov Chain starts
  for (t in 2:N){
    
    #Auxilliary Uniform variable u
    y[t] = runif(1, 0, get_fx(x[t-1]))
    
    #X sampled from markov chain using a symmetric random 
    x[t] = x[t-1] + rnorm(1, mean = 0, sd= sigma)
    
    #Using augmented distribution f(x,u) 
    if(y[t] > get_fx(x[t])) {          #Rejection criterion
      x[t] = x[t-1]
      count_accept = count_accept - 1 
    }
  }
  print(count_accept/N) #Acceptance rate
  x
}

Xs_mh = slice_samp_metropolis(N = 100000, sigma = 5, x0 = 1)
hist(Xs_mh, 100)

plot.ts(Xs_mh, 100) #Note x mix aswell
hist(Xs_mh)

#Time it
N = 5000
start_time = Sys.time()
X_ss = slice_samp_metropolis(N, sigma = 51.2, x0 = 2) 
end_time = Sys.time()
time_elapsed = end_time - start_time
time_elapsed

#Scaled
newMax = 0.04
newMin = 0.0 
range02 <- function(x){ (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin }
scale_fx = range02(fx)

#Best Hist
#par(mfrow = c(2,3))
#X_rw = rwm(N, sigma = 51.2, x0 = 2)
hist(X_ss, 50, main = '', prob = TRUE, col = 'lightcyan2') #400 --> 100
lines(x, scale_fx, type = 'l', lwd=1.5, col = 'red')

#SAME FOR ALL COPY AND PASTE BELOW!

#**********************************************
#Tuning Parameters

#i. Trace plot
plot_ts_samp_ss = function(list_sd, N_samp){
  par(mfrow=c(2,3))
  for (sd_i in list_sd){
    print(sd_i)
    X = slice_samp_metropolis(N = N_samp, sigma = sd_i, x0 = 2)
    plot.ts(X,  xlab = 't', ylab = 'Xt',
            cex.lab=1.5, col = 'black') # main = 'Trace Plot',, main = substitute(paste('Trace Plot, sigma = ', sd_i), sd_i))
  }
  
}

#Apply: sd
N = 100000
list_sd = list(1.6, 3.2, 6.4, 12.8, 25.6, 51.2) #, 102.4, 153.6) #list_sd = list(1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4, 153.6, 204.8)
plot_ts_samp_ss(list_sd, N)
plot_ts_samp_ss(c(102.4, 153.6), N)

#***********************************************
#Diagnostics
list_sd = list(25.6, 51.2, 102.4)
  
ss_get_acf1 = function(list_sd, N_samp){
  list_acf1 <- c() #vector("list", length(list_sd))
  count = 1
  for (sd_i in list_sd){
    X = slice_samp_metropolis(N = N_samp, sigma = sd_i, x0 = 2)
    list_acf1[count] = acf(X, plot = FALSE)[1][1]
    count = count + 1
    
  }
  list_acf1
} 


list_acf = ss_get_acf1(list_sd, N)  
list_acf

#Acceptance Rate
ar_slice_samp_metropolis = function(N, sigma, x0 = 2){
  
  #x- markov chain
  x = vector('numeric', N)
  y = vector('numeric', N)
  x[1] = x0
  count_accept = Nfa
  
  #Loop
  for (t in 2:N){
    
    #y Uniform
    y[t] = runif(1, 0, get_fx(x[t-1]))
    
    #X sampled from markov chain 
    x[t] = x[t-1] + rnorm(1, mean = 0, sd= sigma)
    
    #Using augmented distribution f(x,y) equivalence **
    if(y[t] > get_fx(x[t])) {
      x[t] = x[t-1]
      count_accept = count_accept - 1
    }
  }
  #print('Sigma')
  #print(sigma)
  #print('accept rate')
  count_accept/N
  
}

#Apply 
ss_apply_ar = function(list_sd, N_samp) {
  list_ar <- c() #vector("list", length(list_sd))
  count = 1
  for (sd_i in list_sd){

    list_ar[count] = ar_slice_samp_metropolis(N_samp, sigma = sd_i, x0 = 2)
    count = count + 1
  }
  list_ar
}

#Apply
N = 100000

list_ar = ss_apply_ar(list_sd, 50000)  
list_ar


#ASJD
ss_apply_asjd = function(list_sd, N_samp) {
  list_asjd <- c() #vector("list", length(list_sd))
  count = 1
  for (sd_i in list_sd){
    X = slice_samp_metropolis(N = N_samp, sigma = sd_i, x0 = 2) 
    list_asjd[count] = asjd(X)
    count = count + 1
  }
  list_asjd
}

list_asjd = ss_apply_asjd(list_sd, N)  
list_asjd



#***********************************************
#Assesment
#Ergodic Mean

ss_erg_mean_samp = function(list_sd, N_samp){
  #Plot
  #par(mfrow=c(2,3))
  for (sd_i in list_sd){
    print(sd_i)
    #X = slice_samp_metropolis(N = N_samp, sigma = sd_i, x0 = 2)
    X = HMC_samples(N_HM, epsilon = 0.5, L = 25, q0 = 0)
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

N = 250000
ss_erg_mean_samp(c(51.2), 50000)

#**************
#Hist of samples 
ss_hist_samps = function(list_sd, N_samp){
  #Plot
  par(mfrow=c(2,3))
  for (sd_i in list_sd){
    print(sd_i)
    Xs = slice_samp_metropolis(N = N_samp, sigma = sd_i, x0 = 2)
    hist(Xs, 400, main = '', prob = TRUE)
    #lines(x, fx, xlab='x', ylab='f(x)', 
    #cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='red') #xlim = range(0,1)
    
  }
  
}

ss_hist_samps(list_sd, N)




#Modification

mod_slice_samp_metropolis = function(N, sigma, x0 = 2){
  
  #x- markov chain
  x = vector('numeric', N)
  y = vector('numeric', N)
  x[1] = x0
  count_accept = N
  count_na = 0
  
  #Markov Chain starts
  for (t in 2:N){
    
    #Auxilliary Uniform variable u
    y[t] = runif(1, 0, get_fx(x[t-1]))
    
    #X sampled from markov chain using a symmetric random 
    x[t] = x[t-1] + rnorm(1, mean = 0, sd= sigma)
    
    alpha_accept_new = get_fx(X[t])/(get_fx(X[t]) +get_fx(X[t-1]))
    if (is.na(alpha_accept_new)){ #is.na
      count_na = count_na + 1
    } else if(y[t] > alpha_accept_new) {     #Using augmented distribution f(x,u)     #Rejection criterion
      x[t] = x[t-1]
      count_accept = count_accept - 1 
    }
  }
  print(count_accept/N) #Acceptance rate
  x
}

Xs_mh = mod_slice_samp_metropolis(N = 10000, sigma = 51.2, x0 = 1)

#Hist
hist(Xs_mh, 50, main = '', prob = TRUE, col = 'lightcyan2') #400 --> 100
lines(x, scale_fx, type = 'l', lwd=1.5, col = 'red')
