#Monte Carlo Assignment II
library(BBmisc)

par(mar=c(1,1,1,1))
par(oma=c(2,2,2,2)) #:D 
#op<-par(no.readonly=TRUE)

#Variables - Global
N = 250000
list_sd = list(1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4, 153.6, 204.8) #76.8, list(0.2, 0.4, 0.8, 1.6..)
#list_sd2 = list(0.2, 0.6, 1.6, 3.2, 9.6, 25.6, 51.2, 102.4, 153.6) #76.8, list(0.2, 0.4, 0.8, 1.6..)

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

#fx - Mode point
fx_max = max(fx)
x_max = x[which.max(fx)]

#Plot
par(mfrow=c(1,1)) 
plot(x, fx, xlab='x', ylab='f(x)', main='Density of f(x)', 
     cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='blue') #xlim = range(0,1)
#Plot Mode point
points(x_max, fx_max, col = 'red', lwd = 2.5)

hist(fx, col = 'lightskyblue1')

#Plot Mean
fx_mean = cumsum(fx)/seq_along(fx) 
plot(seq_along(fx_mean), fx_mean,
     xlab='x', ylab='f(x)', main='Ergodic mean of f(x)', 
     cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='blue')


#Normalised
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
scale_fx = range01(fx)
plot(x, scale_fx, type = 'l')

#Hist
hist(fx, prob = TRUE, col = 'lightskyblue1')

#Scaled
newMax = 0.05
newMin = 0.0 
range02 <- function(x){ (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin }
scale_fx = range02(fx)

#Hist
hist(fx, prob = TRUE, col = 'lightcyan2')
lines(x, scale_fx, type = 'l', col = 'red')

#******************************************************
#Part 2: MCMC


#*******************************************************
#Part 1: Random Walk Metropolis

rwm = function(N, sigma, x0 = 2){
  
  #Vectors to store samples
  X = vector('numeric', N)
  X[1] = x0 #Initialise 1st sample in chain in order to start the chain running 
  #count_accept = N #Count of number of accepted values 
  #Loop
  for (t in 2:N){
    X[t] = X[t-1] + rnorm(1, mean = 0, sd= sigma) #Symmetric random walk proposal
    u = runif(1)
    
    if(u > (get_fx(X[t])/get_fx(X[t-1]))) { #Criterion for rejection 
      X[t] = X[t-1] #If criterion for acceptance is not met the next sample in chain is set to current sample
      #count_accept = count_accept - 1
    }
  }
  X #list(X, count_accept) #Return samples and count of accepted values 
}

#Apply
X_rw = rwm(N, sigma = 6.4, x0 = 2)
#i. Trace Plots Plot trace
plot.ts(X_rw)
#Hist samples
hist(X_rw)

#Best
#Time it
N = 5000
start_time = Sys.time()
X_rw = rwm(N, sigma = 51.2, x0 = 2)
end_time = Sys.time()
time_elapsed = end_time - start_time
time_elapsed


#Scaled
newMax = 0.038
newMin = 0.0 
range02 <- function(x){ (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin }
scale_fx = range02(fx)


#Hist
#
#par(mfrow = c(2,3))
#X_rw = rwm(N, sigma = 51.2, x0 = 2)
hist(X_rw, 50, main = '', prob = TRUE, col = 'lightcyan2') #400 --> 100
lines(x, scale_fx, type = 'l', lwd=1.5, col = 'red')

#**************************************************
#Tuning Params + Diagnosis: Apply MCMC - varying sigma 

#i. Trace plot
plot_ts_samp = function(list_sd, N_samp){
  par(mfrow=c(2,3))
  for (sd_i in list_sd){
    print(sd_i)
    X = rwm(N = N_samp, sigma = sd_i, x0 = 2)
    plot.ts(X,  xlab = 't', ylab = 'Xt',
            cex.lab=1.5, col = 'black') # main = 'Trace Plot',, main = substitute(paste('Trace Plot, sigma = ', sd_i), sd_i))
  }
  
}

#Apply: sd
N = 100000
list_sd = list(1.6, 3.2, 6.4, 12.8, 25.6, 51.2) #, 102.4, 153.6) #list_sd = list(1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4, 153.6, 204.8)
plot_ts_samp(list_sd, N)

#Plot small sigma
plot_ts_samp(c(102.4, 153.6), N)

#*************************
#ii. Ergodic Mean

erg_mean_samp = function(list_sd, N_samp){
  #Plot
  par(mfrow=c(2,3))
  for (sd_i in list_sd){
    print(sd_i)
    X = rwm(N = N_samp, sigma = sd_i, x0 = 2)
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
erg_mean_samp(list_sd, N)
erg_mean_samp(c(51.2), 50000)
#erg_mean_samp(c(1.6, 3.2, 6.4),N)

#**************
#Hist of samples 
hist_samps = function(x, scale_fx, list_sd, N_samp){
  #Plot
  par(mfrow=c(2,3))
  for (sd_i in list_sd){
    print(sd_i)
    Xs = rwm(N = N_samp, sigma = sd_i, x0 = 2)
    #hist(Xs, 400, main = '', prob = TRUE)
    hist(Xs, 50, main = '', prob = TRUE, col = 'lightcyan2') #400 --> 100
    lines(x, scale_fx, type = 'l', lwd=1.5, col = 'red')
    #lines(x, fx, xlab='x', ylab='f(x)', 
         #cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='red') #xlim = range(0,1)
    
  }
  
}

hist_samps(list_sd, N)
hist_samps(c(102.4, 153.6), N)
plot(x, fx, xlab='x', ylab='f(x)', 
cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='red') #xlim = range(0,1)

#Diagnostics
list_sd = list(25.6, 51.2, 102.4)
  
#*************************
#iii. Autocorrelation
get_acf1 = function(list_sd, N_samp){
  list_acf1 <- c() #vector("list", length(list_sd))
  count = 1
  for (sd_i in list_sd){
    X = rwm(N = N_samp, sigma = sd_i, x0 = 2)
    list_acf1[count] = acf(X, plot = FALSE)[1][1]
    count = count + 1
    
  }
  list_acf1
} 


list_acf = get_acf1(list_sd, N)  
list_acf

#*************************
#iv. ASJD

#ii. ASJD
asjd = function(X1){
  T2 = length(X1)
  asjd_x = sum((X1[2:T2] - X1[1:T2-1])^2)/(T2-1)
  asjd_x
}

apply_asjd = function(list_sd, N_samp) {
  list_asjd <- c() #vector("list", length(list_sd))
  count = 1
  for (sd_i in list_sd){
    X = rwm(N = N_samp, sigma = sd_i, x0 = 2) 
    list_asjd[count] = asjd(X)
    count = count + 1
  }
  list_asjd
}

list_asjd = apply_asjd(list_sd, N)  
list_asjd


#*************************
#iv. Riemann Sum

riemann_sum = function(X1){
  T = length(X1)
  Xs = sort(X1[1:T])
  rsum = sum( (Xs[2:T]-Xs[1:T-1])*get_fx(Xs[2:T]) )
}

apply_rsum = function(list_sd, N_samp) {
  list_rsum <- c() #vector("list", length(list_sd))
  count = 1
  for (sd_i in list_sd){
    X1 = rwm(N = N_samp, sigma = sd_i, x0 = 2) 
    list_rsum[count] = riemann_sum(X1)
    count = count + 1
  }
  list_rsum
}

list_rsum = apply_rsum(list_sd, N)  
list_rsum

#********************************
#v. Acceptance Rate
rwm_ar = function(N, sigma, x0 = 2){
  
  #Vectors to store samples
  X = vector('numeric', N)
  X[1] = x0 #Initialise 1st sample in chain in order to start the chain running 
  count_accept = N #Count of number of accepted values 
  #Loop
  for (t in 2:N){
    X[t] = X[t-1] + rnorm(1, mean = 0, sd= sigma) #Symmetric random walk proposal
    u = runif(1)
    
    if(u > (get_fx(X[t])/get_fx(X[t-1]))) { #Criterion for rejection 
      X[t] = X[t-1] #If criterion for acceptance is not met the next sample in chain is set to current sample
      count_accept = count_accept - 1
    }
  }
  list(X, count_accept) #Return samples and count of accepted values 
}

#Apply 
apply_ar = function(list_sd, N_samp) {
  list_ar <- c() #vector("list", length(list_sd))
  count = 1
  for (sd_i in list_sd){
    vec_rw = rwm_ar(N = N_samp, sigma = sd_i, x0 = 2) 
    list_ar[count] = vec_rw[[2]]/N_samp
    count = count + 1
  }
  list_ar
}

list_sd2 = list(1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4, 153.6)
list_ar = apply_ar(list_sd2, N)  
list_ar









#*************************
#Draft

#i. Trace plot Apply: Plot
plot_samp = function(list_sd){
  par(mfrow=c(3,3))
  i_row = 1
  j_col = 1
  for (sd_i in list_sd){
    print(sd_i)
    X = rwm(N = 5000, sigma = sd_i, x0 = 2)
    #par(mfrow=c(i_row,j_col)) ##Plot all trace plots together row column
    plot.ts(X) #, main = substitute(paste('Trace Plot, sigma = ', sd_i), sd_i))
    #i_row = i_row + 1
    j_col = j_col + 1
    if(j_col == 3){
      j_col = 1
      i_row = i_row + 1
    }
  }
  
}