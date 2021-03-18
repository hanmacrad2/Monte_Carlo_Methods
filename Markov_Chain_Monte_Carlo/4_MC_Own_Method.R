#4_MC_Own_Method: 
#Note
#Is it -log of the density or exp(density) i.e the probability?? 
library(rootSolve)

#
#epsilon - stepsize for leapfrog steps
#L - the number of leapfrog steps in the trajectory, L

#Function (- log)
get_fx_log = function(x){
  
  #Student number
  sn = c(1,9,8,3,1,6,2) #Other check sn = c(1,9,9,0,6,2,7)
  #Density
  a = sn[7]*exp(-sin((sn[1]*x^2)/(15-sn[1])) - ((x-3-(sn[2]*pi))^2)/(2*(5+sn[3])^2) )
  b = 2*(1+sn[7])*exp(-((x^2)/32))
  c = (10-sn[7])*exp( -cos((sn[4]*x^2)/(15+sn[4])) - ((x+3+(sn[5]*pi))^2)/(2*(5+sn[6])^2))
  fx = a+b+c
  -log(fx)
}

#Gradient
x = seq(-50,80,length=100) 
#x = 0.5
fx_log = get_fx_log(x)
grad_U <- gradient(f = get_fx_log, x)


#****************************************************************
#Samples
#x- markov chain

HMC_samples = function(N, epsilon, L, q0){
  
  #Initialise Parameters
  q_samp = vector('numeric', N)
  current_q = q0
  count_na = 0 
  count_accept = 0
  
  #Start Markov Chain
  for (t in 1:N){
    
    #Current position 
    q = current_q
    #Momentum
    p = rnorm(length(q),0,1) # independent standard normal variates
    #Current momentum
    current_p = p
    
    #Leap frog method
    #U_gradient -> Gradient of Potential energy - Negative log probability of the function 
    grad_U_q = gradient(f = get_fx_log, q)
    
    #A half step for momentum is made at the beginning
    p = p - epsilon * grad_U_q / 2
    
    #Alternate full steps for the position and momentum
    for (i in 1:L){
      
      #A full step for the position is made
      q = q + epsilon * p
      #U_gradient
      grad_U_q = gradient(f = get_fx_log, q)
      
      #A full step for the momentum is made, except at end of trajectory
      if (i!=L){
        p = p - epsilon*grad_U_q
      } 
    }
    
    #U_gradient
    grad_U_q = gradient(f = get_fx_log, q)
    #A half step for momentum at the end.
    p = p - epsilon*grad_U_q / 2
    
    #Momentum is negated at end of trajectory in order to make the proposal symmetric
    p = -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_U = get_fx_log(current_q)# I.e U(current_q)
    current_K = sum(current_p^2) / 2
    proposed_U = get_fx_log(q)
    proposed_K = sum(p^2) / 2
    
    #Accept/Rejection step 
    # Accept or reject the state at end of trajectory.
    #This returns either the position at the end of the trajectory or the initial position
    
    prob_accept = current_U-proposed_U+current_K-proposed_K #Acceptance probability 
    if (is.na(prob_accept)){
      count_na = count_na + 1
    } else if (runif(1) < exp(prob_accept)){
      current_q = q
      count_accept = count_accept + 1
    }
    q_samp[t] = current_q #Store current_q as sample
    
  }
  #print('Accept rate:')
  #print(count_accept/N)
  
  q_samp #Return samples 
}

#Apply
N_HM = 100000
epsilon = 0.3
#current_q = 1 #Or Random np.random.randn(2)
#X_hmc = HMC(fx_log, grad_U, epsilon, length(x), current_q)
q0 = 0
Xs_hmc = HMC_samples(N_HM, epsilon, L = 10, q0)

#Apply
hist(Xs_hmc, prob = TRUE, 400)
lines(x, scale_fx, col = 'red') # type = 'l')

#Trace plot
plot.ts(Xs_hmc)

#Time it
N = 5000
start_time = Sys.time()
X_hmc  = HMC_samples(N, epsilon = 0.5, L = 25, q0 = 0)
end_time = Sys.time()
time_elapsed = end_time - start_time
time_elapsed

#Best
#Scaled
newMax = 0.04
newMin = 0.0 
range02 <- function(x){ (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin }
scale_fx = range02(fx)

#Best Hist
par(mfrow = c(2,3))
#X_rw = rwm(N, sigma = 51.2, x0 = 2)
hist(X_hmc, 50, main = '', prob = TRUE, col = 'lightcyan2') #400 --> 100
lines(x, scale_fx, type = 'l', lwd=1.5, col = 'red')

#**********************************************
#TUNING RESULTS - ALL COMBOS

hamiltonian_results = function(N_HM, q0, list_eps, list_L, x, scale_fx){
  #Plot
  par(mfcol = c(3, 3))
  
  #Params
  count = 0
  list_acf1 <- c()
  #list_ar <- c()
  list_asjd <- c()
  
  for (e in list_eps){
    print('e:')
    print(e)
    for (l in list_L){
      print('l:')
      print(l)
      X = HMC_samples(N_HM, epsilon = e, L = l, q0)
      #Trace
      plot.ts(X,  xlab = 't', ylab = 'Xt',
              cex.lab=1.5, col = 'black')
      #Ergodic Mean
      x_mean = cumsum(X)/seq_along(X)
      print('Mean:')
      print(mean(x_mean[(length(x_mean)-1000):length(x_mean)]))
      print(mean(x_mean[(length(x_mean)-10):length(x_mean)]))
      print('...')
      plot(seq_along(x_mean), x_mean,
           xlab = 't', ylab = 'Xt', #main = 'Ergodic Mean',
           type = 'l', col = 'blue',lwd = '2', cex = 1.5)
      #Samples
      hist(X, 100, main = '', prob = TRUE, col = 'lightcyan2') #400 --> 100
      lines(x, scale_fx, type = 'l', lwd=1.5, col = 'red')
      #Diag
      print('Diagnostics')
      #ACF
      print('ACF:')
      print(acf(X, plot = FALSE)[1][1])
      list_acf1[count] = acf(X, plot = FALSE)[1][1]
      #ASJD
      print('ASJD:')
      print(asjd(X))
      list_asjd[count] = asjd(X)
      
      count = count + 1
    }
    
  }
}

#Apply
q0 = 0
N_HM = 50000
list_eps = c(0.2, 0.3, 0.5)
list_L = c(15, 25, 40)
hamiltonian_results(N_HM, q0, list_eps, list_L, x, scale_fx)



#******************
#Assesemnt: Hist
#Erg mean
hm_hist = function(list_eps, list_L, N_HM, x, scale_fx){
  #Plot
  par(mfrow=c(2,3))
  for (e in list_eps){
    print(e)
    for (l in list_L){
      #Samples
      X = HMC_samples(N_HM, epsilon = e, L = l, q0)
      hist(X, 100, main = '', prob = TRUE, col = 'lightcyan2') #400 --> 100
      lines(x, scale_fx, type = 'l', lwd=1.5, col = 'red')
      
    }
    
  }
} 

hm_hist(list_eps, list_L, N_HM, x, scale_fx)














#**********************************************
#Draft
#Tuning Parameters

#i. Trace plot
plot_ts_samp_hm = function(list_sd, N_samp){
  par(mfrow=c(2,3))
  for (sd_i in list_sd){
    print(sd_i)
    X = slice_samp_metropolis(N = N_samp, sigma = sd_i, x0 = 2)
    plot.ts(X,  xlab = 't', ylab = 'Xt',
            cex.lab=1.5, col = 'black') # main = 'Trace Plot',, main = substitute(paste('Trace Plot, sigma = ', sd_i), sd_i))
  }
  
}



#Erg mean
hm_erg_mean_samp = function(list_eps, list_L, N_HM){
  #Plot
  par(mfrow=c(2,3))
  for (e in list_eps){
    print(e)
  for (l in list_L){
    X = HMC_samples(N_HM, epsilon, L = l, q0)
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
} 

#Apply
list_eps = list(0.3)
list_L = list(10)
hm_erg_mean_samp(list_eps, list_L, N_HM)
erg_mean_samp(list_sd, N)

list_eps = list(0.3)
list_L = list(25)
hm_erg_mean_samp(list_eps, list_L, N_HM)



#Samples
#Erg mean
hm_erg_mean_samp = function(list_eps, list_L, N_HM){
  #Plot
  par(mfrow=c(2,3))
  for (e in list_eps){
    print(e)
    for (l in list_L){
      X = HMC_samples(N_HM, epsilon, L = l, q0)
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
} 


#**********************
#ACCEPTANCE RATE
AR_HMC_samples = function(N, epsilon, L, q0){
  
  q_samp = vector('numeric', N)
  current_q = q0
  count_accept = 0 
  count_na = 0
  for (t in 1:N){
    
    q = current_q
    p = rnorm(length(q),0,1) # independent standard normal variates
    current_p = p
    
    
    #U_gradient
    grad_U_q = gradient(f = get_fx_log, q)
    
    # Make a half step for momentum at the beginning
    p = p - epsilon * grad_U_q / 2
    # Alternate full steps for position and momentum
    for (i in 1:L){
      
      # Make a full step for the position
      q = q + epsilon * p
      #U_gradient
      grad_U_q = gradient(f = get_fx_log, q)
      # Make a full step for the momentum, except at end of trajectory
      if (i!=L){
        p = p - epsilon*grad_U_q
      } 
    }
    #U_gradient
    grad_U_q = gradient(f = get_fx_log, q)
    # Make a half step for momentum at the end.
    p = p - epsilon*grad_U_q / 2
    # Negate momentum at end of trajectory to make the proposal symmetric
    p = -p
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_U = get_fx_log(current_q)# U(current_q)
    current_K = sum(current_p^2) / 2
    proposed_U = get_fx_log(q)
    proposed_K = sum(p^2) / 2
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    prob_accept = current_U-proposed_U+current_K-proposed_K
    if (is.na(prob_accept)){
      count_na = count_na + 1
    } else if (runif(1) < exp(prob_accept)){
      current_q = q
      count_accept = count_accept + 1
    }
    q_samp[t] = current_q
    
  }
  print('Accept rate:')
  print(count_accept/N)
  ar= count_accept/N
  
  list(q_samp, ar)
}

#Apply
arr_apply = function(list_eps, list_L, N_HM){
  list_ar = c()
  count = 0 
  #Plot
  par(mfrow=c(2,3))
  for (e in list_eps){
    print('e')
    print(e)
    for (l in list_L){
      print(l)
      vec = AR_HMC_samples(N_HM, epsilon = e, L = l, q0)
      print(vec[2])
      list_ar[count] = vec[2]
      
    }
    
  }
  list_ar
}

arr_apply(list_eps, list_L, N_HM)
