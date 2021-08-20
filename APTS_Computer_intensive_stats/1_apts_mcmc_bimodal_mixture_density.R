#Computer Intensive Stats - Q2 Bimodal Density - MCMC
library(EnvStats)
#Bimodal mixture density
N <- 100000

components <- sample(1:2,prob=c(0.5, 0.5), size=N, replace=TRUE)
mus <- c(0,10)
sds <- sqrt(c(1,1))

samples <- rnorm(n=N, mean=mus[components], sd=sds[components])

#Plot
plot(density(samples), main="Density of the Bimodal Mixture Model")

#Other 
N <- 100000
components <- sample(1:3,prob=c(0.3,0.5,0.2),size=N,replace=TRUE)
mus <- c(0,10,3)
sds <- sqrt(c(1,1,0.1))
samples <- rnorm(n=N,mean=mus[components],sd=sds[components])

#Hist
hist(samples)

hist2 <- hist(samples, breaks = 80)
hist2$counts <- hist2$counts/sum(hist2$counts)
hist3 = plot(hist2, ylab = 'Density', 
             main = 'Empirical density of Mixture model')

#Mixture density
get_multi_mode_density = function(x){
  
  density = dnormMix(x, mean1 = 0, sd1 = 1, mean2 = 10, sd2 = 1, p.mix = 0.5)
}


#*****************
#*MCMC

#*******************************************************
#Part 1: Random Walk Metropolis

random_walk_metropolis = function(sigma, N = 100000, x0 = 1){
  
  #Vectors to store samples
  X = vector('numeric', N)
  X[1] = x0 #Initialise 1st sample in chain in order to start the chain running 
  count_accept = 0 #Count of number of accepted values 
  #Loop
  for (t in 2:N){
    X[t] = X[t-1] + rnorm(1, mean = 0, sd= sigma) #Symmetric random walk proposal
    u = runif(1)
    
    if(u > (get_multi_mode_density(X[t])/get_multi_mode_density(X[t-1]))) { #Criterion for rejection 
      X[t] = X[t-1] #If criterion not met the next sample <- current sample
      count_accept = count_accept + 1
    }
  }
  print('count_accept')
  print(count_accept/N)
  X  
}

#list(X, count_accept) #Return samples and count of accepted values

#Apply
start_time = Sys.time()
X_rw = random_walk_metropolis(sigma = 5.0)
end_time = Sys.time()
time_elapsed = end_time - start_time
time_elapsed

#Plots
plot.ts(X_rw)

#Hist
hist_rw <- hist(X_rw, breaks = 80)
hist_rw$counts <- hist_rw$counts/sum(hist_rw$counts)
hist_rw = plot(hist_rw, ylab = 'Density', 
             main = 'Empirical density of RW MCMC chain')
print(hist_rw)

#*******************************************************
#Part 1: Slice Sampler 

slice_sampler_metropolised = function(N = 100000, sigma = 2.0, x0 = 0){
  
  #x- markov chain
  x = vector('numeric', N)
  y = vector('numeric', N)
  x[1] = x0
  count_accept = N
  
  #Markov Chain starts
  for (t in 2:N){
    
    #Auxilliary Uniform variable u
    u[t] = runif(1, 0, get_multi_mode_density(x[t-1]))
    
    #X sampled from markov chain using a symmetric random 
    x[t] = x[t-1] + rnorm(1, mean = 0, sd= sigma)
    
    #Using augmented distribution f(x,u) 
    if(u[t] > get_multi_mode_density(x[t])) { #Rejection criterion
      x[t] = x[t-1]
      count_accept = count_accept - 1 
    }
  }
  print(count_accept/N) #Acceptance rate
  x
}

#Time it
start_time = Sys.time()
X_ss = slice_sampler_metropolised()
end_time = Sys.time()
time_elapsed = end_time - start_time
time_elapsed

#Plots
plot.ts(X_ss) #, 100) #Note x mix aswell

#Hist 
hist_ss <- hist(X_ss, breaks = 80)
hist_ss$counts <- hist_ss$counts/sum(hist_ss$counts)
hist_ss = plot(hist_ss, ylab = 'Density', 
               main = 'Empirical density of Metropolised Slice Sampler chain')
print(hist_ss)


#**********************************
#Part 3 Metropolised 

slice_sampler = function(Ns = 5000, w = 0.5, x0 = 1.0){
  
  #Initialise variables
  xs = vector('numeric', Ns)
  x_current = x0; count_accept = 0 
  
  #Iterations
  for (i in 1: Ns) {
    y = runif(1, 0, get_multi_mode_density(x_current)) #uniform sample
    lb = x_current #Left bound
    rb = x_current #Right bound
    
    #Interval w is randomly positioned around x0
    #Expanded in steps of size w until both ends are outside the slice
    while (y < get_multi_mode_density(lb)) {
      lb = lb - w
    }
    while (y < get_multi_mode_density(rb)) {
      rb = rb + w
    }
    
    #X _new - uniformly pick point from interval until a point inside the slice is found
    x_new = runif(1, lb, rb)
    #Points picked that are outside the slice are used to shrink the interval.
    if (y > get_multi_mode_density(x_new)) {
      if (abs(x_new - lb) < abs(x_new - rb)) {
        lb = x_new
      } else {
        lb = y
      }
    }
    else {
      x_current = x_new #Set new value of x to current value
      count_accept = count_accept + 1
    }
    xs[i] = x_current    #Include current x as a sample
  }
  print('count_accept')
  print(count_accept)
  xs
}

#Time it
start_time = Sys.time()
X_ssr = slice_sampler()
end_time = Sys.time()
time_elapsed = end_time - start_time
time_elapsed

#Plots
plot.ts(X_ssr) #, 100) #Note x mix aswell

#Hist 
hist_ssr <- hist(X_ssr, breaks = 80)
hist_ssr$counts <- hist_ssr$counts/sum(hist_ssr$counts)
hist_ssr = plot(hist_ssr, ylab = 'Density', 
               main = 'Empirical density of RW MCMC chain')
print(hist_ssr)

