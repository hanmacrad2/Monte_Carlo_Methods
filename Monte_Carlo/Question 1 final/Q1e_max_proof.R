#Question 1 part e
f_density <- function(x, a, b){
  constant = a*b
  y = constant*(x^(a - 1)*(1-(x^a))^(b - 1))
  
  y
}


##############

inversion_sampler = function(a, b, num_samps){
  
  #Setup
  u = runif(num_samps)
  fu = (1 - (1-u)^(1/b))^(1/a)
  
  fu
  
}

h_sampler = function(g, num_samps){
  
  x_samp = c()
  
  for(i in 1:num_samps) {
    u = 1
    w = 1/2
    iter=0
    threshold=(max(f_density(w, 2, 1/g), f_density(w, 2, g)))/(f_density(w, 2, 1/g) +f_density(w, 2, g))
    
    while(u > threshold) {
      
      X <- inversion_sampler(2, 1/g, 1)
      Y <- inversion_sampler(2, g, 1)
      
      u <- runif(1)
      v <- runif(1)
      
      w <- X*(v>1/2)+Y*(v<1/2)
      iter=iter+1
     
      
  
      threshold=(max(f_density(w, 2, 1/g), f_density(w, 2, g)))/(f_density(w, 2, 1/g) +f_density(w, 2, g))
      if(is.nan(threshold)){threshold=1}
        
      #print(c(threshold,w,iter))
    
    }
    x_samp [i] <- w
  }
  return(x_samp)
  
}
  
#Evaluate
num_samps = 10
g = 5
x_samp = h_sampler(g, num_samps)

#Plot
hist(x_samp, prob = TRUE)

#Target
x = seq(0,1, length = num_samps)
hx1 = fx(x, 2, 5)
hx2 = fx(x, 2, 1/5)
scaling_factor = 1.8
hx1_scaled = hx1/scaling_factor
hx2_scaled = hx2/scaling_factor

#Plots
hist(x_samp, prob = TRUE)
lines(x, hx1_scaled, 
      type = 'l', lwd = 2, col = 'red')
lines(x, hx2_scaled, col = 'green', lwd = 2)

###############
#Function no loop
h_sampler2 = function(g, num_samps){
  
  #1. Get two densities
  h1 = inversion_sampler(2, 1/g, num_samps)
  h2 = inversion_sampler(2, g, num_samps)
  
  #Simulate uniform variables 
  u = runif(num_samps)
  i = runif(num_samps)
  
  #Coin toss between two distributions (Bernoulli)
  w = h1*(i>1/2) + h2*(i<1/2)
  
  #Which density is max 
  threshold=(max(f_density(w, 2, 1/g), f_density(w, 2, g)))/(f_density(w, 2, 1/g) +f_density(w, 2, g))
  accept_w <- u <= threshold
  
  #Accepted samples
  x_samp = w[accept_w]
  
  x_samp
  
}

#Evaluate
num_samps = 1000000
g = 8
start_time = Sys.time()
h_samp = h_sampler2(g, num_samps)
end_time = Sys.time()
time_hs = end_time - start_time
time_hs

#Target
x = seq(0,1, length = num_samps)
hx1 = fx(x, 2, b)
hx2 = fx(x, 2, 1/b)
hx1_scaled = hx1/scaling_factor
hx2_scaled = hx2/scaling_factor



#Plots

hist(h_samp, prob = TRUE, xlim = c(0,1),  
     xlab = 'x', ylab = 'h(x, g)', cex.lab=1.5,
     main = 'Empirical Histogram of h(x,g) (n = 1 million) & Target Densities f(x,g)')
lines(x, hx1_scaled, 
      type = 'l', lwd = 2, col = 'red')
lines(x, hx2_scaled, col = 'green', lwd = 2)
legend('topleft', legend = c('f(x, 2, 2)', 'f(x, 2, 1/2)'), pch = 16, col = c('red', 'green'))



lines(x, hx1_scaled, 
      type = 'l', lwd = 2, col = 'red')
lines(x, hx2_scaled, col = 'green', lwd = 2)

#Plot
par(mfrow=c(1,1))
hist(x_samp, prob = TRUE)

#Draft
h_sampler = function(g, num_samps){
  
  x_samp = c()
  
  for(i in 1:num_samps) {
    u = 1
    w = 1/2
    iter=0
    threshold=(max(f_density(w, 2, 1/g), f_density(w, 2, g)))/(f_density(w, 2, 1/g) +f_density(w, 2, g))
    
    while(u > threshold) {
      
      X <- inversion_sampler(2, 1/g, 1)
      Y <- inversion_sampler(2, g, 1)
      
      u <- runif(1)
      v <- runif(1)
      
      w <- X*(v>1/2)+Y*(v<1/2)
      print(i)
      o = sprintf('u = %f, v = %f, w:', u, v)
      print(o)
      o = sprintf('X = %f, Y = %f, w:', X, Y)
      print(o)
      print(w)
      print('')
      iter=iter+1
      
      
      
      threshold=(max(f_density(w, 2, 1/g), f_density(w, 2, g)))/(f_density(w, 2, 1/g) +f_density(w, 2, g))
      if(is.nan(threshold)){threshold=1}
      
      #print(c(threshold,w,iter))
      
    }
    x_samp [i] <- w
  }
  return(x_samp)
  
}