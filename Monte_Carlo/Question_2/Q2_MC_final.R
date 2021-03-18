#Question 2
library(scales)

#Time it function

#Part c i: Plot function

f = function(y){
  fy = exp(-(y-2)^2) + exp(-(y+2)^2)
  fy
}

fy_eval = function(y){
  fy = exp(-(y-2)^2) + exp(-(y+2)^2)
  fy
}


#Evaluate
y = seq(-5,5, length = 5000)
fy = fy_eval(y)

#Plots
plot(y, fy, type = 'l', ylab ='f(y)', lwd=2, cex.lab=1.5, main = 'Density of f(y)')


#Time it
start_time <- Sys.time()
fy = f_eval(y)
end_time <- Sys.time()
end_time - start_time


#Norm constant
norm_constant = 7.174476 # (1/0.139383) Using wolfram alpha
plot(y, norm_constant*fy, type='l')


#Part c ii.
#Vars
v = seq(-10, 10, length = 1000) #v = seq(-5,5, length = 200)
u = seq(0, 10, length = 1000) #u = seq(0, 1.5, length = 200)

#f
f <- sqrt(f(v / u))
#Evaluate
data <- within(data, f <- sqrt(f(v / u)))



#********************************************************
#Part c ii
A_eval <- function(ubnds, vbnds, dimension1){
  
  #Setup variables
  v = seq(ubnds[1],ubnds[2], length = dimension1) 
  u = seq(vbnds[1], vbnds[2], length = dimension1) 
  data <- expand.grid(v=v, u=u)
  
  #Evaluate function at y = v/u
  data <- within(data, f <- sqrt(fy_eval(v / u))) #Evaluate function for 
  #Check where boundary constraints met i.e  u <= fy(v/u)
  data$check <- data$u <= data$f
  
  #Plot v vs u and colour according to whether the boundary check satisfied or not
  plot(v ~ u, data = data, col = data$check + 1, 
       pch = 16, cex.lab=1.5, main = 'Visualisation of region A')
  legend('topleft', legend = c('Not A', 'A'), col = 1:2, pch = 16)
  
}
  
#Evaluate function
ubnds = c(-2.5, 2.5)
vbnds = c(0.001, 1.02)
A_eval(ubnds, vbnds, 500)

#Plot
plot(v ~ u, data = data, col = data$check + 1, 
       pch = 16, cex.lab=1.5, main = 'Visualisation of region A')
legend('topleft', legend = c('Not A', 'A'), col = 1:2, pch = 16)

  
#Method 2
#Vars
v = seq(-2.5,2.5, length = 500) #v = seq(-5,5, length = 200)
u = seq(0.001, 1.02, length = 500) #u = seq(0, 1.5, length = 200)
data <- expand.grid(v=v, u=u)

#Evaluate
data <- within(data, f <- sqrt(f(v / u)))
#Bounds
data$check <- data$u <= data$f

#Plots
start_time <- Sys.time()
plot(v ~ u, data = data, col = data$check + 1, 
     pch = 16, cex.lab=1.5, main = 'Visualisation of region A')
legend('topleft', legend = c('Not A', 'A'), col = 1:2, pch = 16)

plot(f ~ u, data = data, col = data$check + 1, pch = 16)
legend('topleft', legend = c(FALSE, TRUE), col = 1:2, pch = 16)

plot(f ~ v, data = data, col = data$check + 1, pch = 16, cex.lab=1.5)
legend('topleft', legend = c(FALSE, TRUE), col = 1:2, pch = 16)

#******************
#True/In Boundary 
data_bounded = data[data$check == TRUE, ]
data_bounded$grad = data_bounded$v/data_bounded$u
summary(data_bounded)

#Plot
plot(v ~ u, data = data_bounded, col = 'red', pch = 16, cex.lab=1.5, main = 'Visualisation of region A')

legend('topleft', legend = c('Not A', 'A'), col = 1:2, pch = 16)

#Find B
summary(data_bounded[data_bounded$u > 0,])


#***********************
#Part d 
#Function
A_sampler_loop = function(num_samps, ubnds, vbnds){
  
  #Variables
  us = c()
  vs = c()
  j = 1
  
  for (i in 1:num_samps){
    
    u = runif(1, ubnds[1], ubnds[2]) #uniform [0, 1]
    v = runif(1, vbnds[1], vbnds[2]) #uniform [-2.5, 2.5]
    
    if (u <= sqrt(fy_eval(v/u))) { #If constraint on the boundary of region A is met
      us[j] = u
      vs[j] = v
      j = j+1
      
    }

  }

  return(list(us = us,vs = vs))
}

#***********************
#Part d 
#Function no loop
A_sampler = function(num_samps, ubnds, vbnds){
  
  #1. Variables
  us = c()
  vs = c()
  j = 1
  num_samps = num_proposals*4
  
  #2. Uniformly sample over Rectangle
  u = runif(num_samps, ubnds[1], ubnds[2]) #uniform [0, 1]
  v = runif(num_samps, vbnds[1], vbnds[2]) #uniform [-2.5, 2.5]
  
  #3. Accept samples based on constraint on the boundary of region A
  accept <- u <= sqrt(fy_eval(v/u))
  us = u[accept]
  vs = v[accept]
  
  #Return region A i.e (u,v)
  return(list(us = us,vs = vs))
}

#Samples Apply function
ubnds = c(0.01, 1.4)
vbnds = c(-2.5, 2.5) #ubnds = c(0.0, 2) vbnds = c(-4, 4)
num_proposals = 100000 #000
#num_samps = num_proposals*4

start_time <- Sys.time()
As_100k =  A_sampler(num_proposals, ubnds, vbnds)
end_time <- Sys.time()
time_As = end_time - start_time
print(time_As)
print(length(As$us))

#Plot
par(mfrow=c(1,1))
plot(As_100k$us, As_100k$vs, xlab = 'u', ylab = 'v',
     main = 'Empirical samples of fy (100,000 samples)',
     cex.lab=1.5, col = 'orange', pch = 16)

#**********
#Comparison Plot
par(mfrow=c(2,1))
plot(As$vs, As$us, xlab = 'v', ylab = 'u',
     main = 'Empirical samples of f(v/u) (1,000,000 samples)',
     cex.lab=1.5, col = 'orange', pch = 16)
plot(u ~ v, data = data_bounded, 
     xlab = 'v', ylab = 'u', main = 'Analytical evaluation of f(v/u)',
     col = 'red', pch = 16, cex.lab=1.5)



plot(As$us, As$vs, xlab = 'u', ylab = 'v',
     main = 'Analytical evaluation of f(v/u)',
       cex.lab=1.5, xlim = c(0.001, 1),
       col = alpha('green', 0.15)) #, ylim = vbnds,
points(v ~ u, data = data_bounded, col = alpha('orange', 0.1), cex.lab=1.5, type = 'o')


points(v ~ u, data = data_bounded, 
     col = alpha('orange', 0.2), cex.lab=1.5, type = 'o')
legend('topleft', legend = c('Empirical', 'Analytical'), pch = 16, col = c('orange', 'green'))


plot(v ~ u, data = data_bounded, 
       col = alpha('orange', 0.2), cex.lab=1.5, type = 'o')


lines(v ~ u, data = data_bounded, 
      xlim = c(0.001, 1), col = alpha('blue', 0.4), 
      pch = 16, cex.lab=1.5)  #, main = 'Visualisation of region A')

#************
#Draft
#Function
A_sampler = function(num_samps, ubnds, vbnds){
  
  #Variables
  us = c()
  vs = c()
  j = 1
  
  for (i in 1:num_samps){
    
    u = runif(1, ubnds[1], ubnds[2]) #uniform from a rectangle
    o = sprintf('u = %f', u)
    #print(o)
    v = runif(1, vbnds[1], vbnds[2])
    o = sprintf('v = %f', u)
    #print(o)
    o = sprintf('f = %f', sqrt(fy_eval(v/u)))
    #print(o)
    
    if (u <= sqrt(fy_eval(v/u))) {
      #As[j] = sqrt(f_eval(v/u))
      us[j] = u
      vs[j] = v
      #sprintf('jth sample = %i', j)
      j = j+1
      
    }
    
  }
  
  return(list(us = us,vs = vs))
}
