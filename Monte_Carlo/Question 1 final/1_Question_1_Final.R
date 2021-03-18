f#Assignment 1 Question 1
library(rootSolve)

#Part a Function *Redo
fx <- function(x, a, b){
  y = (x^(a - 1)*(1-(x^a))^(b - 1)) #y = constant*(x^(a - 1)*(1-(x^a))^(b - 1))
  y
}

f_density <- function(x, a, b){
  constant = a*b
  y = constant*(x^(a - 1)*(1-(x^a))^(b - 1)) #y = constant*(x^(a - 1)*(1-(x^a))^(b - 1))
  y
}

#Plots
#Params #1
x = seq(0,1,length=1000)
a = 1/5
b = 1/5
fx1 = fx(x, a, b) #:D 
#Plot
plot(x, fx1, xlab='x', ylab='f(x)', main='Density of f(x, 0.2, 0.2)', xlim = range(0,1), 
     cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='red') 


#Params #2
a = 1
b = 1
fx2 = fx(x, a, b) 
#Plot
plot(x, fx2, xlab='x', ylab='f(x)', main='Density of f(x, 1, 1)', xlim = range(0,1),
     cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='blue') 

#Params #3
a = 5
b = 5
fx3 = fx(x, a, b) 
#Plot
plot(x, fx3, xlab='x', ylab='f(x)', main='Density of f(x, 5, 5)', cex.lab=1.5, cex.main=1.5, type='l', xlim = range(0,1), lwd=2, col='orange') 

#Plot altogether ylim = range(0,3)
plot(x, fx1, xlab='x', ylab='f(x)', main='Pdfs of f(x,a,b)', cex.lab=1.5, type='l', lwd=2, col='red') 
lines(x, fx2, xlab='x', ylab='f(x)', main='f(x, 1, 1) pdf', cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='blue') 
lines(x, fx3, xlab='x', ylab='f(x)', main='f(x, 5, 5) pdf', cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='orange')

legend(0.79, 2.9, legend=c("f(x, 0.2, 0.2)", "f(x,1,1)", "f(x, 5, 5)"),
       col=c("red", "blue", 'orange'), lty=1:3, cex=0.65)


#*********************************************************
#Part c Prelim

#Plot target
a = 2
b = 5
z = seq(0,1, length = 1000)
plot(z, f_density(z, a, b), type = 'l', lwd = 2, col = 'red')
#Point
points(root_fx, f_density(root_fx, a, b), col = 'blue')

#Add bound

#*********************************************************
#Part c: Rejection Sampling

#Step 1: Get bound
evaluate_fx_deriv = function(a, b){
  function (x) (a-1)*x^(a-2)*(1-x^a)^(b-1) -a*(b-1)*(x^(2*a-2)*(1-x^a)^(b-2)) 
}


#Rejection Sampler
rejection_sampler = function(a, b, num_req_samps){
  
  #Step 1: Get bound
  fx_deriv = evaluate_fx_deriv(a,b)
  #Roots
  roots_all <- uniroot.all(fx_deriv, c(0, 1))
  root_fx = roots_all[-1]
  print(root_fx)
  bound = f_density(root_fx, a, b)
  #print(bound)
  
  #Num propo
  num_proposals = num_req_samps + round(num_req_samps*(bound)) #Question - Is that ok?
  #print(num_proposals)
  
  #Rejection algorithm
  x = runif(num_proposals)
  u = runif(num_proposals) 
  
  accept <- u <= (f_density(x, a, b)/bound) #Bound ok as less than or equal to 
  
  
  x_accept = x[accept]
  x_accept = x_accept[1:num_req_samps]
  x_accept
  
}

#Apply rejection sampler
num_req_samps = 1000000
a = 2
b = 1/5
start_time = Sys.time()
fx_rs = rejection_sampler(a, b, num_req_samps)
end_time = Sys.time()
time_rs = end_time - start_time
time_rs

#Plot
z_length = 10000
z = seq(0,1, length = z_length)
hist(fx_rs) #, prob= TRUE)#, 
lines(z, f_density(z, a, b), lwd = 2, col = 'red')

     main = 'Empirical Histogram using Rejection Sampling and Density for f(x,2,5)', 
     xlab = 'x', ylab = 'f(x)', xlim = c(0,1),
     cex.lab=1.5)
lines(z, f_density(z, a, b), lwd = 2, col = 'red')

#*******************************************************
#Part b Inversion Sampling 

inversion_sampler = function(a_is, b_is, num_samps){

  u = runif(num_samps)
  fu = (1 - (1-u)^(1/b_is))^(1/a_is)
  
  fu
}

#Params
a_is = 2
b_is = 1/5
start_time = Sys.time()
fu = inversion_sampler(a_is, 5, 1000000)
end_time = Sys.time()
time_is = end_time - start_time
time_is

#Plot
hist(fu, prob= TRUE, 
     main = 'Empirical Histogram using Inversion Sampling and Density for f(x,2,1/5)', 
     xlab = 'x', ylab = 'f(x)', xlim = c(0,1),
     cex.lab=1.5)
lines(v, f_density(v, a_is, 5), lwd = 2, col = 'red')





#Old
#Target
v = seq(0,1,length=10000) 
plot(v, f_density(v, a_is, b_is), type = 'l', col = 'red')

#Inversion
u = runif(1000000) #seq(0,1,length=1000000) 
fu = (1 - (1-u)^(1/b_is))^(1/a_is)
#hist(fu, prob = TRUE)

