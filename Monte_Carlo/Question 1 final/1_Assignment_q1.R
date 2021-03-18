#Assignment 1
#To do - If time: add (comments)
#Questions
#Question: Inverse
#Q2. *Ask does num samples have to be exact

#Question 1 Part a

#Function * Comment
# a, b > 0 
f_density <- function(x, a, b){
  
  #beta_function = (factorial(alpha - 1)*factorial(beta-1))/factorial(alpha + beta - 1)
  constant = a*b
  y = constant*(x^(a - 1)*(1-(x^a))^(b - 1))
  
  y
}


#Implement
#Params #1
x = seq(0,1,length=1000)
a = 1/5
b = 1/5
fx1 = f_density(x, a, b) #:D 

#Plot
plot(x, fx1, xlab='x', ylab='f(x)', main='f(x, 0.2, 0.2) pdf', xlim = range(0,1), 
     cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='red') 
#dev.off

#Params #2
a = 1
b = 1
fx2 = f_density(x, a, b) 

#Plot
plot(x, fx2, xlab='x', ylab='f(x)', main='f(x, 1, 1) pdf', xlim = range(0,1),
     cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='blue') 

#Params #2
a = 5
b = 5
fx3 = f_density(x, a, b) 

#Plot
#pdf('beta1.pdf')
plot(x, fx3, xlab='x', ylab='f(x)', main='f(x, 5, 5) pdf', cex.lab=1.5, cex.main=1.5, type='l', xlim = range(0,1), lwd=2, col='orange') 

#Plot altogether
plot(x, fx1, xlab='x', ylab='f(x)', main='Beta pdfs', cex.lab=1.5, cex.main=1.5, ylim = range(0,3), type='l', lwd=2, col='red') 
lines(x, fx2, xlab='x', ylab='f(x)', main='f(x, 1, 1) pdf', cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='blue') 
lines(x, fx3, xlab='x', ylab='f(x)', main='f(x, 5, 5) pdf', cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='orange')
legend(0.65, 2.9, legend=c("f(x, 0.2, 0.2)", "f(x,1,1)", "f(x, 5, 5)"),
       col=c("red", "blue", 'orange'), lty=1:3, cex=0.65)


#Check actual function
plot(x,dbeta(x,5,5), type= 'l', lwd=2)


#*********************************************************
#Part c: Rejection Sampling

rejection_sampler = function(a, b, num_req_samps){
  
  #Setup + get bound
  z = seq(0,1, length = 1000)
  fx = f_density(z, a, b)
  bound = max(fx[is.finite(fx)])
  print(bound)
  num_proposals = num_req_samps + round(num_req_samps*(bound)) #Question - Is that ok?
  print(num_proposals)
  
  #Accept/Reject samples
  #runif(num_proposals, min = 0, max=bound)
  #accept <- u <= (f_density(x, a, b)) 
  
  #Accept/Reject samples v2
  x = runif(num_proposals)
  u = runif(num_proposals) 
  
  #Try
  accept <- u <= (f_density(x, a, b)/bound) 
  
  
  x_accept = x[accept]
  x_accept = x_accept[1:num_req_samps]
  x_accept
  
}


#Evaluate
num_req_samps = 1000000
a = 2
b = 5
fx_rs = rejection_sampler(a, b, num_req_samps)

#Plot
z = seq(0,1, length = num_req_samps)
hist(fx_rs, prob= TRUE, 
     main = 'Empirical Histogram and Density for f(x,2,5)', 
     xlab = 'x', ylab = 'f(x)', xlim = c(0,1),
     cex.lab=1.5, cex.main=1.5)
lines(z, f_density(z, a, b), lwd = 2, col = 'red')

#Plot Target density
#Plot
f_target = f_density(z, 2, 5)
plot(z, f_target, xlab='x', ylab='f(x)', main='f(x, 2, 5) pdf', xlim = range(0,1),
     cex.lab=1.5, cex.main=1.5, type='l', lwd=2, col='blue') 
f_density(runif(1), 2,5)

#Evaluate Test/check
num_req_samps = 100000
a = 10
b = 3
fx_beta_rs = rejection_sampler(a, b, num_req_samps)

#Plot
z = seq(0,1, length = num_req_samps)
hist(fx_beta_rs, prob= TRUE, 
     main = 'Empirical Histogram and Density for Beta(2,5)', 
     xlab = 'x', ylab = 'f(x)', xlim = c(0,1),
     cex.lab=1.5, cex.main=1.5)
lines(z, dbeta(z, a, b), lwd = 2, col = 'red')

#Rejection sampler + loop (x use)
rejection_sampler_loop = function(a, b, num_req_samps){

  #Setup
  xs <- c()
  count <- 0
  z = seq(0,1, length = 1000)
  fx = beta_density(z, a, b)
  bound = max(fx[is.finite(fx)])
  print(bound)
  
  for(i in 1:num_req_samps) {
    U <- 2;
    X <- 0
    while(U > beta_density(X,a,b)/bound){
      X <- runif(1)
      U <- runif(1)
      count <- count + 1
    }
    xs[i] <- X
  }
  
  xs
}

#Evaluate + second function :) 

#Evaluate
num_req_samps = 5000
a = 10
b = 3
fx_beta_rs = rejection_sampler_loop(a, b, num_req_samps)

#Plot
z = seq(0,1, length = num_req_samps)
hist(fx_beta_rs, prob= TRUE, 
     main = 'Empirical Histogram and Density for Beta(10,3)', 
     xlab = 'x', ylab = 'f(x)', xlim = c(0,1),
     cex.lab=1.5, cex.main=1.5)
lines(z, dbeta(z, a, b), lwd = 2, col = 'red')


#*********************************************************
#Part e
get_max_function = function(g, num_req_samps){
  
  
  #Setup + get bound (Mabye put in another function)
  z = seq(0,1, length = 1000)
  #f1 - dx
  dx = beta_density(z, 2, g)
  bound_dx = max(dx[is.finite(dx)])
  print(bound_dx)
  #f2 - gx
  gx = beta_density(z, 2, 1/g)
  bound_gx = max(gx[is.finite(gx)])
  print(bound_gx)
  #Use --> num proposals
  max_bound = max(dx[is.finite(dx)], gx[is.finite(gx)])
  num_proposals = num_req_samps*(max_bound + max_bound*0.1) #Question - Is that ok?
  print(num_proposals)
  
  #Vars Setup
  count = 0 
  hs = c()
  
  for (i in 1: num_req_samps){
    
    u = runif(1)
    x = runif(1)
    
    #Find max
    if (beta_density(x, 2, g) >= beta_density(x, 2, 1/g)) {
      fmax = beta_density(x, 2, g)
      bound_fmax = bound_dx
      
    } else if (beta_density(x, 2, 1/g) > beta_density(x, 2, g)) {
      fmax = beta_density(x, 2, 1/g)
      bound_fmax = bound_gx
    }
    
    #Accept/Reject
    if (u <= fmax/bound_fmax) {
      hs[i] = x
    }
    
    #Continue until num_req_samps
    if (length(hs) >= num_req_samps){
      break
    }
    
    #Print progress
    count = count + 1
    if (count%%100000 == 0){
      sprintf('count = %i', count)
    }

    
  }
  hs
}


#Evaluate function
g = 5
num_req_samps = 1000000
h_samp = get_max_function(g, num_req_samps)

#Plot
hist(h_samp, prob= TRUE)
#Plot densities
z = seq(0,1, length = 5000)
lines(z, dbeta(z, 2, 5),type = 'l', col = 'blue')
lines(z, dbeta(z, 2, 1/5), col = 'red')


#***********
#Evaluate
#Plot densities
z = seq(0,1, length = 5000)
plot(z, dbeta(z, 2, 5),type = 'l', col = 'blue')
lines(z, dbeta(z, 2, 1/5), col = 'red')

#Bounds
gx = beta_density(z, 2, 5)
bound_gx = max(gx[is.finite(gx)])
bound_gx

dx = beta_density(z, 2, 1/5)
bound_dx = max(dx[is.finite(dx)])
bound_dx

#Compare
v = runif(1)
fmax2 = max(beta_density(v, 2, 1/5), beta_density(v, 2, 5))

max(fmax2[is.finite(fmax2)])

#Other
plot(z, dbeta(z, 2, 1/5),type = 'l', col = 'blue')
lines(z, dbeta(z, 2, 5), col = 'red')


plot(z, dbeta(z, 10, 3),type = 'l', col = 'blue')
lines(z, dbeta(z, 2, 1/5),type = 'l', col = 'blue')




#*********************************************************
#Part c: Inversion sampling


#Inversion
a_is = 2
b_is = 1/5
u = seq(0,1,length=1000000)
fu = (1 - (1-u)^(1/b))^(1/a)
hist(fu, prob = TRUE)

fx1 = f_density(x, a, b) #:D 
