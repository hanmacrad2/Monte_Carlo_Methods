#Assignment 1 part d: Rejection Sampling

#Function 
f_density <- function(x, a, b){
  
  #beta_function = (factorial(alpha - 1)*factorial(beta-1))/factorial(alpha + beta - 1)
  constant = a*b
  y = constant*(x^(a - 1)*(1-(x^a))^(b - 1))
  
  y
  
}

#Rejection sampling function
rejection_sampler = function(a, b, num_req_samps, bound){
  #Setup
  hs = c()
  count = 0
  j = 1
  for (i in 1: num_req_samps){
    count = count + 1
    x = runif(1)
    u = runif(1, min = 0.0, max = bound)
    
    fx = f_density(x, a, b)
    #print(count)
    #print(x)
    #print('fx')
    #print(fx)
    #print('u')
    #print(u)
    #sprintf('fx = %f', fx)
    #sprintf('u = %f', u)
    
    #Accept/Reject
    if (u <= fx) {
      hs[j] = x
      j = j+1
      #print('j')
      #print(j)
    }
    
  }
  
  hs

}

#Evaluate
num_req_samps = 10000
a = 2
b = 5
h_samp = rejection_sampler(a, b, num_req_samps, bound)

#Plot
z = seq(0,1, length = num_req_samps)
hist(h_samp, prob= TRUE, 
     main = 'Empirical Histogram and Density for f(x,2,5)', 
     xlab = 'x', ylab = 'f(x)', xlim = c(0,1),
     cex.lab=1.5, cex.main=1.5)
lines(z, f_density(z, a, b), lwd = 2, col = 'red')



#Continue until num_req_samps
#if (length(hs) >= num_req_samps){
#  break
#}


#Draft
#Try
a = 2
b = 5
v = runif(1)
f1 = f_density(v, a, b)
f1

#Bound
z = seq(0,1, length = 1000)
fx = f_density(z, a, b)
bound = max(fx[is.finite(fx)])
print(bound)

u = runif(1, min = 0, max = bound)
u

f1

#Plot
plot(z, fx, type = 'l')
lines(z, bound*dunif(z))