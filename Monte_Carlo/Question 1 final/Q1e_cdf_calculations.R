#Question 1e Max --> cdf

#Function - density
fx <- function(x, a, b){
  norm_constant = a*b
  y = norm_constant*(x^(a - 1)*(1-(x^a))^(b - 1))
  y
}

#Function - cdf 
cdf_F = function(x, a, b){
  F = 1 - (1 - x^a)^b
  F
}

#Get xs F1(Xs) and F2(Xs)
b = 5
b =2 
xs = sqrt(1 - b^((2*b/(1-b^2))))
xs

#Get F1
F1_xs = cdf_F(xs, 2, b)
F1_xs

#Get F1
F2_xs = cdf_F(xs, 2, 1/b)
F2_xs

#Constant
C = F1_xs - F2_xs
C

#Proportion of samples
prop = F1_xs/(F1_xs + (1 - F2_xs))

scaling_factor = F1_xs + (1 - F2_xs)
scaling_factor
  
#**********
#Inversion Sampling
num_samps = 5000
  
inversion_sampler = function(a, b, n_samps){ #range_unif
  
  u = runif(n_samps) #min = range_unif[1], max = range_unif[2]
  fu = (1 - (1-u)^(1/b))^(1/a)
  
  fu
}

#**********
#Inversion
num_samps = 1000000

#Target densities
x = seq(0,1, length = num_samps)
hx1 = fx(x, 2, 5)
hx2 = fx(x, 2, 1/5)

#F1 Inversion
n_samps1 = round(prop*num_samps) #range_u = c(0, xs)
#n_samps1 = round(n_samps1*(1/xs))
n_samps1
b = 5

#Method f[]
u = runif(n_samps1) 
f1_samps = (1 - (1-u)^(1/b))^(1/a)
accept_f <- f1_samps <= xs
f1 = f1_samps[accept_f]
hist(f1)

#Plot
hist(h1_samp) #, prob = TRUE)
#lines(x, hx1, 
#      type = 'l', lwd = 2, col = 'red')
#lines(x, hx2, col = 'green', lwd = 2)

#******* F2xs
n_samps2 = num_samps - n_samps1
n_samps2
b =1/5
u = runif(n_samps2,0, 1)
f2_samps = (1 - (1-u)^(1/b))^(1/a)
accept_f2 <- f2_samps > xs
f2 = f2_samps[accept_f2]
hist(f2)


#Combine
h_samp = c(f1, f2)
#h_samp = scaling_factor*h_samp
hist(h_samp) #, prob = TRUE) #, prob = TRUE)

#Plot
hx1_scaled = hx1/scaling_factor
hx2_scaled = hx2/scaling_factor
hist(h_samp, prob = TRUE, xlim = c(0,1),  
     xlab = 'x', ylab = 'h(x, g)', 
     main = 'Empirical Histogram of h(x,g) & Densities f(x, 2, 5), f(x, 2, 1/5)')
lines(x, hx1_scaled, 
      type = 'l', lwd = 2, col = 'red')
lines(x, hx2_scaled, col = 'green', lwd = 2)


#Plot
hist(f1_samps_keep) #, prob = TRUE)
plot(x, hx1, 
      type = 'l', lwd = 2, col = 'red')
lines(x, hx2, col = 'green', lwd = 2)

#Draft
#Method u[0,0,625]
b = 5
u = runif(n_samps1, 0, xs)
h1_samp = (1 - (1-u)^(1/b))^(1/a)