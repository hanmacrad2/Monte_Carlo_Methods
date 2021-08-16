#APTS

#Question 1: Density
get_fx = function(len_x = 100, X){
  
  #Matrix
  count = 1 
  
  while (count <= len_x) {
    for (x_1 in X[1,]) {
      for (x_2 in X[2,]){
        fx[1, count] = x_1
        fx[2, count] = x_2
        fx[3, count] = (4 -2.1*x_1^2 + (x_1^4)/3)*x_1^2 + x_1*x_2 + 4*(x_2^2 -1)*x_2^2
        
        count = count + 1
        
      }
    } 
  }
  
  
  fx
}

#apply to get fx
len_x = 100
X1 = seq(-3,3,length=len_x) 
X2 = seq(-2,2,length=len_x)
X <- matrix(0, nrow=2, ncol= len_x)#
X[1,] = X1
X[2,] = X2
fx = get_fx(X)

#Function v2
get_fx = function(len_x = 100){
  
  #Matrix
  X1 = seq(-3,3,length=len_x) 
  X2 = seq(-2,2,length=len_x) 
  fx = matrix(0, nrow = 3, ncol = len_x^2)
  count = 1 
  
  while (count <= len_x) {
    for (x_1 in X1) {
      for (x_2 in X2){
        fx[1, count] = x_1
        fx[2, count] = x_2
        fx[3, count] = (4 -2.1*x_1^2 + (x_1^4)/3)*x_1^2 + x_1*x_2 + 4*(x_2^2 -1)*x_2^2
        
        count = count + 1
        
      }
    } 
  }

  
  fx
}

#apply to get fx
fx = get_fx()
#Plot
scatterplot3js(fx[1,], fx[2,], fx[3,], phi = 40, theta = 20,
               color=rainbow(length(z)),
               colkey = FALSE,
               cex = .3,
               main = "f(x1, x2)")

#fx - Mode point
fx_max = max(fx[3,])
fx_max
x1_max = X1[which.max(fx[3,])]
x1_max
x2_max = X2[which.max(fx[3,])]
x2_max



#************************************************************************************
#Simulated Annelaing
#Need to make bivariate!
simulated_annealing_multi = function(N, x0 = 0, beta1 = 1, alpha = 1.001, sigma = 1){
  
  X <- matrix(0, nrow=2, ncol= N)
  X[,1] <- x0
  fx_current = get_fx(X)
  beta_current = beta1
  count_na = 0 
  
  
  U <- runif(n)
  for(i in 2:n) {
    Y <- X[,i-1] + mvrnorm(mu=as.vector(c(0,0)), Sigma=Sigma)
    alpha <- exp(sum(abs(X[,i -1])) - sum(abs(Y)))
    if(U[i] < alpha) {
      X[,i] <- Y
    } else {
      X[,i] <- X[,i-1]
    }
  }
  
  #Vectors for storage
  x = vector('numeric', N)
  fn = numeric(N)
  count_na = 0 
  
  #Variables - current values
  x_current = x0
  fx_current = get_fx(x_current)
  beta_current = beta1
  
  #Implement Markov Chain in loop 
  for (i in 1: N){
    x_new = x_current + rnorm(1)*sigma #Normal random walk proposal  #Multinormal ***
    fx_new = get_fx(x_new)
    alpha_accept = (fx_new/fx_current)^beta_current #alpha -acceptance probability
    if (is.na(alpha_accept)){
      count_na = count_na + 1
    }
    else if (alpha_accept > runif(1)){ #Criterion for acceptance -
      x_current <- x_new #If holds set next sample in the chain as the proposed sample
      fx_current <- fx_new
    }
    x[i] = x_current #If not the next sample in the chain is again set to as the sample value at previous step
    fn[i] = fx_current
    beta_current = beta_current*alpha #Update temperature/value of beta. 
    
    
  }
  print('count_na')
  print(count_na)
  
  list(x, fn)
}



#Plots

x2 = seq(-2,1,length=5) 

for (i in x2) {
  print(i)
}

#Test
u <- seq(-5, 5, by = .1)
v <- seq(-5, 5, by = .1)
M <- expand.grid(u,v)

x <- M$Var1
y <- M$Var2

sigma <- matrix(c(1, .5, .5, 1), nrow = 2, byrow = TRUE)
z <- dmvnorm(x = M, sigma = sigma)

scatterplot3js(x, y, z, phi = 40, theta = 20,
               color=rainbow(length(z)),
               colkey = FALSE,
               cex = .3,
               main = "Bivariate Normal")

scatterplot3js(X1, X2, z, phi = 40, theta = 20,
               color=rainbow(length(z)),
               colkey = FALSE,
               cex = .3,
               main = "f(x1, x2)")

#Plots
plot
library(plotly)

# Data: volcano is provided by plotly

# Plot
persp(X1, X2, fx[3,], 
      xlab = "x1", ylab = "x2",
      main = "f(x1, x2)"
)

library(rgl)
library(plot3D)
library(threejs)
library(mvtnorm)

surf3D(x = fx[1,], y = fx[2 ,], z = fx[3 ,], type = "surface")
rglwidget(elementId = "plot3drgl")

z = fx[3 ,]
z = matrix(fx[3 ,], nrow = 1, ncol = len_x^2)
fig <- plot_ly(x = fx[1,], y = fx[2 ,], z = z) %>% add_surface()

fig


#Other
scatterplot3js(X1, X2, z, phi = 40, theta = 20,
               color=rainbow(length(z)),
               colkey = FALSE,
               cex = .3,
               axisLabels=c("x1", "x2", "f(x1, x2"),
               main = "f(x1, x2)")

