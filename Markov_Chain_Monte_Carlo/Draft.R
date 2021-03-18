#MC Drafts

#Own Method
#Gradient
library(tensorflow)
par(mfrow=c(1,1))
tensorflow::au
tf$autograd

#Hamiltonian 

#Hamiltonian
HMC = function(U, grad_U, epsilon, L, current_q)
{
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
    if (i!=L) p = p - epsilon*grad_U_q
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
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  {
    return (q) # accept
  }
  else
  {
    return (current_q) # reject
  }
}

#

#Funtion 
epsilon = 0.3
current_q = 1 #Or Random np.random.randn(2)
X_hmc = HMC(fx_log, grad_U, epsilon, length(x), current_q)

HMC = function (U, grad_U, epsilon, L, current_q)
{
  q = current_q
  p = rnorm(length(q),0,1) # independent standard normal variates
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * grad_U(q)
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  {
    return (q) # accept
  }
  else
  {
    return (current_q) # reject
  }
}

































#*********** Draft *********************
#Proposal: uniform
rwmh_uni = function(N, sigma, x0 = 2){
  
  #Params
  X = vector('numeric', N)
  X[1] = x0
  
  #Loop
  for (t in 2:N){
    X[t] = X[t-1] + runif(1, X[t-1] - sigma, X[t-1] - sigma) 
    u = runif(1)
    
    if(u > (get_fx(X[t])/get_fx(X[t-1]))) {
      X[t] = X[t-1]
    }
  }
  X
}

#Apply
X_rw_uni = rwmh_uni(N2, sigma = 6.4, x0 = 2)
#i. Trace Plots Plot trace
plot.ts(X_rw)
#Hist samples
hist(X_rw)











#f - invariant distribution
Rejection sample + t-distribution

def p(x):
  return st.norm.pdf(x, loc=30, scale=10) + st.norm.pdf(x, loc=80, scale=20)


def q(x):
  return st.norm.pdf(x, loc=50, scale=30)


x = np.arange(-50, 151)
M = max(p(x)/q(x))


def rejection_sampling(iter=1000):
  samples = []

for i in range(iter):
  z = np.random.normal(50, 30)
u = np.random.uniform(0, k*q(z))

if u <= p(z):
  samples.append(z)

return np.array(samples)


#Rejection Sampling
num_proposals = 50000