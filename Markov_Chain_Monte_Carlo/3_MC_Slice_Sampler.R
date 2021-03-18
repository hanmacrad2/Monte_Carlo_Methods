##MC Assign II Part c: Slice sampler

#Params
Ns = 5000
x0 = 0 #Start chain or should be uniform --> check paper 
w = 0.5 # Correct? 
#Function
xs = vector('numeric', Ns)
x_current = x0

#Algorithm II 
for (i in 1: Ns) {
  y = runif(1, 0, get_fx(x_current))
  lb = x_current
  rb = x_current
  
  #interval w is randomly positioned around x0
  #Expanded in steps of size w until both ends are outside the slice
  while (y < get_fx(lb)) {
    lb = lb - w
  }
  while (y < get_fx(rb)) {
    rb = rb + w
  }
  
  #x _new - uniformly pick point from interval until a point inside the slice is found
  x_new = runif(1, lb, rb)
  #Points picked that are outside the slice are used to shrink the interval.
  if (y > get_fx(x_new)) {
    if (abs(x_new - lb) < abs(x_new - rb)) {
      lb = x_new
    } else {
      lb = y
    }
  }
  else {
    x_current = x_new
    
  }
  xs[i] = x_current
}

#Plot/Hist :D 
hist(xs, 100, main = '', prob = TRUE, col = 'lightcyan2') #400 --> 100
lines(x, scale_fx, type = 'l', lwd=1.5, col = 'red')

hist(xs, 50, main = '', prob = TRUE, col = 'lightcyan2') #400 --> 100
lines(x, scale_fx, type = 'l', lwd=1.5, col = 'red')

hist(xs)
plot.ts(xs)



#Algorithm
while(length(xs) < Ns) {
  y = runif(1, 0, get_fx(x_current))
  lb = x_current
  rb = x_current
  
  #interval w is randomly positioned around x0
  #Expanded in steps of size w until both ends are outside the slice
  while (y < get_fx(lb)) {
    lb = lb - w
  }
  while (y < get_fx(rb)) {
    rb = rb + w
  }
  
  #x _new - uniformly pick point from interval until a point inside the slice is found
  x = runif(1, lb, rb)
  #Points picked that are outside the slice are used to shrink the interval.
  if (y > get_fx(x)) {
    if (abs(x - lb) < abs(x - rb)) {
      lb = x
    } else {
      lb = y
    }
  }
  else {
    append(xs, x)
    
  }
  
}

#Plot
hist(xs)

dist = stats.norm(5, 3)
w = 0.5
x = dist.rvs()
niters = 1000
xs = []

while len(xs) < niters:
  y = np.random.uniform(0, dist.pdf(x))
lb = x
rb = x
while y < dist.pdf(lb):
  lb -= w
while y < dist.pdf(rb):
  rb += w
#x _new - uniformly pick point from interval until a point inside the slice is found
x = np.random.uniform(lb, rb)
#Points picked that are outside the slice are used to shrink the interval.
if y > dist.pdf(x):
  if np.abs(x-lb) < np.abs(x-rb):
  lb = x
else:
  lb = y
else:
  xs.append(x)




#Function
slice_sampler2 = function(Ns, sigma, w = 0.5 x0 = 2){
  
  #Initialise variables
  xs = vector('numeric', Ns)
  x_current = x0
  
  #Iterations
  for (i in 1: Ns) {
    y = runif(1, 0, get_fx(x_current)) #uniform sample
    lb = x_current #Left bound
    rb = x_current #Right bound
    
    #interval w is randomly positioned around x0
    #Expanded in steps of size w until both ends are outside the slice
    while (y < get_fx(lb)) {
      lb = lb - w
    }
    while (y < get_fx(rb)) {
      rb = rb + w
    }
    
    #x _new - uniformly pick point from interval until a point inside the slice is found
    x_new = runif(1, lb, rb)
    #Points picked that are outside the slice are used to shrink the interval.
    if (y > get_fx(x_new)) {
      if (abs(x_new - lb) < abs(x_new - rb)) {
        lb = x_new
      } else {
        lb = y
      }
    }
    else {
      x_current = x_new #Set new value of x to current value
      
    }
    xs[i] = x_current    #Include current x as a sample
  }
  
  xs
}
#Params
Ns = 5000
x0 = 0 #Start chain or should be uniform --> check paper 
w = 0.5 # Correct? 
#Function
xs = vector('numeric', Ns)
x_current = x0

#Algorithm II 
for (i in 1: Ns) {
  y = runif(1, 0, get_fx(x_current))
  lb = x_current
  rb = x_current
  
  #interval w is randomly positioned around x0
  #Expanded in steps of size w until both ends are outside the slice
  while (y < get_fx(lb)) {
    lb = lb - w
  }
  while (y < get_fx(rb)) {
    rb = rb + w
  }
  
  #x _new - uniformly pick point from interval until a point inside the slice is found
  x_new = runif(1, lb, rb)
  #Points picked that are outside the slice are used to shrink the interval.
  if (y > get_fx(x_new)) {
    if (abs(x_new - lb) < abs(x_new - rb)) {
      lb = x_new
    } else {
      lb = y
    }
  }
  else {
    x_current = x_new
    
  }
  xs[i] = x_current
}