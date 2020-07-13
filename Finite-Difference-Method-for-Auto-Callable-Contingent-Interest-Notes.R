#Project 2
#input
S0 <- 195.09
K <- 0.78*S0
n <- 28
T <- 459/365
sigma <- 0.24
r <-  0.0254562 
r1 <- 0.0249189
r2 <- 0.0254476 
div <- 0.004
c <- 1000*0.020375

jmax = 300
imax = 459*n

Smin = 0
Smax = 3*S0

delt = T/imax
dels = (Smax - Smin)/jmax

print(delt<1/(sigma*jmax)^2)

S = matrix(0, nrow = jmax+1, ncol = imax+1)

divdate <- c(50*n+1, 141*n+1, 232*n+1, 323*n+1, 414*n+1)

# stock grid
k <- 0
for (i in 1:(imax+1)) {
  if (is.element(i, divdate)) {k<-k+1}
  for (j in 1:(jmax+1)) {
    S[j,i]=dels*(j-1)*(1-div)^k
  }
}

#Trigger grid
VT = matrix(0, nrow = jmax+1, ncol = imax+1)
A = vector(length = jmax+1)
B = vector(length = jmax+1)
C = vector(length = jmax+1)

#payoff at maturity
for (j in 1:(jmax+1)) {
  if(S[j,imax+1] >= S0){
    VT[j, imax+1]=(1000+c)*exp(-r2*3/365)
  } else {
    if(S[j,imax+1] >= K){
      VT[j, imax+1]=((1000*S[j,imax+1]/S0)+c)*exp(-r2*3/365)
    } else{
      VT[j, imax+1]=((1000*S[j,imax+1]/S0))*exp(-r2*3/365)
    }
  }
  A[j] = (0.5*sigma^2*(j-1)^2+0.5*(r-div)*(j-1))*delt
  B[j] = 1-r*delt-sigma^2*(j-1)^2*delt
  C[j] = (0.5*sigma^2*(j-1)^2-0.5*(r-div)*(j-1))*delt
}


#lower boundry
for (i in 1:(imax+1)) {
  VT[1, i]=0
}

#upper boundry
reviewdates <- c(92*n+1, 186*n+1, 277*n+1, 368*n+1, 459*n+1)
paymentdates <- c(97*n+1, 189*n+1, 281*n+1, 371*n+1, 462*n+1)

for(i in (imax+1):1){
  if(i > reviewdates[4]) {
    VT[jmax+1, i]=(1000+c)*exp(-r*(imax+1-i)*delt)*exp(-r2*(3/365))
  } else if (i > reviewdates[3]) {
    VT[jmax+1, i]=(1000+c)*exp(-r*(reviewdates[4]-i)*delt)*exp(-r*(3/365))
  } else if (i > reviewdates[2]) {
    VT[jmax+1, i]=(1000+c)*exp(-r*(reviewdates[3]-i)*delt)*exp(-r*(4/365))
  } else if (i > reviewdates[1]) {
    VT[jmax+1, i]=(1000+c)*exp(-r*(reviewdates[2]-i)*delt)*exp(-r*(3/365))
  } else if (i <= reviewdates[1]) {
    VT[jmax+1, i]=(1000+c)*exp(-r*(reviewdates[1]-i)*delt)*exp(-r*(5/365))
  }
}  


#elsewhere
l <- 0

for (i in imax:1) {
    if (is.element(i, reviewdates)){
      l <- which(reviewdates == i)
      for (j in 2:jmax){
      if(S[j,i] > S0) {
        VT[j,i]=(1000+c)*exp(-r*(paymentdates[l]-reviewdates[l])/365)
      } else {
        if(S[j,i] >= K){
          VT[j,i] = (A[j]*VT[j+1,i+1]+B[j]*VT[j,i+1]+C[j]*VT[j-1,i+1])+c*exp(-r*(paymentdates[l]-reviewdates[l])/365)
        } else {
          VT[j,i] = A[j]*VT[j+1,i+1]+B[j]*VT[j,i+1]+C[j]*VT[j-1,i+1]
        }
      }
    } 
    }  else {
      for (j in 2:jmax){
      VT[j,i] = A[j]*VT[j+1,i+1]+B[j]*VT[j,i+1]+C[j]*VT[j-1,i+1]
      }
    }
}

##No Trigger grid

#Maturity
V = matrix(0, nrow = jmax+1, ncol = imax+1)

for (j in 1:(jmax+1)){
  if (S[j,imax+1] >= K){
    V[j,imax+1] = (1000+c)*exp(-r2*3/365)
  } else {
    V[j,imax+1] = VT[j,imax+1]
  }
}

##Upper Boundary
for (i in (imax+1):1){
  V[jmax+1,i]=VT[jmax+1,i]
}

##elsewhere

a <- 0

for (i in imax:1) {
  for (j in 2:jmax){
    if (S[j,i] >= K) {
      if (is.element(i, reviewdates)){
        a <- which(reviewdates == i)
        for (j in 2:jmax){
          if(S[j,i] >= S0) {
            V[j,i]=(1000+c)*exp(-r*(paymentdates[a]-reviewdates[a])/365)
          } else {
            V[j,i] = (A[j]*V[j+1,i+1]+B[j]*V[j,i+1]+C[j]*V[j-1,i+1])+c*exp(-r*(paymentdates[a]-reviewdates[a])/365)
          }
        }
      } else {
        V[j,i] = A[j]*V[j+1,i+1]+B[j]*V[j,i+1]+C[j]*V[j-1,i+1]
      }
    } else {
      V[j,i] = VT[j,i]
    }
  }
}

V[101,1]
