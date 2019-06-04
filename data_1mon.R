
## data generation: data with 1-month time interval ##

library('mnormt')

gendata.study1 <-function(i,a,b,C,N) {

# for generating initial values for x,m,y
rho1 <- 0
sigma <- matrix(c(1,rho1,rho1,rho1,1,rho1,rho1,rho1,1), 3, 3)
mu <- c(0,0,0)
dset1 <- rmnorm(N, mu, sigma)

# assign saving space for x,m,y
x <- matrix(NA,nrow = N,ncol = 1151)
m <- matrix(NA,nrow = N,ncol = 1151)
y <- matrix(NA,nrow = N,ncol = 1151)
r.xm <- r.xy <- r.my <- NULL
x[,1] <- dset1[,1]
m[,1] <- dset1[,2]
y[,1] <- dset1[,3]

for (j in 1:1150) {
   
   ex <- rnorm(N, mean = 0, sd = sqrt(phix))
   em <- rnorm(N, mean = 0, sd = sqrt(phim))
   ey <- rnorm(N, mean = 0, sd = sqrt(phiy))
   x[,j+1] <- x0 * x[,j] + ex
   m[,j+1] <- a * x[,j] + m0 * m[,j] + em
   y[,j+1] <- b * m[,j] + y0 * y[,j] + C*x[,j] + ey
   r.xm <- c(r.xm,cor(x[,j+1], m[,j+1]))
   r.xy <- c(r.xy,cor(x[,j+1], y[,j+1]))
   r.my <- c(r.my,cor(m[,j+1], y[,j+1]))

}

# this is the generated data matrix with 20 time points; discard the first 900 iterations
#dset.1 <- round(cbind(1:N,x[,900:(900+t-1)],m[,900:(900+t-1)],y[,900:(900+t-1)]),5)
dset <- round(cbind(x[,900:(900+t-1)],m[,900:(900+t-1)],y[,900:(900+t-1)]),5)
return(dset)

}

t <- 250 # number of time points for the data with 1 month time interval;

# autoregressive parameters;
x0 <- .7
m0 <- .6
y0 <- .5

# residual variances;
phix <- .9
phim <- .4
phiy <- .1

# parameters a, b, and c'
a <- c(0,.5)
b <- c(0,.4)
C <- c(0,-.1)
N <- c(100,200,400,1000)
cond <- expand.grid(a,b,C,N)
cond.idx <- as.numeric(commandArgs(trailingOnly=TRUE))
a <- cond[cond.idx,1]
b <- cond[cond.idx,2]
C <- cond[cond.idx,3]
N <- cond[cond.idx,4]
# replications to generate data
for (i in 1001:1500){
  dset = gendata.study1(i,a,b,C,N)
  filename <- paste('/.../cond',cond.idx,'data',i,'.txt', sep='')
  write.table(dset,filename, col.names=F,row.names=F)			
}
