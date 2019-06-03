
#####################################################################
#continuous model: Obtaining drift matrix with 1-month interval;
#####################################################################

library(ctsem)
library(parallel)
t <- 250
ntime <- 20
intval_1 <- c(1,3,6)
cond.idx <- as.numeric(commandArgs(trailingOnly=TRUE))
intval <- intval_1[cond.idx]
index.x <- seq(1,ntime*intval,by=intval)
index.m <- seq((t+1),(t+ntime*intval),by=intval)
index.y <- seq((2*t+1),(2*t+ntime*intval),by=intval)

# save to results files
N <- 1000
cond <- 32 #n=1000
B <- 1000 #bootstrap replications
resfile.est <- paste0('/.../',N,'cestimate',intval,'monci/') # save parameters estimates;
dir.create(resfile.est)


runone <- function(rep) {
  filename <- paste('/.../',cond,'data',rep,'.txt', sep='')
  dset <- read.table(filename)[,-1]
  x <- dset[,index.x]
  m <- dset[,index.m]
  y <- dset[,index.y]
  Data <- Name.var <- NULL
  for (i in 1:ntime) {
    Data <- cbind(Data,x[,i],m[,i],y[,i])
    Name.var <- c(Name.var,paste0('Y1_T',i-1),paste0('Y2_T',i-1),paste0('Y3_T',i-1))
  }
  tintval <- matrix(intval,N,ntime-1)
  Name.tintval <- paste0('dT',1:(ntime-1))
  Data <- data.frame(cbind(Data,tintval))
  names(Data) <- c(Name.var,Name.tintval)
  model <- ctModel(n.latent = 3, n.manifest = 3, Tpoints = 20,MANIFESTVAR=diag(0,3),
                   DRIFT = matrix(c("x0",0,0,"a","m0",0,"c","b","y0"),3,3,byrow = T),
                   LAMBDA = diag(3))
  b <- 1
  lp <- 1
  specdir.est <- specind.est <- NULL
  oind1.est <- oind2.est <- oind6.est <- oind12.est <- NULL
  odir1.est <- odir2.est <- odir6.est <- odir12.est <- NULL
  #bootstrap estimates;
  while (b <= B) {
    Nb <- sample(1:N,replace = T)
    xb <- dset[Nb,index.x]
    mb <- dset[Nb,index.m]
    yb <- dset[Nb,index.y]
    Datab <- NULL
    for (i in 1:ntime) {
      Datab <- cbind(Datab,xb[,i],mb[,i],yb[,i])
    }
    Datab <- data.frame(cbind(Datab,tintval))
    names(Datab) <- c(Name.var,Name.tintval)
    fitb <- try(suppressMessages(ctFit(Datab,dataform = "wide", ctmodelobj = model)))
    if (class(fitb)=="try-error") next
    if (class(fitb)!="try-error") {
      Driftb <- try(summary(fitb)$DRIFT)
      if (class(Driftb)!="try-error") {
        Driftb.ind <- Driftb
        Driftb.dir <- Driftb
        a.est <- expm(Driftb*intval)[2,1]
        b.est <- expm(Driftb*intval)[3,2]
        specdir.est <- c(specdir.est,expm(Driftb*intval)[3,1])
        specind.est <- c(specind.est,a.est*b.est)
        Driftb.ind[3,1] <- 0 
        Driftb.dir[2,1] <- 0
        Driftb.dir[3,2] <- 0
        oind1.est <- c(oind1.est,expm(Driftb.ind)[3,1])
        oind2.est <- c(oind2.est,expm(Driftb.ind*2)[3,1])
        oind6.est <- c(oind6.est,expm(Driftb.ind*6)[3,1])
        oind12.est <- c(oind12.est,expm(Driftb.ind*12)[3,1])
        odir1.est <- c(odir1.est,expm(Driftb.dir*.5)[3,1])
        odir2.est <- c(odir2.est,expm(Driftb.dir)[3,1])
        odir6.est <- c(odir6.est,expm(Driftb.dir*3)[3,1])
        odir12.est <- c(odir12.est,expm(Driftb.dir*6)[3,1])
        b <- b+1
      }
    }
    lp <- lp+1
    if (lp > 5000) break
  }
  ci.specdir <- quantile(specdir.est,probs = c(.025,.975))
  ci.specind <- quantile(specind.est,probs = c(.025,.975))
  ci.oind1 <- quantile(oind1.est,probs = c(.025,.975))
  ci.oind2 <- quantile(oind2.est,probs = c(.025,.975))
  ci.oind6 <- quantile(oind6.est,probs = c(.025,.975))
  ci.oind12 <- quantile(oind12.est,probs = c(.025,.975))
  ci.odir1 <- quantile(odir1.est,probs = c(.025,.975))
  ci.odir2 <- quantile(odir2.est,probs = c(.025,.975))
  ci.odir6 <- quantile(odir6.est,probs = c(.025,.975))
  ci.odir12 <- quantile(odir12.est,probs = c(.025,.975))
  ci.all <- rbind(ci.specind,ci.oind1,ci.oind2,ci.oind6,ci.oind12,ci.specdir,ci.odir1,ci.odir2,ci.odir6,ci.odir12)
  filename <- paste0(resfile.est,'CI95_',rep,'.Rdata')
  save(ci.all,file= filename)
}
out <- mclapply(1001:1500, runone, mc.cores = 12)
