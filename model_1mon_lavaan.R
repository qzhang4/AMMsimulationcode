
#####################################################################
#Obtaining parameter estimates, overall indirect effects, and fit indexes
#for the model with 1-month interval;
#####################################################################
library(lavaan)
library(parallel)

# get the condition for the current estimation;
a <- c(0,.5)
b <- c(0,.4)
C <- c(0,-.1)
N <- c(100,200,400,1000)
cond <- expand.grid(a,b,C,N)
cond.idx <- as.numeric(commandArgs(trailingOnly=TRUE))

# save to results files
resfile.est <- paste0('/.../estimate1mon/','cond',1000+cond.idx, '/') # save parameters estimates;
dir.create(resfile.est)
resfile.fit <- paste0('/.../fit1mon/','cond',1000+cond.idx, '/') # save fit indexes;
dir.create(resfile.fit)

## model estimation; ##
runone <- function(rep) {
  filename <- paste('/.../cond',cond.idx,'data',rep,'.txt', sep='')
  dset <- read.table(filename)
  prefixx <- "x"
  prefixm <- "m"
  prefixy <- "y"
  suffix <- seq(1:250)
  namesx <- paste(prefixx,suffix,sep = "")
  namesm <- paste(prefixm,suffix,sep = "")
  namesy <- paste(prefixy,suffix,sep = "")
  colnames(dset) <- c(namesx,namesm,namesy)
  model_1mon <- '
  x2 ~ dx*x1
  x3 ~ dx*x2
  x4 ~ dx*x3
  x5 ~ dx*x4
  x6 ~ dx*x5
  x7 ~ dx*x6
  x8 ~ dx*x7
  x9 ~ dx*x8
  x10 ~ dx*x9
  x11 ~ dx*x10
  x12 ~ dx*x11
  x13 ~ dx*x12
  x14 ~ dx*x13
  x15 ~ dx*x14
  x16 ~ dx*x15
  x17 ~ dx*x16
  x18 ~ dx*x17
  x19 ~ dx*x18
  x20 ~ dx*x19
  m2 ~ dm*m1+a*x1
  m3 ~ dm*m2+a*x2
  m4 ~ dm*m3+a*x3
  m5 ~ dm*m4+a*x4
  m6 ~ dm*m5+a*x5
  m7 ~ dm*m6+a*x6
  m8 ~ dm*m7+a*x7
  m9 ~ dm*m8+a*x8
  m10 ~ dm*m9+a*x9
  m11 ~ dm*m10+a*x10
  m12 ~ dm*m11+a*x11
  m13 ~ dm*m12+a*x12
  m14 ~ dm*m13+a*x13
  m15 ~ dm*m14+a*x14
  m16 ~ dm*m15+a*x15
  m17 ~ dm*m16+a*x16
  m18 ~ dm*m17+a*x17
  m19 ~ dm*m18+a*x18
  m20 ~ dm*m19+a*x19
  y2 ~ dy*y1+b*m1+c*x1
  y3 ~ dy*y2+b*m2+c*x2
  y4 ~ dy*y3+b*m3+c*x3
  y5 ~ dy*y4+b*m4+c*x4
  y6 ~ dy*y5+b*m5+c*x5
  y7 ~ dy*y6+b*m6+c*x6
  y8 ~ dy*y7+b*m7+c*x7
  y9 ~ dy*y8+b*m8+c*x8
  y10 ~ dy*y9+b*m9+c*x9
  y11 ~ dy*y10+b*m10+c*x10
  y12 ~ dy*y11+b*m11+c*x11
  y13 ~ dy*y12+b*m12+c*x12
  y14 ~ dy*y13+b*m13+c*x13
  y15 ~ dy*y14+b*m14+c*x14
  y16 ~ dy*y15+b*m15+c*x15
  y17 ~ dy*y16+b*m16+c*x16
  y18 ~ dy*y17+b*m17+c*x17
  y19 ~ dy*y18+b*m18+c*x18
  y20 ~ dy*y19+b*m19+c*x19
  x2 ~~ ex*x2
  x3 ~~ ex*x3
  x4 ~~ ex*x4
  x5 ~~ ex*x5
  x6 ~~ ex*x6
  x7 ~~ ex*x7
  x8 ~~ ex*x8
  x9 ~~ ex*x9
  x10 ~~ ex*x10
  x11 ~~ ex*x11
  x12 ~~ ex*x12
  x13 ~~ ex*x13
  x14 ~~ ex*x14
  x15 ~~ ex*x15
  x16 ~~ ex*x16
  x17 ~~ ex*x17
  x18 ~~ ex*x18
  x19 ~~ ex*x19
  x20 ~~ ex*x20
  m2 ~~ em*m2
  m3 ~~ em*m3
  m4 ~~ em*m4
  m5 ~~ em*m5
  m6 ~~ em*m6
  m7 ~~ em*m7
  m8 ~~ em*m8
  m9 ~~ em*m9
  m10 ~~ em*m10
  m11 ~~ em*m11
  m12 ~~ em*m12
  m13 ~~ em*m13
  m14 ~~ em*m14
  m15 ~~ em*m15
  m16 ~~ em*m16
  m17 ~~ em*m17
  m18 ~~ em*m18
  m19 ~~ em*m19
  m20 ~~ em*m20
  y2 ~~ ey*y2
  y3 ~~ ey*y3
  y4 ~~ ey*y4
  y5 ~~ ey*y5
  y6 ~~ ey*y6
  y7 ~~ ey*y7
  y8 ~~ ey*y8
  y9 ~~ ey*y9
  y10 ~~ ey*y10
  y11 ~~ ey*y11
  y12 ~~ ey*y12
  y13 ~~ ey*y13
  y14 ~~ ey*y14
  y15 ~~ ey*y15
  y16 ~~ ey*y16
  y17 ~~ ey*y17
  y18 ~~ ey*y18
  y19 ~~ ey*y19
  y20 ~~ ey*y20
  x1 ~~ covxm*m1
  m1 ~~ covmy*y1
  x1 ~~ covxy*y1
  specind := a*b
  specdir := c
  dx1 := dx^(.01)
  dm1 := dm^(.01)
  dy1 := dy^(.01)
  a1 := a*(dm^(.01)-dx^(.01))/(dm-dx)
  b1 := b*(dy^(.01)-dm^(.01))/(dy-dm)
  c1 := c*(dy^(.01)-dx^(.01))/(dy-dx)+a*b/(dy-dx)*((dy^(.01)-dm^(.01))/(dy-dm)-(dm^(.01)-dx^(.01))/(dm-dx))
  overallind11 := a1*b1/(dy1-dx1)*((dy1^(1/.01)-dm1^(1/.01))/(dy1-dm1)-(dm1^(1/.01)-dx1^(1/.01))/(dm1-dx1))
  overallind12 := a1*b1/(dy1-dx1)*((dy1^(2/.01)-dm1^(2/.01))/(dy1-dm1)-(dm1^(2/.01)-dx1^(2/.01))/(dm1-dx1))
  overallind13 := a1*b1/(dy1-dx1)*((dy1^(6/.01)-dm1^(6/.01))/(dy1-dm1)-(dm1^(6/.01)-dx1^(6/.01))/(dm1-dx1))
  overallind14 := a1*b1/(dy1-dx1)*((dy1^(12/.01)-dm1^(12/.01))/(dy1-dm1)-(dm1^(12/.01)-dx1^(12/.01))/(dm1-dx1))
  overalldir11 := c1*(dy1^(.5/.01)-dx1^(.5/.01))/(dy1-dx1)
  overalldir12 := c1*(dy1^(1/.01)-dx1^(1/.01))/(dy1-dx1)
  overalldir13 := c1*(dy1^(3/.01)-dx1^(3/.01))/(dy1-dx1)
  overalldir14 := c1*(dy1^(6/.01)-dx1^(6/.01))/(dy1-dx1)  
  dx2 := dx^(.034)
  dm2 := dm^(.034)
  dy2 := dy^(.034)
  a2 := a*(dm^(.034)-dx^(.034))/(dm-dx)
  b2 := b*(dy^(.034)-dm^(.034))/(dy-dm)
  c2 := c*(dy^(.034)-dx^(.034))/(dy-dx)+a*b/(dy-dx)*((dy^(.034)-dm^(.034))/(dy-dm)-(dm^(.034)-dx^(.034))/(dm-dx))
  overallind21 := a2*b2/(dy2-dx2)*((dy2^(1/.034)-dm2^(1/.034))/(dy2-dm2)-(dm2^(1/.034)-dx2^(1/.034))/(dm2-dx2))
  overallind22 := a2*b2/(dy2-dx2)*((dy2^(2/.034)-dm2^(2/.034))/(dy2-dm2)-(dm2^(2/.034)-dx2^(2/.034))/(dm2-dx2))
  overallind23 := a2*b2/(dy2-dx2)*((dy2^(6/.034)-dm2^(6/.034))/(dy2-dm2)-(dm2^(6/.034)-dx2^(6/.034))/(dm2-dx2))
  overallind24 := a2*b2/(dy2-dx2)*((dy2^(12/.034)-dm2^(12/.034))/(dy2-dm2)-(dm2^(12/.034)-dx2^(12/.034))/(dm2-dx2))
  overalldir21 := c2*(dy2^(.5/.034)-dx2^(.5/.034))/(dy2-dx2)
  overalldir22 := c2*(dy2^(1/.034)-dx2^(1/.034))/(dy2-dx2)
  overalldir23 := c2*(dy2^(3/.034)-dx2^(3/.034))/(dy2-dx2)
  overalldir24 := c2*(dy2^(6/.034)-dx2^(6/.034))/(dy2-dx2)
  dx3 := dx^(.1)
  dm3 := dm^(.1)
  dy3 := dy^(.1)
  a3 := a*(dm^(.1)-dx^(.1))/(dm-dx)
  b3 := b*(dy^(.1)-dm^(.1))/(dy-dm)
  c3 := c*(dy^(.1)-dx^(.1))/(dy-dx)+a*b/(dy-dx)*((dy^(.1)-dm^(.1))/(dy-dm)-(dm^(.1)-dx^(.1))/(dm-dx))
  overallind31 := a3*b3/(dy3-dx3)*((dy3^(1/.1)-dm3^(1/.1))/(dy3-dm3)-(dm3^(1/.1)-dx3^(1/.1))/(dm3-dx3))
  overallind32 := a3*b3/(dy3-dx3)*((dy3^(2/.1)-dm3^(2/.1))/(dy3-dm3)-(dm3^(2/.1)-dx3^(2/.1))/(dm3-dx3))
  overallind33 := a3*b3/(dy3-dx3)*((dy3^(6/.1)-dm3^(6/.1))/(dy3-dm3)-(dm3^(6/.1)-dx3^(6/.1))/(dm3-dx3))
  overallind34 := a3*b3/(dy3-dx3)*((dy3^(12/.1)-dm3^(12/.1))/(dy3-dm3)-(dm3^(12/.1)-dx3^(12/.1))/(dm3-dx3))
  overalldir31 := c3*(dy3^(.5/.1)-dx3^(.5/.1))/(dy3-dx3)
  overalldir32 := c3*(dy3^(1/.1)-dx3^(1/.1))/(dy3-dx3)
  overalldir33 := c3*(dy3^(3/.1)-dx3^(3/.1))/(dy3-dx3)
  overalldir34 := c3*(dy3^(6/.1)-dx3^(6/.1))/(dy3-dx3)  
'
  fit.1mon <- sem(model_1mon,data=dset,se="boot",bootstrap=1000)
  boot.fit.1mon <- parameterEstimates(fit.1mon, boot.ci.type="perc") # obtain percentile bootstrap intervals;
  estdup <- cbind(boot.fit.1mon$label,round(data.frame(cbind(boot.fit.1mon$est,boot.fit.1mon$ci.lower,boot.fit.1mon$ci.upper)),3))
  est <- estdup[!duplicated(estdup[,1]), 2:4] # point estimate and percentile bootstrap invervals for model parameters;
  indexes <- c(3,4,9,10,19,20,23,29) # chisq, df, cfi, tli,aic,bic,rmsea,srmr;
  fitindex <- fitmeasures(fit.1mon)[indexes] # obtain fit indexes of the model;
  est.file <- paste0(resfile.est, 'est',rep,'.Rdata') 
  save(est,file= est.file)  # save parameter estimates of the current replication to a file;
  fit.file <- paste0(resfile.fit, 'fit',rep,'.Rdata')  
  save(fitindex,file= fit.file) # save fit measures of the current replication to a file;
}

out <- mclapply(1001:1500, runone, mc.cores = 20)
