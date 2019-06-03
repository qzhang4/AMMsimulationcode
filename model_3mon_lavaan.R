
#####################################################################
#Obtaining paramter estimates, overall indirect effects, and fit indexes
#for the model with 3-month interval;
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
resfile.est <- paste0('/.../estimate3mon/','cond',1000+cond.idx, '/') # save parameters estimates;
dir.create(resfile.est)
resfile.fit <- paste0('/.../fit3mon/','cond',1000+cond.idx, '/') # save fit indexes;
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
  colnames(dset) <- c("N",namesx,namesm,namesy)
  model_3mon <- '
  x4 ~ dx*x1
  x7 ~ dx*x4
  x10 ~ dx*x7
  x13 ~ dx*x10
  x16 ~ dx*x13
  x19 ~ dx*x16
  x22 ~ dx*x19
  x25 ~ dx*x22
  x28 ~ dx*x25
  x31 ~ dx*x28
  x34 ~ dx*x31
  x37 ~ dx*x34
  x40 ~ dx*x37
  x43 ~ dx*x40
  x46 ~ dx*x43
  x49 ~ dx*x46
  x52 ~ dx*x49
  x55 ~ dx*x52
  x58 ~ dx*x55
  m4 ~ dm*m1+a*x1
  m7 ~ dm*m4+a*x4
  m10 ~ dm*m7+a*x7
  m13 ~ dm*m10+a*x10
  m16 ~ dm*m13+a*x13
  m19 ~ dm*m16+a*x16
  m22 ~ dm*m19+a*x19
  m25 ~ dm*m22+a*x22
  m28 ~ dm*m25+a*x25
  m31 ~ dm*m28+a*x28
  m34 ~ dm*m31+a*x31
  m37 ~ dm*m34+a*x34
  m40 ~ dm*m37+a*x37
  m43 ~ dm*m40+a*x40
  m46 ~ dm*m43+a*x43
  m49 ~ dm*m46+a*x46
  m52 ~ dm*m49+a*x49
  m55 ~ dm*m52+a*x52
  m58 ~ dm*m55+a*x55
  y4 ~ dy*y1+b*m1+c*x1
  y7 ~ dy*y4+b*m4+c*x4
  y10 ~ dy*y7+b*m7+c*x7
  y13 ~ dy*y10+b*m10+c*x10
  y16 ~ dy*y13+b*m13+c*x13
  y19 ~ dy*y16+b*m16+c*x16
  y22 ~ dy*y19+b*m19+c*x19
  y25 ~ dy*y22+b*m22+c*x22
  y28 ~ dy*y25+b*m25+c*x25
  y31 ~ dy*y28+b*m28+c*x28
  y34 ~ dy*y31+b*m31+c*x31
  y37 ~ dy*y34+b*m34+c*x34
  y40 ~ dy*y37+b*m37+c*x37
  y43 ~ dy*y40+b*m40+c*x40
  y46 ~ dy*y43+b*m43+c*x43
  y49 ~ dy*y46+b*m46+c*x46
  y52 ~ dy*y49+b*m49+c*x49
  y55 ~ dy*y52+b*m52+c*x52
  y58 ~ dy*y55+b*m55+c*x55
  x4 ~~ ex*x4
  x7 ~~ ex*x7
  x10 ~~ ex*x10
  x13 ~~ ex*x13
  x16 ~~ ex*x16
  x19 ~~ ex*x19
  x22 ~~ ex*x22
  x25 ~~ ex*x25
  x28 ~~ ex*x28
  x31 ~~ ex*x31
  x34 ~~ ex*x34
  x37 ~~ ex*x37
  x40 ~~ ex*x40
  x43 ~~ ex*x43
  x46 ~~ ex*x46
  x49 ~~ ex*x49
  x52 ~~ ex*x52
  x55 ~~ ex*x55
  x58 ~~ ex*x58
  m4 ~~ em*m4
  m7 ~~ em*m7
  m10 ~~ em*m10
  m13 ~~ em*m13
  m16 ~~ em*m16
  m19 ~~ em*m19
  m22 ~~ em*m22
  m25 ~~ em*m25
  m28 ~~ em*m28
  m31 ~~ em*m31
  m34 ~~ em*m34
  m37 ~~ em*m37
  m40 ~~ em*m40
  m43 ~~ em*m43
  m46 ~~ em*m46
  m49 ~~ em*m49
  m52 ~~ em*m52
  m55 ~~ em*m55
  m58 ~~ em*m58
  y4 ~~ ey*y4
  y7 ~~ ey*y7
  y10 ~~ ey*y10
  y13 ~~ ey*y13
  y16 ~~ ey*y16
  y19 ~~ ey*y19
  y22 ~~ ey*y22
  y25 ~~ ey*y25
  y28 ~~ ey*y28
  y31 ~~ ey*y31
  y34 ~~ ey*y34
  y37 ~~ ey*y37
  y40 ~~ ey*y40
  y43 ~~ ey*y43
  y46 ~~ ey*y46
  y49 ~~ ey*y49
  y52 ~~ ey*y52
  y55 ~~ ey*y55
  y58 ~~ ey*y58
  x1 ~~ covxm*m1
  m1 ~~ covmy*y1
  x1 ~~ covxy*y1
  specind := a*b
  specdir := c
  dx1 := dx^(.01/3)
  dm1 := dm^(.01/3)
  dy1 := dy^(.01/3)
  a1 := a*(dm^(.01/3)-dx^(.01/3))/(dm-dx)
  b1 := b*(dy^(.01/3)-dm^(.01/3))/(dy-dm)
  c1 := c*(dy^(.01/3)-dx^(.01/3))/(dy-dx)+a*b/(dy-dx)*((dy^(.01/3)-dm^(.01/3))/(dy-dm)-(dm^(.01/3)-dx^(.01/3))/(dm-dx))  
  overallind11 := a1*b1/(dy1-dx1)*((dy1^(1/.01)-dm1^(1/.01))/(dy1-dm1)-(dm1^(1/.01)-dx1^(1/.01))/(dm1-dx1))
  overallind12 := a1*b1/(dy1-dx1)*((dy1^(2/.01)-dm1^(2/.01))/(dy1-dm1)-(dm1^(2/.01)-dx1^(2/.01))/(dm1-dx1))
  overallind13 := a1*b1/(dy1-dx1)*((dy1^(6/.01)-dm1^(6/.01))/(dy1-dm1)-(dm1^(6/.01)-dx1^(6/.01))/(dm1-dx1))
  overallind14 := a1*b1/(dy1-dx1)*((dy1^(12/.01)-dm1^(12/.01))/(dy1-dm1)-(dm1^(12/.01)-dx1^(12/.01))/(dm1-dx1))
  overalldir11 := c1*(dy1^(.5/.01)-dx1^(.5/.01))/(dy1-dx1)
  overalldir12 := c1*(dy1^(1/.01)-dx1^(1/.01))/(dy1-dx1)
  overalldir13 := c1*(dy1^(3/.01)-dx1^(3/.01))/(dy1-dx1)
  overalldir14 := c1*(dy1^(6/.01)-dx1^(6/.01))/(dy1-dx1)  
  dx2 := dx^(.034/3)
  dm2 := dm^(.034/3)
  dy2 := dy^(.034/3)
  a2 := a*(dm^(.034/3)-dx^(.034/3))/(dm-dx)
  b2 := b*(dy^(.034/3)-dm^(.034/3))/(dy-dm)
  c2 := c*(dy^(.034/3)-dx^(.034/3))/(dy-dx)+a*b/(dy-dx)*((dy^(.034/3)-dm^(.034/3))/(dy-dm)-(dm^(.034/3)-dx^(.034/3))/(dm-dx))
  overallind21 := a2*b2/(dy2-dx2)*((dy2^(1/.034)-dm2^(1/.034))/(dy2-dm2)-(dm2^(1/.034)-dx2^(1/.034))/(dm2-dx2))
  overallind22 := a2*b2/(dy2-dx2)*((dy2^(2/.034)-dm2^(2/.034))/(dy2-dm2)-(dm2^(2/.034)-dx2^(2/.034))/(dm2-dx2))
  overallind23 := a2*b2/(dy2-dx2)*((dy2^(6/.034)-dm2^(6/.034))/(dy2-dm2)-(dm2^(6/.034)-dx2^(6/.034))/(dm2-dx2))
  overallind24 := a2*b2/(dy2-dx2)*((dy2^(12/.034)-dm2^(12/.034))/(dy2-dm2)-(dm2^(12/.034)-dx2^(12/.034))/(dm2-dx2))
  overalldir21 := c2*(dy2^(.5/.034)-dx2^(.5/.034))/(dy2-dx2)
  overalldir22 := c2*(dy2^(1/.034)-dx2^(1/.034))/(dy2-dx2)
  overalldir23 := c2*(dy2^(3/.034)-dx2^(3/.034))/(dy2-dx2)
  overalldir24 := c2*(dy2^(6/.034)-dx2^(6/.034))/(dy2-dx2)
  dx3 := dx^(.1/3)
  dm3 := dm^(.1/3)
  dy3 := dy^(.1/3)
  a3 := a*(dm^(.1/3)-dx^(.1/3))/(dm-dx)
  b3 := b*(dy^(.1/3)-dm^(.1/3))/(dy-dm)
  c3 := c*(dy^(.1/3)-dx^(.1/3))/(dy-dx)+a*b/(dy-dx)*((dy^(.1/3)-dm^(.1/3))/(dy-dm)-(dm^(.1/3)-dx^(.1/3))/(dm-dx))
  overallind31 := a3*b3/(dy3-dx3)*((dy3^(1/.1)-dm3^(1/.1))/(dy3-dm3)-(dm3^(1/.1)-dx3^(1/.1))/(dm3-dx3))
  overallind32 := a3*b3/(dy3-dx3)*((dy3^(2/.1)-dm3^(2/.1))/(dy3-dm3)-(dm3^(2/.1)-dx3^(2/.1))/(dm3-dx3))
  overallind33 := a3*b3/(dy3-dx3)*((dy3^(6/.1)-dm3^(6/.1))/(dy3-dm3)-(dm3^(6/.1)-dx3^(6/.1))/(dm3-dx3))
  overallind34 := a3*b3/(dy3-dx3)*((dy3^(12/.1)-dm3^(12/.1))/(dy3-dm3)-(dm3^(12/.1)-dx3^(12/.1))/(dm3-dx3))
  overalldir31 := c3*(dy3^(.5/.1)-dx3^(.5/.1))/(dy3-dx3)
  overalldir32 := c3*(dy3^(1/.1)-dx3^(1/.1))/(dy3-dx3)
  overalldir33 := c3*(dy3^(3/.1)-dx3^(3/.1))/(dy3-dx3)
  overalldir34 := c3*(dy3^(6/.1)-dx3^(6/.1))/(dy3-dx3)  
'
  fit.3mon <- sem(model_3mon,data=dset,se="boot",bootstrap=1000)
  boot.fit.3mon <- parameterEstimates(fit.3mon, boot.ci.type="perc") # obtain percentile bootstrap intervals;
  estdup <- cbind(boot.fit.3mon$label,round(data.frame(cbind(boot.fit.3mon$est,boot.fit.3mon$ci.lower,boot.fit.3mon$ci.upper)),3))
  est <- estdup[!duplicated(estdup[,1]), 2:4] # point estimate and percentile bootstrap invervals for model parameters;
  indexes <- c(3,4,9,10,19,20,23,29) # chisq, df, cfi, tli,aic,bic,rmsea,srmr;
  fitindex <- fitmeasures(fit.3mon)[indexes] # obtain fit indexes of the model;
  est.file <- paste0(resfile.est, 'est',rep,'.Rdata') 
  save(est,file= est.file)  # save parameter estimates of the current replication to a file;
  fit.file <- paste0(resfile.fit, 'fit',rep,'.Rdata')  
  save(fitindex,file= fit.file) # save fit measures of the current replication to a file;
}

out <- mclapply(1001:1500, runone, mc.cores = 20)
