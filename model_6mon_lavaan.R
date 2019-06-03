
#####################################################################
#Obtaining paramter estimates, overall indirect effects, and fit indexes
#for the model with 6-month interval;
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
resfile.est <- paste0('/.../estimate6mon/','cond',1000+cond.idx, '/') # save parameters estimates;
dir.create(resfile.est)
resfile.fit <- paste0('/.../fit6mon/','cond',1000+cond.idx, '/') # save fit indexes;
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
  model_6mon <- '
  x7 ~ dx*x1
  x13 ~ dx*x7
  x19 ~ dx*x13
  x25 ~ dx*x19
  x31 ~ dx*x25
  x37 ~ dx*x31
  x43 ~ dx*x37
  x49 ~ dx*x43
  x55 ~ dx*x49
  x61 ~ dx*x55
  x67 ~ dx*x61
  x73 ~ dx*x67
  x79 ~ dx*x73
  x85 ~ dx*x79
  x91 ~ dx*x85
  x97 ~ dx*x91
  x103 ~ dx*x97
  x109 ~ dx*x103
  x115 ~ dx*x109
  m7 ~ dm*m1+a*x1
  m13 ~ dm*m7+a*x7
  m19 ~ dm*m13+a*x13
  m25 ~ dm*m19+a*x19
  m31 ~ dm*m25+a*x25
  m37 ~ dm*m31+a*x31
  m43 ~ dm*m37+a*x37
  m49 ~ dm*m43+a*x43
  m55 ~ dm*m49+a*x49
  m61 ~ dm*m55+a*x55
  m67 ~ dm*m61+a*x61
  m73 ~ dm*m67+a*x67
  m79 ~ dm*m73+a*x73
  m85 ~ dm*m79+a*x79
  m91 ~ dm*m85+a*x85
  m97 ~ dm*m91+a*x91
  m103 ~ dm*m97+a*x97
  m109 ~ dm*m103+a*x103
  m115 ~ dm*m109+a*x109
  y7 ~ dy*y1+b*m1+c*x1
  y13 ~ dy*y7+b*m7+c*x7
  y19 ~ dy*y13+b*m13+c*x13
  y25 ~ dy*y19+b*m19+c*x19
  y31 ~ dy*y25+b*m25+c*x25
  y37 ~ dy*y31+b*m31+c*x31
  y43 ~ dy*y37+b*m37+c*x37
  y49 ~ dy*y43+b*m43+c*x43
  y55 ~ dy*y49+b*m49+c*x49
  y61 ~ dy*y55+b*m55+c*x55
  y67 ~ dy*y61+b*m61+c*x61
  y73 ~ dy*y67+b*m67+c*x67
  y79 ~ dy*y73+b*m73+c*x73
  y85 ~ dy*y79+b*m79+c*x79
  y91 ~ dy*y85+b*m85+c*x85
  y97 ~ dy*y91+b*m91+c*x91
  y103 ~ dy*y97+b*m97+c*x97
  y109 ~ dy*y103+b*m103+c*x103
  y115 ~ dy*y109+b*m109+c*x109
  x7 ~~ ex*x7
  x13 ~~ ex*x13
  x19 ~~ ex*x19
  x25 ~~ ex*x25
  x31 ~~ ex*x31
  x37 ~~ ex*x37
  x43 ~~ ex*x43
  x49 ~~ ex*x49
  x55 ~~ ex*x55
  x61 ~~ ex*x61
  x67 ~~ ex*x67
  x73 ~~ ex*x73
  x79 ~~ ex*x79
  x85 ~~ ex*x85
  x91 ~~ ex*x91
  x97 ~~ ex*x97
  x103 ~~ ex*x103
  x109 ~~ ex*x109
  x115 ~~ ex*x115
  m7 ~~ em*m7
  m13 ~~ em*m13
  m19 ~~ em*m19
  m25 ~~ em*m25
  m31 ~~ em*m31
  m37 ~~ em*m37
  m43 ~~ em*m43
  m49 ~~ em*m49
  m55 ~~ em*m55
  m61 ~~ em*m61
  m67 ~~ em*m67
  m73 ~~ em*m73
  m79 ~~ em*m79
  m85 ~~ em*m85
  m91 ~~ em*m91
  m97 ~~ em*m97
  m103 ~~ em*m103
  m109 ~~ em*m109
  m115 ~~ em*m115
  y7 ~~ ey*y7
  y13 ~~ ey*y13
  y19 ~~ ey*y19
  y25 ~~ ey*y25
  y31 ~~ ey*y31
  y37 ~~ ey*y37
  y43 ~~ ey*y43
  y49 ~~ ey*y49
  y55 ~~ ey*y55
  y61 ~~ ey*y61
  y67 ~~ ey*y67
  y73 ~~ ey*y73
  y79 ~~ ey*y79
  y85 ~~ ey*y85
  y91 ~~ ey*y91
  y97 ~~ ey*y97
  y103 ~~ ey*y103
  y109 ~~ ey*y109
  y115 ~~ ey*y115
  x1 ~~ covxm*m1
  m1 ~~ covmy*y1
  x1 ~~ covxy*y1
  specind := a*b
  specdir := c
  dx1 := dx^(.01/6)
  dm1 := dm^(.01/6)
  dy1 := dy^(.01/6)
  a1 := a*(dm^(.01/6)-dx^(.01/6))/(dm-dx)
  b1 := b*(dy^(.01/6)-dm^(.01/6))/(dy-dm)
  c1 := c*(dy^(.01/6)-dx^(.01/6))/(dy-dx)+a*b/(dy-dx)*((dy^(.01/6)-dm^(.01/6))/(dy-dm)-(dm^(.01/6)-dx^(.01/6))/(dm-dx))  
  overallind11 := a1*b1/(dy1-dx1)*((dy1^(1/.01)-dm1^(1/.01))/(dy1-dm1)-(dm1^(1/.01)-dx1^(1/.01))/(dm1-dx1))
  overallind12 := a1*b1/(dy1-dx1)*((dy1^(2/.01)-dm1^(2/.01))/(dy1-dm1)-(dm1^(2/.01)-dx1^(2/.01))/(dm1-dx1))
  overallind13 := a1*b1/(dy1-dx1)*((dy1^(6/.01)-dm1^(6/.01))/(dy1-dm1)-(dm1^(6/.01)-dx1^(6/.01))/(dm1-dx1))
  overallind14 := a1*b1/(dy1-dx1)*((dy1^(12/.01)-dm1^(12/.01))/(dy1-dm1)-(dm1^(12/.01)-dx1^(12/.01))/(dm1-dx1))
  overalldir11 := c1*(dy1^(.5/.01)-dx1^(.5/.01))/(dy1-dx1)
  overalldir12 := c1*(dy1^(1/.01)-dx1^(1/.01))/(dy1-dx1)
  overalldir13 := c1*(dy1^(3/.01)-dx1^(3/.01))/(dy1-dx1)
  overalldir14 := c1*(dy1^(6/.01)-dx1^(6/.01))/(dy1-dx1)  
  dx2 := dx^(.034/6)
  dm2 := dm^(.034/6)
  dy2 := dy^(.034/6)
  a2 := a*(dm^(.034/6)-dx^(.034/6))/(dm-dx)
  b2 := b*(dy^(.034/6)-dm^(.034/6))/(dy-dm)
  c2 := c*(dy^(.034/6)-dx^(.034/6))/(dy-dx)+a*b/(dy-dx)*((dy^(.034/6)-dm^(.034/6))/(dy-dm)-(dm^(.034/6)-dx^(.034/6))/(dm-dx))
  overallind21 := a2*b2/(dy2-dx2)*((dy2^(1/.034)-dm2^(1/.034))/(dy2-dm2)-(dm2^(1/.034)-dx2^(1/.034))/(dm2-dx2))
  overallind22 := a2*b2/(dy2-dx2)*((dy2^(2/.034)-dm2^(2/.034))/(dy2-dm2)-(dm2^(2/.034)-dx2^(2/.034))/(dm2-dx2))
  overallind23 := a2*b2/(dy2-dx2)*((dy2^(6/.034)-dm2^(6/.034))/(dy2-dm2)-(dm2^(6/.034)-dx2^(6/.034))/(dm2-dx2))
  overallind24 := a2*b2/(dy2-dx2)*((dy2^(12/.034)-dm2^(12/.034))/(dy2-dm2)-(dm2^(12/.034)-dx2^(12/.034))/(dm2-dx2))
  overalldir21 := c2*(dy2^(.5/.034)-dx2^(.5/.034))/(dy2-dx2)
  overalldir22 := c2*(dy2^(1/.034)-dx2^(1/.034))/(dy2-dx2)
  overalldir23 := c2*(dy2^(3/.034)-dx2^(3/.034))/(dy2-dx2)
  overalldir24 := c2*(dy2^(6/.034)-dx2^(6/.034))/(dy2-dx2)
  dx3 := dx^(.1/6)
  dm3 := dm^(.1/6)
  dy3 := dy^(.1/6)
  a3 := a*(dm^(.1/6)-dx^(.1/6))/(dm-dx)
  b3 := b*(dy^(.1/6)-dm^(.1/6))/(dy-dm)
  c3 := c*(dy^(.1/6)-dx^(.1/6))/(dy-dx)+a*b/(dy-dx)*((dy^(.1/6)-dm^(.1/6))/(dy-dm)-(dm^(.1/6)-dx^(.1/6))/(dm-dx))
  overallind31 := a3*b3/(dy3-dx3)*((dy3^(1/.1)-dm3^(1/.1))/(dy3-dm3)-(dm3^(1/.1)-dx3^(1/.1))/(dm3-dx3))
  overallind32 := a3*b3/(dy3-dx3)*((dy3^(2/.1)-dm3^(2/.1))/(dy3-dm3)-(dm3^(2/.1)-dx3^(2/.1))/(dm3-dx3))
  overallind33 := a3*b3/(dy3-dx3)*((dy3^(6/.1)-dm3^(6/.1))/(dy3-dm3)-(dm3^(6/.1)-dx3^(6/.1))/(dm3-dx3))
  overallind34 := a3*b3/(dy3-dx3)*((dy3^(12/.1)-dm3^(12/.1))/(dy3-dm3)-(dm3^(12/.1)-dx3^(12/.1))/(dm3-dx3))
  overalldir31 := c3*(dy3^(.5/.1)-dx3^(.5/.1))/(dy3-dx3)
  overalldir32 := c3*(dy3^(1/.1)-dx3^(1/.1))/(dy3-dx3)
  overalldir33 := c3*(dy3^(3/.1)-dx3^(3/.1))/(dy3-dx3)
  overalldir34 := c3*(dy3^(6/.1)-dx3^(6/.1))/(dy3-dx3)  
'
  fit.6mon <- sem(model_6mon,data=dset,se="boot",bootstrap=1000)
  boot.fit.6mon <- parameterEstimates(fit.6mon, boot.ci.type="perc") # obtain percentile bootstrap intervals;
  estdup <- cbind(boot.fit.6mon$label,round(data.frame(cbind(boot.fit.6mon$est,boot.fit.6mon$ci.lower,boot.fit.6mon$ci.upper)),3))
  est <- estdup[!duplicated(estdup[,1]), 2:4] # point estimate and percentile bootstrap invervals for model parameters;
  indexes <- c(3,4,9,10,19,20,23,29) # chisq, df, cfi, tli,aic,bic,rmsea,srmr;
  fitindex <- fitmeasures(fit.6mon)[indexes] # obtain fit indexes of the model;
  est.file <- paste0(resfile.est, 'est',rep,'.Rdata') 
  save(est,file= est.file)  # save parameter estimates of the current replication to a file;
  fit.file <- paste0(resfile.fit, 'fit',rep,'.Rdata')  
  save(fitindex,file= fit.file) # save fit measures of the current replication to a file;
}

out <- mclapply(1001:1500, runone, mc.cores = 20)
