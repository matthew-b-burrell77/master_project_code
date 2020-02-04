install.packages("pracma")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("Hmisc")
install.packages("dplyr")
library(dplyr)
library("lubridate")
library(ggplot2)
library(pracma)
library(Hmisc)

library(readxl)
data1 <- read_excel("~/EC 6691/Cleanned Data/data1.xlsx")
View(data1)

expecgap <- data1$`S&W`- data1$CPCYOY
inflaggap <- data1$CPI - data1$CPCYOY

TT = nrow(data1)
sigu = runif(1, 0, 10)
sig0 = runif(1, 0, 10)
sig1 = runif(1, 0, 10)
sig2 = runif(1, 0, 10)

beta0t = matrix(data = 0, nrow = TT, ncol = 1)
beta1t = matrix(data = 0, nrow = TT, ncol = 1)
beta2t = matrix(data = 0, nrow = TT, ncol = 1)

y = matrix(inflaggap, nrow = TT, ncol = 1)
x1 = matrix(expecgap, nrow = TT, ncol = 1)
x2 = matrix(data1$GAP, nrow = TT, ncol = 1)
x3 = matrix(data1$IPD, nrow = TT, ncol = 1)

for (i in 1:TT){
  if (i == 1){
    beta0t[i] = sqrt(sig0)+rnorm(1, mean = 0, sd = 1)
    beta1t[i] = sqrt(sig1)+rnorm(1, mean = 0, sd = 1)
    beta2t[i] = sqrt(sig2)+rnorm(1, mean = 0, sd = 1)
  }
  else{
    beta0t[i] = beta0t[i-1] + sqrt(sig0)*rnorm(1, mean = 0, sd = 1)
    beta1t[i] = beta1t[i-1] + sqrt(sig1)*rnorm(1, mean = 0, sd = 1)
    beta2t[i] = beta2t[i-1] + sqrt(sig2)*rnorm(1, mean = 0, sd = 1)
  }
}

V_PROB <- function(ERR, VAR){
  z=sqrt(VAR^2)
  VAL = (1/sqrt(2*pi*det(z)))*exp(-0.5*t(ERR)/(VAR)*ERR)
  return(VAL)
}

trans <- function(prmtr1){
  prmtr_new = matrix(data = 0, nrow = 1, ncol = 4)
  prmtr_new[1] = exp(-prmtr1[1])
  prmtr_new[2] = exp(-prmtr1[2])
  prmtr_new[3] = exp(-prmtr1[3])
  prmtr_new[4] = exp(-prmtr1[4])
  return(prmtr_new)
}

loglik <- function(prmtr1){
  LOGLIK = matrix(data = 0, nrow = TT, ncol = 1)
  
  PRMTR = trans(prmtr1)
  
  sigu_e = PRMTR[1,1]
  sig0_e = PRMTR[1,2]
  sig1_e = PRMTR[1,3]
  sig2_e = PRMTR[1,4]
 
  
  F_a = matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol = 3)
  
  Q = matrix(c(sig0_e,0,0,0,sig1_e,0,0,0,sig2_e), nrow = 3, ncol = 3)
  
  R = sigu_e
  
  B_d = nrow(Q)
  
  B_LL0 = matrix(data = 0, nrow = B_d, ncol = 1)
  
  plk = matrix(data = 0, nrow = B_d, ncol = B_d)
  plk[1,1] = 10000000
  plk[2,2] = 10000000
  plk[3,3] = 10000000
  
  P_LL0 = plk
  
  for (i in 1:TT){
    Q0 = matrix(c(sig0_e,0,0,0,sig1_e,0,0,0,sig2_e), nrow = 3, ncol = 3)
    
    F0 = matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol =3)
    
    H0 = matrix(c(x1[i,1], x2[i,1], x3[i,1]), nrow = 1, ncol = 3)
    
    # Kalman Filter
    
    B_TL00 = F0 %*% B_LL0
    P_TL00 = F0 %*% P_LL0 %*% t(F0) + Q0
    
    F_CAST00 = y[i] - H0 %*% B_TL00
    SS00 = H0 %*% P_TL00 %*% t(H0) + R
    
    A = P_TL00 %*% t(H0)
    B_TT00 = B_TL00 + mrdivide(A, SS00, pinv =TRUE) %*% F_CAST00
    P_TT00 = (diag(B_d) - mrdivide(A, SS00, pinv =TRUE) %*% H0 %*% P_TL00)
    
    PR_VL00 = V_PROB(F_CAST00,SS00)
    
    B_LL0 = B_TT00
    P_LL0 = P_TT00
    
    LOGLIK[i,1]= -log(PR_VL00)
  }
  los = 5
  
  LIKV = sum(LOGLIK[6:TT,1])
  return(LIKV)
}

prmtr0 = matrix(c(sigu,sig0,sig1,sig2), nrow = 1, ncol = 4)
B = optim(prmtr0,loglik, method = "Nelder-Mead" )

results = trans(B$par)

F_a = matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol = 3)

sigu_e = results[1]
sig0_e = results[2]
sig1_e = results[3]
sig2_e = results[4]

Q = matrix(c(sig0_e,0,0,0,sig1_e,0,0,0,sig2_e), nrow = 3, ncol = 3)

R = sigu_e

B_d = nrow(Q)

B_LL0 = matrix(data=0, nrow = B_d, ncol = 1)

plk = matrix(data=0, nrow = B_d, ncol = B_d)
plk[1,1] = 10000000
plk[2,2] = 10000000
plk[3,3] = 10000000

P_LL0 = plk

beta0t_e = matrix(data=0, nrow = TT, ncol = 1)
beta1t_e = matrix(data=0, nrow = TT, ncol = 1)
beta2t_e = matrix(data=0, nrow = TT, ncol = 1)

LOGLIK = matrix(data = 0, nrow = TT, ncol = 1)

for(i in 1:TT){
  
  Q0 = matrix(c(sig0_e,0,0,0,sig1_e,0,0,0,sig2_e), nrow = 3, ncol = 3)
  
  F0 = matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol =3)
  
  H0 = matrix(c(x1[i,1], x2[i,1], x3[i,1]), nrow = 1, ncol = 3)
  
  #Kalam Filter
  
  B_TL00 = F0 %*% B_LL0
  P_TL00 = F0 %*% P_LL0 %*% t(F0) + Q0
  
  F_CAST00 = y[i,] - H0 %*% B_TL00
  SS00 = H0 %*% P_TL00 %*% t(H0)+ R
  
  A = P_TL00 %*% t(H0)
  B_TT00 =  B_TL00 + mrdivide(A, SS00, pinv =TRUE) %*% F_CAST00
  P_TT00 = (diag(B_d) - mrdivide(A, SS00, pinv =TRUE) %*% H0 %*% P_TL00)
  
  PR_VL00 = V_PROB(F_CAST00,SS00)
  
  B_LL0 = B_TT00
  P_ll0 = P_TT00
  
  beta0t_e[i,1] = B_TT00[1,1]
  beta1t_e[i,1] = B_TT00[2,1]
  beta2t_e[i,1] = B_TT00[3,1]
  
  LOGLIK[i,1]= -log(PR_VL00)
}

plot(beta0t_e ~ data1$observation_date, col="blue" , las=2, type ="l", ann = FALSE, ylim = c(-0.5,2))
title(ylab = expression(theta))
title(main = expression(paste(theta," Over-Time")))
title(xlab = "Quarterly Frequency")

plot(beta1t_e ~ data1$observation_date, col="blue" , las = 2, type ="l", ann = FALSE, ylim = c(-1,1))
title(main="k Over-Time")
title(xlab = "Quarterly Frequency")
title(ylab = "k")

plot(beta2t_e ~ data1$observation_date, col="blue" , las=2, type ="l", ann = FALSE, ylim = c(-1,1))
title(main=expression(paste(gamma," Over-Time")))
title(xlab = "Quarterly Frequency")
title(ylab = expression(gamma))

plot(data1$CPI~data1$observation_date, type="l", ylab = "Inflation", xlab = "Quarterly Frequency")

betas.df = data.frame(data1$observation_date, beta0t_e, beta1t_e,beta2t_e)
betas.df <- betas.df %>%
  mutate(data1.observation_date = as.Date(data1.observation_date, format= "%y.%m.%d"))
       
recessions.df = read.table(textConnection(
  "Peak, Trough
  1960-04-01, 1961-04-01
  1969-10-01, 1971-01-01
  1973-10-01, 1975-01-01
  1980-01-01, 1980-07-01
  1981-07-01, 1982-10-01
  1990-07-01, 1991-01-01
  2001-01-01, 2001-10-01
  2007-10-01, 2009-04-01"), sep=',',
  colClasses=c('Date', 'Date'), header=TRUE)

datap = rbind(betas.df,recessions.df)

g=ggplot(data=betas.df) +
  geom_line(aes(data1.observation_date, beta0t_e),color="blue")+
  labs(title = expression(paste(theta," Over-Time")) , x ="Quarterly Frequency", y = expression(theta))+  theme_bw() 
g = g +  geom_rect(data=recessions.df, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='red', alpha=0.2)
plot(g)

c=ggplot(data=betas.df) +
  geom_line(aes(data1.observation_date, beta1t_e),color="blue")+
  labs(title ="k Over-Time" , x ="Quarterly Frequency", y = "k")+  theme_bw() 
c = c +  geom_rect(data=recessions.df, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='red', alpha=0.2)
plot(c)

j=ggplot(data=betas.df) +
  geom_line(aes(data1.observation_date, beta2t_e),color="blue")+
  labs(title = expression(paste(gamma," Over-Time")) , x ="Quarterly Frequency", y = expression(gamma))+  theme_bw() 
j = j +  geom_rect(data=recessions.df, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='red', alpha=0.2)
plot(j)

plot(x2~y)
