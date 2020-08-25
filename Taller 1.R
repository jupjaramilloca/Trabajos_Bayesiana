## Ejercicio número 1. 


library(BayesDA)
data(cow)


ybar <- mean(cow$protein)
s <- sd(cow$protein)


##a) ----------------------------------------------------------------------------

library(BayesDA)
data(cow)
y_bar <- mean(cow$protein)
s2 <- var(cow$protein)
n <- 50
s0 <- 5
k0 <- 5
mu0 <- 3
vn <- 50.2
kn <- 55
s2_0 <- 2
v0 <- 0.2 
sigma2_n <- ((v0*s2_0)+((n-1)*s2)+((n*k0)*(y_bar-mu0)^2)/kn)/vn
sigma2_n  
  


alpha <- vn/2
beta <- sigma2_n*vn/2
sigma2_n
alpha
beta


library(pscl)
LS <- qinvgamma(p = 0.025,shape = alpha,rate = beta)
LI <- qinvgamma(p = 0.975,shape = alpha,rate = beta)


intervalo_sigma1 <- c(LI,LS) #Intervalo de credibilidad para sigma
intervalo_sigma1


mun <- ((mu0*k0)+(n*y_bar))/kn
param2 <- sigma2_n/(k0+n)
param2

LI <- y_bar-qt(0.975,55.2)*sqrt(param2)
LS <- y_bar+qt(0.975,55.2)*sqrt(param2)

intervalo_theta1 <- c(LI,LS)
intervalo_theta1              ##Interalo de credibilida para theta

## b)--------------------------------------------------------------------------
#intervalo de credibilidad para theta

LS=ybar+qt(0.975,length(cow$protein)-1)*(s/sqrt(length(cow$protein)))
LI=ybar-qt(0.975,length(cow$protein)-1)*s/sqrt(length(cow$protein))

intervalo_theta <- c(LI,LS) 
intervalo_theta

#Intervalo de credibilidad para sigma^2
library(invgamma)
li  <-qinvgamma(p = 0.05,shape = (length(cow$protein)-1)/2,rate = ((length(cow$protein)-1)*s^2)/2)
ls <- qinvgamma(p = 0.975,shape = (length(cow$protein)-1)/2,rate = ((length(cow$protein)-1)*s^2)/2)
intervalos_sigma <-c(li,ls)
intervalos_sigma



## C)---------------------------------------------------------------------------
#Estimar la varianza con media conocida, Implica que la distribución del muestreo
# es una distribución normal con media conocida y Varianza desconocida. 


library(pscl)
prote <- cow$protein    # Variable con distribución normal 
mu=3                    # Se conoce la media 
n=length(cow$protein)    
nu=sum((prote-mu)^2)/n  # Es como una especie de varianza pero con "mu" conocida.



    # Recordemos que la distribución a priori **conjugada** para \sigma^2 es la gamma inversa
alpha=0.1   # se conoce 
beta=0.2    # Se conoce

nu0=2*alpha       
sigma0=2*beta/nu0

a2=(n+nu0)/2
b2=(n*nu+nu0*sigma0)/2

## Intervalo de credibilidad, ya conociendo los parametros de la gamma inversa a2 y b2



qigamma(0.025,a2,b2)## intervalo de credibilidad superior 

qigamma(0.975,a2,b2)## intervalo de credibilidad inferior




## d)---------------------------------------------------------------------------

## Estimar a theta con sigma conocida, Implica que la distribucion muestral provenga 
## de una distribución normal, con una distribución a priori con parametros conocidos

n <- length(cow$protein)

s <- 0.05 # Este paremetro se conoce previamente

Mu0 <- 3   # Este parametro se conoce y es para la distribución a priori
T02 <- 0.02 # Este parametro se conoce y es para la distribución a priori


Mu_n <- ((ybar*n)/s + (Mu0/T02))/((n/s)+(1/T02)) ## Hiperparametros para la posterior
Mu_n
Tn2 <- 1/((n/s)+(1/T02))                         ## Hiperparametros para la posterior

     # entonces la distribución posterior tiene como parametros Mu_n y Tn2


## Un intervalo de credibilidad para theta 

LI <- qnorm(p = 0.05,mean = Mu_n,sd = Tn2) ## intervalo de credibilidad para theta
LS <- qnorm(p= 0.975,mean = Mu_n,sd = Tn2) ## intervalo de credibilidad para theta

LI
LS



### e) ------------------------------------------------------------------------
## Graficos de densidad posterior

set.seed(123)
library(invgamma)

sigmas1 <- rinvgamma(n = 3000,shape = 25.1,rate = 1.70664)
sigmas2 <- rinvgamma(n =3000,shape=(length(cow$protein)-1)/2,rate = ((length(cow$protein)-1)*s^2)/2 )
sigmas3 <- rinvgamma(n = 3000,shape = a2,rate = b2)


#Graficos para densidades de sigma.


plot(density(sigmas1),xlab = expression(sigma^2),xlim=c(0.035,0.30),ylim=c(0,50),main = "Densidades para Sigma")
lines(density(sigmas2),col="red",lty= 2)
lines(density(sigmas3),col="green",lty= 2)
grid()
legend("topright", legend=c("sigma literal A", 
                            "Sigma literal B", 
                            "Sigma literal D"),
       col=c(1,2,3), lty=1,cex=0.8)

#Graficos para densidades de theta.


thetas1 <- rt(3000,mun)*sqrt(param2)+mun
thetas2 <- ybar+rt(3000,length(cow$protein)-1)*(s/sqrt(length(cow$protein)))

thetas3 <- rnorm(n = 3000,mean = Mu_n,sd = sqrt(Tn2))

plot(density(thetas1),ylim=c(0,15),xlim=c(3.1,3.5),xlab = expression(theta),main = "Densidades para la Media")
lines(density(thetas2),col="red")
lines(density(thetas3),col="green")
grid()

legend("topright", legend=c("Media literal A", 
                            "Media literal B", 
                            "Meida literal D"),
       col=c(1,2,3), lty=1,cex=0.8)

##-----------------------------------------------------------------------------
##////////////////////////////////////////////////////////////////////////////
##----------------------------------------------------------------------------


## Ejercicio Número 2. Intervalo de credibilidad para theta



xi<-c(24.42, 1.51, 15.38, 7.50, 9.09, 5.98,15.04,  57.99,  1.40,  10.27)
n<-10
alpha0<-11
beta0<-2.5
alpha<- n+ alpha0
beta<- sum(1/xi)+beta0

LS<-qgamma(0.975,shape= alpha,scale=1/ beta)
LI <- qgamma(0.025,shape= alpha,scale=1/ beta)



intervalo_theta2 <- c(LI,LS)
intervalo_theta2



##-----------------------------------------------------------------------------
##////////////////////////////////////////////////////////////////////////////
##----------------------------------------------------------------------------



## Ejercicio Número 3 , Multinomial


library(MCMCpack)
set.seed(123)
alpha<-c(212,321,344,1302)
thetas<-rdirichlet(3000, alpha)

dim(thetas)

mat <- rbind(c(1,1,1,2,2,2),c(3,3,3,4,4,4),c(5,5,5,6,6,6))
np <-layout(mat = mat,widths = rep.int(x = 1,times = ncol(mat)),heights = rep.int(x = 1,times = nrow(mat)),respect = FALSE)

plot(density(thetas[,1]-thetas[,2]),main="",xlab = expression(theta[1]-theta[2]),
     ylab="Densidad") ## Theta_1 - Theta_2
grid()

plot(density(thetas[,1]-thetas[,3]),main="",xlab = expression(theta[1]-theta[3]),
     ylab="Densidad") ## Theta_1 - Theta_2
grid()

plot(density(thetas[,1]-thetas[,4]),main="",xlab = expression(theta[1]-theta[4]),
     ylab="Densidad") ## Theta_1 - Theta_2

grid()
plot(density(thetas[,2]-thetas[,3]),main="",xlab = expression(theta[2]-theta[3]),
     ylab="Densidad") ## Theta_1 - Theta_2
grid()

plot(density(thetas[,2]-thetas[,4]),main="",xlab = expression(theta[2]-theta[4]),
     ylab="Densidad") ## Theta_1 - Theta_2
grid()

plot(density(thetas[,3]-thetas[,4]),main="",xlab = expression(theta[3]-theta[4]),
     ylab="Densidad") ## Theta_1 - Theta_2
grid()


## literal B 


medias <- apply(X = thetas,MARGIN = 2,FUN = "mean")
medias
Varianza <- apply(X = thetas,MARGIN = 2,FUN = "sd")
Varianza

boxplot(thetas)
lines(medias,col="red")
grid()
