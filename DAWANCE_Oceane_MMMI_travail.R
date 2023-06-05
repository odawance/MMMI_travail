### Package necessaire

library(deSolve) # pour la fonction ode

### Donnees

## Jeu de donnees sur les hospitalisations

setwd("D:/Téléchargements/DAWANCE_Oceane_MMMI_travail")
hospitalisation <- read.table("covid-hospitalizations.csv", header=TRUE, sep=",")
hospitalisation <- hospitalisation[hospitalisation$entity == "Italy",]
hospitalisation <- hospitalisation[substr(hospitalisation$date,1,4) == "2020",]
View(hospitalisation)

## Selection des donnees d'admissions a l'hopital et d'occupations des hopitaux

A <- hospitalisation[hospitalisation$indicator == "Weekly new hospital admissions",]
A$value <- round(A$value/7) # Approximation des admissions journalieres
# Approximation de la repartition moins de 18 ans - adultes
A = list(round(c(A[,5])*(2/100)),round(c(A[,5])*(98/100)))

O <- hospitalisation[hospitalisation$indicator == "Daily hospital occupancy",]
O = list(round(c(O[,5])*(2/100)),round(c(O[,5])*(98/100)))

## Visualisation des donnees d'occupation journaliere des hopitaux

nbre_jours = length(O[[1]])

par(mfrow = c(2,3))

plot(1:nbre_jours, O[[2]],xlab="Nombre de jours depuis le 24 fevrier 2020", ylab="Nombre de personnes hospitalisees",
     col="black",ylim=c(0, max(O[[2]])), main="Occupation journaliere des hopitaux")
lines(1:nbre_jours, O[[1]],type="p", col="blue")
legend("top", legend = c("Donnees (moins de 18 ans)", "Donnees (adultes)"), col = c("blue", "black"), pch = 1)

dev.new()
par(mfrow = c(2,3))

### Modele

## Repartition en differentes periodes

periodes <- list(
  c("26 fevrier", "2020-02-26", "2020-03-31"),
  c("01 avril", "2020-04-01", "2020-05-31"),
  c("01 juin", "2020-06-01", "2020-08-31"),
  c("01 septembre", "2020-09-01", "2020-10-31"),
  c("01 novembre", "2020-11-01", "2020-12-31")
)

## Selection des donnees d'admissions a l'hopital et d'occupations des hopitaux pour chaque periode

A_p <- lapply(periodes, function(periode) {
  A <- hospitalisation[hospitalisation$indicator == "Weekly new hospital admissions",]
  A$value <- round(A$value/7)
  A <- subset(A, date >= periode[2] & date <= periode[3])
  A = list(round(c(A[,5])*(2/100)),round(c(A[,5])*(98/100)))
  return(A)
})

O_p <- lapply(periodes, function(periode) {
  O <- hospitalisation[hospitalisation$indicator == "Daily hospital occupancy",]
  O <- subset(O, date >= periode[2] & date <= periode[3])
  O = list(round(c(O[,5])*(2/100)),round(c(O[,5])*(98/100)))
  return(O)
})

for (i in seq_along(periodes)) { # boucle sur les periodes
  
  periode <- periodes[[i]][1]
  A <- A_p[[i]]
  O <- O_p[[i]]
  
  nbre_jours = length(O[[1]])
  
## Initialisations
  
# Approximation de la taille de la population (moins de 18 ans - adultes) (https://knoema.fr/atlas/Italie/topics/Donn%C3%A9es-d%C3%A9mographiques)
  N = c(59438851-50282820, 50282820) 
  
# Nombre de reproduction de base
  R_0 = c(3,1,1,2,1)
  R_0 = R_0[i]
  
# Parametres pour chaque classe d'age
  gamma = c(1/10, 1/10) # Inverse du temps moyen avant guerison ou deces
  beta = R_0*gamma      # Taux de transmission
  p = c(0.002, 0.1)     # Probabilite d'être hospitalise suite a l'infection
  alpha = p*gamma/(1-p) # Taux d'hospitalisation
  delta = c(1/10, 1/14) # Inverse du temps moyen de l'hospitalisation
  param = list(beta, alpha, delta)
  
# Compartiments initiaux pour chaque classe d'age
  I0 = c(A[[1]][1]/alpha[[1]],A[[2]][1]/alpha[[2]])
  H0 = c(O[[1]][1],O[[2]][1]) 
  estimation = c(0, 1/100, 3/100, 5/100, 10/100)
  R0 = c((estimation[i])*N[[1]],(estimation[i])*N[[2]]) 
  S0 = c(N[[1]] - I0[[1]] - H0[[1]] - R0[[1]],N[[2]] - I0[[2]] - H0[[2]] - R0[[2]])  
  CI = list(S0, I0, H0)
  
## Modele SIHR
  
  t = 0:nbre_jours
  
  SIHR = function(t, X, P){
    
    S = X[1]
    I = X[2]
    H = X[3] 
    
    beta = P[1] 
    alpha = P[2] 
    delta = P[3]
    
    dS = -beta*S*I/(N[[classe]]) 
    dI = beta*S*I/(N[[classe]])-(alpha+gamma[[classe]])*I 
    dH = alpha*I -delta*H 
    dX=c(dS,dI,dH) 
    
    return(list(dX))
  }
  
## Fonction de log-vraisemblance
  
  log_likelihood=function(theta){
    
    CI = theta[1:3]
    param = theta[4:6]
    
    X = ode(CI,t,SIHR,param) # Resolution du systeme d'EDO
    
    o = X[,4]          # Hospitalisation theoriques : H(t)
    a = param[2]*X[,3] # Admissions theoriques : alpha*I(t)
    
    logL_O = dpois(O[[classe]],o,log=T)
    logL_A = dpois(A[[classe]],a,log=T)
    logL = sum(c(logL_O,logL_A))
    
    return(logL)
  }
  
## Optimisation de la log-vraisemblance
  
  theta0 = list(c(sapply(CI, "[[", 1),sapply(param, "[[", 1)),c(sapply(CI, "[[", 2),sapply(param, "[[", 2)))
  opt = list()
  classe=1
  opt = append(opt,list(optim(theta0[[1]],log_likelihood,control=list(fnscale=-1))))
  classe=2
  opt = append(opt,list(optim(theta0[[2]],log_likelihood,control=list(fnscale=-1))))
  
## Resultats de l'optimisation
  
  S0 = c(opt[[1]]$par[1],opt[[2]]$par[1])
  I0 = c(opt[[1]]$par[2],opt[[2]]$par[2])
  H0 = c(opt[[1]]$par[3],opt[[2]]$par[3])
  CI_opt = list(S0,I0,H0)
  
  beta = c(opt[[1]]$par[4],opt[[2]]$par[4])
  alpha = c(opt[[1]]$par[5],opt[[2]]$par[5])
  delta = c(opt[[1]]$par[6],opt[[2]]$par[6])
  param_opt = list(beta,alpha,delta) 
  
## Simulation du modele pour les conditions initiales et parametres estimes
  
  T = nbre_jours
  t = 1:T
  X_opt = lapply(seq_along(param_opt[[1]]), function(i) ode(sapply(CI_opt, "[[", i), t, SIHR, sapply(param_opt, "[[", i)))
  
## Visualisation des resultats optimaux
  
  dev.set(2)
  
  plot(t, O[[1]], type = "p", ylim = c(0, max(O[[2]])), xlab = paste("Nombre de jours depuis le", periode,"2020"), ylab = "Nombre de personnes hospitalisees", main = "Modele SIHR", col = "blue")
  lines(X_opt[[1]][,1],X_opt[[1]][,4], col = "red")
  lines(t, O[[2]], col="black", type="p")
  lines(X_opt[[2]][,1],X_opt[[2]][,4], col = "green")
  position = c("topleft","topright","topright","topleft","center")
  legend(position[i], legend = c("Donnees (moins de 18 ans)", "Modele (moins de 18 ans)", "Donnees (adultes)", "Modele (adultes)"), col = c("blue", "red", "black", "green"), lty = c(0,1,0,1), pch = c(1,NA,1,NA))
  
### Extrapolation
  
  T_extrapolation = T+31
  t_extrapolation = seq_len(T_extrapolation)
  X_extrapolation = lapply(seq_along(param_opt[[1]]), function(i) ode(sapply(CI_opt, "[[", i), t_extrapolation, SIHR, sapply(param_opt, "[[", i)))
  
## Visualisation des resultats de l'extrapolation
  
  dev.set(3)
  
  plot(t_extrapolation, X_extrapolation[[2]][, 4], type = "l", ylim = c(0, max(X_extrapolation[[2]][, 4])), xlab = paste("Nombre de jours depuis le", periode,"2020"), ylab = "Nombre de personnes hospitalisees", 
       main = "Extrapolation du modele SIHR (adultes)", col = "orange", lty=2)
  lines(t, O[[2]], type = "p", col = "black")
  lines(X_opt[[2]][,1],X_opt[[2]][,4], col = "green")
  legend(position[i], legend = c("Donnees", "Modele", "Extrapolation"), col = c("black", "green", "orange"), lty = c(0,1,2), pch = c(1,NA,NA))
}

## Visualisation des donnees d'occupation journaliere des hopitaux jusque fin mars 2021

hospitalisation2 <- read.table("covid-hospitalizations.csv", header=TRUE, sep=",")
hospitalisation2 <- hospitalisation2[hospitalisation2$entity == "Italy",]
hospitalisation2 <- subset(hospitalisation2, date >= "2020-02-24" & date <= "2021-03-31")
View(hospitalisation2)

O2 <- hospitalisation2[hospitalisation2$indicator == "Daily hospital occupancy",]
O2 = list(round(c(O2[,5])*(2/100)),round(c(O2[,5])*(98/100)))

nbre_jours2 = length(O2[[1]])

plot(1:nbre_jours2, O2[[2]],xlab="Nombre de jours depuis le 24 fevrier 2020", ylab="Nombre de personnes hospitalisees",
     col="black",ylim=c(0, max(O2[[2]])), main="Occupation journaliere des hopitaux")
