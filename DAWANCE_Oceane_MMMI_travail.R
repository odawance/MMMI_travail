### Packages necessaires

library(deSolve) # pour la fonction ode
library(mcmc)

### Donnees

## Jeu de donnees sur les hospitalisations

setwd("D:/Téléchargements/DAWANCE_Oceane_MMMI_travail_aout")
hospitalisation <- read.table("covid-hospitalizations.csv", header=TRUE, sep=",")
hospitalisation <- hospitalisation[hospitalisation$entity == "Italy",]
hospitalisation <- hospitalisation[substr(hospitalisation$date,1,4) == "2020",]
hospitalisation <- hospitalisation[-(1:8), ]
#View(hospitalisation)

## Selection des donnees d'admissions a l'hopital et d'occupations des hopitaux

A <- hospitalisation[hospitalisation$indicator == "Weekly new hospital admissions",]
A$value <- round(A$value/7) # Approximation des admissions journalieres
# Approximation de la repartition moins de 18 ans - adultes
# A = list(c(round((A[1:99,5])*(0.7/100)),round((A[100:310,5])*(1.5/100))),
#          c(round((A[1:99,5])*(99.3/100)),round((A[100:310,5])*(98.5/100))))
A = list(c(round((A[1:98,5])*(1.5/100)),round((A[99:310,5])*(1.5/100))),
         c(round((A[1:98,5])*(98.5/100)),round((A[99:310,5])*(98.5/100))))

O <- hospitalisation[hospitalisation$indicator == "Daily hospital occupancy",]
# O = list(c(round((O[1:98,5])*(0.7/100)),round((O[99:310,5])*(1.5/100))),
#          c(round((O[1:98,5])*(99.3/100)),round((O[99:310,5])*(98.5/100))))
O = list(c(round((O[1:98,5])*(1.5/100)),round((O[99:310,5])*(1.5/100))),
         c(round((O[1:98,5])*(98.5/100)),round((O[99:310,5])*(98.5/100))))

## Visualisation des donnees d'occupation journaliere des hopitaux

nbre_jours = length(O[[1]])

par(mfrow = c(1,2))

plot(1:nbre_jours, O[[1]],xlab="Nombre de jours depuis le 26 fevrier 2020", ylab="Nombre de personnes hospitalisees",
     col="blue")#, main="Occupation journaliere des hopitaux\n (moins de 18 ans)")
legend("topleft", legend = "Donnees (moins de 18 ans)", col = "blue", pch = 1)
plot(1:nbre_jours, O[[2]],xlab="Nombre de jours depuis le 26 fevrier 2020", ylab="Nombre de personnes hospitalisees",
     col="blue")#, main="Occupation journaliere des hopitaux \n (adultes)")
legend("topleft", legend = "Donnees (adultes)", col = "blue", pch = 1)

### Modele

## Initialisations

# Approximation de la taille de la population (moins de 18 ans - adultes)
N = c(59438851-50282820, 50282820) 

# Nombre de reproduction de base
R_0 = c(3,1,1,2,1)

# Parametres pour chaque classe d'age
gamma = c(1/10, 1/15) # Inverse du temps moyen avant guerison ou deces

contact_matrix = matrix(c(12.75, 10.85, 2.24, 14.54), nrow = 2, ncol = 2)
vp = eigen(contact_matrix/gamma)
vp_principale = max(vp$values)

beta1 = R_0[1]/vp_principale*c(1,1) # Taux de transmission
beta2 = R_0[2]/vp_principale*c(1,1)
beta3 = R_0[3]/vp_principale*c(1,1)
beta4 = R_0[4]/vp_principale*c(1,1)
beta5 = R_0[5]/vp_principale*c(1,1)
p = c(0.000333, 0.0048) # Taux d'hospitalisation
delta = c(1/10, 1/15) # Inverse du temps moyen de l'hospitalisation
param = list(beta1, beta2, beta3, beta4, beta5, p, delta,gamma)

# Compartiments initiaux pour chaque classe d'age
I0 = c(A[[1]][1]/(p[1]*gamma[1]),A[[2]][1]/(p[2]*gamma[2]))
H0 = c(O[[1]][1],O[[2]][1]) 
R0 = c(0,0)
S0 = c(N[1] - I0[1] - H0[1] - R0[1],N[2] - I0[2] - H0[2] - R0[2])  
CI = list(S0, I0, H0)

## Modele SIHR

t = 0:nbre_jours

SIHR = function(t, X, P){
  
  S1 = X[1]
  S2 = X[2]
  I1 = X[3]
  I2 = X[4]
  H1 = X[5]
  H2 = X[6]
  
  beta <- if (t >= 0 && t <= 12) {P[[1]]}
  else if (t > 12 && t <= 98) {P[[2]]}
  else if (t > 98 && t <= 201) {P[[3]]}
  else if (t > 201 && t <= 252) {P[[4]]}
  else {P[[5]]}
  
  p = P[[6]] 
  delta = P[[7]]
  gamma = P[[8]]

  dS1 = -(beta[1])*S1*(12.75*I1/(N[1])+10.85*I2/(N[2]))
  dS2 = -(beta[2])*S2*(2.24*I1/(N[1])+14.54*I2/(N[2]))
  dI1 = (beta[1])*S1*(12.75*I1/(N[1])+10.85*I2/(N[2]))-(gamma[1])*I1
  dI2 = (beta[2])*S2*(2.24*I1/(N[1])+14.54*I2/(N[2]))-(gamma[2])*I2
  dH1 = p[1]*(gamma[1])*I1-(delta[1])*H1
  dH2 = p[2]*(gamma[2])*I2-(delta[2])*H2 
  dX=c(dS1,dS2,dI1,dI2,dH1,dH2) 
  
  return(list(dX))
}

theta0 = c(CI[[2]][1],CI[[2]][2],CI[[3]][1],CI[[3]][2],
           param[[1]][1], param[[1]][2], param[[2]][1], param[[2]][2], param[[3]][1], param[[3]][2], param[[4]][1],
           param[[4]][2], param[[5]][1], param[[5]][2], param[[6]][1], param[[6]][2], param[[7]][1], param[[7]][2], param[[8]][1], param[[8]][2])
CI_X = c(N[1]-theta0[1]-theta0[3],N[2]-theta0[2]-theta0[4],theta0[1],theta0[2],theta0[3],theta0[4])
param_X = list(c(theta0[5],theta0[6]),c(theta0[7],theta0[8]),c(theta0[9],theta0[10]),
               c(theta0[11],theta0[12]),c(theta0[13],theta0[14]),c(theta0[15],theta0[16]),c(theta0[17],theta0[18]),c(theta0[19],theta0[20]))
X = ode(CI_X,t,SIHR,param_X, method="bdf")
plot(X)

## Fonction de log-vraisemblance

log_likelihood=function(log_theta){
  
  theta = exp(log_theta)
  
  if(theta[1] > N[1]/100 || theta[2] > N[2]/90 || theta[3] > 100 || theta[4] > 1000) {return(-Inf)}
  
  if (theta[15] > 1 || theta[16] > 1 || theta[17] > 1 || theta[18] > 1 || theta[19] > 1 || theta[20] > 1) {return(-Inf)}
  
  CI = c(N[1]-theta[1]-theta[3],N[2]-theta[2]-theta[4],theta[1],theta[2],theta[3],theta[4])
  param = list(c(theta[5],theta[6]),c(theta[7],theta[8]),c(theta[9],theta[10]),
               c(theta[11],theta[12]),c(theta[13],theta[14]),c(theta[15],theta[16]),c(theta[17],theta[18]),c(theta[19],theta[20]))
  
  X = ode(CI,t,SIHR,param, method="bdf") # Resolution du systeme d'EDO
  
  o1 = X[,6] # Charge hospitaliere theorique : H(t)
  o2 = X[,7]
  a1 = (param[[6]][1])*(param[[8]][1])*X[,4] # Admissions theoriques : p*gamma*I(t)
  a2 = (param[[6]][2])*(param[[8]][2])*X[,5]
  
  logL_O1 = dpois(O[[1]],o1,log=T)
  logL_O2 = dpois(O[[2]],o2,log=T)
  logL_A1 = dpois(A[[1]],a1,log=T)
  logL_A2 = dpois(A[[2]],a2,log=T)
  logL = sum(c(logL_O1,logL_O2,logL_A1,logL_A2))
  
  return(logL)
}

theta0 = log(c(CI[[2]][1],CI[[2]][2],CI[[3]][1],CI[[3]][2],
           param[[1]][1], param[[1]][2], param[[2]][1], param[[2]][2], param[[3]][1], param[[3]][2], param[[4]][1],
           param[[4]][2], param[[5]][1], param[[5]][2], param[[6]][1], param[[6]][2], param[[7]][1], param[[7]][2], param[[8]][1], param[[8]][2]))
#log_likelihood(theta0)

## Optimisation de la log-vraisemblance

time_opt <- system.time(opt <- list(optim(theta0,log_likelihood,control=list(fnscale=-1))))
opt[[1]]$par = exp(opt[[1]]$par)

## Resultats de l'optimisation

S01 = N[1]-opt[[1]]$par[1]-opt[[1]]$par[3];
S02 = N[2]-opt[[1]]$par[2]-opt[[1]]$par[3];
I01 = opt[[1]]$par[1];
I02 = opt[[1]]$par[2];
H01 = opt[[1]]$par[3];
H02 = opt[[1]]$par[4];
CI_opt = c(S01,S02,I01,I02,H01,H02);

beta1 = c(opt[[1]]$par[5],opt[[1]]$par[6]);
beta2 = c(opt[[1]]$par[7],opt[[1]]$par[8]);
beta3 = c(opt[[1]]$par[9],opt[[1]]$par[10]);
beta4 = c(opt[[1]]$par[11],opt[[1]]$par[12]);
beta5 = c(opt[[1]]$par[13],opt[[1]]$par[14]);
p = c(opt[[1]]$par[15],opt[[1]]$par[16]);
delta = c(opt[[1]]$par[17],opt[[1]]$par[18]);
gamma = c(opt[[1]]$par[19],opt[[1]]$par[20]);
param_opt = list(beta1,beta2,beta3,beta4,beta5,p,delta,gamma) 

## Simulation du modele pour les conditions initiales et parametres estimes

T = nbre_jours
t = 1:T
X_opt = ode(CI_opt,t,SIHR,param_opt, method="bdf")

## Visualisation des resultats optimaux

par(mfrow = c(1,2))

plot(t, O[[1]], type = "p", ylim = c(0, 800), xlab = paste("Nombre de jours depuis le 26 fevrier 2020"), ylab = "Nombre de personnes hospitalisees", col = "blue")#, main = "Modele SIHR (moins de 18 ans)"
lines(X_opt[,1],X_opt[,6], col = "red", lwd = 2)
legend("topleft", legend = c("Donnees (moins de 18 ans)", "Modele (moins de 18 ans)"), col = c("blue", "red"), lty = c(0,1), lwd = c(1,2), pch = c(1,NA))

plot(t, O[[2]], type = "p", ylim = c(0, 50000), xlab = paste("Nombre de jours depuis le 26 fevrier 2020"), ylab = "Nombre de personnes hospitalisees", col = "blue")#, main = "Modele SIHR (adultes)"
lines(X_opt[,1],X_opt[,7], col = "red", lwd = 2)
legend("topleft", legend = c("Donnees (adultes)", "Modele (adultes)"), col = c("blue", "red"), lty = c(0,1), lwd = c(1,2), pch = c(1,NA))

## MCMC

theta_opt = log(c(CI_opt[3],CI_opt[4],CI_opt[5],CI_opt[6],
             param_opt[[1]][1],param_opt[[1]][2],param_opt[[2]][1],param_opt[[2]][2],param_opt[[3]][1],param_opt[[3]][2],param_opt[[4]][1],param_opt[[4]][2],
             param_opt[[5]][1],param_opt[[5]][2],param_opt[[6]][1],param_opt[[6]][2],param_opt[[7]][1],param_opt[[7]][2],param_opt[[8]][1],param_opt[[8]][2]))
#time_out <- system.time(out <- metrop(log_likelihood, theta_opt, 100000, scale=0.0005))
time_out <- system.time(out <- metrop(log_likelihood, theta_opt, nbatch=100000, blen=10, scale=0.0005))
out$accept
out$batch = exp(out$batch)

#theta_out = log(out$batch[nrow(out$batch),])
#time_out <- system.time(out <- metrop(log_likelihood, theta_out, nbatch=100000, blen=10, scale=0.0005))
#out$batch = exp(out$batch)

#indices=(nrow(out$batch) - 4999):nrow(out$batch)
#indices=indices[seq(1, length(indices), by = 50)]
indices=(nrow(out$batch) - 99999):nrow(out$batch)
indices=indices[seq(1, length(indices), by = 100)]
out_select <- out$batch[indices, ]

param_finaux <- apply(out_select[, 5:20], 2, function(col) {
  c(moyenne = mean(col), quantile_025 = quantile(col, probs = 0.025), quantile_975 = quantile(col, probs = 0.975))
})

results1 = function(out){
  
  CI_out = c(N[1]-out[1]-out[3],N[2]-out[2]-out[1],out[1],out[2],out[3],out[4]);
  param_out = list(c(out[5],out[6]),c(out[7],out[8]),c(out[9],out[10]),c(out[11],out[12]),
                   c(out[13],out[14]),c(out[15],out[16]),c(out[17],out[18]),c(out[19],out[20]));
  
  # param_out = list(c(out[5],out[6]),c(out[5],out[6]),c(out[5],out[6]),c(out[5],out[6]),
  #                  c(out[5],out[6]),c(out[15],out[16]),c(out[17],out[18]),c(out[19],out[20]));
  
  X_out = ode(CI_out,t,SIHR,param_out, method="bdf");
  
  lines(X_out[,1],X_out[,6], col = "red")

  return(X_out[,6])
}

results2 = function(out){
  
  t = 1:nbre_jours
  
  CI_out = c(N[1]-out[1]-out[3],N[2]-out[2]-out[1],out[1],out[2],out[3],out[4]);
  param_out = list(c(out[5],out[6]),c(out[7],out[8]),c(out[9],out[10]),c(out[11],out[12]),
                   c(out[13],out[14]),c(out[15],out[16]),c(out[17],out[18]),c(out[19],out[20]));
  
  # param_out = list(c(out[5],out[6]),c(out[5],out[6]),c(out[5],out[6]),c(out[5],out[6]),
  #                  c(out[5],out[6]),c(out[15],out[16]),c(out[17],out[18]),c(out[19],out[20]));
  
  X_out = ode(CI_out,t,SIHR,param_out, method="bdf");

  lines(X_out[,1],X_out[,7], col = "red")
  
  return(X_out[,7])
}

resume <- function(matrice) {
  resultats <- apply(matrice, 1, function(ligne) {
    moyenne <- mean(ligne)
    quantile_025 <- quantile(ligne, probs = 0.025)
    quantile_975 <- quantile(ligne, probs = 0.975)
    
    c(moyenne = moyenne, quantile_025 = quantile_025, quantile_975 = quantile_975)
  })
  return(resultats)
}

par(mfrow = c(1,2))

plot(t, O[[1]], type = "p", ylim = c(0, 800), xlab = paste("Nombre de jours depuis le 26 fevrier 2020"), ylab = "Nombre de personnes hospitalisees", col = "blue")#, main = "Modele SIHR (moins de 18 ans)"
#plot(t, O[[1]], type = "p", ylim = c(0, 600), xlab = paste("Nombre de jours depuis le 26 fevrier 2020"), ylab = "Nombre de personnes hospitalisees", col = "blue")#, main = "Modele SIHR (moins de 18 ans)"
#plot(t, O[[1]], type = "p", ylim = c(0, 5000), xlab = paste("Nombre de jours depuis le 26 fevrier 2020"), ylab = "Nombre de personnes hospitalisees", col = "blue")#, main = "Modele SIHR (moins de 18 ans)"
#plot(t, O[[1]], type = "p", ylim = c(0, 6000), xlab = paste("Nombre de jours depuis le 26 fevrier 2020"), ylab = "Nombre de personnes hospitalisees", col = "blue")#, main = "Modele SIHR (moins de 18 ans)"

H1_out <- c()
for (i in 1:nrow(out_select)) {
  H1_out <- cbind(H1_out, results1(out_select[i, ])) 
}

legend("topleft", legend = c("Donnees (moins de 18 ans)", "Modele (moins de 18 ans)"), col = c("blue", "red"), lty = c(0,1), lwd = c(1,2), pch = c(1,NA))
#legend("topleft", legend = c("Donnees (moins de 18 ans)", "Modele : moyenne (moins de 18 ans)", "Modele : quantiles (0.025 ; 0.975) (moins de 18 ans)"), col = c("blue", "red", "black"), lty = c(0,1,2), lwd = c(1,2,2), pch = c(1,NA,NA))
#legend("topright", legend = c("Donnees (moins de 18 ans)", "Modele (moins de 18 ans)"), col = c("blue", "red"), lty = c(0,1), lwd = c(1,2), pch = c(1,NA))
#legend("topright", legend = c("Donnees (moins de 18 ans)", "Modele : moyenne (moins de 18 ans)", "Modele : quantiles (0.025 ; 0.975) (moins de 18 ans)"), col = c("blue", "red", "black"), lty = c(0,1,2), lwd = c(1,2,2), pch = c(1,NA,NA))

#lines(resume(H1_out)[1,], col = "red", lwd = 2)
#lines(resume(H1_out)[2,], col = "black", lty = 2, lwd = 2)
#lines(resume(H1_out)[3,], col = "black", lty = 2, lwd = 2)

plot(t, O[[2]], type = "p", ylim = c(0, 50000), xlab = paste("Nombre de jours depuis le 26 fevrier 2020"), ylab = "Nombre de personnes hospitalisees", col = "blue")#, main = "Modele SIHR (adultes)", col = "blue")
#plot(t, O[[2]], type = "p", ylim = c(0, 40000), xlab = paste("Nombre de jours depuis le 26 fevrier 2020"), ylab = "Nombre de personnes hospitalisees", col = "blue")#, main = "Modele SIHR (adultes)", col = "blue")
#plot(t, O[[2]], type = "p", ylim = c(0, 200000), xlab = paste("Nombre de jours depuis le 26 fevrier 2020"), ylab = "Nombre de personnes hospitalisees", col = "blue")#, main = "Modele SIHR (adultes)", col = "blue")
#plot(t, O[[2]], type = "p", ylim = c(0, 150000), xlab = paste("Nombre de jours depuis le 26 fevrier 2020"), ylab = "Nombre de personnes hospitalisees", col = "blue")#, main = "Modele SIHR (adultes)", col = "blue")

H2_out <- c()
for (i in 1:nrow(out_select)) {
  H2_out <- cbind(H2_out, results2(out_select[i, ])) 
}

legend("topleft", legend = c("Donnees (adultes)", "Modele (adultes)"), col = c("blue", "red"), lty = c(0,1), lwd = c(1,2), pch = c(1,NA))
#legend("topleft", legend = c("Donnees (adultes)", "Modele : moyenne (adultes)", "Modele : quantiles (0.025 ; 0.975) (adultes)"), col = c("blue", "red", "black"), lty = c(0,1,2), lwd = c(1,2,2), pch = c(1,NA,NA))
#legend("topright", legend = c("Donnees (adultes)", "Modele (adultes)"), col = c("blue", "red"), lty = c(0,1), lwd = c(1,2), pch = c(1,NA))
#legend("topright", legend = c("Donnees (adultes)", "Modele : moyenne (adultes)", "Modele : quantiles (0.025 ; 0.975) (adultes)"), col = c("blue", "red", "black"), lty = c(0,1,2), lwd = c(1,2,2), pch = c(1,NA,NA))

#lines(resume(H2_out)[1,], col = "red", lwd = 2)
#lines(resume(H2_out)[2,], col = "black", lty = 2, lwd = 2)
#lines(resume(H2_out)[3,], col = "black", lty = 2, lwd = 2)
