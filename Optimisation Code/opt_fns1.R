# function to install and require packages
# https://gist.github.com/stevenworthington/3178163
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

##Function to generate balances given max balance at terminal age
bal_list <- function(age_term_bal, factor){
  ages <- rev(age_start:(age_retire - 1))
  start_bal_list <- matrix(nrow = length(ages), ncol = 2)
  start_bal_list [ , 1] <- ages
  colnames(start_bal_list) <- c("Age", "Start_Bal")
  for(idx in 1:40){
    start_bal_list [idx , 2] <- round_any(age_term_bal / (factor ** (idx - 1)), 10000, f = ceiling)
  }
  as.data.frame(start_bal_list)
}

##Function to generate balances given max balance
fnbalance <- function(itBalStartMax, spacing) {
  fnbalance <- rev(seq(0, itBalStartMax, by = spacing)) 
  fnbalance <- fnbalance[fnbalance >= 1]
  fnbalance
}

##Function to do schumaker interpolation 
Interpol.Fun1 <- function(Bals){   #for ordered data only
  Interpol.consum <- list(as.numeric(schu_spline$Spline(Bals)))
  Interpol.consum <- c(do.call("cbind",Interpol.consum))
  return(Interpol.consum)
}

Interpol.Fun2 <- function(Bals, LB, UB){   #for ordered data only
  Interpol.consum <- list(as.numeric(schu_spline$Spline(Bals)))
  Interpol.consum <- c(do.call("cbind",Interpol.consum))
  Interpol.consum <- ifelse(Interpol.consum < LB, LB, Interpol.consum)
  Interpol.consum <- ifelse(Interpol.consum > UB, UB, Interpol.consum)
  return(Interpol.consum)
}


#Taxation parameters
K1 = 10275; K2 = 41775; K3 = 89075; K4 = 170050; K5 = 215950; K6 = 539900 # income tax brackets
S1 = 0.1; S2 = 0.12; S3 = 0.22; S4 = 0.24; S5 = 0.32; S6 = 0.35; S7 = 0.37 # marginal income tax rates

#Starting point for each tax bracket
b1 <- 0
b2 <- S1 * K1
b3 <- S1 * K1 + S2 * (K2 - K1)
b4 <- S1 * K1 + S2 * (K2 - K1) + S3 * (K3 - K2) 
b5 <- S1 * K1 + S2 * (K2 - K1) + S3 * (K3 - K2) + S4 * (K4 - K3)
b6 <- S1 * K1 + S2 * (K2 - K1) + S3 * (K3 - K2) + S4 * (K4 - K3) + S5 * (K5 - K4)
b7 <- S1 * K1 + S2 * (K2 - K1) + S3 * (K3 - K2) + S4 * (K4 - K3) + S5 * (K5 - K4) + S6 * (K6 - K5)

#US federal income tax calculation
income_tax <- function(income){    
  Tax_br_1 <- (b1 + S1 * income) * (income >= 0 & income <= K1)
  Tax_br_2 <- (b2 + S2 * (income - K1)) * (income > K1 & income <= K2)
  Tax_br_3 <- (b3 + S3 * (income - K2)) * (income > K2 & income <= K3)
  Tax_br_4 <- (b4 + S4 * (income - K3)) * (income > K3 & income <= K4)
  Tax_br_5 <- (b5 + S5 * (income - K4)) * (income > K4 & income <= K5)
  Tax_br_6 <- (b6 + S6 * (income - K5)) * (income > K5 & income <= K6)
  Tax_br_7 <- (b7 + S7 * (income - K6)) * (income > K6)
  Tax <- (Tax_br_1 + Tax_br_2 + Tax_br_3 + Tax_br_4 + Tax_br_5 + Tax_br_6 + Tax_br_7)
  Tax
}

#function to calculate tax rate
tax_rate <- function(income){    
  Tax_rate_1 <- S1 * (income >= 0 & income <= K1)
  Tax_rate_2 <- S2 * (income > K1 & income <= K2)
  Tax_rate_3 <- S3 * (income > K2 & income <= K3)
  Tax_rate_4 <- S4 * (income > K3 & income <= K4)
  Tax_rate_5 <- S5 * (income > K4 & income <= K5)
  Tax_rate_6 <- S6 * (income > K5 & income <= K6)
  Tax_rate_7 <- S7 * (income > K6)
  rate <- (Tax_rate_1 + Tax_rate_2 + Tax_rate_3 + Tax_rate_4 + Tax_rate_5 + Tax_rate_6 + Tax_rate_7)
  rate
}

#CRRA utility calculation
utility_cal1 <- function(cons) { #w/o tax
  utility <- ((cons ^ (1 - rho)) / (1 - rho))
  utility
}

utility_cal2 <- function(cons) { #with tax
  posttax_cons <- cons - income_tax(cons)
  utility <- ((posttax_cons ^ (1 - rho)) / (1 - rho))
  utility
}


##Function to optimise pi_t & Wt at terminal age - 1
temp_fun1_B <- function(w) {
  cons <- salary * w
  M_t1 <- (M_t + salary * (1 - w)) * (1 + (1 - y) * itRetB + y * itRetA)
  - utility_cal2(cons) - beta * annuity2 * sum(itRetProb * (utility_cal2(M_t1 / annuity1)))
}

temp_fun1_S <- function(param) {
  cons <- salary * param[1]
  M_t1 <- (M_t + salary * (1 - param[1])) * (1 + (1 - param[2]) * itRetB + param[2] * itRetA)
  - utility_cal2(cons) - beta * annuity2 * sum(itRetProb * (utility_cal2(M_t1 / annuity1)))
}

##Function to optimise pi_t & Wt at terminal age - 2 onwards
temp_fun2_B <- function(w) {
  cons <- salary * w
  
  interp_endbal <- (M_t + salary - cons) * (1 + (1 - y) * itRetB + y * itRetA)  
  interpolated_consum <- Interpol.Fun1(Bals = interp_endbal[order(interp_endbal)]) #using ordered interp_endbal can reduce running time slightly
  interpolated_utility <- (interpolated_consum ^ (1 - rho)) / (1 - rho)
  
  utility1 <- utility_cal2(cons) 
  utility2 <- beta * sum(interpolated_utility * itRetProb) 
  utility <- -(utility1 + utility2)
  utility
}

temp_fun2_S <- function(param) {
  cons <- salary * param[1]
  
  interp_endbal <- (M_t + salary - cons) * (1 + (1 - param[2]) * itRetB + param[2] * itRetA)
  interpolated_consum <- Interpol.Fun1(Bals = interp_endbal[order(interp_endbal)]) #using ordered interp_endbal can reduce running time slightly 
  interpolated_utility <- (interpolated_consum ^ (1 - rho)) / (1 - rho)
  
  utility1 <- utility_cal2(cons) 
  utility2 <- beta * sum(interpolated_utility * itRetProb) 
  utility <- -(utility1 + utility2)
  utility
}


