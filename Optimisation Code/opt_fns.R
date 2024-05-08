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
  fnbalance
}


##Functions to do interpolation 
Interpol.Fun1 <- function(X.Data, Y.Data, X, LB, UB){   #for value_m (linear interpolation)
  Interpol.Obj <- approxExtrap(x = X.Data, y = Y.Data, xout = X, rule = 2)
  Interpol.Obj <- c(do.call("cbind", Interpol.Obj[2]))
  Interpol.Obj <- ifelse(Interpol.Obj < LB, LB, Interpol.Obj)
  Interpol.Obj <- ifelse(Interpol.Obj > UB, UB, Interpol.Obj)
  return(Interpol.Obj)
}

Interpol.Fun2 <- function(schu_spline, Bals){   #for utility  (schumaker interpolation; relatively more monotonic and smooth)
  Interpol.obj <- list(as.numeric(schu_spline$Spline(Bals)))
  Interpol.obj <- c(do.call("cbind",Interpol.obj))
  Interpol.obj <- ifelse(Interpol.obj < 0, 0, Interpol.obj)
  return(Interpol.obj)
}

Interpol.Fun3 <- function(Consum.Data, Allocation.Data, Cons){  # for allocation (linear interpolation)
  Interpol.Obj <- approxExtrap(x = Consum.Data, y = Allocation.Data, xout = Cons, rule = 2)
  Interpol.Obj <- c(do.call("cbind", Interpol.Obj[2]))
  Interpol.Obj <- ifelse(Interpol.Obj < L_bound, L_bound, Interpol.Obj)
  Interpol.Obj <- ifelse(Interpol.Obj > U_bound, U_bound, Interpol.Obj)
  return(Interpol.Obj)
}


#Function to find out the highest eqprop which leads to the highest prob
u_max <- function (x) {  
  max(which(x == max(x)))
}

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

#function to calculate utility at terminal age - 1
utility_pre_terminal <- function(param) {
  consumption <- salary * param[1] 
  prop1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = consumption)
  itBalEnd <- (itBalStart - consumption + salary) * (1 + itRetA * prop1 + itRetB * (1 - prop1))
  utility1 <- utility_cal2(consumption) 
  utility2 <- sum(utility_cal2(itBalEnd / annuity1) * itRetProb) * beta * annuity2
  utility <- - (utility1 + utility2)
  utility
}

#function to calculate utility at terminal age - 2 onwards
utility_xx <- function(param) {
  consumption <- salary * param[1] 
  prop1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = consumption)
  itBalEnd <- (itBalStart - consumption + salary) * (1 + itRetA * prop1 + itRetB * (1 - prop1))
  utility1 <- utility_cal2(consumption) 
  interpolated_consum <-  Interpol.Fun2(schu_spline = schu_spline1, Bals = itBalEnd)
  interpolated_utility <- utility_cal1(interpolated_consum)
  utility2 <- sum(interpolated_utility * itRetProb) * beta
  utility <- - (utility1 + utility2)
  utility
}
