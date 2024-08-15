#CRRA Projections (updated with two risk measures included)

## Initialisation  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

source("Initial_proj_functions_CRRA_updated.r")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr","matrixStats")

ipak(libs)

# Loading returns and optimisation results
## optimisation results
aa_decision_CRRA_s = read.csv("worker_S_case_rho4_salary70000_base_v5.csv")
aa_decision_CRRA_b = read.csv("worker_B_case_rho4_salary70000_base_v5.csv")
aa_decision_CRRA_h1 = read.csv("mgmt_H1_rho_4salary70000_base1_v5.1.csv")
aa_decision_CRRA_h2 = read.csv("mgmt_H2_rho_4salary70000_base1_v5.1.csv")
colnames(aa_decision_CRRA_s)[7] <- "eqprop"

## returns
return_2assets <- read.csv("Ret_2_asset.csv", header = T)
itRetA <- as.numeric(return_2assets$DE)
itRetB <- as.numeric(return_2assets$Cash)
len <- nrow(return_2assets) 

return_2assets_new <- read.csv("Ret_2_asset_alpha1.csv", header = T)
itRetA_new <- as.numeric(return_2assets_new$DE)
itRetB_new <- as.numeric(return_2assets_new$Cash)

## parameters
rf <- mean(itRetB)               
real_discount <- 0.96        # subjective discount rate
age_retire <- 65
age_start <- 25
rho <- 4                     # utility parameter
size_proj <- 10000           # approximate size of the return vector for each year
salary <- 70000
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
terminal_target <- 857868.51 # corresponding targeted terminal wealth for investment target B (in H1&H2 cases)

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

##setting upper bound & lower bound for risky asset allocation proportion (worker/manager)
U_bound1_1 <- 1
U_bound1_2 <- 1.16
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
UB_pi_s <- min(U_bound1_1, U_bound2) 
UB_pi_m<- min(U_bound1_2, U_bound2)

L_bound1_1 <- 0 
L_bound1_2 <- -0.16
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
LB_pi_s <- max(L_bound1_1, L_bound2) 
LB_pi_m<- max(L_bound1_2, L_bound2) 

# Simulation scenario variables 
start_bal <- c(0, 9500, 19000)

# Set the size of projection matrix
sim_year_total <- age_retire - age_start 
sampling_set <- ceiling(size_proj / len) 
vertical_size <- len * sampling_set


# Settings for Util_mat   
Util_mat <- matrix(nrow = 1, ncol = 35) ##
colnames(Util_mat) <- c("Start.Bal", "Average.Final.Bal.s", "Average.Pretax.CEC.s", "Average.Posttax.CEC.s","Average.Final.Bal.h1", "Average.Pretax.CEC.h1", "Average.Posttax.CEC.h1", "Average.Final.Bal.b", "Average.Pretax.CEC.b", "Average.Posttax.CEC.b",
                        "Average.Final.Bal.h2", "Average.Pretax.CEC.h2", "Average.Posttax.CEC.h2", "Extra.Fee.Rate.h1.vs.s_cec", "Extra.Fee.Rate.h2.vs.s_cec", "Extra.Fee.Rate.s.vs.b_cec", "Extra.Fee.Rate.s.vs.b_cec", "Extra.Fee.Rate.h1.vs.b_cec",
                        "Extra.Fee.Rate.h2.vs.b_cec", "VaR.Final.Bal.s.95percent", "VaR.Final.Bal.b.95percent", "VaR.Final.Bal.h1.95percent", "VaR.Final.Bal.h2.95percent", "CVaR.Final.Bal.s.95percent", "CVaR.Final.Bal.b.95percent", "CVaR.Final.Bal.h1.95percent", "CVaR.Final.Bal.h2.95percent",
                        "VaR.Final.Bal.s.99percent", "VaR.Final.Bal.b.99percent", "VaR.Final.Bal.h1.99percent", "VaR.Final.Bal.h2.99percent", "CVaR.Final.Bal.s.99percent", "CVaR.Final.Bal.b.99percent", "CVaR.Final.Bal.h1.99percent", "CVaR.Final.Bal.h2.99percent")
Util_mat <- as.data.frame(Util_mat)

Outlier <- 2
###########################################################################################################################
# Generate Util_mat (compare among all four cases)
Util <- foreach(indx = 1:length(start_bal), .combine = rbind) %dopar% {
  
  bal_indx <- start_bal[indx]
  cat(paste0("bal", bal_indx, "in progress\n"))
  
  proj_CRRA_s(bal_start = bal_indx, aa_decision = aa_decision_CRRA_s, 
              returnData = return_2assets, sim_year_total = sim_year_total, vertical_size = vertical_size, alloboundL = LB_pi_s, alloboundU = UB_pi_s)
  proj_CRRA_h1(bal_start = bal_indx, aa_decision = aa_decision_CRRA_h1, 
              returnData = return_2assets, sim_year_total = sim_year_total, vertical_size = vertical_size, alloboundL = LB_pi_m, alloboundU = UB_pi_m)
  proj_CRRA_h2(bal_start = bal_indx, aa_decision = aa_decision_CRRA_h2,  #diff. risky return due to alpha
               returnData = return_2assets_new, sim_year_total = sim_year_total, vertical_size = vertical_size, alloboundL = LB_pi_s, alloboundU = UB_pi_s)
  proj_CRRA_b(bal_start = bal_indx, aa_decision = aa_decision_CRRA_b, 
              returnData = return_2assets, sim_year_total = sim_year_total, vertical_size = vertical_size, alloboundL = LB_pi_s, alloboundU = UB_pi_s)
  
  Util_mat[1, 1] <- bal_indx 
  
  ave_bal1 <- Util_DF_CRRA_s[ ,3]
  ave_bal2 <- Util_DF_CRRA_h1[ ,3]
  ave_bal3 <- Util_DF_CRRA_b[ ,3]
  ave_bal4 <- Util_DF_CRRA_h2[ ,3]
  
  Util_mat[1, 2] <- ave_bal1
  Util_mat[1, 5] <- ave_bal2
  Util_mat[1, 8] <- ave_bal3
  Util_mat[1, 11] <- ave_bal4
  
  cec1 <- Util_DF_CRRA_s[ ,4]
  cec2 <- Util_DF_CRRA_h1[ ,4]
  cec3 <- Util_DF_CRRA_b[ ,4]
  cec4 <- Util_DF_CRRA_h2[ ,4]
  
  Util_mat[1,3] <- cec_fun(cec1)
  Util_mat[1,6] <- cec_fun(cec2)
  Util_mat[1,9] <- cec_fun(cec3)
  Util_mat[1,12] <- cec_fun(cec4)
  
  Util_mat[1,4] <- Util_mat[1,3] - income_tax(Util_mat[1,3])
  Util_mat[1,7] <- Util_mat[1,6] - income_tax(Util_mat[1,6])
  Util_mat[1,10] <- Util_mat[1,9] - income_tax(Util_mat[1,9])
  Util_mat[1,13] <- Util_mat[1,12] - income_tax(Util_mat[1,12])
  
  load("ending_bal_s.RData"); load("ending_bal_h1.RData"); load("ending_bal_b.RData"); load("ending_bal_h2.RData")
  
  #compute extra management fee rate
  Util_mat[1, 14] <- x_mgt_fee1_1(base_cec = cec1, f_cec = cec2, bal_start = bal_indx, aa_decision = aa_decision_CRRA_h1,
                                returnData = return_2assets, sim_year_total = sim_year_total, vertical_size = vertical_size, alloboundL = LB_pi_m, alloboundU = UB_pi_m)
  Util_mat[1, 15] <- x_mgt_fee1_1(base_cec = cec1, f_cec = cec4, bal_start = bal_indx, aa_decision = aa_decision_CRRA_h2,
                                returnData = return_2assets_new, sim_year_total = sim_year_total, vertical_size = vertical_size, alloboundL = LB_pi_s, alloboundU = UB_pi_s)
  Util_mat[1, 16] <- x_mgt_fee1_2(base_cec = cec3, f_cec = cec1, bal_start = bal_indx, aa_decision = aa_decision_CRRA_s,
                                returnData = return_2assets, sim_year_total = sim_year_total, vertical_size = vertical_size, alloboundL = LB_pi_s, alloboundU = UB_pi_s)
  
  Util_mat[1, 17] <- x_mgt_fee1_2(base_cec = cec3, f_cec = cec1, bal_start = bal_indx, aa_decision = aa_decision_CRRA_s,
                                  returnData = return_2assets, sim_year_total = sim_year_total, vertical_size = vertical_size, alloboundL = LB_pi_s, alloboundU = UB_pi_s)
  Util_mat[1, 18] <- x_mgt_fee1_1(base_cec = cec3, f_cec = cec2, bal_start = bal_indx, aa_decision = aa_decision_CRRA_h1,
                                  returnData = return_2assets, sim_year_total = sim_year_total, vertical_size = vertical_size, alloboundL = LB_pi_m, alloboundU = UB_pi_m)
  Util_mat[1, 19] <- x_mgt_fee1_1(base_cec = cec3, f_cec = cec4, bal_start = bal_indx, aa_decision = aa_decision_CRRA_h2,
                                  returnData = return_2assets_new, sim_year_total = sim_year_total, vertical_size = vertical_size, alloboundL = LB_pi_s, alloboundU = UB_pi_s)
  
  #recall 95% VaR results for terminal balances
  VaR_bal1 <- Util_DF_CRRA_s[ ,5]
  VaR_bal2 <- Util_DF_CRRA_b[ ,5]
  VaR_bal3 <- Util_DF_CRRA_h1[ ,5]
  VaR_bal4 <- Util_DF_CRRA_h2[ ,5]
  
  Util_mat[1, 20] <- VaR_bal1
  Util_mat[1, 21] <- VaR_bal2
  Util_mat[1, 22] <- VaR_bal3
  Util_mat[1, 23] <- VaR_bal4
  
  #recall 95% CVaR results for terminal balances
  CVaR_bal1 <- Util_DF_CRRA_s[ ,6]
  CVaR_bal2 <- Util_DF_CRRA_b[ ,6]
  CVaR_bal3 <- Util_DF_CRRA_h1[ ,6]
  CVaR_bal4 <- Util_DF_CRRA_h2[ ,6]
  
  Util_mat[1, 24] <- CVaR_bal1
  Util_mat[1, 25] <- CVaR_bal2
  Util_mat[1, 26] <- CVaR_bal3
  Util_mat[1, 27] <- CVaR_bal4
  
  #recall 99% VaR results for terminal balances
  VaR_bal5 <- Util_DF_CRRA_s[ ,7]
  VaR_bal6 <- Util_DF_CRRA_b[ ,7]
  VaR_bal7 <- Util_DF_CRRA_h1[ ,7]
  VaR_bal8 <- Util_DF_CRRA_h2[ ,7]
  
  Util_mat[1, 28] <- VaR_bal5
  Util_mat[1, 29] <- VaR_bal6
  Util_mat[1, 30] <- VaR_bal7
  Util_mat[1, 31] <- VaR_bal8
  
  #recall 99% CVaR results for terminal balances
  CVaR_bal5 <- Util_DF_CRRA_s[ ,8]
  CVaR_bal6 <- Util_DF_CRRA_b[ ,8]
  CVaR_bal7 <- Util_DF_CRRA_h1[ ,8]
  CVaR_bal8 <- Util_DF_CRRA_h2[ ,8]
  
  Util_mat[1, 32] <- CVaR_bal5
  Util_mat[1, 33] <- CVaR_bal6
  Util_mat[1, 34] <- CVaR_bal7
  Util_mat[1, 35] <- CVaR_bal8
  
  Util_mat
}

write.csv(Util, "Projection_summary_crra1.csv")
