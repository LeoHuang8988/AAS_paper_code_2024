## Projection functions

## function to install and require packages
# https://gist.github.com/stevenworthington/3178163
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

## function to generate matrices of risky asset returns - for projection
gen_mats <- function(data, return_size = len, n = sampling_set, years = sim_year_total, ret_seed = 54){
  sampling_index <- 1:return_size
  vertical_size <- return_size * n
  dummy_index <- 1:(n * years)
  set.seed(ret_seed)
  sample_list <- lapply(dummy_index, function(x) base::sample(sampling_index, size = length(sampling_index), replace = F))
  sample_vector <- do.call(c, sample_list)
  sample_matrix <- matrix(sample_vector, ncol = sim_year_total, byrow = FALSE)
  
  asset1_vector <- data[sample_vector, 2]
  asset1_matrix <- matrix(asset1_vector, ncol = sim_year_total, byrow = FALSE)
  
  return(list(asset1_matrix = asset1_matrix))
}

## function to do linear interpolation - for projection (applies to eqprop(worker/manager) & consprop(worker/manager))
Interpol.Fun <- function(Balance.Data, Allocation.Data, Bals, Rule, LB, UB){   
  Interpol.Obj <- approxExtrap(x = Balance.Data, y = Allocation.Data, xout = Bals, rule = Rule)
  Interpol.Obj <- c(do.call("cbind", Interpol.Obj[2]))
  Interpol.Obj <- ifelse(Interpol.Obj < LB, LB, Interpol.Obj)
  Interpol.Obj <- ifelse(Interpol.Obj > UB, UB, Interpol.Obj)
  return(Interpol.Obj)
}

## Settings for income tax calculations
K1 = 10275; K2 = 41775; K3 = 89075; K4 = 170050; K5 = 215950; K6 = 539900 # income tax brackets
S1 = 0.1; S2 = 0.12; S3 = 0.22; S4 = 0.24; S5 = 0.32; S6 = 0.35; S7 = 0.37 # marginal income tax rates

## Starting point for each tax bracket
b1 <- 0
b2 <- S1 * K1
b3 <- S1 * K1 + S2 * (K2 - K1)
b4 <- S1 * K1 + S2 * (K2 - K1) + S3 * (K3 - K2) 
b5 <- S1 * K1 + S2 * (K2 - K1) + S3 * (K3 - K2) + S4 * (K4 - K3)
b6 <- S1 * K1 + S2 * (K2 - K1) + S3 * (K3 - K2) + S4 * (K4 - K3) + S5 * (K5 - K4)
b7 <- S1 * K1 + S2 * (K2 - K1) + S3 * (K3 - K2) + S4 * (K4 - K3) + S5 * (K5 - K4) + S6 * (K6 - K5)

## US federal income tax calculation
income_tax <- function(income){    #true tax function 
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

## function to calculate certainty equivalent consumption (CEC)
cec_fun <- function(cec){
  dis_fac <- seq(1:sim_year_total)
  dis_fac1 <- c(1, real_discount ^ dis_fac)
  
  temp_fun1 <- function(x1){ # adjust for buying annuity at age 65
    util <- ((x1 - income_tax(x1)) ^ (1 - rho)) / (1 - rho)
    util_vec <- rep(util, sim_year_total)
    sim_util1 <- sum(util_vec * dis_fac1[-(sim_year_total + 1)]) + util * dis_fac1[(sim_year_total + 1)] * annuity2
    sim_util <- (sim_util1 * (1 - rho)) ^ (1 / (1 - rho))
    sim_util - cec
  }
  
  y <- uniroot(f = temp_fun1, interval = c(0, 100000), extendInt = "yes")$root
  y
}

# Projection function for worker
proj_CRRA_s <- function(bal_start, aa_decision, returnData, sim_year_total, vertical_size, alloboundL, alloboundU){
  
  rets <- gen_mats(returnData)
  asset1_matrix <- rets$asset1_matrix
  
  aa1_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year
  aa2_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year
  con_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year
  ret_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # end of year
  bal_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # end of year
  consum_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
  posttax_consum_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
  utility_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
  
  for (sim_year_idx in 1:sim_year_total) {
    
    age_idx <- age_start + sim_year_idx - 1
    
    cat(paste0("CRRA projection(S) - age ", age_idx, " in progress\n"))
    
    LB_C <- ifelse(age_idx < age_tr, LB_C1_rate, LB_C2_rate)
    
    aa_temp <- aa_decision[rev(which(aa_decision$age == age_idx)), c(2,6,7)]
    
    ## assign asset allocations/consprops based on interpolated value
    if (sim_year_idx > 1) {
      aa1_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$eqprop, Bals = bal_matrix[, sim_year_idx - 1], Outlier, LB = alloboundL, UB = alloboundU)
      aa2_matrix[, sim_year_idx] <- 1 - aa1_matrix[, sim_year_idx]
      con_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_matrix[, sim_year_idx - 1], Outlier, LB = LB_C, UB = 1)
    } else {
      aa1_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$eqprop, Bals = bal_start, Outlier, LB = alloboundL, UB = alloboundU)
      aa2_matrix[, sim_year_idx] <- 1 - aa1_matrix[, sim_year_idx]
      con_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_start, Outlier , LB = LB_C, UB = 1)
    }
    
    ## calculating investment return
    ret_matrix[, sim_year_idx] <- 1 + aa1_matrix[, sim_year_idx] * asset1_matrix[, sim_year_idx] + aa2_matrix[, sim_year_idx] * rf 
    
    ## calculating balance matrix
    if (sim_year_idx > 1) {
      bal_matrix[, sim_year_idx] <- (bal_matrix[, sim_year_idx - 1] + salary * (1 - con_matrix[, sim_year_idx])) * ret_matrix[, sim_year_idx]
    } else {
      bal_matrix[, sim_year_idx] <- (bal_start + salary * (1 - con_matrix[, sim_year_idx])) * ret_matrix[, sim_year_idx]
    }
  }
  
  ## bal_matrix is end of year matrix; adjust, append a column in the beginning (i.e. end of year 0)
  bal_matrix <- cbind(rep(bal_start, vertical_size), bal_matrix)
  
  ## con_matrix is start of year matrix; adjust, append a column in the end (i.e. age_terminal, consumes everything)
  con_matrix <- cbind(con_matrix, rep(1, vertical_size))
  
  ## calculate pre-tax consumption amount
  consum_matrix <- salary * con_matrix 
  consum_matrix[, sim_year_total + 1] <- bal_matrix[, sim_year_total + 1] / annuity1 # adjusting annuity payment in column for age 65. Needs to be treated accordingly
  
  ## calculate post-tax consumption amount
  posttax_consum_matrix <- consum_matrix - income_tax(consum_matrix)
  
  ## discount factor (for utility)
  dis_fac <- seq(1:sim_year_total)
  dis_fac1 <- c(1, real_discount ^ dis_fac)
  dis_fac_matrix <- matrix(rep(dis_fac1, vertical_size), byrow = TRUE, ncol = sim_year_total + 1) # start of year
  
  ## calculating total utility 
  utility_matrix <- (posttax_consum_matrix ^ (1 - rho)) / (1 - rho)
  utility_matrix[, sim_year_total + 1] <- utility_matrix[, sim_year_total + 1] * annuity2
  
  utility_matrix <- utility_matrix * dis_fac_matrix # adjust for discounting
  
  ## calculating median and mean of ending balances
  end_bal_s <- bal_matrix[,sim_year_total + 1]
  median_end_bal <- quantile(end_bal_s, 0.5)
  mean_end_bal <- mean(end_bal_s)
  
  ## Calculating averages
  age <- seq(25, 65, 1)
  mean_aa1 <- c(apply(aa1_matrix, 2, mean),0)
  mean_aa2 <- c(apply(aa2_matrix, 2, mean),0)
  mean_bal <- apply(bal_matrix, 2, mean)
  mean_con <- apply(con_matrix, 2, mean)
  mean_consum <- apply(consum_matrix, 2, mean)
  mean_posttax_consum <- apply(posttax_consum_matrix, 2, mean)
  mean_decisions_s <- cbind(age, mean_bal, mean_aa1, mean_aa2, mean_con, mean_consum, mean_posttax_consum)
  colnames(mean_decisions_s) <- c("Age", "Avg_Bal", "Avg_E", "Avg_AC",  "Avg_conprop", "Avg_pretax_cons", "Avg_posttax_cons")
  
  ## Calculating medians
  median_aa1 <- c(apply(aa1_matrix, 2, median),0)
  median_aa2 <- c(apply(aa2_matrix, 2, median),0)
  median_bal <- apply(bal_matrix, 2, median)
  median_con <- apply(con_matrix, 2, median)
  median_consum <- apply(consum_matrix, 2, median)
  median_posttax_consum <- apply(posttax_consum_matrix, 2, median)
  median_decisions_s <- cbind(age, median_bal, median_aa1, median_aa2, median_con, median_consum, median_posttax_consum)
  colnames(median_decisions_s) <- c("Age", "Med_Bal", "Med_E", "Med_AC", "Med_conprop", "Med_pretax_cons", "Med_posttax_cons")
  
  ## Calculating 25th percentile (lower quartile)
  lowquar_aa1 <- c(colQuantiles(aa1_matrix, probs = 0.05),0) 
  lowquar_aa2 <- c(colQuantiles(aa2_matrix, probs = 0.05),0) 
  lowquar_bal <- colQuantiles(bal_matrix, probs = 0.05)
  lowquar_con <- colQuantiles(con_matrix, probs = 0.05)
  lowquar_consum <- colQuantiles(consum_matrix, probs = 0.05)
  lowquar_posttax_consum <- colQuantiles(posttax_consum_matrix, probs = 0.05)
  lowquar_decisions_s <- cbind(age, lowquar_bal, lowquar_aa1, lowquar_aa2, lowquar_con, lowquar_consum, lowquar_posttax_consum)
  colnames(lowquar_decisions_s) <- c("Age", "lowquar_Bal", "lowquar_E", "lowquar_AC", "lowquar_conprop", "lowquar_pretax_cons", "lowquar_posttax_cons")
  
  ## Calculating 75th percentile (upper quartile)
  upquar_aa1 <- c(colQuantiles(aa1_matrix, probs = 0.95),0) 
  upquar_aa2 <- c(colQuantiles(aa2_matrix, probs = 0.95),0) 
  upquar_bal <- colQuantiles(bal_matrix, probs = 0.95)
  upquar_con <- colQuantiles(con_matrix, probs = 0.95)
  upquar_consum <- colQuantiles(consum_matrix, probs = 0.95)
  upquar_posttax_consum <- colQuantiles(posttax_consum_matrix, probs = 0.95)
  upquar_decisions_s <- cbind(age, upquar_bal, upquar_aa1, upquar_aa2, upquar_con, upquar_consum, upquar_posttax_consum)
  colnames(upquar_decisions_s) <- c("Age", "upquar_Bal", "upquar_E", "upquar_AC", "upquar_conprop", "upquar_pretax_cons", "upquar_posttax_cons")
  
  ## Utility and CEC
  sim_utility <- apply(utility_matrix, 1, sum)
  sim_pseudoConsum <- (sim_utility * (1 - rho)) ^ (1 / (1 - rho))

  ## Construct a projection vector
  Util.DF <- data.frame(Start.Bal = bal_start, Median.Final.Bal = median_end_bal, Mean.Final.Bal = mean_end_bal, Average.CEC = mean(sim_pseudoConsum), VaR.Final.Bal1 =  as.numeric(quantile(end_bal_s, 0.05)), CVaR.Final.Bal1 =  as.numeric(mean(end_bal_s[end_bal_s <= quantile(end_bal_s, 0.05)])), VaR.Final.Bal2 =  as.numeric(quantile(end_bal_s, 0.01)), CVaR.Final.Bal2 =  as.numeric(mean(end_bal_s[end_bal_s <= quantile(end_bal_s, 0.01)])))
  
  ## Save final fund balances
  save(end_bal_s, file = "ending_bal_s.RData")
  
  # ## Need to export interesting results
  # str_prefix <- paste("projection_CRRA(S)_start_age", age_start, "bal", bal_start, "_", sep = "_") # prefix (user to specify)
  # interest_outputs_s <- c("mean_decisions_s", "median_decisions_s","lowquar_decisions_s","upquar_decisions_s" )
  # l.df <- lapply(interest_outputs_s, function(x) write.csv(get(x), paste0(str_prefix, x, ".csv"), row.names = F))
  # save(median_decisions_s, file = paste0("worker_med_projection_bal_", bal_start, ".RData"))##
  # save(lowquar_decisions_s, file = paste0("worker_lowquar_projection_bal_", bal_start, ".RData"))##
  # save(upquar_decisions_s, file = paste0("worker_upquar_projection_bal_", bal_start, ".RData"))##
  
  assign(paste0("Util_DF_CRRA_s"), Util.DF, envir = .GlobalEnv)
}


# projection function for passive fund manager
proj_CRRA_h1 <- function(bal_start, aa_decision, returnData, sim_year_total, vertical_size, alloboundL, alloboundU){ 
  
  rets <- gen_mats(returnData)
  asset1_matrix <- rets$asset1_matrix
  
  aa1_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year 
  aa2_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year 
  con_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year
  ret_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # end of year
  bal_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # end of year
  consum_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
  posttax_consum_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
  utility_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year

  for (sim_year_idx in 1:sim_year_total) {

    age_idx <- age_start + sim_year_idx - 1
    
    cat(paste0("CRRA projection(H1) - age ", age_idx, " in progress\n"))
    
    LB_C <- ifelse(age_idx < age_tr, LB_C1_rate, LB_C2_rate)

    aa_temp <- aa_decision[rev(which(aa_decision$age == age_idx)), c(2,4,6)]
    
    ## assign asset weights/consprop/value_m based on interpolated value
    if (sim_year_idx > 1) {
      aa1_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$wA, Bals = bal_matrix[, sim_year_idx - 1], Outlier, LB = alloboundL, UB = alloboundU)
      aa2_matrix[, sim_year_idx] <- 1 - aa1_matrix[, sim_year_idx]
      con_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_matrix[, sim_year_idx - 1], Outlier, LB = LB_C, UB = 1)
    } else {
      aa1_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$wA, Bals = bal_start, Outlier, LB = alloboundL, UB = alloboundU)
      aa2_matrix[, sim_year_idx] <- 1 - aa1_matrix[, sim_year_idx]
      con_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_start, Outlier, LB = LB_C, UB = 1)
      }
    
    ## calculating investment return
    ret_matrix[, sim_year_idx] <- 1 + aa1_matrix[, sim_year_idx] * asset1_matrix[, sim_year_idx] + aa2_matrix[, sim_year_idx] * rf
    
    ## calculating balance matrix
    if (sim_year_idx > 1) {
      bal_matrix[, sim_year_idx] <- (bal_matrix[, sim_year_idx - 1] + salary * (1 - con_matrix[, sim_year_idx])) * ret_matrix[, sim_year_idx]
    } else {
      bal_matrix[, sim_year_idx] <- (bal_start + salary * (1 - con_matrix[, sim_year_idx])) * ret_matrix[, sim_year_idx]
    }
  }
  
  ## bal_matrix is end of year matrix; adjust, append a column in the beginning (i.e. end of year 0)
  bal_matrix <- cbind(rep(bal_start, vertical_size), bal_matrix)
  
  ## con_matrix is start of year matrix; adjust, append a column in the end (i.e. age_terminal, consumes everything)
  con_matrix <- cbind(con_matrix, rep(1, vertical_size)) 
  
  ## pre-tax consumption amount
  consum_matrix <- salary * con_matrix 
  consum_matrix[, sim_year_total + 1] <- bal_matrix[, sim_year_total + 1] / annuity1 # adjusting annuity payment in column for age 65. Needs to be treated accordingly
  
  ## calculate post-tax consumption amount
  posttax_consum_matrix <- consum_matrix - income_tax(consum_matrix)
  
  ## discount factor (for utility)
  dis_fac <- seq(1:sim_year_total)
  dis_fac1 <- c(1, real_discount ^ dis_fac)
  dis_fac_matrix <- matrix(rep(dis_fac1, vertical_size), byrow = TRUE, ncol = sim_year_total + 1) # start of year
  
  ## calculating total utility 
  utility_matrix <- (posttax_consum_matrix ^ (1 - rho)) / (1 - rho)
  utility_matrix[, sim_year_total + 1] <- utility_matrix[, sim_year_total + 1] * annuity2
  
  utility_matrix <- utility_matrix * dis_fac_matrix # adjust for discounting
  
  ## calculating median and mean of ending balances
  end_bal_h1 <- bal_matrix[,sim_year_total + 1]
  median_end_bal1 <- quantile(end_bal_h1, 0.5)
  mean_end_bal1 <- mean(end_bal_h1)
  
  ## Calculating averages
  age <- seq(25, 65, 1)
  mean_aa1 <- c(apply(aa1_matrix, 2, mean),0)
  mean_aa2 <- c(apply(aa2_matrix, 2, mean),0)
  mean_bal <- apply(bal_matrix, 2, mean)
  mean_con <- apply(con_matrix, 2, mean)
  mean_consum <- apply(consum_matrix, 2, mean)
  mean_posttax_consum <- apply(posttax_consum_matrix, 2, mean)
  mean_decisions_h1 <- cbind(age, mean_bal, mean_aa1, mean_aa2, mean_con, mean_consum, mean_posttax_consum)
  colnames(mean_decisions_h1) <- c("Age", "Avg_Bal", "Avg_E", "Avg_AC",  "Avg_conprop", "Avg_pretax_cons", "Avg_posttax_cons")
  
  ## Calculating medians
  median_aa1 <- c(apply(aa1_matrix, 2, median),0)
  median_aa2 <- c(apply(aa2_matrix, 2, median),0)
  median_bal <- apply(bal_matrix, 2, median)
  median_con <- apply(con_matrix, 2, median)
  median_consum <- apply(consum_matrix, 2, median)
  median_posttax_consum <- apply(posttax_consum_matrix, 2, median)
  median_decisions_h1 <- cbind(age, median_bal, median_aa1, median_aa2, median_con, median_consum, median_posttax_consum)
  colnames(median_decisions_h1) <- c("Age", "Med_Bal", "Med_E", "Med_AC", "Med_conprop", "Med_pretax_cons", "Med_posttax_cons")
  
  ## Calculating 25th percentile (lower quartile)
  lowquar_aa1 <- c(colQuantiles(aa1_matrix, probs = 0.05),0) 
  lowquar_aa2 <- c(colQuantiles(aa2_matrix, probs = 0.05),0) 
  lowquar_bal <- colQuantiles(bal_matrix, probs = 0.05)
  lowquar_con <- colQuantiles(con_matrix, probs = 0.05)
  lowquar_consum <- colQuantiles(consum_matrix, probs = 0.05)
  lowquar_posttax_consum <- colQuantiles(posttax_consum_matrix, probs = 0.05)
  lowquar_decisions_h1 <- cbind(age, lowquar_bal, lowquar_aa1, lowquar_aa2, lowquar_con, lowquar_consum, lowquar_posttax_consum)
  colnames(lowquar_decisions_h1) <- c("Age", "lowquar_Bal", "lowquar_E", "lowquar_AC", "lowquar_conprop", "lowquar_pretax_cons", "lowquar_posttax_cons")
  
  ## Calculating 75th percentile (upper quartile)
  upquar_aa1 <- c(colQuantiles(aa1_matrix, probs = 0.95),0) 
  upquar_aa2 <- c(colQuantiles(aa2_matrix, probs = 0.95),0) 
  upquar_bal <- colQuantiles(bal_matrix, probs = 0.95)
  upquar_con <- colQuantiles(con_matrix, probs = 0.95)
  upquar_consum <- colQuantiles(consum_matrix, probs = 0.95)
  upquar_posttax_consum <- colQuantiles(posttax_consum_matrix, probs = 0.95)
  upquar_decisions_h1 <- cbind(age, upquar_bal, upquar_aa1, upquar_aa2, upquar_con, upquar_consum, upquar_posttax_consum)
  colnames(upquar_decisions_h1) <- c("Age", "upquar_Bal", "upquar_E", "upquar_AC", "upquar_conprop", "upquar_pretax_cons", "upquar_posttax_cons")
  
  ## Utility and CEC
  sim_utility <- apply(utility_matrix, 1, sum)
  sim_pseudoConsum <- (sim_utility * (1 - rho)) ^ (1 / (1 - rho))
  
  ## Construct a projection vector
  Util.DF <- data.frame(Start.Bal = bal_start, Median.Final.Bal = median_end_bal1, Mean.Final.Bal = mean_end_bal1, Average.CEC = mean(sim_pseudoConsum), VaR.Final.Bal1 =  as.numeric(quantile(end_bal_h1, 0.05)), CVaR.Final.Bal1 =  as.numeric(mean(end_bal_h1[end_bal_h1 <= quantile(end_bal_h1, 0.05)])), VaR.Final.Bal2 =  as.numeric(quantile(end_bal_h1, 0.01)), CVaR.Final.Bal2 =  as.numeric(mean(end_bal_h1[end_bal_h1 <= quantile(end_bal_h1, 0.01)])))
  
  ## Save final fund balances
  save(end_bal_h1, file = "ending_bal_h1.RData")
  
  # ## Need to export interesting results
  # str_prefix <- paste("projection_CRRA(H1)_start_age", age_start, "bal", bal_start, "_", sep = "_") # prefix (user to specify)
  # interest_outputs_h1 <- c("mean_decisions_h1", "median_decisions_h1", "lowquar_decisions_h1", "upquar_decisions_h1")
  # l.df <- lapply(interest_outputs_h1, function(x) write.csv(get(x), paste0(str_prefix, x, ".csv"), row.names = F))
  # save(median_decisions_h1, file = paste0("manager1_med_projection_bal_", bal_start, ".RData"))##
  # save(lowquar_decisions_h1, file = paste0("manager1_lowquar_projection_bal_", bal_start, ".RData"))##
  # save(upquar_decisions_h1, file = paste0("manager1_upquar_projection_bal_", bal_start, ".RData"))##
  
  assign(paste0("Util_DF_CRRA_h1"), Util.DF, envir = .GlobalEnv)
}

# projection function for active fund manager
proj_CRRA_h2 <- function(bal_start, aa_decision, returnData, sim_year_total, vertical_size, alloboundL, alloboundU){ 
  
  rets <- gen_mats(returnData)
  asset1_matrix <- rets$asset1_matrix
  
  aa1_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year 
  aa2_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year 
  con_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year
  ret_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # end of year
  bal_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # end of year
  consum_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
  posttax_consum_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
  utility_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
  
  for (sim_year_idx in 1:sim_year_total) {
    
    
    age_idx <- age_start + sim_year_idx - 1
    
    cat(paste0("CRRA projection(H2) - age ", age_idx, " in progress\n"))
    
    LB_C <- ifelse(age_idx < age_tr, LB_C1_rate, LB_C2_rate)
    
    aa_temp <- aa_decision[rev(which(aa_decision$age == age_idx)), c(2,4,6)]
    
    ## assign asset weights/consprop/value_m based on interpolated value
    if (sim_year_idx > 1) {
      aa1_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$wA, Bals = bal_matrix[, sim_year_idx - 1], Outlier, LB = alloboundL, UB = alloboundU)
      aa2_matrix[, sim_year_idx] <- 1 - aa1_matrix[, sim_year_idx]
      con_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_matrix[, sim_year_idx - 1], Outlier, LB = LB_C, UB = 1)
    } else {
      aa1_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$wA, Bals = bal_start, Outlier, LB = alloboundL, UB = alloboundU)
      aa2_matrix[, sim_year_idx] <- 1 - aa1_matrix[, sim_year_idx]
      con_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_start, Outlier, LB = LB_C, UB = 1)
    }
    
    ## calculating investment return
    ret_matrix[, sim_year_idx] <- 1 + aa1_matrix[, sim_year_idx] * asset1_matrix[, sim_year_idx] + aa2_matrix[, sim_year_idx] * rf
    
    ## calculating balance matrix
    if (sim_year_idx > 1) {
      bal_matrix[, sim_year_idx] <- (bal_matrix[, sim_year_idx - 1] + salary * (1 - con_matrix[, sim_year_idx])) * ret_matrix[, sim_year_idx]
    } else {
      bal_matrix[, sim_year_idx] <- (bal_start + salary * (1 - con_matrix[, sim_year_idx])) * ret_matrix[, sim_year_idx]
    }
  }
  
  ## bal_matrix is end of year matrix; adjust, append a column in the beginning (i.e. end of year 0)
  bal_matrix <- cbind(rep(bal_start, vertical_size), bal_matrix)
  
  ## con_matrix is start of year matrix; adjust, append a column in the end (i.e. age_terminal, consumes everything)
  con_matrix <- cbind(con_matrix, rep(1, vertical_size)) 
  
  ## pre-tax consumption amount
  consum_matrix <- salary * con_matrix 
  consum_matrix[, sim_year_total + 1] <- bal_matrix[, sim_year_total + 1] / annuity1 # adjusting annuity payment in column for age 65. Needs to be treated accordingly
  
  ## calculate post-tax consumption amount
  posttax_consum_matrix <- consum_matrix - income_tax(consum_matrix)
  
  ## discount factor (for utility)
  dis_fac <- seq(1:sim_year_total)
  dis_fac1 <- c(1, real_discount ^ dis_fac)
  dis_fac_matrix <- matrix(rep(dis_fac1, vertical_size), byrow = TRUE, ncol = sim_year_total + 1) # start of year
  
  ## calculating total utility 
  utility_matrix <- (posttax_consum_matrix ^ (1 - rho)) / (1 - rho)
  utility_matrix[, sim_year_total + 1] <- utility_matrix[, sim_year_total + 1] * annuity2
  
  utility_matrix <- utility_matrix * dis_fac_matrix # adjust for discounting
  
  ## calculating median and mean of ending balances
  end_bal_h2 <- bal_matrix[,sim_year_total + 1]
  median_end_bal <- quantile(end_bal_h2, 0.5)
  mean_end_bal <- mean(end_bal_h2)
  
  ## Calculating averages
  age <- seq(25, 65, 1)
  mean_aa1 <- c(apply(aa1_matrix, 2, mean),0)
  mean_aa2 <- c(apply(aa2_matrix, 2, mean),0)
  mean_bal <- apply(bal_matrix, 2, mean)
  mean_con <- apply(con_matrix, 2, mean)
  mean_consum <- apply(consum_matrix, 2, mean)
  mean_posttax_consum <- apply(posttax_consum_matrix, 2, mean)
  mean_decisions_h2 <- cbind(age, mean_bal, mean_aa1, mean_aa2, mean_con, mean_consum, mean_posttax_consum)
  colnames(mean_decisions_h2) <- c("Age", "Avg_Bal", "Avg_E", "Avg_AC",  "Avg_conprop", "Avg_pretax_cons", "Avg_posttax_cons")
  
  ## Calculating medians
  median_aa1 <- c(apply(aa1_matrix, 2, median),0)
  median_aa2 <- c(apply(aa2_matrix, 2, median),0)
  median_bal <- apply(bal_matrix, 2, median)
  median_con <- apply(con_matrix, 2, median)
  median_consum <- apply(consum_matrix, 2, median)
  median_posttax_consum <- apply(posttax_consum_matrix, 2, median)
  median_decisions_h2 <- cbind(age, median_bal, median_aa1, median_aa2, median_con, median_consum, median_posttax_consum)
  colnames(median_decisions_h2) <- c("Age", "Med_Bal", "Med_E", "Med_AC", "Med_conprop", "Med_pretax_cons", "Med_posttax_cons")
  
  
  ## Calculating 25th percentile (lower quartile)
  lowquar_aa1 <- c(colQuantiles(aa1_matrix, probs = 0.05),0) 
  lowquar_aa2 <- c(colQuantiles(aa2_matrix, probs = 0.05),0) 
  lowquar_bal <- colQuantiles(bal_matrix, probs = 0.05)
  lowquar_con <- colQuantiles(con_matrix, probs = 0.05)
  lowquar_consum <- colQuantiles(consum_matrix, probs = 0.05)
  lowquar_posttax_consum <- colQuantiles(posttax_consum_matrix, probs = 0.05)
  lowquar_decisions_h2 <- cbind(age, lowquar_bal, lowquar_aa1, lowquar_aa2, lowquar_con, lowquar_consum, lowquar_posttax_consum)
  colnames(lowquar_decisions_h2) <- c("Age", "lowquar_Bal", "lowquar_E", "lowquar_AC", "lowquar_conprop", "lowquar_pretax_cons", "lowquar_posttax_cons")
  
  ## Calculating 75th percentile (upper quartile)
  upquar_aa1 <- c(colQuantiles(aa1_matrix, probs = 0.95),0) 
  upquar_aa2 <- c(colQuantiles(aa2_matrix, probs = 0.95),0) 
  upquar_bal <- colQuantiles(bal_matrix, probs = 0.95)
  upquar_con <- colQuantiles(con_matrix, probs = 0.95)
  upquar_consum <- colQuantiles(consum_matrix, probs = 0.95)
  upquar_posttax_consum <- colQuantiles(posttax_consum_matrix, probs = 0.95)
  upquar_decisions_h2 <- cbind(age, upquar_bal, upquar_aa1, upquar_aa2, upquar_con, upquar_consum, upquar_posttax_consum)
  colnames(upquar_decisions_h2) <- c("Age", "upquar_Bal", "upquar_E", "upquar_AC", "upquar_conprop", "upquar_pretax_cons", "upquar_posttax_cons")
  
  ## Utility and CEC
  sim_utility <- apply(utility_matrix, 1, sum)
  sim_pseudoConsum <- (sim_utility * (1 - rho)) ^ (1 / (1 - rho))
  
  ## Construct a projection vector
  Util.DF <- data.frame(Start.Bal = bal_start, Median.Final.Bal = median_end_bal, Mean.Final.Bal = mean_end_bal, Average.CEC = mean(sim_pseudoConsum), VaR.Final.Bal1 =  as.numeric(quantile(end_bal_h2, 0.05)), CVaR.Final.Bal1 =  as.numeric(mean(end_bal_h2[end_bal_h2 <= quantile(end_bal_h2, 0.05)])), VaR.Final.Bal2 =  as.numeric(quantile(end_bal_h2, 0.01)), CVaR.Final.Bal2 =  as.numeric(mean(end_bal_h2[end_bal_h2 <= quantile(end_bal_h2, 0.01)])))
  
  ## Save final fund balances
  save(end_bal_h2, file = "ending_bal_h2.RData")
  
  # ## Need to export interesting results
  # str_prefix <- paste("projection_CRRA(H2)_start_age", age_start, "bal", bal_start, "_", sep = "_") # prefix (user to specify)
  # interest_outputs_h2 <- c("mean_decisions_h2", "median_decisions_h2", "lowquar_decisions_h2", "upquar_decisions_h2")
  # l.df <- lapply(interest_outputs_h2, function(x) write.csv(get(x), paste0(str_prefix, x, ".csv"), row.names = F))
  # save(median_decisions_h2, file = paste0("manager2_med_projection_bal_", bal_start, ".RData"))##
  # save(lowquar_decisions_h2, file = paste0("manager2_lowquar_projection_bal_", bal_start, ".RData"))##
  # save(upquar_decisions_h2, file = paste0("manager2_upquar_projection_bal_", bal_start, ".RData"))##
  
  assign(paste0("Util_DF_CRRA_h2"), Util.DF, envir = .GlobalEnv)
}


# Projection function for benchmark investment
proj_CRRA_b <- function(bal_start, aa_decision, returnData, sim_year_total, vertical_size, alloboundL, alloboundU){
  
  rets <- gen_mats(returnData)
  asset1_matrix <- rets$asset1_matrix
  
  aa1_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year
  aa2_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year
  con_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year
  ret_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # end of year
  bal_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # end of year
  consum_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
  posttax_consum_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
  utility_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
  
  for (sim_year_idx in 1:sim_year_total) {

    age_idx <- age_start + sim_year_idx - 1
    
    cat(paste0("CRRA projection(B) - age ", age_idx, " in progress\n"))
    
    LB_C <- ifelse(age_idx < age_tr, LB_C1_rate, LB_C2_rate)
    
    aa_temp <- aa_decision[rev(which(aa_decision$age == age_idx)), c(2,6,7)]
    
    ## assign asset weights/consprop based on interpolated value
    if (sim_year_idx > 1) {
      aa1_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$eqprop, Bals = bal_matrix[, sim_year_idx - 1], Outlier, LB = alloboundL, UB = alloboundU)
      aa2_matrix[, sim_year_idx] <- 1 - aa1_matrix[, sim_year_idx]
      con_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_matrix[, sim_year_idx - 1], Outlier, LB = LB_C, UB = 1)
      
    } else {
      aa1_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$eqprop, Bals = bal_start, Outlier, LB = alloboundL, UB = alloboundU)
      aa2_matrix[, sim_year_idx] <- 1 - aa1_matrix[, sim_year_idx]
      con_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_start, Outlier, LB = LB_C, UB = 1)
    }
    
    ## calculating investment return
    ret_matrix[, sim_year_idx] <- 1 + aa1_matrix[, sim_year_idx] * asset1_matrix[, sim_year_idx] + aa2_matrix[, sim_year_idx] * rf 
    
    ## calculating balance matrix
    if (sim_year_idx > 1) {
      bal_matrix[, sim_year_idx] <- (bal_matrix[, sim_year_idx - 1] + salary * (1 - con_matrix[, sim_year_idx])) * ret_matrix[, sim_year_idx]
    } else {
      bal_matrix[, sim_year_idx] <- (bal_start + salary * (1 - con_matrix[, sim_year_idx])) * ret_matrix[, sim_year_idx]
    }
  }
  
  ## bal_matrix is end of year matrix; adjust, append a column in the beginning (i.e. end of year 0)
  bal_matrix <- cbind(rep(bal_start, vertical_size), bal_matrix)
  
  ## con_matrix is start of year matrix; adjust, append a column in the end (i.e. age_terminal, consumes everything)
  con_matrix <- cbind(con_matrix, rep(1, vertical_size))
  
  ## pre-tax consumption amount
  consum_matrix <- salary * con_matrix 
  consum_matrix[, sim_year_total + 1] <- bal_matrix[, sim_year_total + 1] / annuity1 # adjusting annuity payment in column for age 65. Needs to be treated accordingly
  
  ## calculate post-tax consumption amount
  posttax_consum_matrix <- consum_matrix - income_tax(consum_matrix)
  
  ## discount factor (for utility)
  dis_fac <- seq(1:sim_year_total)
  dis_fac1 <- c(1, real_discount ^ dis_fac)
  dis_fac_matrix <- matrix(rep(dis_fac1, vertical_size), byrow = TRUE, ncol = sim_year_total + 1) # start of year
  
  ## calculating total utility 
  utility_matrix <- (posttax_consum_matrix ^ (1 - rho)) / (1 - rho)
  utility_matrix[, sim_year_total + 1] <- utility_matrix[, sim_year_total + 1] * annuity2
  
  utility_matrix <- utility_matrix * dis_fac_matrix # adjust for discounting
  
  ## calculating median and mean of ending balances
  end_bal_b <- bal_matrix[,sim_year_total + 1]
  median_end_bal <- quantile(end_bal_b, 0.5)
  mean_end_bal <- mean(end_bal_b)
  
  ## Calculating averages
  age <- seq(25, 65, 1)
  mean_aa1 <- c(apply(aa1_matrix, 2, mean),0)
  mean_aa2 <- c(apply(aa2_matrix, 2, mean),0)
  mean_bal <- apply(bal_matrix, 2, mean)
  mean_con <- apply(con_matrix, 2, mean)
  mean_consum <- apply(consum_matrix, 2, mean)
  mean_posttax_consum <- apply(posttax_consum_matrix, 2, mean)
  mean_decisions_b <- cbind(age, mean_bal, mean_aa1, mean_aa2, mean_con, mean_consum, mean_posttax_consum)
  colnames(mean_decisions_b) <- c("Age", "Avg_Bal", "Avg_E", "Avg_AC",  "Avg_conprop", "Avg_pretax_cons", "Avg_posttax_cons")
  
  ## Calculating medians
  median_aa1 <- c(apply(aa1_matrix, 2, median),0)
  median_aa2 <- c(apply(aa2_matrix, 2, median),0)
  median_bal <- apply(bal_matrix, 2, median)
  median_con <- apply(con_matrix, 2, median)
  median_consum <- apply(consum_matrix, 2, median)
  median_posttax_consum <- apply(posttax_consum_matrix, 2, median)
  median_decisions_b <- cbind(age, median_bal, median_aa1, median_aa2, median_con, median_consum, median_posttax_consum)
  colnames(median_decisions_b) <- c("Age", "Med_Bal", "Med_E", "Med_AC", "Med_conprop", "Med_pretax_cons", "Med_posttax_cons")
  
  ## Calculating 25th percentile (lower quartile)
  lowquar_aa1 <- c(colQuantiles(aa1_matrix, probs = 0.05),0) 
  lowquar_aa2 <- c(colQuantiles(aa2_matrix, probs = 0.05),0) 
  lowquar_bal <- colQuantiles(bal_matrix, probs = 0.05)
  lowquar_con <- colQuantiles(con_matrix, probs = 0.05)
  lowquar_consum <- colQuantiles(consum_matrix, probs = 0.05)
  lowquar_posttax_consum <- colQuantiles(posttax_consum_matrix, probs = 0.05)
  lowquar_decisions_b <- cbind(age, lowquar_bal, lowquar_aa1, lowquar_aa2, lowquar_con, lowquar_consum, lowquar_posttax_consum)
  colnames(lowquar_decisions_b) <- c("Age", "lowquar_Bal", "lowquar_E", "lowquar_AC", "lowquar_conprop", "lowquar_pretax_cons", "lowquar_posttax_cons")
  
  ## Calculating 75th percentile (upper quartile)
  upquar_aa1 <- c(colQuantiles(aa1_matrix, probs = 0.95),0) 
  upquar_aa2 <- c(colQuantiles(aa2_matrix, probs = 0.95),0) 
  upquar_bal <- colQuantiles(bal_matrix, probs = 0.95)
  upquar_con <- colQuantiles(con_matrix, probs = 0.95)
  upquar_consum <- colQuantiles(consum_matrix, probs = 0.95)
  upquar_posttax_consum <- colQuantiles(posttax_consum_matrix, probs = 0.95)
  upquar_decisions_b <- cbind(age, upquar_bal, upquar_aa1, upquar_aa2, upquar_con, upquar_consum, upquar_posttax_consum)
  colnames(upquar_decisions_b) <- c("Age", "upquar_Bal", "upquar_E", "upquar_AC", "upquar_conprop", "upquar_pretax_cons", "upquar_posttax_cons")
  
  ## Utility and CEC
  sim_utility <- apply(utility_matrix, 1, sum)
  sim_pseudoConsum <- (sim_utility * (1 - rho)) ^ (1 / (1 - rho))
  
  ## Construct a projection vector
  Util.DF <- data.frame(Start.Bal = bal_start, Median.Final.Bal = median_end_bal, Mean.Final.Bal = mean_end_bal, Average.CEC = mean(sim_pseudoConsum), VaR.Final.Bal1 =  as.numeric(quantile(end_bal_b, 0.05)), CVaR.Final.Bal1 =  as.numeric(mean(end_bal_b[end_bal_b <= quantile(end_bal_b, 0.05)])), VaR.Final.Bal2 =  as.numeric(quantile(end_bal_b, 0.01)), CVaR.Final.Bal2 =  as.numeric(mean(end_bal_b[end_bal_b <= quantile(end_bal_b, 0.01)])))
  
  ## Save final fund balances
  save(end_bal_b, file = "ending_bal_b.RData")
  
  # ## Need to export interesting results
  # str_prefix <- paste("projection_CRRA(Benchmark)_start_age", age_start, "bal", bal_start, "_", sep = "_") # prefix (user to specify)
  # interest_outputs_b <- c("mean_decisions_b", "median_decisions_b", "lowquar_decisions_b", "upquar_decisions_b")
  # l.df <- lapply(interest_outputs_b, function(x) write.csv(get(x), paste0(str_prefix, x, ".csv"), row.names = F))
  # save(median_decisions_b, file = paste0("benchmark_med_projection_bal_", bal_start, ".RData"))##
  # save(lowquar_decisions_b, file = paste0("benchmark_lowquar_projection_bal_", bal_start, ".RData"))##
  # save(upquar_decisions_b, file = paste0("benchmark_upquar_projection_bal_", bal_start, ".RData"))##
  
  assign(paste0("Util_DF_CRRA_b"), Util.DF, envir = .GlobalEnv)
}


## calculate extra management fee (given the same CEC value)
#(1) for passive fund manager
x_mgt_fee1_1 <- function(base_cec, f_cec, bal_start, aa_decision, returnData, sim_year_total, vertical_size, alloboundL, alloboundU){
  
  rets <- gen_mats(returnData)
  asset1_matrix <- rets$asset1_matrix
  
  j <- 0
  
  opt_func <- function (y){
    j <- j + 1
    cat(paste0("mgt_fee iter ", j, " in progress\n"))
    
    aa1_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year 
    aa2_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year 
    con_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year 
    ret_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # end of year
    bal_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # end of year
    consum_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
    posttax_consum_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
    utility_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
    
    for (sim_year_idx in 1:sim_year_total) {
      
      
      age_idx <- age_start + sim_year_idx - 1
      
      LB_C <- ifelse(age_idx < age_tr, LB_C1_rate, LB_C2_rate)
      
      aa_temp <- aa_decision[rev(which(aa_decision$age == age_idx)), c(2,4,6)]
      
      # assign asset weights/consprop based on interpolated value
      if (sim_year_idx > 1) {
        aa1_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$wA, Bals = bal_matrix[, sim_year_idx - 1], Outlier, LB = alloboundL, UB = alloboundU)
        aa2_matrix[, sim_year_idx] <- 1 - aa1_matrix[, sim_year_idx]
        con_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_matrix[, sim_year_idx - 1], Outlier, LB = LB_C, UB = 1)
      } else {
        aa1_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$wA, Bals = bal_start, Outlier, LB = alloboundL, UB = alloboundU)
        aa2_matrix[, sim_year_idx] <- 1 - aa1_matrix[, sim_year_idx]
        con_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_start, Outlier, LB = LB_C, UB = 1)
      }
      
      # calculating investment return
      ret_matrix[, sim_year_idx] <- 1 + aa1_matrix[, sim_year_idx] * asset1_matrix[, sim_year_idx] + aa2_matrix[, sim_year_idx] * rf
      
      # calculating balance matrix (including extra management fee y (%))
      if (sim_year_idx > 1) {
        bal_matrix[, sim_year_idx] <- (bal_matrix[, sim_year_idx - 1] + salary * (1 - con_matrix[, sim_year_idx])) * ret_matrix[, sim_year_idx] * (1 - y)
      } else {
        bal_matrix[, sim_year_idx] <- (bal_start + salary * (1 - con_matrix[, sim_year_idx])) * ret_matrix[, sim_year_idx] * (1 - y)
      }
    }
    
    # adjust bal_matrix & con_matrix
    bal_matrix <- cbind(rep(bal_start, vertical_size), bal_matrix)
    con_matrix <- cbind(con_matrix, rep(1, vertical_size)) 
    
    # calculate pre-tax consumption amount
    consum_matrix <- salary * con_matrix 
    consum_matrix[, sim_year_total + 1] <- bal_matrix[, sim_year_total + 1] / annuity1 # adjusting annuity payment in column for age 65. Needs to be treated accordingly
    
    # calculate post-tax consumption amount
    posttax_consum_matrix <- consum_matrix - income_tax(consum_matrix)
    
    # discount factor (for utility)
    dis_fac <- seq(1:sim_year_total)
    dis_fac1 <- c(1, real_discount ^ dis_fac)
    dis_fac_matrix <- matrix(rep(dis_fac1, vertical_size), byrow = TRUE, ncol = sim_year_total + 1) # start of year
    
    # calculating total utility 
    utility_matrix <- (posttax_consum_matrix ^ (1 - rho)) / (1 - rho)
    utility_matrix[, sim_year_total + 1] <- utility_matrix[, sim_year_total + 1] * annuity2
    
    utility_matrix <- utility_matrix * dis_fac_matrix # adjust for discounting
    
    sim_utility <- apply(utility_matrix, 1, sum)
    sim_pseudoConsum <- (sim_utility * (1 - rho)) ^ (1 / (1 - rho))
    
    diff <- mean(sim_pseudoConsum) - base_cec
    diff
  }
  
  set.seed(1)
  y <- uniroot(f = opt_func, interval = c(0, 0.5), f.lower = (f_cec - base_cec), extendInt = "yes", tol = 1e-8)$root
  y
}

#(2) for worker/active fund manager
x_mgt_fee1_2 <- function(base_cec, f_cec, bal_start, aa_decision, returnData, sim_year_total, vertical_size, alloboundL, alloboundU){
  
  rets <- gen_mats(returnData)
  asset1_matrix <- rets$asset1_matrix
  
  j <- 0
  
  opt_func <- function (y){
    j <- j + 1
    cat(paste0("mgt_fee iter ", j, " in progress\n"))
    
    aa1_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year 
    aa2_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year 
    con_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # start of year 
    ret_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # end of year
    bal_matrix <- matrix(ncol = sim_year_total, nrow = vertical_size) # end of year
    consum_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
    posttax_consum_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
    utility_matrix <- matrix(ncol = sim_year_total + 1, nrow = vertical_size) # start of year
    
    for (sim_year_idx in 1:sim_year_total) {
      
      
      age_idx <- age_start + sim_year_idx - 1
      
      LB_C <- ifelse(age_idx < age_tr, LB_C1_rate, LB_C2_rate)
      
      aa_temp <- aa_decision[rev(which(aa_decision$age == age_idx)), c(2,6,7)]
      
      # assign asset weights/consprop based on interpolated value
      if (sim_year_idx > 1) {
        aa1_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$eqprop, Bals = bal_matrix[, sim_year_idx - 1], Outlier, LB = alloboundL, UB = alloboundU)
        aa2_matrix[, sim_year_idx] <- 1 - aa1_matrix[, sim_year_idx]
        con_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_matrix[, sim_year_idx - 1], Outlier, LB = LB_C, UB = 1)
      } else {
        aa1_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$eqprop, Bals = bal_start, Outlier, LB = alloboundL, UB = alloboundU)
        aa2_matrix[, sim_year_idx] <- 1 - aa1_matrix[, sim_year_idx]
        con_matrix[, sim_year_idx] <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_start, Outlier, LB = LB_C, UB = 1)
      }
      
      # calculating investment return
      ret_matrix[, sim_year_idx] <- 1 + aa1_matrix[, sim_year_idx] * asset1_matrix[, sim_year_idx] + aa2_matrix[, sim_year_idx] * rf
      
      # calculating balance matrix (including extra management fee y (%))
      if (sim_year_idx > 1) {
        bal_matrix[, sim_year_idx] <- (bal_matrix[, sim_year_idx - 1] + salary * (1 - con_matrix[, sim_year_idx])) * ret_matrix[, sim_year_idx] * (1 - y)
      } else {
        bal_matrix[, sim_year_idx] <- (bal_start + salary * (1 - con_matrix[, sim_year_idx])) * ret_matrix[, sim_year_idx] * (1 - y)
      }
    }
    
    # adjust bal_matrix & con_matrix
    bal_matrix <- cbind(rep(bal_start, vertical_size), bal_matrix)
    con_matrix <- cbind(con_matrix, rep(1, vertical_size)) 
    
    # calculate pre-tax consumption amount
    consum_matrix <- salary * con_matrix 
    consum_matrix[, sim_year_total + 1] <- bal_matrix[, sim_year_total + 1] / annuity1 # adjusting annuity payment in column for age 65. Needs to be treated accordingly
    
    # calculate post-tax consumption amount
    posttax_consum_matrix <- consum_matrix - income_tax(consum_matrix)
    
    # discount factor (for utility)
    dis_fac <- seq(1:sim_year_total)
    dis_fac1 <- c(1, real_discount ^ dis_fac)
    dis_fac_matrix <- matrix(rep(dis_fac1, vertical_size), byrow = TRUE, ncol = sim_year_total + 1) # start of year
    
    # calculating total utility 
    utility_matrix <- (posttax_consum_matrix ^ (1 - rho)) / (1 - rho)
    utility_matrix[, sim_year_total + 1] <- utility_matrix[, sim_year_total + 1] * annuity2
    
    utility_matrix <- utility_matrix * dis_fac_matrix # adjust for discounting
    
    sim_utility <- apply(utility_matrix, 1, sum)
    sim_pseudoConsum <- (sim_utility * (1 - rho)) ^ (1 / (1 - rho))
    
    diff <- mean(sim_pseudoConsum) - base_cec
    diff
  }
  
  set.seed(1)
  y <- uniroot(f = opt_func, interval = c(0, 0.5), f.lower = (f_cec - base_cec), extendInt = "yes", tol = 1e-8)$root
  y
}

