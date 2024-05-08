#================H1 Case BASE1 (allocation range is [-0.16,1.16])=========================#

## Initialisation  
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("opt_fns.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores()-1,outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))

#importing target allocation data
alloData <- read.csv("Target_asset_allocation.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 4 
l <- 0.7 # replacement ratio of post-retirement yearly consumption level
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
Wt <- 0.5 # weight assigned to benchmark A
WT <- 1 # weight assigned to target B
C_star <- salary - income_tax(salary) # consumption target (pre-retirement)
B_0 <- l * C_star # post-tax post-retirement consumption target

#starting balance at each age
factor <- 1.04
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
max_start_bal <- 3000000 ##maximum starting balance, as balances above these have same outcomes
start_bal_list$Start_Bal <- ifelse(start_bal_list$Start_Bal > max_start_bal, max_start_bal, start_bal_list$Start_Bal)
bal_spacing <- 5000 ## how fine the balance grid should be (max = 10000)

##Adjustments for Resolution
numC <- 51       # resolution for consumption vector 
numC1 <- 101     # resolution for consumption vector (for age 64)
numPi <- 6601   # resolution for equity allocation % vector (incremental value=(1.16+0.16)/6600=0.0002)
numPi1 <- 13201  # resolution for equity allocation % vector (incremental value=(1.16+0.16)/13200=0.0001) (for age 64)

width1 <- 1000 # first reoptimisation on consumption between +/- width1
width2 <- 100  # second reoptimisation on consumption between +/- width2
width3 <- 50   # third reoptimisation on consumption between +/- width3

#setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1.16 
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- -0.16 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation controls
controlDE1 <- list(reltol = 1e-10, itermax = 300)
controlDE2 <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()

## at terminal age - 1 (age 64)
  age_x <- age_retire - 1
  cat(paste("age_x", age_x, "in progress \n"))
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  A <- itRetA * TAA_1[nrow(alloData)] + itRetB * TAA_2[nrow(alloData)] #set target return
  A_mat <- matrix(rep(A, numPi1), byrow = TRUE, ncol = length(A), nrow =  numPi1)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons 
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[1] <- length(balances)
  
  # set consumption vector
  consum_vec <- seq(LB_C, salary, length.out = numC1)
  
  # compute prob 1 matrix
  eqprop_vec = seq(L_bound, U_bound, length.out = numPi1)
  ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB
  vec <- rep(1,length(consum_vec))
  count1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(count1_mat) / length(itRetA)

  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
   
    itBalStart <- balances[y]
    set.seed(1)
    
    #optimisation for asset allocation
    eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
      inc_mat <- itBalEnd_mat / annuity1
      posttax_inc_mat <- inc_mat - income_tax(inc_mat)
      count2_mat <- (posttax_inc_mat - B_0) > 0
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- WT * prob2
      value <- value1 + value2
      
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val[j] <- value[max(which(value == max(value)))]
      val1[j] <- value1[max(which(value == max(value)))]
      val2[j] <- value2[max(which(value == max(value)))]
    }
    
    inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_pre_terminal, 
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE1)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]
    
    ################################
    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC1)
    
    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
      inc_mat <- itBalEnd_mat / annuity1
      posttax_inc_mat <- inc_mat - income_tax(inc_mat)
      count2_mat <- (posttax_inc_mat - B_0) > 0
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- WT * prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val[j] <- value[max(which(value == max(value)))]
      val1[j] <- value1[max(which(value == max(value)))]
      val2[j] <- value2[max(which(value == max(value)))]
    }
    
    inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_pre_terminal, 
                         lower = B1_p,##
                         upper = B2_p,
                         control = controlDE1)
    optParam <- listOptim$optim$bestmem
    cons2 <- salary * optParam[1]
    
    # re-optimisation parameter [2]
    B3 <- cons2 - width2
    B3 <- ifelse(B3 < LB_C, LB_C, B3)
    B4 <- cons2 + width2
    B4 <- ifelse(B4 > salary, salary, B4)
    B3_p = B3 / salary
    B4_p = B4 / salary
    consum_vec <- seq(B3, B4, length.out = numC1)
    
    #re-optimisation for asset allocation [2]
    eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
      inc_mat <- itBalEnd_mat / annuity1
      posttax_inc_mat <- inc_mat - income_tax(inc_mat)
      count2_mat <- (posttax_inc_mat - B_0) > 0
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- WT * prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val[j] <- value[max(which(value == max(value)))]
      val1[j] <- value1[max(which(value == max(value)))]
      val2[j] <- value2[max(which(value == max(value)))]
    }
    
    inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #re-optimisation for consumption [2]
    listOptim <- DEoptim(fn = utility_pre_terminal, 
                         lower = B3_p,##
                         upper = B4_p,
                         control = controlDE1)
    optParam <- listOptim$optim$bestmem
    cons3 <- salary * optParam[1]
    
    ################################
    # re-optimisation parameter [3]
    B5 <- cons3 - width3
    B5 <- ifelse(B5 < LB_C, LB_C, B5)
    B6 <- cons3 + width3
    B6 <- ifelse(B6 > salary, salary, B6)
    B5_p = B5 / salary
    B6_p = B6 / salary
    consum_vec <- seq(B5, B6, length.out = numC1)
    
    #re-optimisation for asset allocation [3]
    eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
      inc_mat <- itBalEnd_mat / annuity1
      posttax_inc_mat <- inc_mat - income_tax(inc_mat)
      count2_mat <- (posttax_inc_mat - B_0) > 0
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- WT * prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val[j] <- value[max(which(value == max(value)))]
      val1[j] <- value1[max(which(value == max(value)))]
      val2[j] <- value2[max(which(value == max(value)))]
    }
    
    inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2
    valfun2_ord = valfun2[order(consum)]  # small to large
    
    #re-optimisation for consumption [3]
    listOptim <- DEoptim(fn = utility_pre_terminal, 
                         lower = B5_p,##
                         upper = B6_p,
                         control = controlDE1)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]
    
    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = WT)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)
    
    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[nrow(alloData)], 
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility, 
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[1]  

  
# pre retirement ages (from age 25 to age 63)

# eqprop vector & return matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB

#age_x = 63
for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))

  start.time <- Sys.time()

  Upper <- WT + (age_retire - age_x - 1) * Wt
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[(age_retire - age_x)] <- length(balances)

  consum_vec <- seq(LB_C, salary, length.out = numC) #consumption vector

  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))

  itPrevBal = resultsPrev$balance
  itPrevBal_ord = itPrevBal[order(itPrevBal)]  # small to large
  itPrevpseudoCons = resultsPrev$pseudoConsum
  itPrevpseudoCons_ord = itPrevpseudoCons[order(itPrevBal)]  # small to large
  itPrevVal = resultsPrev$value_manager
  itPrevVal_ord = itPrevVal[order(itPrevBal)]  # small to large

  schu_spline1 <- Schumaker(x = itPrevBal_ord, y = itPrevpseudoCons_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline for pseudoconsum at age x+1

  A <- TAA_1[(age_x - 24)] * itRetA + TAA_2[(age_x - 24)] * itRetB  # set target return
  A_mat = matrix(rep(A, numPi), byrow = TRUE, ncol = length(A), nrow =  numPi)

  # compute prob1 matrix
  vec <- rep(1,length(consum_vec))
  cal1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(cal1_mat) / length(itRetA)

  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    set.seed(1)
    itBalStart <- balances[y]

    #optimisation for asset allocation
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL

    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)

    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large

    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_xx,
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]

    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC)

    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL

    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }  
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)

    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1_vec
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2_vec
    valfun2_ord = valfun2[order(consum)]  # small to large

    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_xx,
                         lower = B1_p,#
                         upper = B2_p,
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]

    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = Upper)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)

    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[(age_x - 24)],
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility,
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)

  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}


all.data.list <- paste0("result_x", age_start:(age_retire - 1))##
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t7 <- do.call(rbind, l.df)

save(data_vfi_t7, file = paste0("mgmt_H1_rho_",rho, "salary", salary, "_base1_v5.1.RData"))##
write.csv(data_vfi_t7, paste0("mgmt_H1_rho_",rho, "salary", salary, "_base1_v5.1.csv"), row.names = F)##

###########################################################################################################################################
#================H1 Case BASE2 (allocation range is [0,1])=========================#

## Initialisation  
rm(list=ls())

source("opt_fns.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores()-1,outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))

#importing target allocation data
alloData <- read.csv("Target_asset_allocation.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 4 
l <- 0.7 # replacement ratio of post-retirement yearly consumption level
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
Wt <- 0.5 # weight assigned to benchmark A
WT <- 1 # weight assigned to target B
C_star <- salary - income_tax(salary) # consumption target (pre-retirement)
B_0 <- l * C_star # post-tax post-retirement consumption target

#starting balance at each age
factor <- 1.04
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
max_start_bal <- 3000000 ##maximum starting balance, as balances above these have same outcomes
start_bal_list$Start_Bal <- ifelse(start_bal_list$Start_Bal > max_start_bal, max_start_bal, start_bal_list$Start_Bal)
bal_spacing <- 5000 ## how fine the balance grid should be (max = 10000)

##Adjustments for Resolution
numC <- 51       # resolution for consumption vector 
numC1 <- 101     # resolution for consumption vector (for age 64)
numPi <- 5001   # resolution for equity allocation % vector (incremental value=(1+0)/5000=0.0002)
numPi1 <- 10001  # resolution for equity allocation % vector (incremental value=(1+0)/10000=0.0001) (for age 64)

width1 <- 1000 # first reoptimisation on consumption between +/- width1
width2 <- 100  # second reoptimisation on consumption between +/- width2
width3 <- 50   # third reoptimisation on consumption between +/- width3

#setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1 
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- 0
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation controls
controlDE1 <- list(reltol = 1e-10, itermax = 300)
controlDE2 <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()
###################################
## at terminal age - 1 (age 64)
age_x <- age_retire - 1
cat(paste("age_x", age_x, "in progress \n"))
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
A <- itRetA * TAA_1[nrow(alloData)] + itRetB * TAA_2[nrow(alloData)] #set target return
A_mat <- matrix(rep(A, numPi1), byrow = TRUE, ncol = length(A), nrow =  numPi1)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons 
LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
len_grid[1] <- length(balances)

# set consumption vector
consum_vec <- seq(LB_C, salary, length.out = numC1)

# compute prob 1 matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi1)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB
vec <- rep(1,length(consum_vec))
count1_mat <- (ret_mat - A_mat) > 0
prob1 <- rowSums(count1_mat) / length(itRetA)


result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
  
  itBalStart <- balances[y]
  set.seed(1)
  
  #optimisation for asset allocation
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #optimisation for consumption
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = c(LB_C_rate),
                       upper = c(1),
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons1 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [1]
  B1 <- cons1 - width1
  B1 <- ifelse(B1 < LB_C, LB_C, B1)
  B2 <- cons1 + width1
  B2 <- ifelse(B2 > salary, salary, B2)
  B1_p = B1 / salary
  B2_p = B2 / salary
  consum_vec <- seq(B1, B2, length.out = numC1)
  
  #re-optimisation for asset allocation [1]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [1]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B1_p,##
                       upper = B2_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons2 <- salary * optParam[1]
  
  # re-optimisation parameter [2]
  B3 <- cons2 - width2
  B3 <- ifelse(B3 < LB_C, LB_C, B3)
  B4 <- cons2 + width2
  B4 <- ifelse(B4 > salary, salary, B4)
  B3_p = B3 / salary
  B4_p = B4 / salary
  consum_vec <- seq(B3, B4, length.out = numC1)
  
  #re-optimisation for asset allocation [2]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [2]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B3_p,##
                       upper = B4_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons3 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [3]
  B5 <- cons3 - width3
  B5 <- ifelse(B5 < LB_C, LB_C, B5)
  B6 <- cons3 + width3
  B6 <- ifelse(B6 > salary, salary, B6)
  B5_p = B5 / salary
  B6_p = B6 / salary
  consum_vec <- seq(B5, B6, length.out = numC1)
  
  #re-optimisation for asset allocation [3]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  valfun1 = inv_result$val1
  valfun1_ord = valfun1[order(consum)]  # small to large
  valfun2 = inv_result$val2
  valfun2_ord = valfun2[order(consum)]  # small to large
  
  #re-optimisation for consumption [3]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B5_p,##
                       upper = B6_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  cons <- salary * optParam[1]
  
  eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
  val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
  val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = WT)
  
  val_1 = val_1_1 + val_1_2
  ending_balance = itBalStart - cons + salary
  aftax_cons = cons - income_tax(cons)
  
  resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[nrow(alloData)], 
                 consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility, 
                 pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
  resultsIt
}
rownames(result_xx) <- NULL
assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1]  
############################################################
# pre retirement ages (from age 25 to age 63)

# eqprop vector & return matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB

#age_x = 63
for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  Upper <- WT + (age_retire - age_x - 1) * Wt
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[(age_retire - age_x)] <- length(balances)
  
  consum_vec <- seq(LB_C, salary, length.out = numC) #consumption vector
  
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  
  itPrevBal = resultsPrev$balance
  itPrevBal_ord = itPrevBal[order(itPrevBal)]  # small to large
  itPrevpseudoCons = resultsPrev$pseudoConsum
  itPrevpseudoCons_ord = itPrevpseudoCons[order(itPrevBal)]  # small to large
  itPrevVal = resultsPrev$value_manager
  itPrevVal_ord = itPrevVal[order(itPrevBal)]  # small to large
  
  schu_spline1 <- Schumaker(x = itPrevBal_ord, y = itPrevpseudoCons_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline for pseudoconsum at age x+1
  
  A <- TAA_1[(age_x - 24)] * itRetA + TAA_2[(age_x - 24)] * itRetB  # set target return
  A_mat = matrix(rep(A, numPi), byrow = TRUE, ncol = length(A), nrow =  numPi)
  
  # compute prob1 matrix
  vec <- rep(1,length(consum_vec))
  cal1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(cal1_mat) / length(itRetA)
  
  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    set.seed(1)
    itBalStart <- balances[y]
    
    #optimisation for asset allocation
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_xx,
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]
    
    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC)
    
    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }  
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1_vec
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2_vec
    valfun2_ord = valfun2[order(consum)]  # small to large
    
    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_xx,
                         lower = B1_p,#
                         upper = B2_p,
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]
    
    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = Upper)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)
    
    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[(age_x - 24)],
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility,
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire - 1))##
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t7 <- do.call(rbind, l.df)

save(data_vfi_t7, file = paste0("mgmt_H1_rho_",rho, "salary", salary, "_base2_v5.1.RData"))##
write.csv(data_vfi_t7, paste0("mgmt_H1_rho_",rho, "salary", salary, "_base2_v5.1.csv"), row.names = F)##

###########################################################################################################################################
#================H1 Case base3 (allocation range is [-0.32,1.32])=========================#

## Initialisation  
rm(list=ls())

source("opt_fns.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores(),outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))

#importing target allocation data
alloData <- read.csv("Target_asset_allocation.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 4 
l <- 0.7 # replacement ratio of post-retirement yearly consumption level
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
Wt <- 0.5 # weight assigned to benchmark A
WT <- 1 # weight assigned to target B
C_star <- salary - income_tax(salary) # consumption target (pre-retirement)
B_0 <- l * C_star # post-tax post-retirement consumption target

#starting balance at each age
factor <- 1.04
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
max_start_bal <- 3000000 ##maximum starting balance, as balances above these have same outcomes
start_bal_list$Start_Bal <- ifelse(start_bal_list$Start_Bal > max_start_bal, max_start_bal, start_bal_list$Start_Bal)
bal_spacing <- 5000 ## how fine the balance grid should be (max = 10000)

##Adjustments for Resolution
numC <- 51       # resolution for consumption vector 
numC1 <- 101     # resolution for consumption vector (for age 64)
numPi <- 8201   # resolution for equity allocation % vector (incremental value=(1.32+0.32)/8200=0.0002)
numPi1 <- 16401  # resolution for equity allocation % vector (incremental value=(1.32+0.32)/16400=0.0001) (for age 64)

width1 <- 1000 # first reoptimisation on consumption between +/- width1
width2 <- 100  # second reoptimisation on consumption between +/- width2
width3 <- 50   # third reoptimisation on consumption between +/- width3

#setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1.32
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- -0.32
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation controls
controlDE1 <- list(reltol = 1e-10, itermax = 300)
controlDE2 <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()
###################################
## at terminal age - 1 (age 64)
age_x <- age_retire - 1
cat(paste("age_x", age_x, "in progress \n"))
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
A <- itRetA * TAA_1[nrow(alloData)] + itRetB * TAA_2[nrow(alloData)] #set target return
A_mat <- matrix(rep(A, numPi1), byrow = TRUE, ncol = length(A), nrow =  numPi1)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons 
LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
len_grid[1] <- length(balances)

# set consumption vector
consum_vec <- seq(LB_C, salary, length.out = numC1)

# compute prob 1 matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi1)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB
vec <- rep(1,length(consum_vec))
count1_mat <- (ret_mat - A_mat) > 0
prob1 <- rowSums(count1_mat) / length(itRetA)


result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
  
  itBalStart <- balances[y]
  set.seed(1)
  
  #optimisation for asset allocation
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #optimisation for consumption
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = c(LB_C_rate),
                       upper = c(1),
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons1 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [1]
  B1 <- cons1 - width1
  B1 <- ifelse(B1 < LB_C, LB_C, B1)
  B2 <- cons1 + width1
  B2 <- ifelse(B2 > salary, salary, B2)
  B1_p = B1 / salary
  B2_p = B2 / salary
  consum_vec <- seq(B1, B2, length.out = numC1)
  
  #re-optimisation for asset allocation [1]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [1]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B1_p,##
                       upper = B2_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons2 <- salary * optParam[1]
  
  # re-optimisation parameter [2]
  B3 <- cons2 - width2
  B3 <- ifelse(B3 < LB_C, LB_C, B3)
  B4 <- cons2 + width2
  B4 <- ifelse(B4 > salary, salary, B4)
  B3_p = B3 / salary
  B4_p = B4 / salary
  consum_vec <- seq(B3, B4, length.out = numC1)
  
  #re-optimisation for asset allocation [2]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [2]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B3_p,##
                       upper = B4_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons3 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [3]
  B5 <- cons3 - width3
  B5 <- ifelse(B5 < LB_C, LB_C, B5)
  B6 <- cons3 + width3
  B6 <- ifelse(B6 > salary, salary, B6)
  B5_p = B5 / salary
  B6_p = B6 / salary
  consum_vec <- seq(B5, B6, length.out = numC1)
  
  #re-optimisation for asset allocation [3]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  valfun1 = inv_result$val1
  valfun1_ord = valfun1[order(consum)]  # small to large
  valfun2 = inv_result$val2
  valfun2_ord = valfun2[order(consum)]  # small to large
  
  #re-optimisation for consumption [3]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B5_p,##
                       upper = B6_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  cons <- salary * optParam[1]
  
  eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
  val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
  val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = WT)
  
  val_1 = val_1_1 + val_1_2
  ending_balance = itBalStart - cons + salary
  aftax_cons = cons - income_tax(cons)
  
  resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[nrow(alloData)], 
                 consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility, 
                 pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
  resultsIt
}
rownames(result_xx) <- NULL
assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1]  
############################################################
# pre retirement ages (from age 25 to age 63)

# eqprop vector & return matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB

#age_x = 63
for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  Upper <- WT + (age_retire - age_x - 1) * Wt
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[(age_retire - age_x)] <- length(balances)
  
  consum_vec <- seq(LB_C, salary, length.out = numC) #consumption vector
  
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  
  itPrevBal = resultsPrev$balance
  itPrevBal_ord = itPrevBal[order(itPrevBal)]  # small to large
  itPrevpseudoCons = resultsPrev$pseudoConsum
  itPrevpseudoCons_ord = itPrevpseudoCons[order(itPrevBal)]  # small to large
  itPrevVal = resultsPrev$value_manager
  itPrevVal_ord = itPrevVal[order(itPrevBal)]  # small to large
  
  schu_spline1 <- Schumaker(x = itPrevBal_ord, y = itPrevpseudoCons_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline for pseudoconsum at age x+1
  
  A <- TAA_1[(age_x - 24)] * itRetA + TAA_2[(age_x - 24)] * itRetB  # set target return
  A_mat = matrix(rep(A, numPi), byrow = TRUE, ncol = length(A), nrow =  numPi)
  
  # compute prob1 matrix
  vec <- rep(1,length(consum_vec))
  cal1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(cal1_mat) / length(itRetA)
  
  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    set.seed(1)
    itBalStart <- balances[y]
    
    #optimisation for asset allocation
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_xx,
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]
    
    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC)
    
    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }  
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1_vec
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2_vec
    valfun2_ord = valfun2[order(consum)]  # small to large
    
    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_xx,
                         lower = B1_p,#
                         upper = B2_p,
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]
    
    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = Upper)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)
    
    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[(age_x - 24)],
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility,
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire - 1))##
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t7 <- do.call(rbind, l.df)

save(data_vfi_t7, file = paste0("mgmt_H1_rho_",rho, "salary", salary, "_base3_v5.1.RData"))##
write.csv(data_vfi_t7, paste0("mgmt_H1_rho_",rho, "salary", salary, "_base3_v5.1.csv"), row.names = F)##

###########################################################################################################################################
#================H1 Case sens1_1 (salary = $55000)=========================#

## Initialisation  
rm(list=ls())

source("opt_fns.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores(),outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))

#importing target allocation data
alloData <- read.csv("Target_asset_allocation.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 55000
rho <- 4 
l <- 0.7 # replacement ratio of post-retirement yearly consumption level
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
Wt <- 0.5 # weight assigned to benchmark A
WT <- 1 # weight assigned to target B
C_star <- salary - income_tax(salary) # consumption target (pre-retirement)
B_0 <- l * C_star # post-tax post-retirment consumption target

#starting balance at each age
factor <- 1.04
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
max_start_bal <- 3000000 ##maximum starting balance, as balances above these have same outcomes
start_bal_list$Start_Bal <- ifelse(start_bal_list$Start_Bal > max_start_bal, max_start_bal, start_bal_list$Start_Bal)
bal_spacing <- 5000 ## how fine the balance grid should be (max = 10000)

##Adjustments for Resolution
numC <- 51       # resolution for consumption vector 
numC1 <- 101     # resolution for consumption vector (for age 64)
numPi <- 6601   # resolution for equity allocation % vector (incremental value=(1.16+0.16)/6600=0.0002)
numPi1 <- 13201  # resolution for equity allocation % vector (incremental value=(1.16+0.16)/13200=0.0001) (for age 64)

width1 <- 1000 # first reoptimisation on consumption between +/- width1
width2 <- 100  # second reoptimisation on consumption between +/- width2
width3 <- 50   # third reoptimisation on consumption between +/- width3

#setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1.16 
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- -0.16 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation controls
controlDE1 <- list(reltol = 1e-10, itermax = 300)
controlDE2 <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()
###################################
## at terminal age - 1 (age 64)
age_x <- age_retire - 1
cat(paste("age_x", age_x, "in progress \n"))
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
A <- itRetA * TAA_1[nrow(alloData)] + itRetB * TAA_2[nrow(alloData)] #set target return
A_mat <- matrix(rep(A, numPi1), byrow = TRUE, ncol = length(A), nrow =  numPi1)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons 
LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
len_grid[1] <- length(balances)

# set consumption vector
consum_vec <- seq(LB_C, salary, length.out = numC1)

# compute prob 1 matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi1)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB
vec <- rep(1,length(consum_vec))
count1_mat <- (ret_mat - A_mat) > 0
prob1 <- rowSums(count1_mat) / length(itRetA)


result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
  
  itBalStart <- balances[y]
  set.seed(1)
  
  #optimisation for asset allocation
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #optimisation for consumption
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = c(LB_C_rate),
                       upper = c(1),
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons1 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [1]
  B1 <- cons1 - width1
  B1 <- ifelse(B1 < LB_C, LB_C, B1)
  B2 <- cons1 + width1
  B2 <- ifelse(B2 > salary, salary, B2)
  B1_p = B1 / salary
  B2_p = B2 / salary
  consum_vec <- seq(B1, B2, length.out = numC1)
  
  #re-optimisation for asset allocation [1]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [1]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B1_p,##
                       upper = B2_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons2 <- salary * optParam[1]
  
  # re-optimisation parameter [2]
  B3 <- cons2 - width2
  B3 <- ifelse(B3 < LB_C, LB_C, B3)
  B4 <- cons2 + width2
  B4 <- ifelse(B4 > salary, salary, B4)
  B3_p = B3 / salary
  B4_p = B4 / salary
  consum_vec <- seq(B3, B4, length.out = numC1)
  
  #re-optimisation for asset allocation [2]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [2]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B3_p,##
                       upper = B4_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons3 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [3]
  B5 <- cons3 - width3
  B5 <- ifelse(B5 < LB_C, LB_C, B5)
  B6 <- cons3 + width3
  B6 <- ifelse(B6 > salary, salary, B6)
  B5_p = B5 / salary
  B6_p = B6 / salary
  consum_vec <- seq(B5, B6, length.out = numC1)
  
  #re-optimisation for asset allocation [3]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  valfun1 = inv_result$val1
  valfun1_ord = valfun1[order(consum)]  # small to large
  valfun2 = inv_result$val2
  valfun2_ord = valfun2[order(consum)]  # small to large
  
  #re-optimisation for consumption [3]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B5_p,##
                       upper = B6_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  cons <- salary * optParam[1]
  
  eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
  val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
  val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = WT)
  
  val_1 = val_1_1 + val_1_2
  ending_balance = itBalStart - cons + salary
  aftax_cons = cons - income_tax(cons)
  
  resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[nrow(alloData)], 
                 consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility, 
                 pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
  resultsIt
}
rownames(result_xx) <- NULL
assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1]  
############################################################
# pre retirement ages (from age 25 to age 63)

# eqprop vector & return matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB

#age_x = 63
for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  Upper <- WT + (age_retire - age_x - 1) * Wt
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[(age_retire - age_x)] <- length(balances)
  
  consum_vec <- seq(LB_C, salary, length.out = numC) #consumption vector
  
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  
  itPrevBal = resultsPrev$balance
  itPrevBal_ord = itPrevBal[order(itPrevBal)]  # small to large
  itPrevpseudoCons = resultsPrev$pseudoConsum
  itPrevpseudoCons_ord = itPrevpseudoCons[order(itPrevBal)]  # small to large
  itPrevVal = resultsPrev$value_manager
  itPrevVal_ord = itPrevVal[order(itPrevBal)]  # small to large
  
  schu_spline1 <- Schumaker(x = itPrevBal_ord, y = itPrevpseudoCons_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline for pseudoconsum at age x+1
  
  A <- TAA_1[(age_x - 24)] * itRetA + TAA_2[(age_x - 24)] * itRetB  # set target return
  A_mat = matrix(rep(A, numPi), byrow = TRUE, ncol = length(A), nrow =  numPi)
  
  # compute prob1 matrix
  vec <- rep(1,length(consum_vec))
  cal1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(cal1_mat) / length(itRetA)
  
  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    set.seed(1)
    itBalStart <- balances[y]
    
    #optimisation for asset allocation
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_xx,
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]
    
    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC)
    
    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }  
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1_vec
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2_vec
    valfun2_ord = valfun2[order(consum)]  # small to large
    
    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_xx,
                         lower = B1_p,#
                         upper = B2_p,
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]
    
    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = Upper)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)
    
    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[(age_x - 24)],
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility,
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire - 1))##
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t7 <- do.call(rbind, l.df)

save(data_vfi_t7, file = paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens1_1_v5.1.RData"))##
write.csv(data_vfi_t7, paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens1_1_v5.1.csv"), row.names = F)##

###########################################################################################################################################

#================H1 Case sens1_2 (salary = $85000)=========================#

## Initialisation  
rm(list=ls())

source("opt_fns.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores(),outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))

#importing target allocation data
alloData <- read.csv("Target_asset_allocation.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 85000
rho <- 4 
l <- 0.7 # replacement ratio of post-retirement yearly consumption level
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
Wt <- 0.5 # weight assigned to benchmark A
WT <- 1 # weight assigned to target B
C_star <- salary - income_tax(salary) # consumption target (pre-retirement)
B_0 <- l * C_star # post-tax post-retirment consumption target

#starting balance at each age
factor <- 1.04
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
max_start_bal <- 3000000 ##maximum starting balance, as balances above these have same outcomes
start_bal_list$Start_Bal <- ifelse(start_bal_list$Start_Bal > max_start_bal, max_start_bal, start_bal_list$Start_Bal)
bal_spacing <- 5000 ## how fine the balance grid should be (max = 10000)

##Adjustments for Resolution
numC <- 51       # resolution for consumption vector 
numC1 <- 101     # resolution for consumption vector (for age 64)
numPi <- 6601   # resolution for equity allocation % vector (incremental value=(1.16+0.16)/6600=0.0002)
numPi1 <- 13201  # resolution for equity allocation % vector (incremental value=(1.16+0.16)/13200=0.0001) (for age 64)

width1 <- 1000 # first reoptimisation on consumption between +/- width1
width2 <- 100  # second reoptimisation on consumption between +/- width2
width3 <- 50   # third reoptimisation on consumption between +/- width3

#setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1.16 
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- -0.16 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation controls
controlDE1 <- list(reltol = 1e-10, itermax = 300)
controlDE2 <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()
###################################
## at terminal age - 1 (age 64)
age_x <- age_retire - 1
cat(paste("age_x", age_x, "in progress \n"))
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
A <- itRetA * TAA_1[nrow(alloData)] + itRetB * TAA_2[nrow(alloData)] #set target return
A_mat <- matrix(rep(A, numPi1), byrow = TRUE, ncol = length(A), nrow =  numPi1)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons 
LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
len_grid[1] <- length(balances)

# set consumption vector
consum_vec <- seq(LB_C, salary, length.out = numC1)

# compute prob 1 matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi1)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB
vec <- rep(1,length(consum_vec))
count1_mat <- (ret_mat - A_mat) > 0
prob1 <- rowSums(count1_mat) / length(itRetA)


result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
  
  itBalStart <- balances[y]
  set.seed(1)
  
  #optimisation for asset allocation
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #optimisation for consumption
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = c(LB_C_rate),
                       upper = c(1),
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons1 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [1]
  B1 <- cons1 - width1
  B1 <- ifelse(B1 < LB_C, LB_C, B1)
  B2 <- cons1 + width1
  B2 <- ifelse(B2 > salary, salary, B2)
  B1_p = B1 / salary
  B2_p = B2 / salary
  consum_vec <- seq(B1, B2, length.out = numC1)
  
  #re-optimisation for asset allocation [1]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [1]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B1_p,##
                       upper = B2_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons2 <- salary * optParam[1]
  
  # re-optimisation parameter [2]
  B3 <- cons2 - width2
  B3 <- ifelse(B3 < LB_C, LB_C, B3)
  B4 <- cons2 + width2
  B4 <- ifelse(B4 > salary, salary, B4)
  B3_p = B3 / salary
  B4_p = B4 / salary
  consum_vec <- seq(B3, B4, length.out = numC1)
  
  #re-optimisation for asset allocation [2]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [2]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B3_p,##
                       upper = B4_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons3 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [3]
  B5 <- cons3 - width3
  B5 <- ifelse(B5 < LB_C, LB_C, B5)
  B6 <- cons3 + width3
  B6 <- ifelse(B6 > salary, salary, B6)
  B5_p = B5 / salary
  B6_p = B6 / salary
  consum_vec <- seq(B5, B6, length.out = numC1)
  
  #re-optimisation for asset allocation [3]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  valfun1 = inv_result$val1
  valfun1_ord = valfun1[order(consum)]  # small to large
  valfun2 = inv_result$val2
  valfun2_ord = valfun2[order(consum)]  # small to large
  
  #re-optimisation for consumption [3]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B5_p,##
                       upper = B6_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  cons <- salary * optParam[1]
  
  eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
  val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
  val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = WT)
  
  val_1 = val_1_1 + val_1_2
  ending_balance = itBalStart - cons + salary
  aftax_cons = cons - income_tax(cons)
  
  resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[nrow(alloData)], 
                 consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility, 
                 pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
  resultsIt
}
rownames(result_xx) <- NULL
assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1]  
############################################################
# pre retirement ages (from age 25 to age 63)

# eqprop vector & return matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB

#age_x = 63
for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  Upper <- WT + (age_retire - age_x - 1) * Wt
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[(age_retire - age_x)] <- length(balances)
  
  consum_vec <- seq(LB_C, salary, length.out = numC) #consumption vector
  
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  
  itPrevBal = resultsPrev$balance
  itPrevBal_ord = itPrevBal[order(itPrevBal)]  # small to large
  itPrevpseudoCons = resultsPrev$pseudoConsum
  itPrevpseudoCons_ord = itPrevpseudoCons[order(itPrevBal)]  # small to large
  itPrevVal = resultsPrev$value_manager
  itPrevVal_ord = itPrevVal[order(itPrevBal)]  # small to large
  
  schu_spline1 <- Schumaker(x = itPrevBal_ord, y = itPrevpseudoCons_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline for pseudoconsum at age x+1
  
  A <- TAA_1[(age_x - 24)] * itRetA + TAA_2[(age_x - 24)] * itRetB  # set target return
  A_mat = matrix(rep(A, numPi), byrow = TRUE, ncol = length(A), nrow =  numPi)
  
  # compute prob1 matrix
  vec <- rep(1,length(consum_vec))
  cal1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(cal1_mat) / length(itRetA)
  
  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    set.seed(1)
    itBalStart <- balances[y]
    
    #optimisation for asset allocation
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_xx,
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]
    
    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC)
    
    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }  
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1_vec
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2_vec
    valfun2_ord = valfun2[order(consum)]  # small to large
    
    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_xx,
                         lower = B1_p,#
                         upper = B2_p,
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]
    
    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = Upper)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)
    
    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[(age_x - 24)],
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility,
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire - 1))##
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t7 <- do.call(rbind, l.df)

save(data_vfi_t7, file = paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens1_2_v5.1.RData"))##
write.csv(data_vfi_t7, paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens1_2_v5.1.csv"), row.names = F)##

###########################################################################################################################################

#================H1 Case sens2_1 (rho = 3)=========================#

## Initialisation  
rm(list=ls())

source("opt_fns.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores(),outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))

#importing target allocation data
alloData <- read.csv("Target_asset_allocation.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 3 
l <- 0.7 # replacement ratio of post-retirement yearly consumption level
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
Wt <- 0.5 # weight assigned to benchmark A
WT <- 1 # weight assigned to target B
C_star <- salary - income_tax(salary) # consumption target (pre-retirement)
B_0 <- l * C_star # post-tax post-retirment consumption target

#starting balance at each age
factor <- 1.04
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
max_start_bal <- 3000000 ##maximum starting balance, as balances above these have same outcomes
start_bal_list$Start_Bal <- ifelse(start_bal_list$Start_Bal > max_start_bal, max_start_bal, start_bal_list$Start_Bal)
bal_spacing <- 5000 ## how fine the balance grid should be (max = 10000)

##Adjustments for Resolution
numC <- 51       # resolution for consumption vector 
numC1 <- 101     # resolution for consumption vector (for age 64)
numPi <- 6601   # resolution for equity allocation % vector (incremental value=(1.16+0.16)/6600=0.0002)
numPi1 <- 13201  # resolution for equity allocation % vector (incremental value=(1.16+0.16)/13200=0.0001) (for age 64)

width1 <- 1000 # first reoptimisation on consumption between +/- width1
width2 <- 100  # second reoptimisation on consumption between +/- width2
width3 <- 50   # third reoptimisation on consumption between +/- width3

#setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1.16 
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- -0.16 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation controls
controlDE1 <- list(reltol = 1e-10, itermax = 300)
controlDE2 <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()
###################################
## at terminal age - 1 (age 64)
age_x <- age_retire - 1
cat(paste("age_x", age_x, "in progress \n"))
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
A <- itRetA * TAA_1[nrow(alloData)] + itRetB * TAA_2[nrow(alloData)] #set target return
A_mat <- matrix(rep(A, numPi1), byrow = TRUE, ncol = length(A), nrow =  numPi1)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons 
LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
len_grid[1] <- length(balances)

# set consumption vector
consum_vec <- seq(LB_C, salary, length.out = numC1)

# compute prob 1 matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi1)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB
vec <- rep(1,length(consum_vec))
count1_mat <- (ret_mat - A_mat) > 0
prob1 <- rowSums(count1_mat) / length(itRetA)


result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
  
  itBalStart <- balances[y]
  set.seed(1)
  
  #optimisation for asset allocation
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #optimisation for consumption
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = c(LB_C_rate),
                       upper = c(1),
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons1 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [1]
  B1 <- cons1 - width1
  B1 <- ifelse(B1 < LB_C, LB_C, B1)
  B2 <- cons1 + width1
  B2 <- ifelse(B2 > salary, salary, B2)
  B1_p = B1 / salary
  B2_p = B2 / salary
  consum_vec <- seq(B1, B2, length.out = numC1)
  
  #re-optimisation for asset allocation [1]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [1]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B1_p,##
                       upper = B2_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons2 <- salary * optParam[1]
  
  # re-optimisation parameter [2]
  B3 <- cons2 - width2
  B3 <- ifelse(B3 < LB_C, LB_C, B3)
  B4 <- cons2 + width2
  B4 <- ifelse(B4 > salary, salary, B4)
  B3_p = B3 / salary
  B4_p = B4 / salary
  consum_vec <- seq(B3, B4, length.out = numC1)
  
  #re-optimisation for asset allocation [2]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [2]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B3_p,##
                       upper = B4_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons3 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [3]
  B5 <- cons3 - width3
  B5 <- ifelse(B5 < LB_C, LB_C, B5)
  B6 <- cons3 + width3
  B6 <- ifelse(B6 > salary, salary, B6)
  B5_p = B5 / salary
  B6_p = B6 / salary
  consum_vec <- seq(B5, B6, length.out = numC1)
  
  #re-optimisation for asset allocation [3]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  valfun1 = inv_result$val1
  valfun1_ord = valfun1[order(consum)]  # small to large
  valfun2 = inv_result$val2
  valfun2_ord = valfun2[order(consum)]  # small to large
  
  #re-optimisation for consumption [3]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B5_p,##
                       upper = B6_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  cons <- salary * optParam[1]
  
  eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
  val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
  val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = WT)
  
  val_1 = val_1_1 + val_1_2
  ending_balance = itBalStart - cons + salary
  aftax_cons = cons - income_tax(cons)
  
  resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[nrow(alloData)], 
                 consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility, 
                 pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
  resultsIt
}
rownames(result_xx) <- NULL
assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1]  
############################################################
# pre retirement ages (from age 25 to age 63)

# eqprop vector & return matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB

#age_x = 63
for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  Upper <- WT + (age_retire - age_x - 1) * Wt
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[(age_retire - age_x)] <- length(balances)
  
  consum_vec <- seq(LB_C, salary, length.out = numC) #consumption vector
  
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  
  itPrevBal = resultsPrev$balance
  itPrevBal_ord = itPrevBal[order(itPrevBal)]  # small to large
  itPrevpseudoCons = resultsPrev$pseudoConsum
  itPrevpseudoCons_ord = itPrevpseudoCons[order(itPrevBal)]  # small to large
  itPrevVal = resultsPrev$value_manager
  itPrevVal_ord = itPrevVal[order(itPrevBal)]  # small to large
  
  schu_spline1 <- Schumaker(x = itPrevBal_ord, y = itPrevpseudoCons_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline for pseudoconsum at age x+1
  
  A <- TAA_1[(age_x - 24)] * itRetA + TAA_2[(age_x - 24)] * itRetB  # set target return
  A_mat = matrix(rep(A, numPi), byrow = TRUE, ncol = length(A), nrow =  numPi)
  
  # compute prob1 matrix
  vec <- rep(1,length(consum_vec))
  cal1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(cal1_mat) / length(itRetA)
  
  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    set.seed(1)
    itBalStart <- balances[y]
    
    #optimisation for asset allocation
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_xx,
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]
    
    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC)
    
    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }  
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1_vec
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2_vec
    valfun2_ord = valfun2[order(consum)]  # small to large
    
    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_xx,
                         lower = B1_p,#
                         upper = B2_p,
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]
    
    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = Upper)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)
    
    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[(age_x - 24)],
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility,
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire - 1))##
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t7 <- do.call(rbind, l.df)

save(data_vfi_t7, file = paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens2_1_v5.1.RData"))##
write.csv(data_vfi_t7, paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens2_1_v5.1.csv"), row.names = F)##

###########################################################################################################################################

#================H1 Case sens2_2 (rho = 5)=========================#

## Initialisation  
rm(list=ls())

source("opt_fns.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores(),outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))

#importing target allocation data
alloData <- read.csv("Target_asset_allocation.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 5 
l <- 0.7 # replacement ratio of post-retirement yearly consumption level
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
Wt <- 0.5 # weight assigned to benchmark A
WT <- 1 # weight assigned to target B
C_star <- salary - income_tax(salary) # consumption target (pre-retirement)
B_0 <- l * C_star # post-tax post-retirment consumption target

#starting balance at each age
factor <- 1.04
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
max_start_bal <- 3000000 ##maximum starting balance, as balances above these have same outcomes
start_bal_list$Start_Bal <- ifelse(start_bal_list$Start_Bal > max_start_bal, max_start_bal, start_bal_list$Start_Bal)
bal_spacing <- 5000 ## how fine the balance grid should be (max = 10000)

##Adjustments for Resolution
numC <- 51       # resolution for consumption vector 
numC1 <- 101     # resolution for consumption vector (for age 64)
numPi <- 6601   # resolution for equity allocation % vector (incremental value=(1.16+0.16)/6600=0.0002)
numPi1 <- 13201  # resolution for equity allocation % vector (incremental value=(1.16+0.16)/13200=0.0001) (for age 64)

width1 <- 1000 # first reoptimisation on consumption between +/- width1
width2 <- 100  # second reoptimisation on consumption between +/- width2
width3 <- 50   # third reoptimisation on consumption between +/- width3

#setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1.16 
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- -0.16 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation controls
controlDE1 <- list(reltol = 1e-10, itermax = 300)
controlDE2 <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()
###################################
## at terminal age - 1 (age 64)
age_x <- age_retire - 1
cat(paste("age_x", age_x, "in progress \n"))
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
A <- itRetA * TAA_1[nrow(alloData)] + itRetB * TAA_2[nrow(alloData)] #set target return
A_mat <- matrix(rep(A, numPi1), byrow = TRUE, ncol = length(A), nrow =  numPi1)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons 
LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
len_grid[1] <- length(balances)

# set consumption vector
consum_vec <- seq(LB_C, salary, length.out = numC1)

# compute prob 1 matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi1)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB
vec <- rep(1,length(consum_vec))
count1_mat <- (ret_mat - A_mat) > 0
prob1 <- rowSums(count1_mat) / length(itRetA)


result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
  
  itBalStart <- balances[y]
  set.seed(1)
  
  #optimisation for asset allocation
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #optimisation for consumption
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = c(LB_C_rate),
                       upper = c(1),
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons1 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [1]
  B1 <- cons1 - width1
  B1 <- ifelse(B1 < LB_C, LB_C, B1)
  B2 <- cons1 + width1
  B2 <- ifelse(B2 > salary, salary, B2)
  B1_p = B1 / salary
  B2_p = B2 / salary
  consum_vec <- seq(B1, B2, length.out = numC1)
  
  #re-optimisation for asset allocation [1]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [1]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B1_p,##
                       upper = B2_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons2 <- salary * optParam[1]
  
  # re-optimisation parameter [2]
  B3 <- cons2 - width2
  B3 <- ifelse(B3 < LB_C, LB_C, B3)
  B4 <- cons2 + width2
  B4 <- ifelse(B4 > salary, salary, B4)
  B3_p = B3 / salary
  B4_p = B4 / salary
  consum_vec <- seq(B3, B4, length.out = numC1)
  
  #re-optimisation for asset allocation [2]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [2]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B3_p,##
                       upper = B4_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons3 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [3]
  B5 <- cons3 - width3
  B5 <- ifelse(B5 < LB_C, LB_C, B5)
  B6 <- cons3 + width3
  B6 <- ifelse(B6 > salary, salary, B6)
  B5_p = B5 / salary
  B6_p = B6 / salary
  consum_vec <- seq(B5, B6, length.out = numC1)
  
  #re-optimisation for asset allocation [3]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  valfun1 = inv_result$val1
  valfun1_ord = valfun1[order(consum)]  # small to large
  valfun2 = inv_result$val2
  valfun2_ord = valfun2[order(consum)]  # small to large
  
  #re-optimisation for consumption [3]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B5_p,##
                       upper = B6_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  cons <- salary * optParam[1]
  
  eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
  val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
  val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = WT)
  
  val_1 = val_1_1 + val_1_2
  ending_balance = itBalStart - cons + salary
  aftax_cons = cons - income_tax(cons)
  
  resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[nrow(alloData)], 
                 consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility, 
                 pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
  resultsIt
}
rownames(result_xx) <- NULL
assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1]  
############################################################
# pre retirement ages (from age 25 to age 63)

# eqprop vector & return matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB

#age_x = 63
for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  Upper <- WT + (age_retire - age_x - 1) * Wt
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[(age_retire - age_x)] <- length(balances)
  
  consum_vec <- seq(LB_C, salary, length.out = numC) #consumption vector
  
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  
  itPrevBal = resultsPrev$balance
  itPrevBal_ord = itPrevBal[order(itPrevBal)]  # small to large
  itPrevpseudoCons = resultsPrev$pseudoConsum
  itPrevpseudoCons_ord = itPrevpseudoCons[order(itPrevBal)]  # small to large
  itPrevVal = resultsPrev$value_manager
  itPrevVal_ord = itPrevVal[order(itPrevBal)]  # small to large
  
  schu_spline1 <- Schumaker(x = itPrevBal_ord, y = itPrevpseudoCons_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline for pseudoconsum at age x+1
  
  A <- TAA_1[(age_x - 24)] * itRetA + TAA_2[(age_x - 24)] * itRetB  # set target return
  A_mat = matrix(rep(A, numPi), byrow = TRUE, ncol = length(A), nrow =  numPi)
  
  # compute prob1 matrix
  vec <- rep(1,length(consum_vec))
  cal1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(cal1_mat) / length(itRetA)
  
  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    set.seed(1)
    itBalStart <- balances[y]
    
    #optimisation for asset allocation
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_xx,
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]
    
    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC)
    
    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }  
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1_vec
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2_vec
    valfun2_ord = valfun2[order(consum)]  # small to large
    
    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_xx,
                         lower = B1_p,#
                         upper = B2_p,
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]
    
    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = Upper)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)
    
    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[(age_x - 24)],
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility,
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire - 1))##
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t7 <- do.call(rbind, l.df)

save(data_vfi_t7, file = paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens2_2_v5.1.RData"))##
write.csv(data_vfi_t7, paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens2_2_v5.1.csv"), row.names = F)##


###########################################################################################################################################
#Note: Since sense 3_1 and sense 3_2 were found to be inappropriate for inclusion in our analysis, the relevant code has been omitted.

#================H1 Case sens4_1 (-1% ERP)=========================#

## Initialisation  
rm(list=ls())

source("opt_fns.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores()-1,outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset_sens1.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))

#importing target allocation data
alloData <- read.csv("Target_asset_allocation.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 4 
l <- 0.7 # replacement ratio of post-retirement yearly consumption level
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
Wt <- 0.5 # weight assigned to benchmark A
WT <- 1 # weight assigned to target B
C_star <- salary - income_tax(salary) # consumption target (pre-retirement)
B_0 <- l * C_star # post-tax post-retirment consumption target

#starting balance at each age
factor <- 1.04
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
max_start_bal <- 3000000 ##maximum starting balance, as balances above these have same outcomes
start_bal_list$Start_Bal <- ifelse(start_bal_list$Start_Bal > max_start_bal, max_start_bal, start_bal_list$Start_Bal)
bal_spacing <- 5000 ## how fine the balance grid should be (max = 10000)

##Adjustments for Resolution
numC <- 51       # resolution for consumption vector 
numC1 <- 101     # resolution for consumption vector (for age 64)
numPi <- 6601   # resolution for equity allocation % vector (incremental value=(1.16+0.16)/6600=0.0002)
numPi1 <- 13201  # resolution for equity allocation % vector (incremental value=(1.16+0.16)/13200=0.0001) (for age 64)

width1 <- 1000 # first reoptimisation on consumption between +/- width1
width2 <- 100  # second reoptimisation on consumption between +/- width2
width3 <- 50   # third reoptimisation on consumption between +/- width3

#setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1.16 
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- -0.16 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation controls
controlDE1 <- list(reltol = 1e-10, itermax = 300)
controlDE2 <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()
###################################
## at terminal age - 1 (age 64)
age_x <- age_retire - 1
cat(paste("age_x", age_x, "in progress \n"))
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
A <- itRetA * TAA_1[nrow(alloData)] + itRetB * TAA_2[nrow(alloData)] #set target return
A_mat <- matrix(rep(A, numPi1), byrow = TRUE, ncol = length(A), nrow =  numPi1)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons 
LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
len_grid[1] <- length(balances)

# set consumption vector
consum_vec <- seq(LB_C, salary, length.out = numC1)

# compute prob 1 matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi1)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB
vec <- rep(1,length(consum_vec))
count1_mat <- (ret_mat - A_mat) > 0
prob1 <- rowSums(count1_mat) / length(itRetA)


result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {

  itBalStart <- balances[y]
  set.seed(1)
  
  #optimisation for asset allocation
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #optimisation for consumption
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = c(LB_C_rate),
                       upper = c(1),
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons1 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [1]
  B1 <- cons1 - width1
  B1 <- ifelse(B1 < LB_C, LB_C, B1)
  B2 <- cons1 + width1
  B2 <- ifelse(B2 > salary, salary, B2)
  B1_p = B1 / salary
  B2_p = B2 / salary
  consum_vec <- seq(B1, B2, length.out = numC1)
  
  #re-optimisation for asset allocation [1]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [1]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B1_p,##
                       upper = B2_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons2 <- salary * optParam[1]
  
  # re-optimisation parameter [2]
  B3 <- cons2 - width2
  B3 <- ifelse(B3 < LB_C, LB_C, B3)
  B4 <- cons2 + width2
  B4 <- ifelse(B4 > salary, salary, B4)
  B3_p = B3 / salary
  B4_p = B4 / salary
  consum_vec <- seq(B3, B4, length.out = numC1)
  
  #re-optimisation for asset allocation [2]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [2]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B3_p,##
                       upper = B4_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons3 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [3]
  B5 <- cons3 - width3
  B5 <- ifelse(B5 < LB_C, LB_C, B5)
  B6 <- cons3 + width3
  B6 <- ifelse(B6 > salary, salary, B6)
  B5_p = B5 / salary
  B6_p = B6 / salary
  consum_vec <- seq(B5, B6, length.out = numC1)
  
  #re-optimisation for asset allocation [3]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  valfun1 = inv_result$val1
  valfun1_ord = valfun1[order(consum)]  # small to large
  valfun2 = inv_result$val2
  valfun2_ord = valfun2[order(consum)]  # small to large
  
  #re-optimisation for consumption [3]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B5_p,##
                       upper = B6_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  cons <- salary * optParam[1]
  
  eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
  val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
  val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = WT)
  
  val_1 = val_1_1 + val_1_2
  ending_balance = itBalStart - cons + salary
  aftax_cons = cons - income_tax(cons)
  
  resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[nrow(alloData)], 
                 consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility, 
                 pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
  resultsIt
}
rownames(result_xx) <- NULL
assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1]  
############################################################
# pre retirement ages (from age 25 to age 63)

# eqprop vector & return matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB

#age_x = 63
for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  Upper <- WT + (age_retire - age_x - 1) * Wt
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[(age_retire - age_x)] <- length(balances)
  
  consum_vec <- seq(LB_C, salary, length.out = numC) #consumption vector
  
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  
  itPrevBal = resultsPrev$balance
  itPrevBal_ord = itPrevBal[order(itPrevBal)]  # small to large
  itPrevpseudoCons = resultsPrev$pseudoConsum
  itPrevpseudoCons_ord = itPrevpseudoCons[order(itPrevBal)]  # small to large
  itPrevVal = resultsPrev$value_manager
  itPrevVal_ord = itPrevVal[order(itPrevBal)]  # small to large
  
  schu_spline1 <- Schumaker(x = itPrevBal_ord, y = itPrevpseudoCons_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline for pseudoconsum at age x+1
  
  A <- TAA_1[(age_x - 24)] * itRetA + TAA_2[(age_x - 24)] * itRetB  # set target return
  A_mat = matrix(rep(A, numPi), byrow = TRUE, ncol = length(A), nrow =  numPi)
  
  # compute prob1 matrix
  vec <- rep(1,length(consum_vec))
  cal1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(cal1_mat) / length(itRetA)
  
  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    set.seed(1)
    itBalStart <- balances[y]
    
    #optimisation for asset allocation
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_xx,
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]
    
    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC)
    
    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }  
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1_vec
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2_vec
    valfun2_ord = valfun2[order(consum)]  # small to large
    
    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_xx,
                         lower = B1_p,#
                         upper = B2_p,
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]
    
    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = Upper)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)
    
    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[(age_x - 24)],
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility,
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire - 1))##
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t7 <- do.call(rbind, l.df)

save(data_vfi_t7, file = paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens4_1_v5.1.RData"))##
write.csv(data_vfi_t7, paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens4_1_v5.1.csv"), row.names = F)##

###########################################################################################################################################

#================H1 Case sens4_2 (+1% ERP)=========================#

## Initialisation  
rm(list=ls())

source("opt_fns.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores(),outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset_sens2.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))

#importing target allocation data
alloData <- read.csv("Target_asset_allocation.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 4 
l <- 0.7 # replacement ratio of post-retirement yearly consumption level
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
Wt <- 0.5 # weight assigned to benchmark A
WT <- 1 # weight assigned to target B
C_star <- salary - income_tax(salary) # consumption target (pre-retirement)
B_0 <- l * C_star # post-tax post-retirment consumption target

#starting balance at each age
factor <- 1.04
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
max_start_bal <- 3000000 ##maximum starting balance, as balances above these have same outcomes
start_bal_list$Start_Bal <- ifelse(start_bal_list$Start_Bal > max_start_bal, max_start_bal, start_bal_list$Start_Bal)
bal_spacing <- 5000 ## how fine the balance grid should be (max = 10000)

##Adjustments for Resolution
numC <- 51       # resolution for consumption vector 
numC1 <- 101     # resolution for consumption vector (for age 64)
numPi <- 6601   # resolution for equity allocation % vector (incremental value=(1.16+0.16)/6600=0.0002)
numPi1 <- 13201  # resolution for equity allocation % vector (incremental value=(1.16+0.16)/13200=0.0001) (for age 64)

width1 <- 1000 # first reoptimisation on consumption between +/- width1
width2 <- 100  # second reoptimisation on consumption between +/- width2
width3 <- 50   # third reoptimisation on consumption between +/- width3

#setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1.16 
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- -0.16 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation controls
controlDE1 <- list(reltol = 1e-10, itermax = 300)
controlDE2 <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()
###################################
## at terminal age - 1 (age 64)
age_x <- age_retire - 1
cat(paste("age_x", age_x, "in progress \n"))
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
A <- itRetA * TAA_1[nrow(alloData)] + itRetB * TAA_2[nrow(alloData)] #set target return
A_mat <- matrix(rep(A, numPi1), byrow = TRUE, ncol = length(A), nrow =  numPi1)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons 
LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
len_grid[1] <- length(balances)

# set consumption vector
consum_vec <- seq(LB_C, salary, length.out = numC1)

# compute prob 1 matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi1)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB
vec <- rep(1,length(consum_vec))
count1_mat <- (ret_mat - A_mat) > 0
prob1 <- rowSums(count1_mat) / length(itRetA)


result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
  
  itBalStart <- balances[y]
  set.seed(1)
  
  #optimisation for asset allocation
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #optimisation for consumption
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = c(LB_C_rate),
                       upper = c(1),
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons1 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [1]
  B1 <- cons1 - width1
  B1 <- ifelse(B1 < LB_C, LB_C, B1)
  B2 <- cons1 + width1
  B2 <- ifelse(B2 > salary, salary, B2)
  B1_p = B1 / salary
  B2_p = B2 / salary
  consum_vec <- seq(B1, B2, length.out = numC1)
  
  #re-optimisation for asset allocation [1]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [1]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B1_p,##
                       upper = B2_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons2 <- salary * optParam[1]
  
  # re-optimisation parameter [2]
  B3 <- cons2 - width2
  B3 <- ifelse(B3 < LB_C, LB_C, B3)
  B4 <- cons2 + width2
  B4 <- ifelse(B4 > salary, salary, B4)
  B3_p = B3 / salary
  B4_p = B4 / salary
  consum_vec <- seq(B3, B4, length.out = numC1)
  
  #re-optimisation for asset allocation [2]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [2]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B3_p,##
                       upper = B4_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons3 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [3]
  B5 <- cons3 - width3
  B5 <- ifelse(B5 < LB_C, LB_C, B5)
  B6 <- cons3 + width3
  B6 <- ifelse(B6 > salary, salary, B6)
  B5_p = B5 / salary
  B6_p = B6 / salary
  consum_vec <- seq(B5, B6, length.out = numC1)
  
  #re-optimisation for asset allocation [3]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  valfun1 = inv_result$val1
  valfun1_ord = valfun1[order(consum)]  # small to large
  valfun2 = inv_result$val2
  valfun2_ord = valfun2[order(consum)]  # small to large
  
  #re-optimisation for consumption [3]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B5_p,##
                       upper = B6_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  cons <- salary * optParam[1]
  
  eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
  val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
  val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = WT)
  
  val_1 = val_1_1 + val_1_2
  ending_balance = itBalStart - cons + salary
  aftax_cons = cons - income_tax(cons)
  
  resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[nrow(alloData)], 
                 consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility, 
                 pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
  resultsIt
}
rownames(result_xx) <- NULL
assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1]  
############################################################
# pre retirement ages (from age 25 to age 63)

# eqprop vector & return matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB

#age_x = 63
for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  Upper <- WT + (age_retire - age_x - 1) * Wt
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[(age_retire - age_x)] <- length(balances)
  
  consum_vec <- seq(LB_C, salary, length.out = numC) #consumption vector
  
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  
  itPrevBal = resultsPrev$balance
  itPrevBal_ord = itPrevBal[order(itPrevBal)]  # small to large
  itPrevpseudoCons = resultsPrev$pseudoConsum
  itPrevpseudoCons_ord = itPrevpseudoCons[order(itPrevBal)]  # small to large
  itPrevVal = resultsPrev$value_manager
  itPrevVal_ord = itPrevVal[order(itPrevBal)]  # small to large
  
  schu_spline1 <- Schumaker(x = itPrevBal_ord, y = itPrevpseudoCons_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline for pseudoconsum at age x+1
  
  A <- TAA_1[(age_x - 24)] * itRetA + TAA_2[(age_x - 24)] * itRetB  # set target return
  A_mat = matrix(rep(A, numPi), byrow = TRUE, ncol = length(A), nrow =  numPi)
  
  # compute prob1 matrix
  vec <- rep(1,length(consum_vec))
  cal1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(cal1_mat) / length(itRetA)
  
  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    set.seed(1)
    itBalStart <- balances[y]
    
    #optimisation for asset allocation
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_xx,
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]
    
    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC)
    
    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }  
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1_vec
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2_vec
    valfun2_ord = valfun2[order(consum)]  # small to large
    
    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_xx,
                         lower = B1_p,#
                         upper = B2_p,
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]
    
    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = Upper)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)
    
    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[(age_x - 24)],
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility,
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire - 1))##
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t7 <- do.call(rbind, l.df)

save(data_vfi_t7, file = paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens4_2_v5.1.RData"))##
write.csv(data_vfi_t7, paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens4_2_v5.1.csv"), row.names = F)##

###########################################################################################################################################

#================H1 Case sens5 (using vanguard benchmark)=========================#

## Initialisation  
rm(list=ls())

source("opt_fns.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores(),outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))

#importing target allocation data
alloData <- read.csv("Target_asset_allocation_sens.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 4 
l <- 0.7 # replacement ratio of post-retirement yearly consumption level
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
Wt <- 0.5 # weight assigned to benchmark A
WT <- 1 # weight assigned to target B
C_star <- salary - income_tax(salary) # consumption target (pre-retirement)
B_0 <- l * C_star # post-tax post-retirement consumption target

#starting balance at each age
factor <- 1.04
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
max_start_bal <- 3000000 ##maximum starting balance, as balances above these have same outcomes
start_bal_list$Start_Bal <- ifelse(start_bal_list$Start_Bal > max_start_bal, max_start_bal, start_bal_list$Start_Bal)
bal_spacing <- 5000 ## how fine the balance grid should be (max = 10000)

##Adjustments for Resolution
numC <- 51       # resolution for consumption vector 
numC1 <- 101     # resolution for consumption vector (for age 64)
numPi <- 6601   # resolution for equity allocation % vector (incremental value=(1.16+0.16)/6600=0.0002)
numPi1 <- 13201  # resolution for equity allocation % vector (incremental value=(1.16+0.16)/13200=0.0001) (for age 64)

width1 <- 1000 # first reoptimisation on consumption between +/- width1
width2 <- 100  # second reoptimisation on consumption between +/- width2
width3 <- 50   # third reoptimisation on consumption between +/- width3

#setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1.16 
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- -0.16 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation controls
controlDE1 <- list(reltol = 1e-10, itermax = 300)
controlDE2 <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()
###################################
## at terminal age - 1 (age 64)
age_x <- age_retire - 1
cat(paste("age_x", age_x, "in progress \n"))
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
A <- itRetA * TAA_1[nrow(alloData)] + itRetB * TAA_2[nrow(alloData)] #set target return
A_mat <- matrix(rep(A, numPi1), byrow = TRUE, ncol = length(A), nrow =  numPi1)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons 
LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
len_grid[1] <- length(balances)

# set consumption vector
consum_vec <- seq(LB_C, salary, length.out = numC1)

# compute prob 1 matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi1)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB
vec <- rep(1,length(consum_vec))
count1_mat <- (ret_mat - A_mat) > 0
prob1 <- rowSums(count1_mat) / length(itRetA)


result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
  
  itBalStart <- balances[y]
  set.seed(1)
  
  #optimisation for asset allocation
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #optimisation for consumption
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = c(LB_C_rate),
                       upper = c(1),
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons1 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [1]
  B1 <- cons1 - width1
  B1 <- ifelse(B1 < LB_C, LB_C, B1)
  B2 <- cons1 + width1
  B2 <- ifelse(B2 > salary, salary, B2)
  B1_p = B1 / salary
  B2_p = B2 / salary
  consum_vec <- seq(B1, B2, length.out = numC1)
  
  #re-optimisation for asset allocation [1]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [1]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B1_p,##
                       upper = B2_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons2 <- salary * optParam[1]
  
  # re-optimisation parameter [2]
  B3 <- cons2 - width2
  B3 <- ifelse(B3 < LB_C, LB_C, B3)
  B4 <- cons2 + width2
  B4 <- ifelse(B4 > salary, salary, B4)
  B3_p = B3 / salary
  B4_p = B4 / salary
  consum_vec <- seq(B3, B4, length.out = numC1)
  
  #re-optimisation for asset allocation [2]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [2]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B3_p,##
                       upper = B4_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons3 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [3]
  B5 <- cons3 - width3
  B5 <- ifelse(B5 < LB_C, LB_C, B5)
  B6 <- cons3 + width3
  B6 <- ifelse(B6 > salary, salary, B6)
  B5_p = B5 / salary
  B6_p = B6 / salary
  consum_vec <- seq(B5, B6, length.out = numC1)
  
  #re-optimisation for asset allocation [3]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  valfun1 = inv_result$val1
  valfun1_ord = valfun1[order(consum)]  # small to large
  valfun2 = inv_result$val2
  valfun2_ord = valfun2[order(consum)]  # small to large
  
  #re-optimisation for consumption [3]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B5_p,##
                       upper = B6_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  cons <- salary * optParam[1]
  
  eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
  val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
  val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = WT)
  
  val_1 = val_1_1 + val_1_2
  ending_balance = itBalStart - cons + salary
  aftax_cons = cons - income_tax(cons)
  
  resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[nrow(alloData)], 
                 consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility, 
                 pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
  resultsIt
}
rownames(result_xx) <- NULL
assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1]  
############################################################
# pre retirement ages (from age 25 to age 63)

# eqprop vector & return matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB

#age_x = 63
for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  Upper <- WT + (age_retire - age_x - 1) * Wt
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[(age_retire - age_x)] <- length(balances)
  
  consum_vec <- seq(LB_C, salary, length.out = numC) #consumption vector
  
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  
  itPrevBal = resultsPrev$balance
  itPrevBal_ord = itPrevBal[order(itPrevBal)]  # small to large
  itPrevpseudoCons = resultsPrev$pseudoConsum
  itPrevpseudoCons_ord = itPrevpseudoCons[order(itPrevBal)]  # small to large
  itPrevVal = resultsPrev$value_manager
  itPrevVal_ord = itPrevVal[order(itPrevBal)]  # small to large
  
  schu_spline1 <- Schumaker(x = itPrevBal_ord, y = itPrevpseudoCons_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline for pseudoconsum at age x+1
  
  A <- TAA_1[(age_x - 24)] * itRetA + TAA_2[(age_x - 24)] * itRetB  # set target return
  A_mat = matrix(rep(A, numPi), byrow = TRUE, ncol = length(A), nrow =  numPi)
  
  # compute prob1 matrix
  vec <- rep(1,length(consum_vec))
  cal1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(cal1_mat) / length(itRetA)
  
  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    set.seed(1)
    itBalStart <- balances[y]
    
    #optimisation for asset allocation
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_xx,
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]
    
    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC)
    
    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }  
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1_vec
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2_vec
    valfun2_ord = valfun2[order(consum)]  # small to large
    
    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_xx,
                         lower = B1_p,#
                         upper = B2_p,
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]
    
    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = Upper)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)
    
    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[(age_x - 24)],
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility,
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire - 1))##
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t7 <- do.call(rbind, l.df)

save(data_vfi_t7, file = paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens5_v5.1.RData"))##
write.csv(data_vfi_t7, paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens5_v5.1.csv"), row.names = F)##


###########################################################################################################################################

#================H1 Case sens6_1 (replacement rate = 0.6)=========================#

## Initialisation  
rm(list=ls())

source("opt_fns.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores(),outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))

#importing target allocation data
alloData <- read.csv("Target_asset_allocation.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 4 
l <- 0.6 # replacement ratio of post-retirement yearly consumption level
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
Wt <- 0.5 # weight assigned to benchmark A
WT <- 1 # weight assigned to target B
C_star <- salary - income_tax(salary) # consumption target (pre-retirement)
B_0 <- l * C_star # post-tax post-retirement consumption target

#starting balance at each age
factor <- 1.04
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
max_start_bal <- 3000000 ##maximum starting balance, as balances above these have same outcomes
start_bal_list$Start_Bal <- ifelse(start_bal_list$Start_Bal > max_start_bal, max_start_bal, start_bal_list$Start_Bal)
bal_spacing <- 5000 ## how fine the balance grid should be (max = 10000)

##Adjustments for Resolution
numC <- 51       # resolution for consumption vector 
numC1 <- 101     # resolution for consumption vector (for age 64)
numPi <- 6601   # resolution for equity allocation % vector (incremental value=(1.16+0.16)/6600=0.0002)
numPi1 <- 13201  # resolution for equity allocation % vector (incremental value=(1.16+0.16)/13200=0.0001) (for age 64)

width1 <- 1000 # first reoptimisation on consumption between +/- width1
width2 <- 100  # second reoptimisation on consumption between +/- width2
width3 <- 50   # third reoptimisation on consumption between +/- width3

#setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1.16 
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- -0.16 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation controls
controlDE1 <- list(reltol = 1e-10, itermax = 300)
controlDE2 <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()
###################################
## at terminal age - 1 (age 64)
age_x <- age_retire - 1
cat(paste("age_x", age_x, "in progress \n"))
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
A <- itRetA * TAA_1[nrow(alloData)] + itRetB * TAA_2[nrow(alloData)] #set target return
A_mat <- matrix(rep(A, numPi1), byrow = TRUE, ncol = length(A), nrow =  numPi1)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons 
LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
len_grid[1] <- length(balances)

# set consumption vector
consum_vec <- seq(LB_C, salary, length.out = numC1)

# compute prob 1 matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi1)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB
vec <- rep(1,length(consum_vec))
count1_mat <- (ret_mat - A_mat) > 0
prob1 <- rowSums(count1_mat) / length(itRetA)


result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
  
  itBalStart <- balances[y]
  set.seed(1)
  
  #optimisation for asset allocation
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #optimisation for consumption
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = c(LB_C_rate),
                       upper = c(1),
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons1 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [1]
  B1 <- cons1 - width1
  B1 <- ifelse(B1 < LB_C, LB_C, B1)
  B2 <- cons1 + width1
  B2 <- ifelse(B2 > salary, salary, B2)
  B1_p = B1 / salary
  B2_p = B2 / salary
  consum_vec <- seq(B1, B2, length.out = numC1)
  
  #re-optimisation for asset allocation [1]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [1]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B1_p,##
                       upper = B2_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons2 <- salary * optParam[1]
  
  # re-optimisation parameter [2]
  B3 <- cons2 - width2
  B3 <- ifelse(B3 < LB_C, LB_C, B3)
  B4 <- cons2 + width2
  B4 <- ifelse(B4 > salary, salary, B4)
  B3_p = B3 / salary
  B4_p = B4 / salary
  consum_vec <- seq(B3, B4, length.out = numC1)
  
  #re-optimisation for asset allocation [2]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [2]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B3_p,##
                       upper = B4_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons3 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [3]
  B5 <- cons3 - width3
  B5 <- ifelse(B5 < LB_C, LB_C, B5)
  B6 <- cons3 + width3
  B6 <- ifelse(B6 > salary, salary, B6)
  B5_p = B5 / salary
  B6_p = B6 / salary
  consum_vec <- seq(B5, B6, length.out = numC1)
  
  #re-optimisation for asset allocation [3]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  valfun1 = inv_result$val1
  valfun1_ord = valfun1[order(consum)]  # small to large
  valfun2 = inv_result$val2
  valfun2_ord = valfun2[order(consum)]  # small to large
  
  #re-optimisation for consumption [3]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B5_p,##
                       upper = B6_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  cons <- salary * optParam[1]
  
  eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
  val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
  val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = WT)
  
  val_1 = val_1_1 + val_1_2
  ending_balance = itBalStart - cons + salary
  aftax_cons = cons - income_tax(cons)
  
  resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[nrow(alloData)], 
                 consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility, 
                 pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
  resultsIt
}
rownames(result_xx) <- NULL
assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1]  
############################################################
# pre retirement ages (from age 25 to age 63)

# eqprop vector & return matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB

#age_x = 63
for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  Upper <- WT + (age_retire - age_x - 1) * Wt
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[(age_retire - age_x)] <- length(balances)
  
  consum_vec <- seq(LB_C, salary, length.out = numC) #consumption vector
  
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  
  itPrevBal = resultsPrev$balance
  itPrevBal_ord = itPrevBal[order(itPrevBal)]  # small to large
  itPrevpseudoCons = resultsPrev$pseudoConsum
  itPrevpseudoCons_ord = itPrevpseudoCons[order(itPrevBal)]  # small to large
  itPrevVal = resultsPrev$value_manager
  itPrevVal_ord = itPrevVal[order(itPrevBal)]  # small to large
  
  schu_spline1 <- Schumaker(x = itPrevBal_ord, y = itPrevpseudoCons_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline for pseudoconsum at age x+1
  
  A <- TAA_1[(age_x - 24)] * itRetA + TAA_2[(age_x - 24)] * itRetB  # set target return
  A_mat = matrix(rep(A, numPi), byrow = TRUE, ncol = length(A), nrow =  numPi)
  
  # compute prob1 matrix
  vec <- rep(1,length(consum_vec))
  cal1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(cal1_mat) / length(itRetA)
  
  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    set.seed(1)
    itBalStart <- balances[y]
    
    #optimisation for asset allocation
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_xx,
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]
    
    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC)
    
    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }  
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1_vec
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2_vec
    valfun2_ord = valfun2[order(consum)]  # small to large
    
    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_xx,
                         lower = B1_p,#
                         upper = B2_p,
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]
    
    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = Upper)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)
    
    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[(age_x - 24)],
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility,
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire - 1))##
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t7 <- do.call(rbind, l.df)

save(data_vfi_t7, file = paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens6_1_v5.1.RData"))##
write.csv(data_vfi_t7, paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens6_1_v5.1.csv"), row.names = F)##

###########################################################################################################################################

#================H1 Case sens6_2 (replacement rate = 0.8)=========================#

## Initialisation  
rm(list=ls())

source("opt_fns.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores(),outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))

#importing target allocation data
alloData <- read.csv("Target_asset_allocation.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 4 
l <- 0.8 # replacement ratio of post-retirement yearly consumption level
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
Wt <- 0.5 # weight assigned to benchmark A
WT <- 1 # weight assigned to target B
C_star <- salary - income_tax(salary) # consumption target (pre-retirement)
B_0 <- l * C_star # post-tax post-retirement consumption target

#starting balance at each age
factor <- 1.04
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
max_start_bal <- 3000000 ##maximum starting balance, as balances above these have same outcomes
start_bal_list$Start_Bal <- ifelse(start_bal_list$Start_Bal > max_start_bal, max_start_bal, start_bal_list$Start_Bal)
bal_spacing <- 5000 ## how fine the balance grid should be (max = 10000)

##Adjustments for Resolution
numC <- 51       # resolution for consumption vector 
numC1 <- 101     # resolution for consumption vector (for age 64)
numPi <- 6601   # resolution for equity allocation % vector (incremental value=(1.16+0.16)/6600=0.0002)
numPi1 <- 13201  # resolution for equity allocation % vector (incremental value=(1.16+0.16)/13200=0.0001) (for age 64)

width1 <- 1000 # first reoptimisation on consumption between +/- width1
width2 <- 100  # second reoptimisation on consumption between +/- width2
width3 <- 50   # third reoptimisation on consumption between +/- width3

#setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1.16 
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- -0.16 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation controls
controlDE1 <- list(reltol = 1e-10, itermax = 300)
controlDE2 <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()
###################################
## at terminal age - 1 (age 64)
age_x <- age_retire - 1
cat(paste("age_x", age_x, "in progress \n"))
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
A <- itRetA * TAA_1[nrow(alloData)] + itRetB * TAA_2[nrow(alloData)] #set target return
A_mat <- matrix(rep(A, numPi1), byrow = TRUE, ncol = length(A), nrow =  numPi1)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons 
LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
len_grid[1] <- length(balances)

# set consumption vector
consum_vec <- seq(LB_C, salary, length.out = numC1)

# compute prob 1 matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi1)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB
vec <- rep(1,length(consum_vec))
count1_mat <- (ret_mat - A_mat) > 0
prob1 <- rowSums(count1_mat) / length(itRetA)


result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
  
  itBalStart <- balances[y]
  set.seed(1)
  
  #optimisation for asset allocation
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #optimisation for consumption
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = c(LB_C_rate),
                       upper = c(1),
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons1 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [1]
  B1 <- cons1 - width1
  B1 <- ifelse(B1 < LB_C, LB_C, B1)
  B2 <- cons1 + width1
  B2 <- ifelse(B2 > salary, salary, B2)
  B1_p = B1 / salary
  B2_p = B2 / salary
  consum_vec <- seq(B1, B2, length.out = numC1)
  
  #re-optimisation for asset allocation [1]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [1]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B1_p,##
                       upper = B2_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons2 <- salary * optParam[1]
  
  # re-optimisation parameter [2]
  B3 <- cons2 - width2
  B3 <- ifelse(B3 < LB_C, LB_C, B3)
  B4 <- cons2 + width2
  B4 <- ifelse(B4 > salary, salary, B4)
  B3_p = B3 / salary
  B4_p = B4 / salary
  consum_vec <- seq(B3, B4, length.out = numC1)
  
  #re-optimisation for asset allocation [2]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [2]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B3_p,##
                       upper = B4_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons3 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [3]
  B5 <- cons3 - width3
  B5 <- ifelse(B5 < LB_C, LB_C, B5)
  B6 <- cons3 + width3
  B6 <- ifelse(B6 > salary, salary, B6)
  B5_p = B5 / salary
  B6_p = B6 / salary
  consum_vec <- seq(B5, B6, length.out = numC1)
  
  #re-optimisation for asset allocation [3]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  valfun1 = inv_result$val1
  valfun1_ord = valfun1[order(consum)]  # small to large
  valfun2 = inv_result$val2
  valfun2_ord = valfun2[order(consum)]  # small to large
  
  #re-optimisation for consumption [3]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B5_p,##
                       upper = B6_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  cons <- salary * optParam[1]
  
  eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
  val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
  val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = WT)
  
  val_1 = val_1_1 + val_1_2
  ending_balance = itBalStart - cons + salary
  aftax_cons = cons - income_tax(cons)
  
  resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[nrow(alloData)], 
                 consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility, 
                 pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
  resultsIt
}
rownames(result_xx) <- NULL
assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1]  
############################################################
# pre retirement ages (from age 25 to age 63)

# eqprop vector & return matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB

#age_x = 63
for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  Upper <- WT + (age_retire - age_x - 1) * Wt
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[(age_retire - age_x)] <- length(balances)
  
  consum_vec <- seq(LB_C, salary, length.out = numC) #consumption vector
  
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  
  itPrevBal = resultsPrev$balance
  itPrevBal_ord = itPrevBal[order(itPrevBal)]  # small to large
  itPrevpseudoCons = resultsPrev$pseudoConsum
  itPrevpseudoCons_ord = itPrevpseudoCons[order(itPrevBal)]  # small to large
  itPrevVal = resultsPrev$value_manager
  itPrevVal_ord = itPrevVal[order(itPrevBal)]  # small to large
  
  schu_spline1 <- Schumaker(x = itPrevBal_ord, y = itPrevpseudoCons_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline for pseudoconsum at age x+1
  
  A <- TAA_1[(age_x - 24)] * itRetA + TAA_2[(age_x - 24)] * itRetB  # set target return
  A_mat = matrix(rep(A, numPi), byrow = TRUE, ncol = length(A), nrow =  numPi)
  
  # compute prob1 matrix
  vec <- rep(1,length(consum_vec))
  cal1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(cal1_mat) / length(itRetA)
  
  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    set.seed(1)
    itBalStart <- balances[y]
    
    #optimisation for asset allocation
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_xx,
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]
    
    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC)
    
    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }  
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1_vec
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2_vec
    valfun2_ord = valfun2[order(consum)]  # small to large
    
    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_xx,
                         lower = B1_p,#
                         upper = B2_p,
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]
    
    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = Upper)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)
    
    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[(age_x - 24)],
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility,
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire - 1))##
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t7 <- do.call(rbind, l.df)

save(data_vfi_t7, file = paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens6_2_v5.1.RData"))##
write.csv(data_vfi_t7, paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens6_2_v5.1.csv"), row.names = F)##

###########################################################################################################################################

#================H1 Case sens7_1 (weight in target B = 0.5)=========================#

## Initialisation  
rm(list=ls())

source("opt_fns.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores(),outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))

#importing target allocation data
alloData <- read.csv("Target_asset_allocation.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 4 
l <- 0.7 # replacement ratio of post-retirement yearly consumption level
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
Wt <- 0.5 # weight assigned to benchmark A
WT <- 0.5 # weight assigned to target B
C_star <- salary - income_tax(salary) # consumption target (pre-retirement)
B_0 <- l * C_star # post-tax post-retirement consumption target

#starting balance at each age
factor <- 1.04
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
max_start_bal <- 3000000 ##maximum starting balance, as balances above these have same outcomes
start_bal_list$Start_Bal <- ifelse(start_bal_list$Start_Bal > max_start_bal, max_start_bal, start_bal_list$Start_Bal)
bal_spacing <- 5000 ## how fine the balance grid should be (max = 10000)

##Adjustments for Resolution
numC <- 51       # resolution for consumption vector 
numC1 <- 101     # resolution for consumption vector (for age 64)
numPi <- 6601   # resolution for equity allocation % vector (incremental value=(1.16+0.16)/6600=0.0002)
numPi1 <- 13201  # resolution for equity allocation % vector (incremental value=(1.16+0.16)/13200=0.0001) (for age 64)

width1 <- 1000 # first reoptimisation on consumption between +/- width1
width2 <- 100  # second reoptimisation on consumption between +/- width2
width3 <- 50   # third reoptimisation on consumption between +/- width3

#setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1.16 
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- -0.16 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation controls
controlDE1 <- list(reltol = 1e-10, itermax = 300)
controlDE2 <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()
###################################
## at terminal age - 1 (age 64)
age_x <- age_retire - 1
cat(paste("age_x", age_x, "in progress \n"))
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
A <- itRetA * TAA_1[nrow(alloData)] + itRetB * TAA_2[nrow(alloData)] #set target return
A_mat <- matrix(rep(A, numPi1), byrow = TRUE, ncol = length(A), nrow =  numPi1)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons 
LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
len_grid[1] <- length(balances)

# set consumption vector
consum_vec <- seq(LB_C, salary, length.out = numC1)

# compute prob 1 matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi1)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB
vec <- rep(1,length(consum_vec))
count1_mat <- (ret_mat - A_mat) > 0
prob1 <- rowSums(count1_mat) / length(itRetA)


result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
  
  itBalStart <- balances[y]
  set.seed(1)
  
  #optimisation for asset allocation
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #optimisation for consumption
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = c(LB_C_rate),
                       upper = c(1),
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons1 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [1]
  B1 <- cons1 - width1
  B1 <- ifelse(B1 < LB_C, LB_C, B1)
  B2 <- cons1 + width1
  B2 <- ifelse(B2 > salary, salary, B2)
  B1_p = B1 / salary
  B2_p = B2 / salary
  consum_vec <- seq(B1, B2, length.out = numC1)
  
  #re-optimisation for asset allocation [1]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [1]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B1_p,##
                       upper = B2_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons2 <- salary * optParam[1]
  
  # re-optimisation parameter [2]
  B3 <- cons2 - width2
  B3 <- ifelse(B3 < LB_C, LB_C, B3)
  B4 <- cons2 + width2
  B4 <- ifelse(B4 > salary, salary, B4)
  B3_p = B3 / salary
  B4_p = B4 / salary
  consum_vec <- seq(B3, B4, length.out = numC1)
  
  #re-optimisation for asset allocation [2]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [2]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B3_p,##
                       upper = B4_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons3 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [3]
  B5 <- cons3 - width3
  B5 <- ifelse(B5 < LB_C, LB_C, B5)
  B6 <- cons3 + width3
  B6 <- ifelse(B6 > salary, salary, B6)
  B5_p = B5 / salary
  B6_p = B6 / salary
  consum_vec <- seq(B5, B6, length.out = numC1)
  
  #re-optimisation for asset allocation [3]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  valfun1 = inv_result$val1
  valfun1_ord = valfun1[order(consum)]  # small to large
  valfun2 = inv_result$val2
  valfun2_ord = valfun2[order(consum)]  # small to large
  
  #re-optimisation for consumption [3]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B5_p,##
                       upper = B6_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  cons <- salary * optParam[1]
  
  eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
  val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
  val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = WT)
  
  val_1 = val_1_1 + val_1_2
  ending_balance = itBalStart - cons + salary
  aftax_cons = cons - income_tax(cons)
  
  resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[nrow(alloData)], 
                 consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility, 
                 pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
  resultsIt
}
rownames(result_xx) <- NULL
assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1]  
############################################################
# pre retirement ages (from age 25 to age 63)

# eqprop vector & return matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB

#age_x = 63
for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  Upper <- WT + (age_retire - age_x - 1) * Wt
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[(age_retire - age_x)] <- length(balances)
  
  consum_vec <- seq(LB_C, salary, length.out = numC) #consumption vector
  
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  
  itPrevBal = resultsPrev$balance
  itPrevBal_ord = itPrevBal[order(itPrevBal)]  # small to large
  itPrevpseudoCons = resultsPrev$pseudoConsum
  itPrevpseudoCons_ord = itPrevpseudoCons[order(itPrevBal)]  # small to large
  itPrevVal = resultsPrev$value_manager
  itPrevVal_ord = itPrevVal[order(itPrevBal)]  # small to large
  
  schu_spline1 <- Schumaker(x = itPrevBal_ord, y = itPrevpseudoCons_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline for pseudoconsum at age x+1
  
  A <- TAA_1[(age_x - 24)] * itRetA + TAA_2[(age_x - 24)] * itRetB  # set target return
  A_mat = matrix(rep(A, numPi), byrow = TRUE, ncol = length(A), nrow =  numPi)
  
  # compute prob1 matrix
  vec <- rep(1,length(consum_vec))
  cal1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(cal1_mat) / length(itRetA)
  
  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    set.seed(1)
    itBalStart <- balances[y]
    
    #optimisation for asset allocation
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_xx,
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]
    
    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC)
    
    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }  
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1_vec
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2_vec
    valfun2_ord = valfun2[order(consum)]  # small to large
    
    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_xx,
                         lower = B1_p,#
                         upper = B2_p,
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]
    
    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = Upper)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)
    
    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[(age_x - 24)],
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility,
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire - 1))##
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t7 <- do.call(rbind, l.df)

save(data_vfi_t7, file = paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens7_1_v5.1.RData"))##
write.csv(data_vfi_t7, paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens7_1_v5.1.csv"), row.names = F)##


###########################################################################################################################################

#================H1 Case sens7_2 (weight in target B = 5)=========================#

## Initialisation  
rm(list=ls())

source("opt_fns.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores(),outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))

#importing target allocation data
alloData <- read.csv("Target_asset_allocation.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 4 
l <- 0.7 # replacement ratio of post-retirement yearly consumption level
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free interest rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before age 50
limit2 <- 27000 # 401(k) contribution limit on or after age 50
age_tr <- 50 # threshold age for changing contribution limit
Wt <- 0.5 # weight assigned to benchmark A
WT <- 5 # weight assigned to target B
C_star <- salary - income_tax(salary) # consumption target (pre-retirement)
B_0 <- l * C_star # post-tax post-retirment consumption target

#starting balance at each age
factor <- 1.04
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
max_start_bal <- 3000000 ##maximum starting balance, as balances above these have same outcomes
start_bal_list$Start_Bal <- ifelse(start_bal_list$Start_Bal > max_start_bal, max_start_bal, start_bal_list$Start_Bal)
bal_spacing <- 5000 ## how fine the balance grid should be (max = 10000)

##Adjustments for Resolution
numC <- 51       # resolution for consumption vector 
numC1 <- 101     # resolution for consumption vector (for age 64)
numPi <- 6601   # resolution for equity allocation % vector (incremental value=(1.16+0.16)/6600=0.0002)
numPi1 <- 13201  # resolution for equity allocation % vector (incremental value=(1.16+0.16)/13200=0.0001) (for age 64)

width1 <- 1000 # first reoptimisation on consumption between +/- width1
width2 <- 100  # second reoptimisation on consumption between +/- width2
width3 <- 50   # third reoptimisation on consumption between +/- width3

#setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1.16 
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- -0.16 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation controls
controlDE1 <- list(reltol = 1e-10, itermax = 300)
controlDE2 <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()
###################################
## at terminal age - 1 (age 64)
age_x <- age_retire - 1
cat(paste("age_x", age_x, "in progress \n"))
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
A <- itRetA * TAA_1[nrow(alloData)] + itRetB * TAA_2[nrow(alloData)] #set target return
A_mat <- matrix(rep(A, numPi1), byrow = TRUE, ncol = length(A), nrow =  numPi1)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons 
LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
len_grid[1] <- length(balances)

# set consumption vector
consum_vec <- seq(LB_C, salary, length.out = numC1)

# compute prob 1 matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi1)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB
vec <- rep(1,length(consum_vec))
count1_mat <- (ret_mat - A_mat) > 0
prob1 <- rowSums(count1_mat) / length(itRetA)

result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
  
  itBalStart <- balances[y]
  set.seed(1)
  
  #optimisation for asset allocation
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #optimisation for consumption
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = c(LB_C_rate),
                       upper = c(1),
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons1 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [1]
  B1 <- cons1 - width1
  B1 <- ifelse(B1 < LB_C, LB_C, B1)
  B2 <- cons1 + width1
  B2 <- ifelse(B2 > salary, salary, B2)
  B1_p = B1 / salary
  B2_p = B2 / salary
  consum_vec <- seq(B1, B2, length.out = numC1)
  
  #re-optimisation for asset allocation [1]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [1]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B1_p,##
                       upper = B2_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons2 <- salary * optParam[1]
  
  # re-optimisation parameter [2]
  B3 <- cons2 - width2
  B3 <- ifelse(B3 < LB_C, LB_C, B3)
  B4 <- cons2 + width2
  B4 <- ifelse(B4 > salary, salary, B4)
  B3_p = B3 / salary
  B4_p = B4 / salary
  consum_vec <- seq(B3, B4, length.out = numC1)
  
  #re-optimisation for asset allocation [2]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  
  #re-optimisation for consumption [2]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B3_p,##
                       upper = B4_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  cons3 <- salary * optParam[1]
  
  ################################
  # re-optimisation parameter [3]
  B5 <- cons3 - width3
  B5 <- ifelse(B5 < LB_C, LB_C, B5)
  B6 <- cons3 + width3
  B6 <- ifelse(B6 > salary, salary, B6)
  B5_p = B5 / salary
  B6_p = B6 / salary
  consum_vec <- seq(B5, B6, length.out = numC1)
  
  #re-optimisation for asset allocation [3]
  eqprop <- NULL; val <- NULL; val1 <- NULL; val2 <- NULL
  
  
  for (j in seq_along(consum_vec)) {
    consumption <- consum_vec[j]
    itBalEnd_mat <- (itBalStart - consumption + salary) * (1 + ret_mat)
    inc_mat <- itBalEnd_mat / annuity1
    posttax_inc_mat <- inc_mat - income_tax(inc_mat)
    count2_mat <- (posttax_inc_mat - B_0) > 0
    prob2 <- rowSums(count2_mat) / length(itRetA)
    value1 <- Wt * prob1
    value2 <- WT * prob2
    value <- value1 + value2
    eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
    val[j] <- value[max(which(value == max(value)))]
    val1[j] <- value1[max(which(value == max(value)))]
    val2[j] <- value2[max(which(value == max(value)))]
  }
  
  inv_result <- data.frame(consum_vec, eqprop, val, val1, val2)
  
  consum = inv_result$consum_vec
  consum_ord = consum[order(consum)]  # small to large
  allo = inv_result$eqprop
  allo_ord = allo[order(consum)]  # small to large
  valfun1 = inv_result$val1
  valfun1_ord = valfun1[order(consum)]  # small to large
  valfun2 = inv_result$val2
  valfun2_ord = valfun2[order(consum)]  # small to large
  
  #re-optimisation for consumption [3]
  listOptim <- DEoptim(fn = utility_pre_terminal, 
                       lower = B5_p,##
                       upper = B6_p,
                       control = controlDE1)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  cons <- salary * optParam[1]
  
  eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
  val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
  val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = WT)
  
  val_1 = val_1_1 + val_1_2
  ending_balance = itBalStart - cons + salary
  aftax_cons = cons - income_tax(cons)
  
  resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[nrow(alloData)], 
                 consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility, 
                 pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
  resultsIt
}
rownames(result_xx) <- NULL
assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1]  
############################################################
# pre retirement ages (from age 25 to age 63)

# eqprop vector & return matrix
eqprop_vec = seq(L_bound, U_bound, length.out = numPi)
ret_mat = eqprop_vec %o% itRetA + (1 - eqprop_vec) %o% itRetB

#age_x = 63
for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  Upper <- WT + (age_retire - age_x - 1) * Wt
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum cons
  LB_C_rate <- ifelse(age_x < age_tr, LB_C1_rate, LB_C2_rate) #set minimum cons %
  len_grid[(age_retire - age_x)] <- length(balances)
  
  consum_vec <- seq(LB_C, salary, length.out = numC) #consumption vector
  
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  
  itPrevBal = resultsPrev$balance
  itPrevBal_ord = itPrevBal[order(itPrevBal)]  # small to large
  itPrevpseudoCons = resultsPrev$pseudoConsum
  itPrevpseudoCons_ord = itPrevpseudoCons[order(itPrevBal)]  # small to large
  itPrevVal = resultsPrev$value_manager
  itPrevVal_ord = itPrevVal[order(itPrevBal)]  # small to large
  
  schu_spline1 <- Schumaker(x = itPrevBal_ord, y = itPrevpseudoCons_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline for pseudoconsum at age x+1
  
  A <- TAA_1[(age_x - 24)] * itRetA + TAA_2[(age_x - 24)] * itRetB  # set target return
  A_mat = matrix(rep(A, numPi), byrow = TRUE, ncol = length(A), nrow =  numPi)
  
  # compute prob1 matrix
  vec <- rep(1,length(consum_vec))
  cal1_mat <- (ret_mat - A_mat) > 0
  prob1 <- rowSums(cal1_mat) / length(itRetA)
  
  result_xx <- foreach(y = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    set.seed(1)
    itBalStart <- balances[y]
    
    #optimisation for asset allocation
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    
    #optimisation for consumption
    listOptim <- DEoptim(fn = utility_xx,
                         lower = c(LB_C_rate),
                         upper = c(1),
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    cons1 <- salary * optParam[1]
    
    # re-optimisation parameter [1]
    B1 <- cons1 - width1
    B1 <- ifelse(B1 < LB_C, LB_C, B1)
    B2 <- cons1 + width1
    B2 <- ifelse(B2 > salary, salary, B2)
    B1_p = B1 / salary
    B2_p = B2 / salary
    consum_vec <- seq(B1, B2, length.out = numC)
    
    #re-optimisation for asset allocation [1]
    eqprop <- NULL; val_vec <- NULL; val1_vec <- NULL; val2_vec <- NULL
    
    for (j in seq_along(consum_vec)) {
      consumption <- consum_vec[j]
      itBalEnd_mat = (itBalStart - consumption + salary) * (1 + ret_mat)
      count2_mat = matrix(rep(0,length(itRetA)*length(eqprop_vec)),nrow=length(eqprop_vec),ncol=length(itRetA))
      for (ret_idx in 1:length(itRetA)) {
        count2_mat[,ret_idx] <- Interpol.Fun1(X.Data = itPrevBal_ord, Y.Data = itPrevVal_ord, X = itBalEnd_mat[,ret_idx], LB = 0, UB = Upper)
      }  
      prob2 <- rowSums(count2_mat) / length(itRetA)
      value1 <- Wt * prob1
      value2 <- prob2
      value <- value1 + value2
      eqprop[j] <- eqprop_vec[max(which(value == max(value)))]
      val_vec[j] <- value[max(which(value == max(value)))]
      val1_vec[j] <- value1[max(which(value == max(value)))]
      val2_vec[j] <- value2[max(which(value == max(value)))]
    }
    inv_result <- data.frame(consum_vec, eqprop, val_vec, val1_vec, val2_vec)
    
    consum = inv_result$consum_vec
    consum_ord = consum[order(consum)]  # small to large
    allo = inv_result$eqprop
    allo_ord = allo[order(consum)]  # small to large
    valfun1 = inv_result$val1_vec
    valfun1_ord = valfun1[order(consum)]  # small to large
    valfun2 = inv_result$val2_vec
    valfun2_ord = valfun2[order(consum)]  # small to large
    
    #re-optimisation for consumption [1]
    listOptim <- DEoptim(fn = utility_xx,
                         lower = B1_p,#
                         upper = B2_p,
                         control = controlDE2)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    cons <- salary * optParam[1]
    
    eqprop_1 <- Interpol.Fun3(Consum.Data = consum_ord, Allocation.Data = allo_ord, Cons = cons)
    val_1_1 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun1_ord, X = cons, LB = 0, UB = Wt)
    val_1_2 <-  Interpol.Fun1(X.Data = consum_ord, Y.Data = valfun2_ord, X = cons, LB = 0, UB = Upper)
    
    val_1 = val_1_1 + val_1_2
    ending_balance = itBalStart - cons + salary
    aftax_cons = cons - income_tax(cons)
    
    resultsIt <- c(age = age_x, balance = itBalStart, ending_balance = ending_balance, wA = eqprop_1, wA_T = TAA_1[(age_x - 24)],
                   consprop = cons/salary, pretax_cons = cons, aftertax_cons = aftax_cons, utility = utility,
                   pseudoConsum = newObj, value_manager = val_1, value_manager1 = val_1_1, value_manager2 = val_1_2)
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire - 1))##
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t7 <- do.call(rbind, l.df)

save(data_vfi_t7, file = paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens7_2_v5.1.RData"))##
write.csv(data_vfi_t7, paste0("mgmt_H1_rho_",rho, "salary", salary, "_sens7_2_v5.1.csv"), row.names = F)##
