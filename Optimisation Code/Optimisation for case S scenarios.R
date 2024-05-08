#================S Case BASE Scenario=========================#
## Initialisation  
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("opt_fns1.R")

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
mean(itRetA);mean(itRetB)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 4
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before threshold age
limit2 <- 27000 # 401(k) contribution limit on or after threshold age
age_tr <- 50 # threshold age for changing contribution limit

##start balance at each age
factor <- 1.05  #lower factor may cause extrapolations in optimisation process; higher factor may cause extrapolations in projection process 
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
bal_spacing <- 1000 

##setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- 0 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

##optimisation control parameters
controlDE <- list(reltol = 1e-10, itermax = 400)


start.time <- Sys.time()
##########################################################################
## age at retirement - 1 (age 64)
age_x <- age_retire - 1
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2)
UB_w <- salary/salary
LB_w <- LB_C/salary
len_grid[1] <- length(balances)

result_xx <- NULL

result_xx <- foreach(iter = 1:length(balances), .combine = rbind) %dopar% {
  M_t <- balances[iter]
  
  #compute W_t &  pi_t
  set.seed(1)
  listOptim <- DEoptim(fn = temp_fun1_S, 
                       lower = c(LB_w,L_bound),
                       upper = c(UB_w,U_bound),
                       control = controlDE)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  
  W_t <- salary * optParam[1]
  
  resultsIt <- c(age = age_x, balance = M_t, pretax_cons = W_t, aftax_cons = W_t - income_tax(W_t), ending_balance =  M_t - W_t + salary, 
                 consprop = optParam[1], eqprop = optParam[2], utility = utility, pseudoConsum = newObj)  
  resultsIt
}

rownames(result_xx) <- NULL
result_x64 <- data.frame(result_xx)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1] 


## pre retirement ages (from age 25 to age 63)

for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  #setting up new set of balances
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  result_xx <- NULL
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2)
  UB_w <- salary/salary
  LB_w <- LB_C/salary
  len_grid[(age_retire - age_x)] <- length(balances)
  
  xdata = resultsPrev$balance
  xdata_ord = xdata[order(xdata)]  # Sall to large
  ydata = resultsPrev$pseudoConsum
  ydata_ord = ydata[order(xdata)]  # Sall to large
  
  schu_spline <- Schumaker(x = xdata_ord, y = ydata_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline
  
  
  result_xx <- foreach(iter = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    M_t <- balances[iter]
    
    #compute W_t & pi_t
    set.seed(1)
    listOptim <- DEoptim(fn = temp_fun2_S, 
                         lower = c(LB_w,L_bound),
                         upper = c(UB_w,U_bound),
                         control = controlDE)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    
    W_t <- salary * optParam[1]
    
    resultsIt <- c(age = age_x, balance = M_t, pretax_cons = W_t, aftax_cons = W_t - income_tax(W_t), ending_balance =  M_t - W_t + salary, 
                   consprop = optParam[1], eqprop = optParam[2], utility = utility, pseudoConsum = newObj)  
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire-1))
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t5 <- do.call(rbind, l.df)
data_vfi_t5 <- data_vfi_t5[data_vfi_t5$balance >= 0,] # Just keeping positive values of starting balance. 

save(data_vfi_t5, file = paste0("worker_S_case_rho",rho, "_salary", salary, "_base_v5.RData"))##
write.csv(data_vfi_t5, paste0("worker_S_case_rho",rho, "_salary", salary, "_base_v5.csv"), row.names = F)##

#########################################################################################################################################

#================S case sens1_1 (salary = $55,000)=========================#
## Initialisation  
rm(list=ls())

source("opt_fns1.R")

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
mean(itRetA);mean(itRetB)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 55000
rho <- 4
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before threshold age
limit2 <- 27000 # 401(k) contribution limit on or after threshold age
age_tr <- 50 # threshold age for changing contribution limit

##start balance at each age
factor <- 1.05  #lower factor may cause extrapolations in optimisation process; higher factor may cause extrapolations in projection process 
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
bal_spacing <- 1000 

##setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- 0 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

##optimisation control parameters
controlDE <- list(reltol = 1e-10, itermax = 300)


start.time <- Sys.time()
##########################################################################
## age at retirement - 1 (age 64)
age_x <- age_retire - 1
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2)
UB_w <- salary/salary
LB_w <- LB_C/salary
len_grid[1] <- length(balances)

result_xx <- NULL


result_xx <- foreach(iter = 1:length(balances), .combine = rbind) %dopar% {
  M_t <- balances[iter]
  
  #compute W_t &  pi_t
  set.seed(1)
  listOptim <- DEoptim(fn = temp_fun1_S, 
                       lower = c(LB_w,L_bound),
                       upper = c(UB_w,U_bound),
                       control = controlDE)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  
  W_t <- salary * optParam[1]
  
  resultsIt <- c(age = age_x, balance = M_t, pretax_cons = W_t, aftax_cons = W_t - income_tax(W_t), ending_balance =  M_t - W_t + salary, 
                 consprop = optParam[1], eqprop = optParam[2], utility = utility, pseudoConsum = newObj)  
  resultsIt
}

rownames(result_xx) <- NULL
result_x64 <- data.frame(result_xx)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1] 


## pre retirement ages (from age 25 to age 63)

for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  #setting up new set of balances
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  result_xx <- NULL
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2)
  UB_w <- salary/salary
  LB_w <- LB_C/salary
  len_grid[(age_retire - age_x)] <- length(balances)
  
  xdata = resultsPrev$balance
  xdata_ord = xdata[order(xdata)]  # Sall to large
  ydata = resultsPrev$pseudoConsum
  ydata_ord = ydata[order(xdata)]  # Sall to large
  
  schu_spline <- Schumaker(x = xdata_ord, y = ydata_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline
  
  
  result_xx <- foreach(iter = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    M_t <- balances[iter]
    
    #compute W_t & pi_t
    set.seed(1)
    listOptim <- DEoptim(fn = temp_fun2_S, 
                         lower = c(LB_w,L_bound),
                         upper = c(UB_w,U_bound),
                         control = controlDE)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    
    W_t <- salary * optParam[1]
    
    resultsIt <- c(age = age_x, balance = M_t, pretax_cons = W_t, aftax_cons = W_t - income_tax(W_t), ending_balance =  M_t - W_t + salary, 
                   consprop = optParam[1], eqprop = optParam[2], utility = utility, pseudoConsum = newObj)  
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire-1))
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t5 <- do.call(rbind, l.df)
data_vfi_t5 <- data_vfi_t5[data_vfi_t5$balance >= 0,] # Just keeping positive values of starting balance. 

save(data_vfi_t5, file = paste0("worker_S_case_rho",rho, "_salary", salary, "_sens1_1_v5.RData"))##
write.csv(data_vfi_t5, paste0("worker_S_case_rho",rho, "_salary", salary, "_sens1_1_v5.csv"), row.names = F)##


#########################################################################################################################################

#================S case sens1_2 (salary = $85,000)=========================#
## Initialisation  
rm(list=ls())

source("opt_fns1.R")

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
mean(itRetA);mean(itRetB)

#importing target allocation data
alloData <- read.csv("Target_asset_allocation.csv", header = T, stringsAsFactors = F)[-1] ###
TAA_1 <- as.numeric(alloData$TAA_1)
TAA_2 <- as.numeric(alloData$TAA_2)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 85000
rho <- 4
annuity1 <- 19.2745 # immediate life annuity at age 65 (based on risk-free rate)
annuity2 <- 13.7469 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 19500 # 401(k) contribution limit before threshold age
limit2 <- 26000 # 401(k) contribution limit on or after threshold age
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
age_tr <- 50 # threshold age for changing contribution limit

##start balance at each age
factor <- 1.09   #lower factor may cause extrapolations in optimisation process; higher factor may cause extrapolations in projection process 
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
bal_spacing <- 1000 

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

#optimisation control pparameters
controlDE <- list(reltol = 1e-10, itermax = 200)

start.time <- Sys.time()

## at retirement - 1 (age 64)
age_x <- age_retire - 1
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2) #set minimum pretax consumption
result_xx <- NULL

UB_w <- salary/salary
LB_w <- LB_C/salary
len_grid[1] <- length(balances)

# set benchmark investment strategy
y <- TAA_1[nrow(alloData)]


result_xx <- foreach(iter = 1:length(balances), .combine = rbind) %dopar% {
  M_t <- balances[iter]
  
  #compute W_t &  pi_t
  set.seed(1)
  listOptim <- DEoptim(fn = temp_fun1_S, 
                       lower = c(LB_w,L_bound),
                       upper = c(UB_w,U_bound),
                       control = controlDE)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  
  W_t <- salary * optParam[1]
  
  resultsIt <- c(age = age_x, balance = M_t, pretax_cons = W_t, aftax_cons = W_t - income_tax(W_t), ending_balance =  M_t - W_t + salary, 
                 consprop = optParam[1], eqprop = y, utility = utility, pseudoConsum = newObj)  
  resultsIt
}

rownames(result_xx) <- NULL
result_x64 <- data.frame(result_xx)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1] 

## pre retirement ages (from age 25 to age 63)

for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  #set up new set of balances
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  result_xx <- NULL
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2)
  UB_w <- salary/salary
  LB_w <- LB_C/salary
  
  #determine benchmark investment strategy
  y <- TAA_1[(age_x - 24)]
  len_grid[(age_retire - age_x)] <- length(balances)
  
  xdata = resultsPrev$balance
  xdata_ord = xdata[order(xdata)]  # Sall to large
  ydata = resultsPrev$pseudoConsum
  ydata_ord = ydata[order(xdata)]  # Sall to large
  
  schu_spline <- Schumaker(x = xdata_ord, y = ydata_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline
  
  
  result_xx <- foreach(iter = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    M_t <- balances[iter]
    
    #compute W_t
    set.seed(1)
    listOptim <- DEoptim(fn = temp_fun2_S, 
                         lower = c(LB_w,L_bound),
                         upper = c(UB_w,L_bound),
                         control = controlDE)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    
    W_t <- salary * optParam[1]
    
    resultsIt <- c(age = age_x, balance = M_t, pretax_cons = W_t, aftax_cons = W_t - income_tax(W_t), ending_balance =  M_t - W_t + salary, 
                   consprop = optParam[1], eqprop = y, utility = utility, pseudoConsum = newObj)  
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire-1))
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t6 <- do.call(rbind, l.df)
data_vfi_t6 <- data_vfi_t6[data_vfi_t6$balance >= 0,] # Just keeping positive values of starting balance. 

save(data_vfi_t6, file = paste0("worker_S_case_rho",rho, "_salary", salary, "_sens1_2_v5.RData"))##
write.csv(data_vfi_t6, paste0("worker_S_case_rho",rho, "_salary", salary, "_sens1_2_v5.csv"), row.names = F)##


#########################################################################################################################################

#================S case sens2_1 (rho = 3)=========================#
## Initialisation  
rm(list=ls())

source("opt_fns1.R")

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
mean(itRetA);mean(itRetB)


##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 3
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before threshold age
limit2 <- 27000 # 401(k) contribution limit on or after threshold age
age_tr <- 50 # threshold age for changing contribution limit

##start balance at each age
factor <- 1.05  #lower factor may cause extrapolations in optimisation process; higher factor may cause extrapolations in projection process 
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
bal_spacing <- 1000 

##setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- 0 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

##optimisation control parameters
controlDE <- list(reltol = 1e-10, itermax = 300)


start.time <- Sys.time()
##########################################################################
## age at retirement - 1 (age 64)
age_x <- age_retire - 1
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2)
UB_w <- salary/salary
LB_w <- LB_C/salary
len_grid[1] <- length(balances)

result_xx <- NULL

result_xx <- foreach(iter = 1:length(balances), .combine = rbind) %dopar% {
  M_t <- balances[iter]
  
  #compute W_t &  pi_t
  set.seed(1)
  listOptim <- DEoptim(fn = temp_fun1_S, 
                       lower = c(LB_w,L_bound),
                       upper = c(UB_w,U_bound),
                       control = controlDE)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  
  W_t <- salary * optParam[1]
  
  resultsIt <- c(age = age_x, balance = M_t, pretax_cons = W_t, aftax_cons = W_t - income_tax(W_t), ending_balance =  M_t - W_t + salary, 
                 consprop = optParam[1], eqprop = optParam[2], utility = utility, pseudoConsum = newObj)  
  resultsIt
}

rownames(result_xx) <- NULL
result_x64 <- data.frame(result_xx)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1] 


## pre retirement ages (from age 25 to age 63)

for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  #setting up new set of balances
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  result_xx <- NULL
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2)
  UB_w <- salary/salary
  LB_w <- LB_C/salary
  len_grid[(age_retire - age_x)] <- length(balances)
  
  xdata = resultsPrev$balance
  xdata_ord = xdata[order(xdata)]  # Sall to large
  ydata = resultsPrev$pseudoConsum
  ydata_ord = ydata[order(xdata)]  # Sall to large
  
  schu_spline <- Schumaker(x = xdata_ord, y = ydata_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline
  
  
  result_xx <- foreach(iter = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    M_t <- balances[iter]
    
    #compute W_t & pi_t
    set.seed(1)
    listOptim <- DEoptim(fn = temp_fun2_S, 
                         lower = c(LB_w,L_bound),
                         upper = c(UB_w,U_bound),
                         control = controlDE)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    
    W_t <- salary * optParam[1]
    
    resultsIt <- c(age = age_x, balance = M_t, pretax_cons = W_t, aftax_cons = W_t - income_tax(W_t), ending_balance =  M_t - W_t + salary, 
                   consprop = optParam[1], eqprop = optParam[2], utility = utility, pseudoConsum = newObj)  
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire-1))
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t5 <- do.call(rbind, l.df)
data_vfi_t5 <- data_vfi_t5[data_vfi_t5$balance >= 0,] # Just keeping positive values of starting balance. 

save(data_vfi_t5, file = paste0("worker_S_case_rho",rho, "_salary", salary, "_sens2_1_v5.RData"))##
write.csv(data_vfi_t5, paste0("worker_S_case_rho",rho, "_salary", salary, "_sens2_1_v5.csv"), row.names = F)##


#########################################################################################################################################

#================S case sens2_2 (rho = 5)=========================#
## Initialisation  
rm(list=ls())

source("opt_fns1.R")

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
mean(itRetA);mean(itRetB)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 5
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before threshold age
limit2 <- 27000 # 401(k) contribution limit on or after threshold age
age_tr <- 50 # threshold age for changing contribution limit

##start balance at each age
factor <- 1.05  #lower factor may cause extrapolations in optimisation process; higher factor may cause extrapolations in projection process 
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
bal_spacing <- 1000 

##setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- 0 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

##optimisation control parameters
controlDE <- list(reltol = 1e-10, itermax = 300)


start.time <- Sys.time()
##########################################################################
## age at retirement - 1 (age 64)
age_x <- age_retire - 1
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2)
UB_w <- salary/salary
LB_w <- LB_C/salary
len_grid[1] <- length(balances)

result_xx <- NULL

result_xx <- foreach(iter = 1:length(balances), .combine = rbind) %dopar% {
  M_t <- balances[iter]
  
  #compute W_t &  pi_t
  set.seed(1)
  listOptim <- DEoptim(fn = temp_fun1_S, 
                       lower = c(LB_w,L_bound),
                       upper = c(UB_w,U_bound),
                       control = controlDE)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  
  W_t <- salary * optParam[1]
  
  resultsIt <- c(age = age_x, balance = M_t, pretax_cons = W_t, aftax_cons = W_t - income_tax(W_t), ending_balance =  M_t - W_t + salary, 
                 consprop = optParam[1], eqprop = optParam[2], utility = utility, pseudoConsum = newObj)  
  resultsIt
}

rownames(result_xx) <- NULL
result_x64 <- data.frame(result_xx)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1] 


## pre retirement ages (from age 25 to age 63)

for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  #setting up new set of balances
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  result_xx <- NULL
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2)
  UB_w <- salary/salary
  LB_w <- LB_C/salary
  len_grid[(age_retire - age_x)] <- length(balances)
  
  xdata = resultsPrev$balance
  xdata_ord = xdata[order(xdata)]  # Sall to large
  ydata = resultsPrev$pseudoConsum
  ydata_ord = ydata[order(xdata)]  # Sall to large
  
  schu_spline <- Schumaker(x = xdata_ord, y = ydata_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline
  
  
  result_xx <- foreach(iter = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    M_t <- balances[iter]
    
    #compute W_t & pi_t
    set.seed(1)
    listOptim <- DEoptim(fn = temp_fun2_S, 
                         lower = c(LB_w,L_bound),
                         upper = c(UB_w,U_bound),
                         control = controlDE)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    
    W_t <- salary * optParam[1]
    
    resultsIt <- c(age = age_x, balance = M_t, pretax_cons = W_t, aftax_cons = W_t - income_tax(W_t), ending_balance =  M_t - W_t + salary, 
                   consprop = optParam[1], eqprop = optParam[2], utility = utility, pseudoConsum = newObj)  
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire-1))
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t5 <- do.call(rbind, l.df)
data_vfi_t5 <- data_vfi_t5[data_vfi_t5$balance >= 0,] # Just keeping positive values of starting balance. 

save(data_vfi_t5, file = paste0("worker_S_case_rho",rho, "_salary", salary, "_sens2_2_v5.RData"))##
write.csv(data_vfi_t5, paste0("worker_S_case_rho",rho, "_salary", salary, "_sens2_2_v5.csv"), row.names = F)##


#########################################################################################################################################
#Note: Since sense 3_1 and sense 3_2 were found to be inappropriate for inclusion in our analysis, the relevant code has been omitted.

#================S case sens4_1 (-1% ERP)=========================#
## Initialisation  
rm(list=ls())

source("opt_fns1.R")

libs <- c("parallel", "doParallel", "foreach", "doSNOW", "DEoptim", "Hmisc", "plyr", "dplyr", "schumaker")
ipak(libs)

#makeCluster utilises parallel processing and makes code print to console
clus <- makeCluster(detectCores(),outfile = "")
#Declare the package(s) to be used at each node
decPackDeopt <- clusterEvalQ(clus,library(DEoptim))
decPackDeopt <- clusterEvalQ(clus,library(Hmisc))
registerDoSNOW(clus)

#importing historical return data
returnData <- read.csv("Ret_2_asset_sens1.csv", header = T, stringsAsFactors = F)[-1]
itRetA <- as.numeric(returnData$DE)
itRetB <- as.numeric(returnData$Cash)
itRetProb <- rep(1 / nrow(returnData), nrow(returnData))
mean(itRetA);mean(itRetB)
#0.05; 0.01

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 4
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before threshold age
limit2 <- 27000 # 401(k) contribution limit on or after threshold age
age_tr <- 50 # threshold age for changing contribution limit

##start balance at each age
factor <- 1.05  #lower factor may cause extrapolations in optimisation process; higher factor may cause extrapolations in projection process 
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
bal_spacing <- 1000 

##setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- 0 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

##optimisation control parameters
controlDE <- list(reltol = 1e-10, itermax = 300)


start.time <- Sys.time()
##########################################################################
## age at retirement - 1 (age 64)
age_x <- age_retire - 1
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2)
UB_w <- salary/salary
LB_w <- LB_C/salary
len_grid[1] <- length(balances)

result_xx <- NULL


result_xx <- foreach(iter = 1:length(balances), .combine = rbind) %dopar% {
  M_t <- balances[iter]
  
  #compute W_t &  pi_t
  set.seed(1)
  listOptim <- DEoptim(fn = temp_fun1_S, 
                       lower = c(LB_w,L_bound),
                       upper = c(UB_w,U_bound),
                       control = controlDE)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  
  W_t <- salary * optParam[1]
  
  resultsIt <- c(age = age_x, balance = M_t, pretax_cons = W_t, aftax_cons = W_t - income_tax(W_t), ending_balance =  M_t - W_t + salary, 
                 consprop = optParam[1], eqprop = optParam[2], utility = utility, pseudoConsum = newObj)  
  resultsIt
}

rownames(result_xx) <- NULL
result_x64 <- data.frame(result_xx)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1] 


## pre retirement ages (from age 25 to age 63)

for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  #setting up new set of balances
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  result_xx <- NULL
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2)
  UB_w <- salary/salary
  LB_w <- LB_C/salary
  len_grid[(age_retire - age_x)] <- length(balances)
  
  xdata = resultsPrev$balance
  xdata_ord = xdata[order(xdata)]  # Sall to large
  ydata = resultsPrev$pseudoConsum
  ydata_ord = ydata[order(xdata)]  # Sall to large
  
  schu_spline <- Schumaker(x = xdata_ord, y = ydata_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline
  
  
  result_xx <- foreach(iter = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    M_t <- balances[iter]
    
    #compute W_t & pi_t
    set.seed(1)
    listOptim <- DEoptim(fn = temp_fun2_S, 
                         lower = c(LB_w,L_bound),
                         upper = c(UB_w,U_bound),
                         control = controlDE)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    
    W_t <- salary * optParam[1]
    
    resultsIt <- c(age = age_x, balance = M_t, pretax_cons = W_t, aftax_cons = W_t - income_tax(W_t), ending_balance =  M_t - W_t + salary, 
                   consprop = optParam[1], eqprop = optParam[2], utility = utility, pseudoConsum = newObj)  
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire-1))
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t5 <- do.call(rbind, l.df)
data_vfi_t5 <- data_vfi_t5[data_vfi_t5$balance >= 0,] # Just keeping positive values of starting balance. 

save(data_vfi_t5, file = paste0("worker_S_case_rho",rho, "_salary", salary, "_sens4_1_v5.RData"))##
write.csv(data_vfi_t5, paste0("worker_S_case_rho",rho, "_salary", salary, "_sens4_1_v5.csv"), row.names = F)##


#########################################################################################################################################

#================S case sens4_2 (+1% ERP)=========================#
## Initialisation  
rm(list=ls())

source("opt_fns1.R")

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
mean(itRetA);mean(itRetB)

##Parameters
age_retire <- 65
age_start <- 25
salary <- 70000
rho <- 4
annuity1 <- 18.1313 # immediate life annuity at age 65 (based on risk-free rate)
annuity2 <- 13.3559 # immediate life annuity at age 65 (based on subjective discount factor)
beta <- 0.96 # subjective discount factor
limit1 <- 20500 # 401(k) contribution limit before threshold age
limit2 <- 27000 # 401(k) contribution limit on or after threshold age
age_tr <- 50 # threshold age for changing contribution limit

##start balance at each age
factor <- 1.05  #lower factor may cause extrapolations in optimisation process; higher factor may cause extrapolations in projection process 
age_term_bal <- 5000000
start_bal_list <- bal_list(age_term_bal, factor)
bal_spacing <- 1000 

##setting upper bound & lower bound for the equity allocation proportion
U_bound1 <- 1
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
U_bound <- min(U_bound1, U_bound2) 
L_bound1 <- 0 
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
L_bound <- max(L_bound1, L_bound2) 

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50

#set vectors counting number of balance grids/ running time each age
len_grid <- rep(0,(age_retire - age_start))
run_time <- rep(0,(age_retire - age_start))

##optimisation control parameters
controlDE <- list(reltol = 1e-10, itermax = 300)


start.time <- Sys.time()
##########################################################################
## age at retirement - 1 (age 64)
age_x <- age_retire - 1
bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
balances <- fnbalance(bal_start, bal_spacing)
LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2)
UB_w <- salary/salary
LB_w <- LB_C/salary
len_grid[1] <- length(balances)

result_xx <- NULL


result_xx <- foreach(iter = 1:length(balances), .combine = rbind) %dopar% {
  M_t <- balances[iter]
  
  #compute W_t &  pi_t
  set.seed(1)
  listOptim <- DEoptim(fn = temp_fun1_S, 
                       lower = c(LB_w,L_bound),
                       upper = c(UB_w,U_bound),
                       control = controlDE)
  optParam <- listOptim$optim$bestmem
  optObjVal <- - listOptim$optim$bestval
  utility <- optObjVal
  newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
  
  W_t <- salary * optParam[1]
  
  resultsIt <- c(age = age_x, balance = M_t, pretax_cons = W_t, aftax_cons = W_t - income_tax(W_t), ending_balance =  M_t - W_t + salary, 
                 consprop = optParam[1], eqprop = optParam[2], utility = utility, pseudoConsum = newObj)  
  resultsIt
}

rownames(result_xx) <- NULL
result_x64 <- data.frame(result_xx)

end.time <- Sys.time()
run_time[1] <- as.numeric(round(end.time - start.time,2), units = "mins")
run_time[1] 


## pre retirement ages (from age 25 to age 63)

for (age_x in rev(age_start:(age_retire - 2))) {
  cat(paste("age_x", age_x, "in progress \n"))
  
  start.time <- Sys.time()
  
  #setting up new set of balances
  resultsPrev <- as.data.frame(get(paste0("result_x", age_x + 1)))
  result_xx <- NULL
  bal_start <- start_bal_list[which(start_bal_list$Age == age_x), 2]
  balances <- fnbalance(bal_start, bal_spacing)
  LB_C <- ifelse(age_x < age_tr, LB_C1, LB_C2)
  UB_w <- salary/salary
  LB_w <- LB_C/salary
  len_grid[(age_retire - age_x)] <- length(balances)
  
  xdata = resultsPrev$balance
  xdata_ord = xdata[order(xdata)]  # Sall to large
  ydata = resultsPrev$pseudoConsum
  ydata_ord = ydata[order(xdata)]  # Sall to large
  
  schu_spline <- Schumaker(x = xdata_ord, y = ydata_ord, Vectorised = TRUE, Extrapolation = c("Linear")) # fit a schumaker spline
  
  
  result_xx <- foreach(iter = 1:length(balances), .combine = rbind, .packages = "schumaker") %dopar% {
    M_t <- balances[iter]
    
    #compute W_t & pi_t
    set.seed(1)
    listOptim <- DEoptim(fn = temp_fun2_S, 
                         lower = c(LB_w,L_bound),
                         upper = c(UB_w,U_bound),
                         control = controlDE)
    optParam <- listOptim$optim$bestmem
    optObjVal <- - listOptim$optim$bestval
    utility <- optObjVal
    newObj <- (utility * (1 - rho)) ^ (1 / (1 - rho))
    
    W_t <- salary * optParam[1]
    
    resultsIt <- c(age = age_x, balance = M_t, pretax_cons = W_t, aftax_cons = W_t - income_tax(W_t), ending_balance =  M_t - W_t + salary, 
                   consprop = optParam[1], eqprop = optParam[2], utility = utility, pseudoConsum = newObj)  
    resultsIt
  }
  rownames(result_xx) <- NULL
  assign(paste0("result_x", age_x), result_xx, envir = .GlobalEnv)
  
  end.time <- Sys.time()
  run_time[(age_retire - age_x)]  <- as.numeric(round(end.time - start.time,2), units = "mins")
  run_time[(age_retire - age_x)]
}

all.data.list <- paste0("result_x", age_start:(age_retire-1))
l.df <- lapply(all.data.list, function(x) as.data.frame(get(x)))
data_vfi_t5 <- do.call(rbind, l.df)
data_vfi_t5 <- data_vfi_t5[data_vfi_t5$balance >= 0,] # Just keeping positive values of starting balance. 

save(data_vfi_t5, file = paste0("worker_S_case_rho",rho, "_salary", salary, "_sens4_2_v5.RData"))##
write.csv(data_vfi_t5, paste0("worker_S_case_rho",rho, "_salary", salary, "_sens4_2_v5.csv"), row.names = F)##

#########################################################################################################################################