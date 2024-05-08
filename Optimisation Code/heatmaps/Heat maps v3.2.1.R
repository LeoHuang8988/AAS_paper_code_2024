## Heat maps 

## Initialisation  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

#########################################################
### (A) Installing and loading required packages
#########################################################

ipak <- function(pkg){
  # function to install and require packages
  # https://gist.github.com/stevenworthington/3178163
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
libs <- c("ggplot2", "RColorBrewer", "scales", "grid", "Hmisc", "parallel", "gridExtra", "gtools", "lattice")
ipak(libs)
#####################################################################
### (B) Reading in data and transforming it into appropriate format
#####################################################################
SM_base = read.csv("worker_S_case_rho4_salary70000_base_v5.csv")
BM_base = read.csv("worker_B_case_rho4_salary70000_base_v5.csv")
HM_base1 = read.csv("mgmt_H1_rho_4salary70000_base1_v5.1.csv")
HM_base2 = read.csv("mgmt_H2_rho_4salary70000_base1_v5.1.csv")

## parameters
age_retire <- 65
age_start <- 25         
salary <- 70000
limit1 <- 20500 # 401(k) contribution limit before threshold age
limit2 <- 27000 # 401(k) contribution limit on or after threshold age
age_tr <- 50 # threshold age for changing contribution limit

## function to do linear interpolation (applies to eqprop(w/m) & consprop(w/m))
Interpol.Fun <- function(Balance.Data, Allocation.Data, Bals, Rule, LB, UB){   
  Interpol.Obj <- approxExtrap(x = Balance.Data, y = Allocation.Data, xout = Bals, rule = Rule)
  Interpol.Obj <- c(do.call("cbind", Interpol.Obj[2]))
  Interpol.Obj <- ifelse(Interpol.Obj < LB, LB, Interpol.Obj)
  Interpol.Obj <- ifelse(Interpol.Obj > UB, UB, Interpol.Obj)
  return(Interpol.Obj)
}

##setting lower bound for consumption
LB_C1 <- salary - limit1 # lower bound for consumption before age 50
LB_C1_rate <- LB_C1 / salary
LB_C2 <- salary - limit2 # lower bound for consumption on or after age 50
LB_C2_rate <- LB_C2 / salary

##setting upper bound & lower bound for risky asset allocation proportion (w/m)
return_2assets <- read.csv("Ret_2_asset.csv", header = T)
itRetA <- as.numeric(return_2assets$DE)
itRetB <- as.numeric(return_2assets$Cash)

U_bound1_1 <- 1
U_bound1_2 <- 1.16
U_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - min(itRetA))  
UB_pi_w <- min(U_bound1_1, U_bound2) 
UB_pi_m<- min(U_bound1_2, U_bound2)

L_bound1_1 <- 0 
L_bound1_2 <- -0.16
L_bound2 <- (1 + unique(itRetB))/(unique(itRetB) - max(itRetA))
LB_pi_w <- max(L_bound1_1, L_bound2) 
LB_pi_m<- max(L_bound1_2, L_bound2) 

##set age sequence and compute length of age sequence
age_seq <- seq(25,64,by=1)
sim_year_total <- length(age_seq)

##set balance sequence
bal_seq <- seq(0,2000000,by=10000)

# (1) worker's case [optimisation based on VFI - case S]
bal_seq_w <- bal_seq  

con_vec <- NULL
pi_vec <- NULL
aa_temp <- NULL
consprop_vec <- NULL
allocation_vec <- NULL
age_vec <- NULL
balance_vec <- NULL

for (i in 1:sim_year_total){
  age_idx <- 25 + i - 1
  LB_C_rate <- ifelse(age_idx < age_tr, LB_C1_rate, LB_C2_rate)
  aa_temp <- SM_base[rev(which(SM_base$age == age_idx)), c(2,6,7)]
  con_vec <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_seq_w, Outlier, LB = LB_C_rate, UB = 1) #consumption rate
  consprop_vec <- c(consprop_vec, con_vec)
  pi_vec <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$eqprop.par2, Bals = bal_seq_w, Outlier, LB = LB_pi_w, UB = UB_pi_w) #allocation rate
  allocation_vec <- c(allocation_vec, pi_vec)
  age_vec <- c(age_vec, rep(age_idx,length(con_vec)))
  balance_vec <- c(balance_vec, bal_seq_w)
}

data_sm <- cbind.data.frame(age_vec,balance_vec,consprop_vec,allocation_vec)
data_sm$balance_vec <- data_sm$balance_vec / 1000
rangeCon <- range(data_sm$consprop_vec)
lowCon <- 0.6
upCon <- 1
rangeAllo <- range(data_sm$allocation_vec)
lowAllo <- 0
upAllo <- 1.2

# heatplot with consumption % (worker)
con_sm <- ggplot(data = data_sm, aes(x=age_vec,y=balance_vec))+
  geom_tile(aes(fill = consprop_vec))+ 
  scale_fill_gradient2(low =  "darkred",mid = "white", high = "darkblue",name = "",midpoint = mean(rangeCon),limits = c(lowCon,upCon),labels = percent(round(seq(lowCon,upCon,(upCon-lowCon)/4),2)),breaks = seq(lowCon,upCon,(upCon-lowCon)/4),guide = "colorbar") + 
  xlab("Age") + 
  ylab("Balance ($'000)") + 
  scale_x_discrete(limits=seq(25,65,by=5)) + 
  scale_y_continuous(limits = c(0,2000)) +
  theme(legend.position="bottom") + 
  theme(text=element_text(size=9), axis.text = element_text(size=9))

jpeg("heatmap_S_con_v5.1.1.jpeg", height = 8, width = 10, units = 'cm', res = 300)
con_sm
dev.off()

# heatplot with risky asset allocation (worker)
allo_sm <- ggplot(data = data_sm, aes(x=age_vec,y=balance_vec))+
  geom_tile(aes(fill = allocation_vec))+ 
  scale_fill_gradient2(low =  "darkred",mid = "white", high = "darkblue",name = "",midpoint = mean(rangeAllo),limits = c(lowAllo,upAllo),labels = percent(round(seq(lowAllo,upAllo,(upAllo-lowAllo)/4),2)),breaks = seq(lowAllo,upAllo,(upAllo-lowAllo)/4),guide = "colorbar") + 
  xlab("Age") + 
  ylab("Balance ($'000)") + 
  scale_x_discrete(limits=seq(25,65,by=5)) + 
  scale_y_continuous(limits = c(0,2000)) +
  theme(legend.position="bottom") + 
  theme(text=element_text(size=9), axis.text = element_text(size=9))

jpeg("heatmap_S_pi_v5.1.1.jpeg", height = 8, width = 10, units = 'cm', res = 300)
allo_sm
dev.off()


# (2) benchmark's case [optimisation based on VFI - case B]
bal_seq_w_b <- bal_seq  

con_vec <- NULL
pi_vec <- NULL
aa_temp <- NULL
consprop_vec <- NULL
allocation_vec <- NULL
age_vec <- NULL
balance_vec <- NULL

for (i in 1:sim_year_total){
  age_idx <- 25 + i - 1
  LB_C_rate <- ifelse(age_idx < age_tr, LB_C1_rate, LB_C2_rate)
  aa_temp <- BM_base[rev(which(BM_base$age == age_idx)), c(2,6,7)]
  con_vec <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_seq_w_b, Outlier, LB = LB_C_rate, UB = 1) #consumption rate
  consprop_vec <- c(consprop_vec, con_vec)
  pi_vec <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$eqprop, Bals = bal_seq_w_b, Outlier, LB = LB_pi_w, UB = UB_pi_w) #allocation rate
  allocation_vec <- c(allocation_vec, pi_vec)
  age_vec <- c(age_vec, rep(age_idx,length(con_vec)))
  balance_vec <- c(balance_vec, bal_seq_w)
}

data_bm <- cbind.data.frame(age_vec,balance_vec,consprop_vec,allocation_vec)
data_bm$balance_vec <- data_bm$balance_vec / 1000
rangeCon <- range(data_bm$consprop_vec)
lowCon <- 0.6
upCon <- 1
rangeAllo <- range(data_bm$allocation_vec)
lowAllo <- 0
upAllo <- 1.2

# heatplot with consumption % (benchmark)
con_bm <- ggplot(data = data_bm, aes(x=age_vec,y=balance_vec))+
  geom_tile(aes(fill = consprop_vec))+ 
  scale_fill_gradient2(low =  "darkred",mid = "white", high = "darkblue",name = "",midpoint = mean(rangeCon),limits = c(lowCon,upCon),labels = percent(round(seq(lowCon,upCon,(upCon-lowCon)/4),2)),breaks = seq(lowCon,upCon,(upCon-lowCon)/4),guide = "colorbar") + 
  xlab("Age") + 
  ylab("Balance ($'000)") + 
  scale_x_discrete(limits=seq(25,65,by=5)) + 
  scale_y_continuous(limits = c(0,2000)) +
  theme(legend.position="bottom") + 
  theme(text=element_text(size=9), axis.text = element_text(size=9))

jpeg("heatmap_B_con_v5.1.1.jpeg", height = 8, width = 10, units = 'cm', res = 300)
con_bm
dev.off()

# heatplot with risky asset allocation (benchmark)
allo_bm <- ggplot(data = data_bm, aes(x=age_vec,y=balance_vec))+
  geom_tile(aes(fill = allocation_vec))+ 
  scale_fill_gradient2(low =  "darkred",mid = "white", high = "darkblue",name = "",midpoint = mean(rangeAllo),limits = c(lowAllo,upAllo),labels = percent(round(seq(lowAllo,upAllo,(upAllo-lowAllo)/4),2)),breaks = seq(lowAllo,upAllo,(upAllo-lowAllo)/4),guide = "colorbar") + 
  xlab("Age") + 
  ylab("Balance ($'000)") + 
  scale_x_discrete(limits=seq(25,65,by=5)) + 
  scale_y_continuous(limits = c(0,2000)) +
  theme(legend.position="bottom") + 
  theme(text=element_text(size=9), axis.text = element_text(size=9))

jpeg("heatmap_B_pi_v5.1.1.jpeg", height = 8, width = 10, units = 'cm', res = 300)
allo_bm
dev.off()


# (3) passive manager's case [optimisation based on VFI - case H1]
bal_seq_m <- bal_seq #pre-consumption balance

con_vec <- NULL
pi_vec <- NULL
aa_temp <- NULL
consprop_vec <- NULL
allocation_vec <- NULL
age_vec <- NULL
balance_vec <- NULL

for (i in 1:sim_year_total){
  age_idx <- 25 + i - 1
  LB_C_rate <- ifelse(age_idx < age_tr, LB_C1_rate, LB_C2_rate)
  aa_temp <- HM_base1[rev(which(HM_base1$age == age_idx)), c(2,4,6)]
  con_vec <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_seq_m, Outlier, LB = LB_C_rate, UB = 1) #consumption rate
  consprop_vec <- c(consprop_vec, con_vec)
  pi_vec <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$wA, Bals = bal_seq_m, Outlier, LB = LB_pi_m, UB = UB_pi_m) #allocation rate
  allocation_vec <- c(allocation_vec, pi_vec)
  age_vec <- c(age_vec, rep(age_idx,length(con_vec)))
  balance_vec <- c(balance_vec, bal_seq_w)
}

data_hm1 <- cbind.data.frame(age_vec,balance_vec,consprop_vec,allocation_vec)
data_hm1$balance_vec <- data_hm1$balance_vec / 1000
rangeCon <- range(data_hm1$consprop_vec)
lowCon <- 0.6
upCon <- 1
rangeAllo <- range(data_hm1$allocation_vec)
lowAllo <- 0
upAllo <- 1.2

# heatplot with consumption % (manager)
con_H1 <- ggplot(data = data_hm1, aes(x=age_vec,y=balance_vec))+
  geom_tile(aes(fill = consprop_vec))+ 
  scale_fill_gradient2(low =  "darkred",mid = "white", high = "darkblue",name = "",midpoint = mean(rangeCon),limits = c(lowCon,upCon),labels = percent(round(seq(lowCon,upCon,(upCon-lowCon)/4),2)),breaks = seq(lowCon,upCon,(upCon-lowCon)/4),guide = "colorbar") + 
  xlab("Age") + 
  ylab("Balance ($'000)") + 
  scale_x_discrete(limits=seq(25,65,by=5)) + 
  scale_y_continuous(limits = c(0,2000)) +
  theme(legend.position="bottom") + 
  theme(text=element_text(size=9), axis.text = element_text(size=9))

jpeg("heatmap_H1_con_v5.1.1.jpeg", height = 8, width = 10, units = 'cm', res = 300)
con_H1
dev.off()

# heatplot with risky asset allocation (manager)
allo_H1 <- ggplot(data = data_hm1, aes(x=age_vec,y=balance_vec))+
  geom_tile(aes(fill = allocation_vec))+ 
  scale_fill_gradient2(low =  "darkred",mid = "white", high = "darkblue",name = "",midpoint = mean(rangeAllo),limits = c(lowAllo,upAllo),labels = percent(round(seq(lowAllo,upAllo,(upAllo-lowAllo)/4),2)),breaks = seq(lowAllo,upAllo,(upAllo-lowAllo)/4),guide = "colorbar") + 
  xlab("Age") + 
  ylab("Balance ($'000)") + 
  scale_x_discrete(limits=seq(25,65,by=5)) + 
  scale_y_continuous(limits = c(0,2000)) +
  theme(legend.position="bottom") + 
  theme(text=element_text(size=9), axis.text = element_text(size=9))

jpeg("heatmap_H1_pi_v5.1.1.jpeg", height = 8, width = 10, units = 'cm', res = 300)
allo_H1
dev.off()

# (4) active manager's case [optimisation based on VFI - case H2]
bal_seq_m <- bal_seq #pre-consumption balance

con_vec <- NULL
pi_vec <- NULL
aa_temp <- NULL
consprop_vec <- NULL
allocation_vec <- NULL
age_vec <- NULL
balance_vec <- NULL

for (i in 1:sim_year_total){
  age_idx <- 25 + i - 1
  LB_C_rate <- ifelse(age_idx < age_tr, LB_C1_rate, LB_C2_rate)
  aa_temp <- HM_base2[rev(which(HM_base2$age == age_idx)), c(2,4,6)]
  con_vec <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$consprop.par1, Bals = bal_seq_m, Outlier, LB = LB_C_rate, UB = 1) #consumption rate
  consprop_vec <- c(consprop_vec, con_vec)
  pi_vec <- Interpol.Fun(Balance.Data = aa_temp$balance, Allocation.Data = aa_temp$wA, Bals = bal_seq_m, Outlier, LB = LB_pi_w, UB = UB_pi_w) #allocation rate
  allocation_vec <- c(allocation_vec, pi_vec)
  age_vec <- c(age_vec, rep(age_idx,length(con_vec)))
  balance_vec <- c(balance_vec, bal_seq_w)
}

data_hm2 <- cbind.data.frame(age_vec,balance_vec,consprop_vec,allocation_vec)
data_hm2$balance_vec <- data_hm2$balance_vec / 1000
rangeCon <- range(data_hm2$consprop_vec)
lowCon <- 0.6
upCon <- 1
rangeAllo <- range(data_hm2$allocation_vec)
lowAllo <- 0
upAllo <- 1.2

# heatplot with consumption % (manager)
con_H2 <- ggplot(data = data_hm2, aes(x=age_vec,y=balance_vec))+
  geom_tile(aes(fill = consprop_vec))+ 
  scale_fill_gradient2(low =  "darkred",mid = "white", high = "darkblue",name = "",midpoint = mean(rangeCon),limits = c(lowCon,upCon),labels = percent(round(seq(lowCon,upCon,(upCon-lowCon)/4),2)),breaks = seq(lowCon,upCon,(upCon-lowCon)/4),guide = "colorbar") + 
  xlab("Age") + 
  ylab("Balance ($'000)") + 
  scale_x_discrete(limits=seq(25,65,by=5)) + 
  scale_y_continuous(limits = c(0,2000)) +
  theme(legend.position="bottom") + 
  theme(text=element_text(size=9), axis.text = element_text(size=9))

jpeg("heatmap_H2_con_v5.1.1.jpeg", height = 8, width = 10, units = 'cm', res = 300)
con_H2
dev.off()

# heatplot with risky asset allocation (manager)
allo_H2 <- ggplot(data = data_hm2, aes(x=age_vec,y=balance_vec))+
  geom_tile(aes(fill = allocation_vec))+ 
  scale_fill_gradient2(low =  "darkred",mid = "white", high = "darkblue",name = "",midpoint = mean(rangeAllo),limits = c(lowAllo,upAllo),labels = percent(round(seq(lowAllo,upAllo,(upAllo-lowAllo)/4),2)),breaks = seq(lowAllo,upAllo,(upAllo-lowAllo)/4),guide = "colorbar") + 
  xlab("Age") + 
  ylab("Balance ($'000)") + 
  scale_x_discrete(limits=seq(25,65,by=5)) + 
  scale_y_continuous(limits = c(0,2000)) +
  theme(legend.position="bottom") + 
  theme(text=element_text(size=9), axis.text = element_text(size=9))

jpeg("heatmap_H2_pi_v5.1.1.jpeg", height = 8, width = 10, units = 'cm', res = 300)
allo_H2
dev.off()

#generate a combined plot with 2 by 2 setting
library(ggpubr)

jpeg("heatmap_pi_v5.1.1.jpeg", height = 20, width = 28, units = 'cm', res = 500)
ggarrange(allo_sm, allo_bm, allo_H1,allo_H2, ncol = 2, nrow = 2,  common.legend = TRUE, legend = "right")
dev.off()

jpeg("heatmap_cons_v5.1.1.jpeg", height = 20, width = 28, units = 'cm', res = 500)
ggarrange(con_sm, con_bm, con_H1,con_H2, ncol = 2, nrow = 2,  common.legend = TRUE, legend = "right")
dev.off()
