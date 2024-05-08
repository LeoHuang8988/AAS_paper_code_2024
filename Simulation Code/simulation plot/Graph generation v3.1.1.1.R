## Construct ggplots
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load projection data (median results at initial balance of $0)
load("worker_med_projection_bal_0.RData");  load("benchmark_med_projection_bal_0.RData"); load("manager1_med_projection_bal_0.RData"); load("manager2_med_projection_bal_0.RData")
med_decisions_s <- as.data.frame(median_decisions_s);  med_decisions_b <- as.data.frame(median_decisions_b)
med_decisions_h1 <- as.data.frame(median_decisions_h1); med_decisions_h2 <- as.data.frame(median_decisions_h2)
library(ggplot2)

# Test loading data
#(a) worker's case
plot(med_decisions_s$Age[1:40],med_decisions_s$Med_E[1:40]) #Risky asset allocation
plot(med_decisions_s$Age[1:40],med_decisions_s$Med_Bal[2:41]) #Fund balance 
plot(med_decisions_s$Age[1:40],med_decisions_s$Med_pretax_cons[1:40]) #Pre-tax consumption
plot(med_decisions_s$Age[1:40],med_decisions_s$Med_posttax_cons[1:40]) #Post-tax consumption
#(b) passive manager's case
plot(med_decisions_h1$Age[1:40],med_decisions_h1$Med_E[1:40])
plot(med_decisions_h1$Age[1:40],med_decisions_h1$Med_Bal[2:41])
plot(med_decisions_h1$Age[1:40],med_decisions_h1$Med_pretax_cons[1:40])
plot(med_decisions_h1$Age[1:40],med_decisions_h1$Med_posttax_cons[1:40])
#(c) benchmark's case
plot(med_decisions_b$Age[1:40],med_decisions_b$Med_E[1:40])
plot(med_decisions_b$Age[1:40],med_decisions_b$Med_Bal[2:41])
plot(med_decisions_b$Age[1:40],med_decisions_b$Med_pretax_cons[1:40])
plot(med_decisions_b$Age[1:40],med_decisions_b$Med_posttax_cons[1:40])
#(b) active manager's case
plot(med_decisions_h2$Age[1:40],med_decisions_h2$Med_E[1:40])
plot(med_decisions_h2$Age[1:40],med_decisions_h2$Med_Bal[2:41])
plot(med_decisions_h2$Age[1:40],med_decisions_h2$Med_pretax_cons[1:40])
plot(med_decisions_h2$Age[1:40],med_decisions_h2$Med_posttax_cons[1:40])

###plot set 1 - median (separate into 4 different variables)###

## Case 1 -- Balance 
# The data have a common independent variable (x)
x <- med_decisions_s$Age[1:40]

# Generate 4 different sets of outputs
y1 <- med_decisions_s$Med_Bal[2:41] / 1000
y2 <- med_decisions_b$Med_Bal[2:41] / 1000
y3 <- med_decisions_h1$Med_Bal[2:41] / 1000
y4 <- med_decisions_h2$Med_Bal[2:41] / 1000
y <- list(y1, y2, y3, y4)

# Make a ggplot
data1 <- data.frame(Cases=rep(c("S", "B", "H1", "H2"), each=40),
                    Age=rep(x,4),
                    Balance=c(y1,y2,y3,y4))
data1$Age <- as.numeric(as.vector(data1$Age))

jpeg("Med_Simu_bal_v5.1.1.1.jpeg", height = 8, width = 11, units = 'cm', res = 300)

data1$Cases <- factor(data1$Cases, levels=c("S", "B", "H1", "H2"))

med_bal = ggplot(data=data1, aes(x=Age, y=Balance, group=Cases, linetype=Cases, color=Cases)) + 
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_linetype_manual(values=c("B" = "solid", "S" = "twodash", "H1" = "dotdash", "H2" = "dotted"),
                        labels=c("S", "B", expression(H[1]), expression(H[2]))) +  # Manually set linetypes
  scale_color_manual(values=c("B" = "blue", "S" = "blue", "H1" = "red", "H2" = "red"),
                     labels=c("S", "B", expression(H[1]), expression(H[2]))) +  # Manually set colors
  geom_line(size=0.8) + 
  theme(legend.position="bottom") + 
  labs(y= "Balance ($'000)", x = "Age") + 
  scale_x_discrete(limits=seq(25,65,by=5)) + 
  ylim(0, 1500) + 
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

med_bal

dev.off()


## Case 2 -- Asset Allocation
# The data have a common independent variable (x)
x <- med_decisions_s$Age[1:40]

# Generate 4 different sets of outputs
y1 <- med_decisions_s$Med_E[1:40] * 100
y2 <- med_decisions_b$Med_E[1:40] * 100
y3 <- med_decisions_h1$Med_E[1:40] * 100
y4 <- med_decisions_h2$Med_E[1:40] * 100
y <- list(y1, y2, y3, y4)

# Make a ggplot
data2 <- data.frame(Cases=rep(c("S", "B", "H1", "H2"), each=40),
                    Age=rep(x,4),
                    Allocation=c(y1,y2,y3,y4))
data2$Age <- as.numeric(as.vector(data2$Age))

jpeg("Med_Simu_allo_v5.1.1.1.jpeg", height = 8, width = 11, units = 'cm', res = 300)

data2$Cases <- factor(data2$Cases, levels=c("S", "B", "H1", "H2"))

med_allo = ggplot(data=data2, aes(x=Age, y=Allocation, group=Cases, linetype=Cases, color=Cases)) + 
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))+  
  scale_linetype_manual(values=c("B" = "solid", "S" = "twodash", "H1" = "dotdash", "H2" = "dotted"),
                        labels=c("S", "B", expression(H[1]), expression(H[2]))) +  # Manually set linetypes
  scale_color_manual(values=c("B" = "blue", "S" = "blue", "H1" = "red", "H2" = "red"),
                     labels=c("S", "B", expression(H[1]), expression(H[2]))) +  # Manually set colors
  geom_line(size=0.8) + 
  theme(legend.position="bottom") + 
  labs(y= "Allocation to Risky Asset (%)", x = "Age") + 
  scale_x_discrete(limits=seq(25,70,by=5)) + 
  ylim(0, 125)+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

med_allo

dev.off()

## Case 3 -- Pre-tax Consumption
# The data have a common independent variable (x)
x <- med_decisions_s$Age[1:40]

# Generate 4 different sets of outputs
y1 <- med_decisions_s$Med_pretax_cons[1:40] / 1000
y2 <- med_decisions_b$Med_pretax_cons[1:40] / 1000
y3 <- med_decisions_h1$Med_pretax_cons[1:40] / 1000
y4 <- med_decisions_h2$Med_pretax_cons[1:40] / 1000
y <- list(y1, y2, y3, y4)

# Make a ggplot
data3 <- data.frame(Cases=rep(c("S", "B", "H1", "H2"), each=40),
                    Age=rep(x,4),
                    Consumption=c(y1,y2,y3,y4))
data3$Age <- as.numeric(as.vector(data3$Age))

jpeg("Med_Simu_Pretax_Consump_v5.1.1.1.jpeg", height = 8, width = 11, units = 'cm', res = 300)

data3$Cases <- factor(data3$Cases, levels=c("S", "B", "H1", "H2"))

med_precon = ggplot(data=data3, aes(x=Age, y=Consumption, group=Cases, linetype=Cases, color=Cases)) + 
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_linetype_manual(values=c("B" = "solid", "S" = "twodash", "H1" = "dotdash", "H2" = "dotted"),
                        labels=c("S", "B", expression(H[1]), expression(H[2]))) +  # Manually set linetypes
  scale_color_manual(values=c("B" = "blue", "S" = "blue", "H1" = "red", "H2" = "red"),
                     labels=c("S", "B", expression(H[1]), expression(H[2]))) +  # Manually set colors
  geom_line(size=0.8) + 
  theme(legend.position="bottom") + 
  labs(y= "Pretax Consumption ($'000)", x = "Age") + 
  scale_x_discrete(limits=seq(25,70,by=5)) + 
  ylim(40, 70)+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

med_precon

dev.off()


## Case 4 -- Post-tax Consumption
# The data have a common independent variable (x)
x <- med_decisions_s$Age[1:40]

# Generate 4 different sets of outputs
y1 <- med_decisions_s$Med_posttax_cons[1:40] / 1000
y2 <- med_decisions_b$Med_posttax_cons[1:40] / 1000
y3 <- med_decisions_h1$Med_posttax_cons[1:40] / 1000
y4 <- med_decisions_h2$Med_posttax_cons[1:40] / 1000
y <- list(y1, y2, y3, y4)

# Make a ggplot 
data4 <- data.frame(Cases=rep(c("S", "B", "H1", "H2"), each=40),
                    Age=rep(x,4),
                    Consumption=c(y1,y2,y3,y4))
data4$Age <- as.numeric(as.vector(data4$Age))

jpeg("Med_Simu_Posttax_Consump_v5.1.1.1.jpeg", height = 8, width = 11, units = 'cm', res = 300)

data4$Cases <- factor(data4$Cases, levels=c("S", "B", "H1", "H2"))

med_postcon <- ggplot(data=data4, aes(x=Age, y=Consumption, group=Cases, linetype=Cases, color=Cases)) + theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_line(size=0.8) +
  scale_linetype_manual(values=c("B" = "solid", "S" = "twodash", "H1" = "dotdash", "H2" = "dotted"),
                        labels=c("S", "B", expression(H[1]), expression(H[2]))) +  # Manually set linetypes
  scale_color_manual(values=c("B" = "blue", "S" = "blue", "H1" = "red", "H2" = "red"),
                     labels=c("S", "B", expression(H[1]), expression(H[2]))) +  # Manually set colors
  theme(legend.position="bottom") +
  labs(y= "Posttax Consumption ($'000)", x = "Age") +
  scale_x_discrete(limits=seq(25,70,by=5)) +
  ylim(40, 60) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

med_postcon


dev.off()
