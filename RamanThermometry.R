library(tidyverse)
library(dplyr)
library(readr)
coef <- 1.43988

## check your working directory and change if needed
getwd()

#setwd()

## what is the expected temperature?
temperature <- 5

## what is the file name? Then read in the data and fill file variable

filename <- paste0("data_1.dat")

col_names <- names(read_delim(filename,"\t", n_max = 0))
data <- read_delim(filename,"\t", col_names = TRUE, skip = 2)


col_names <- c("W", "I")
data <- data[1:2] %>% setNames(col_names) 


plot(data$W, data$I, ylab = "Intensity", xlab = "wavenumber")


## we only take data into account above the Rayleigh peak and below a treshold T_cutoff: 3 * k_B * T:

Tcut <- 3* coef * temperature

ind1 <- which(data[1] >= -Tcut)[1]
ind2 <- which(data[1] >= -7)[1]
ind3 <- which(data[1] >= 5)[1]
ind4 <- ifelse(is.na(which(data[1] >= Tcut)[1]), which(data[1] == max(data[1])), which(data[1] >= Tcut)[1])

## we use the frequency axis of the Anti Stokes response also for the Stokes part  

dat <- data.frame("W" = abs(data[ind2:ind1,1]), "I_AS"= data[ind2:ind1,2], "I_S" = approx(data$W[ind3:ind4], data$I[ind3:ind4], abs(data$W[ind2:ind1]), method = "linear")[2]) 

head(dat)
col_names2 <- c("w1", "I_AS", "I_S")

calc1 <-dat %>% na.omit() %>% setNames(col_names2)

filename3 = paste0("IS_IAS_", temperature, "K.pdf")
p <- calc1 %>% ggplot() + geom_line(aes(x = w1, y = I_AS), color = 'blue') + geom_line(aes(x = w1, y = I_S), color = 'red') + ylab("Raman Intensity") + xlab("wave number (cm$^{-1}$)") + ylim(0,5500) + xlim(0,Tcut+10)
p + theme_bw()
ggsave(filename3, dpi = 300)


## the background is the minimal value of the Anti-Stokes part - 1 
background <- min(calc1$I_AS)
B <- background-1

## now find the temperature
ln <- log((calc1$I_S - B)/(calc1$I_AS-B))
temp <- data.frame("W" =calc1$w1,  "temp"  = coef * calc1$w1 * 1/ln)


T_uncorrect <- temp %>% summarize(mean = mean(temp), sd = sd(temp))

p <- temp %>% ggplot() + geom_point(aes(x = W, y = temp), color = 'blue') + ylab("Temperature (Kelvin)") + xlab("wave number (cm$^{-1}$)") + ylim(0,temperature + 10) + xlim(0,Tcut+10)
p + theme_bw()
filename4 = paste0("temperature_wavenumber_", temperature, "K.pdf")
ggsave(filename4, dpi = 300)


## weighing function based on statical error propagation

er <- ln*ln/(coef* abs(calc1$w1) * sqrt((calc1$I_S/(B-calc1$I_S)^2 + calc1$I_AS/(calc1$I_AS-B)^2)))
Z <- 1/sum(er)
WeightT <- data.frame("W" =calc1$w1, "temp"  = temp$temp, "error" = er,  "P"  = Z * er)
WeightT <- WeightT %>% mutate("PT"  = P*temp)
T_av <- sum(PT)

p <- WeightT %>% ggplot(aes(x = W, y = temp))+ geom_point(color = 'blue')+ geom_errorbar(aes(ymin = temp-error/2, ymax = temp+error/2), color = 'blue') + ylab("Temperature (Kelvin)") + xlab("wave number (cm$^{-1}$)") + ylim(0,temperature + 10) + xlim(0,Tcut+10)
p + theme_bw() + geom_text(x=15, y=2, label=paste("Average weighted temperature =", round(T_av,digits=2), "K"))
filename5 = paste0("temperature_wavenumber_errorbars", temperature, "K.pdf")
ggsave(filename5, dpi = 300)


p <- WeightT %>% ggplot(aes(x = W, y = P))+ geom_point(color = 'blue')+ ylab("Weight") + xlab("wave number (cm$^{-1}$)") + ylim(0,0.1) + xlim(0,Tcut+10)
p + theme_bw()
filename6 = paste0("Weights_wavenumber_errorbars", temperature, "K.pdf")
ggsave(filename6, dpi = 300)


dfT <- data.frame("temp_analytical" = T_av, "temp_measured" = temperature)

filename2 = paste0("T_analytical_all.txt")
write.table(dfT, file = filename2, row.names = FALSE, col.names = TRUE)


