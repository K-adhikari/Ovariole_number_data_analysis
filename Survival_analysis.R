library(survival)
library(ggplot2)
library(survminer)
library(splitstackshape)
library (ggpubr)
library(data.table)
library(coxme)
library(emmeans)
library(multcomp)



# Read the data and assign 'high' or 'low' categories to each genotype based on their ovariole number

setwd("~/Documents/Data/Ovariole_number/Results/Survivorship")
data<- read.csv("Survivorship_data.csv", header = T)

data$ovariole_number<- ifelse(data$Strain=="RAL370", "High",
                              ifelse(data$Strain=="RAL129", "High",
                                     ifelse(data$Strain=="RAL486","High",
                                            ifelse(data$Strain=="RAL737", "High",
                                                   ifelse(data$Strain=="RAL799", "High",
                                                          ifelse(data$Strain=="RAL443", "High",
                                                                 ifelse(data$Strain=="Sham", "None", "Low")))))))




# Prepare data to be used for survival package

data$Strain<- as.factor(data$Strain)
DT<- as.data.table(data)
DT1<- expandRows(DT, "Flies")
DT2<- as.data.frame(DT1) 
DT2$Combo<- paste(DT2$Strain, DT2$Mating_status)



# Fit survival curves

# Here Time means the observation time
# Censor: 1 - died ; 0 - Survived at the end of the study

km_fit <- survfit(Surv(Time, Censor) ~ Combo, data= DT2)

survp<- ggsurvplot(km_fit, data = DT2, risk.table = TRUE, legend="right", risk.table.col="strata",
                   risk.table.y.text = FALSE,
                   risk.table.height=0.35, xlim=c(0,120), break.x.by=24)
survp$plot
survp$table



# Saving the survival plot

pdf("survplot_all_updated.pdf", width = 8, height = 8)
print(survp$plot, newpage = FALSE)
dev.off()



# Fitting a cox mixed effect model

fit_m<- coxme(formula= Surv(Time, Censor) ~ Strain + (1|Block), data=DT1)
summary(fit_m)



# Compare survival between genotypes

A<- emmeans(fit_m, pairwise~Strain)
A
A_means<- emmeans(fit_m, "Strain")
A_means_cld<- cld(A_means, Letters= letters, alpha = 0.05)
A_means_cld
