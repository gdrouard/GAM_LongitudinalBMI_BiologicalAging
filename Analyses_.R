# data object stores BMI values (long format)
data <- read.csv("XXX")
data$PERSON_NB <- as.character(data$PERSON_NB)

#proteomic clock: PAC & HPS
PAC <- read.csv("XXX/PAC.txt")
rownames(PAC) <- PAC$person_id
PAC$person_id <- as.character(PAC$person_id)

HPS <- read.csv("XXX/HPS.txt")
rownames(HPS) <- HPS$person_id
HPS$person_id <- as.character(HPS$person_id)

#ProtAge
ProtAgeDat <- read.csv("XXX/Finn_Twin_ProtAge_Dec_02_2024.csv")
KEY_IDs <- read.csv("~XXX/KEY_IDs.csv")
rownames(ProtAgeDat) <- ProtAgeDat$A_IDs
rownames(KEY_IDs) <- KEY_IDs$A_IDs
rownames(ProtAgeDat) <- KEY_IDs[rownames(ProtAgeDat),"personid"] #math IDs

## keep only overlaping individuals
inter <- intersect(data$PERSON_NB,rownames(PAC))
data <- data[data$PERSON_NB %in% inter,]
PAC <- PAC[inter,]
HPS <- HPS[inter,]
ProtAgeDat <- ProtAgeDat[inter,]


#### plot data

head(data)

##################################################
# Plot the individual trajectories

require(ggplot2)

ggplot(data, aes(x = AGE, y = BMI)) +
  geom_point() +
  geom_line(aes(group = PERSON_NB), alpha = 0.7) +
  labs(title = "Individual Trajectories of BMI", 
       x = "Age of participants", 
       y = "Body mass index") +
  theme_minimal() +
  theme(legend.position = "none")

#ggsave("BMItrajectories.png", plot = last_plot(), dpi = 600, width = 18, height = 14, units = "cm")


############################################################
##### plot measurement frequencies

require(dplyr)
# Step 1: Calculate the number of measurements per individual
measurement_counts <- data %>%
  group_by(PERSON_NB) %>%
  summarise(measurement_count = n())

# Step 2: Create a frequency table
frequency_table <- measurement_counts %>%
  group_by(measurement_count) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

# Step 3: Plot the frequency table with ggplot2
ggplot(frequency_table, aes(x = measurement_count, y = percentage)) +
  geom_bar(stat = "identity", fill = "black") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            vjust = -0.1, size = 3.8,hjust = -0.5) +
  labs(x = "Number of Measurements per Individual", 
       y = "Percentage of Individuals") +
  theme_minimal() + coord_flip() + scale_x_continuous(breaks = 1:9) +
  ylim(0,85) + theme(axis.text.y = element_text(face = "bold"),
                     axis.text.x = element_text(face = "bold"))

#ggsave("BMIfrequencies.png", plot = last_plot(), dpi = 600, width = 10, height = 18, units = "cm")


###########################################
### calculate changes in BMI


library(lme4)  # For linear mixed-effects models
library(lmerTest)  # For p-values in linear mixed-effects models
library(broom)  # For tidying model outputs

data$decades <- data$AGE/10 - 1.8 # decades since the individuals turned 18 (decades works better than years for model fitting)

# Fit the linear mixed-effects model with linear and quadratic terms
model_quad <- lmer(BMI ~ decades + I(decades^2) + (decades + I(decades^2) | PERSON_NB), data = data)

# Display the summary of the model
summary(model_quad)

# Fit the linear mixed-effects model without the quadratic term
model_linear <- lmer(BMI ~ decades + (decades | PERSON_NB), data = data)

# Compare the two models using ANOVA
anova(model_linear, model_quad)

# compare R2 of models
library(performance)

# R² for the linear model
r2_linear <- r2(model_linear)
print(r2_linear)

r2_nonlinear <- r2(model_quad)
print(r2_nonlinear)


# Extract individual-level coefficients (random effects)
individual_coeffs.NL <- coef(model_quad)$PERSON_NB
individual_coeffs.L <- coef(model_linear)$PERSON_NB

# View the first few rows of the coefficients
head(individual_coeffs.NL)

head(individual_coeffs.L)

# Load ggplot2 for plotting linear model growth factors
library(ggplot2)

# Plot the intercepts
ggplot(individual_coeffs.L, aes(x = `(Intercept)`)) +
  geom_histogram(binwidth = 0.5, fill = "blue", alpha = 0.7) +
  labs(title = "Distribution of Intercepts", x = "Intercept", y = "Count")

# Plot the linear slopes
ggplot(individual_coeffs.L, aes(x = `decades`)) +
  geom_histogram(binwidth = 0.01, fill = "green", alpha = 0.7) +
  labs(title = "Distribution of Linear Slopes", x = "Linear Slope", y = "Count")


###### modeling with EAA

# Set seed for reproducibility
set.seed(1871)

#linear model
dim(individual_coeffs.L)
colnames(individual_coeffs.L) <- c("Intercept","LinearSlope")

### we also add BMI at blood sampling
library(readxl)
old_bp_12042021 <- read_excel("XXX/old_bp_12042021.xlsx") # Covariate data
old_bp_12042021 <- old_bp_12042021[!duplicated(old_bp_12042021$PERSON_NB),]
old_bp_12042021 <- as.data.frame(old_bp_12042021)
rownames(old_bp_12042021) <- old_bp_12042021$PERSON_NB
# Calculate BMI at BS from weight/height
individual_coeffs.L$BMI_BS <- 
  old_bp_12042021[rownames(individual_coeffs.L),"PAI_79"]/((old_bp_12042021[rownames(individual_coeffs.L),"PIT_78"]/100)**2)


#### Add clock values

individual_coeffs.L$PAC <- PAC[rownames(individual_coeffs.L),"PAC"]
individual_coeffs.L$EAA_PAC <- PAC[rownames(individual_coeffs.L),"EAA_PAC"]

individual_coeffs.L$HPS <- HPS[rownames(individual_coeffs.L),"HPS"]
individual_coeffs.L$EAA_HPS <- HPS[rownames(individual_coeffs.L),"EAA_HPS"]

ProtAgeDat$ProtAge_res <- scale(residuals(lm(ProtAgeDat$ProtAge~HPS$ChronoAge)))
individual_coeffs.L$EAA_ProtAge <- ProtAgeDat[rownames(individual_coeffs.L),"ProtAge_res"]
individual_coeffs.L$EAA_ProtAge <- ProtAgeDat[rownames(individual_coeffs.L),"ProtAge_res"]


#epigenetic clocks

EpiAge_EAA_VP_epi <- read.csv("XXX/EpiAge_EAA_VP_epi.csv") #Epigenetic clocks
EpiAge_EAA_VP_epi <- as.data.frame(EpiAge_EAA_VP_epi)
rownames(EpiAge_EAA_VP_epi) <- EpiAge_EAA_VP_epi$PersonID

EpiClocks <- EpiAge_EAA_VP_epi

individual_coeffs.L$EAA_Horvath <- EpiClocks[rownames(individual_coeffs.L),"AgeAccel_PCHorvath1"]

individual_coeffs.L$EAA_Dunedin <- EpiClocks[rownames(individual_coeffs.L),"DunedinPACE"]

individual_coeffs.L$EAA_Hannum <- EpiClocks[rownames(individual_coeffs.L),"AgeAccel_PCHannum"]

individual_coeffs.L$EAA_GrimAge <- EpiClocks[rownames(individual_coeffs.L),"AgeAccel_PCGrimAge"]

individual_coeffs.L$EAA_PhenoAge <- EpiClocks[rownames(individual_coeffs.L),"AgeAccel_PCPhenoAge"]

individual_coeffs.L$EAA_GrimAge2 <- EpiClocks[rownames(individual_coeffs.L),"AgeAccelGrim2"]

### add adipose age
OrganClocks_EHEPI_1stGen <- read.csv("XXX/OrganClocks_EHEPI.csv", row.names=1)
individual_coeffs.L$Adipose_unadjusted <- OrganClocks_EHEPI_1stGen[rownames(individual_coeffs.L),"Adipose"]
individual_coeffs.L$Conventional_unadjusted <- OrganClocks_EHEPI_1stGen[rownames(individual_coeffs.L),"Conventional"]
table(rownames(PAC)==rownames(individual_coeffs.L))
individual_coeffs.L$Adipose <- residuals(lm(individual_coeffs.L$Adipose_unadjusted~PAC$ChronoAge))
individual_coeffs.L$Conventional <- residuals(lm(individual_coeffs.L$Conventional_unadjusted~PAC$ChronoAge))


# rescale all clocks with rank based transformation

rint <- function(x){return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))}

colnames(individual_coeffs.L)
for(i in 4:ncol(individual_coeffs.L)){
  individual_coeffs.L[,i] <- rint(individual_coeffs.L[,i])
}


### changes in BMI need to be in year, not decade
individual_coeffs.L$LinearSlope <- individual_coeffs.L$LinearSlope/10

# We add family number! 

individual_coeffs.L$family_nb <- old_bp_12042021[rownames(individual_coeffs.L),"FAMILY_NB"]

####################################################################
## Function to check GAM assumptions

save_gam_check_plots <- function(m, predictor_name, clock_name) {
  # Create the filename
  filename <- paste(predictor_name, clock_name, ".png", sep = "_")
  
  # Open a PNG device
  png(filename, width = 800, height = 800)
  
  # Set up the plotting area to 2x2 for four plots
  par(mfrow = c(2, 2))
  
  # Generate the gam.check plots
  gam.check(m)
  
  # Close the graphics device to save the file
  dev.off()
}

# Example usage:
# save_gam_check_plots(m, "Baseline", "C1")

setwd("XXX\\GAM assumptions plots")


######################################################################
### analyses with GAM
######################################################################

# we create two datasets (Table 2/Table 3) that we gonna fill in
sumstatsT2 <- as.data.frame(matrix(NA,nrow=3*9,ncol=8))
sumstatsT3 <- as.data.frame(matrix(NA,nrow=2*9,ncol=10))

colnames(sumstatsT2) <- c("Investigated","Aging clock","edf","F","p","Dev. explained","R2","Delta_AIC")
colnames(sumstatsT3) <- c("Investigated","Aging clock","edf","F","p","F_ti","p_ti","Dev. explained","R2","Delta_AIC")

library(mgcv)

###########################################################
# intercept only

# PAC
m <- bam(EAA_PAC ~ s(Intercept) + s(family_nb,bs="re",k = 1), 
         family = gaussian, data = individual_coeffs.L, method = "REML")
summary(m)
#plot(m)
m.linear <- bam(EAA_PAC ~ Intercept + s(family_nb,bs="re",k=1), 
                family = gaussian, data = individual_coeffs.L, method = "REML")
anova(m,m.linear)
temp <- summary(m)
sumstatsT2[1,] <- c("Intercept","PAC",temp$s.table[1,1],temp$s.table[1,3],temp$s.table[1,4],
                    temp$dev.expl,temp$r.sq,AIC(m)-AIC(m.linear))
save_gam_check_plots(m, "Baseline", "PAC")

# repeat this script with other clocks to fill in lines 1:9 of sumstatsT2


###########################################################
# Slope only

### PAC
m <- bam(EAA_PAC ~ s(LinearSlope) + Intercept + s(family_nb, bs = "re", k = 1), 
         family = gaussian, data = individual_coeffs.L, method = "REML")
summary(m)

m.linear <- bam(EAA_PAC ~ LinearSlope + Intercept + s(family_nb, bs = "re", k = 1), 
                family = gaussian, data = individual_coeffs.L, method = "REML")

m_no_intercept <- bam(EAA_PAC ~ s(LinearSlope) + s(family_nb, bs = "re", k = 1), 
                      family = gaussian, data = individual_coeffs.L, method = "REML")

null_model <- bam(EAA_PAC ~ s(family_nb, bs = "re", k = 1), family = gaussian, data = individual_coeffs.L, method = "REML")

dev_expl_no_intercept <- (deviance(null_model) - deviance(m_no_intercept)) / deviance(null_model)
r2_no_intercept <- 1 - (deviance(m_no_intercept) / deviance(null_model))

temp <- summary(m)

sumstatsT2[10, ] <- c(
  "Slope", "PAC",temp$s.table[1, 1], temp$s.table[1, 3],temp$s.table[1, 4],  
  dev_expl_no_intercept, r2_no_intercept, AIC(m) - AIC(m.linear)  
)

save_gam_check_plots(m, "Change", "PAC")

# repeat this script with other clocks to fill in lines 9 additional rows in sumstatsT2


###########################################################
# BMI at blood sampling

# PAC
m <- bam(EAA_PAC ~ s(BMI_BS) + s(family_nb,bs="re",k = 1), 
         family = gaussian, data = individual_coeffs.L, method = "REML")
summary(m)
#plot(m)
m.linear <- bam(EAA_PAC ~ BMI_BS + s(family_nb,bs="re",k=1), 
                family = gaussian, data = individual_coeffs.L, method = "REML")
anova(m,m.linear)
temp <- summary(m)
sumstatsT2[19,] <- c("BMI_BS","PAC",temp$s.table[1,1],temp$s.table[1,3],temp$s.table[1,4],
                     temp$dev.expl,temp$r.sq,AIC(m)-AIC(m.linear))
save_gam_check_plots(m, "BMI", "PAC")

# repeat this script with other clocks to fill in lines 9 additional rows in sumstatsT2

####################################################################
##### Interaction models: Intercept x Slope

### Repeat scripts below to fill in sumsStats3 table, changing the clocks

# PAC
m <- bam(EAA_PAC ~ te(Intercept,LinearSlope) + s(family_nb,bs="re",k = 1), 
         family = gaussian, data = individual_coeffs.L, method = "REML")
summary(m)

m_ti <- bam(EAA_PAC ~ s(Intercept) + s(LinearSlope) + ti(Intercept,LinearSlope) + s(family_nb,bs="re",k = 1), 
            family = gaussian, data = individual_coeffs.L, method = "REML")
summary(m_ti)

#plot(m)
vis.gam(m, view = c("Intercept", "LinearSlope"), plot.type = "persp",theta=-30,phi=-5)
vis.gam(m, view = c("Intercept", "LinearSlope"), plot.type = "contour",
        xlab="Fitted baseline BMI",ylab="Change in BMI") #4.51*4.38 dim

m.linear <- gam(EAA_PAC ~ Intercept + LinearSlope + Intercept:LinearSlope + s(family_nb,bs="re",k = 1), 
                family = gaussian, data = individual_coeffs.L, method = "REML")

print(anova(m.linear,m))
temp <- summary(m)
sumstatsT3[1,] <- c("Intercept:Slope","PAC",temp$s.table[1,1],temp$s.table[1,3],temp$s.table[1,4],
                    summary(m_ti)$s.table["ti(Intercept,LinearSlope)",c("F","p-value")],
                    temp$dev.expl,temp$r.sq,AIC(m)-AIC(m.linear))
save_gam_check_plots(m, "inter_BaselinexChange", "PAC")



####################################################################
##### Interaction models: Intercept x BMI_BS


# PAC
m <- bam(EAA_PAC ~ te(Intercept,BMI_BS) + s(family_nb,bs="re",k = 1), 
         family = gaussian, data = individual_coeffs.L, method = "REML")
summary(m)

m_ti <- bam(EAA_PAC ~ s(Intercept) + s(BMI_BS) + ti(Intercept,BMI_BS) + s(family_nb,bs="re",k = 1), 
            family = gaussian, data = individual_coeffs.L, method = "REML")
summary(m_ti)

#plot(m)
vis.gam(m, view = c("Intercept", "BMI_BS"), plot.type = "persp",theta=-30,phi=-5)
vis.gam(m, view = c("Intercept", "BMI_BS"), plot.type = "contour",
        xlab="Fitted baseline BMI",ylab="Change in BMI") #4.51*4.38 dim

m.linear <- gam(EAA_PAC ~ Intercept + BMI_BS + Intercept:BMI_BS + s(family_nb,bs="re",k = 1), 
                family = gaussian, data = individual_coeffs.L, method = "REML")

print(anova(m.linear,m))
temp <- summary(m)
sumstatsT3[10,] <- c("Intercept:BMI","PAC",temp$s.table[1,1],temp$s.table[1,3],temp$s.table[1,4],
                     summary(m_ti)$s.table["ti(Intercept,BMI_BS)",c("F","p-value")],
                     temp$dev.expl,temp$r.sq,AIC(m)-AIC(m.linear))
save_gam_check_plots(m, "inter_BaselinexBMI", "PAC")



#############################################################################

################### SAVE TABLES

setwd("XXX\\Analyses")
#write.csv(sumstatsT2,"Table2_raw.csv",row.names = F)
#write.csv(sumstatsT3,"Table3_raw.csv",row.names = F)


##########################################################################

## checking GAMS:
gam.check(m)

##########################################################
### Get better visualizations?? examples below

library(mgcViz)

m <- bam(EAA_PAC ~ s(Intercept) + s(family_nb,bs="re",k=1), 
         family = gaussian, data = individual_coeffs.L, method = "REML")
summary(m)

b <- getViz(m)

# Extract the edf from the model summary
edf_value <- summary(m)$s.table[1, "edf"]

# Plot the smooth term and modify the y-axis label
plot(b, select = 1) +
  labs(y = paste0("Smooth effect of intercept (edf=", round(edf_value, 2), ")")) +
  ggtitle("PAC")



### make 3D plots
m <- bam(EAA_Dunedin ~ te(Intercept,BMI_BS) + s(family_nb,bs="re",k = 1), 
         family = gaussian, data = individual_coeffs.L, method = "REML")
vis.gam(m, view = c("Intercept", "BMI_BS"), plot.type = "persp",
        theta=-50,phi=-5,
        xlab="Baseline BMI",ylab="BMI (last measurement)")
m <- bam(EAA_ProtAge ~ te(Intercept,BMI_BS) + s(family_nb,bs="re",k = 1), 
         family = gaussian, data = individual_coeffs.L, method = "REML")
vis.gam(m, view = c("Intercept", "BMI_BS"), plot.type = "persp",
        theta=-50,phi=-5,
        xlab="Baseline BMI",ylab="BMI (last measurement)")
m <- bam(EAA_HPS ~ te(Intercept,BMI_BS) + s(family_nb,bs="re",k = 1), 
         family = gaussian, data = individual_coeffs.L, method = "REML")
vis.gam(m, view = c("Intercept", "BMI_BS"), plot.type = "persp",
        theta=50,phi=10,
        xlab="Baseline BMI",ylab="BMI (last measurement)")


m <- bam(EAA_Dunedin ~ te(Intercept,LinearSlope) + s(family_nb,bs="re",k = 1), 
         family = gaussian, data = individual_coeffs.L, method = "REML")
vis.gam(m, view = c("Intercept", "LinearSlope"), plot.type = "persp",
        theta=-50,phi=-5,
        xlab="Baseline BMI",ylab="Change in BMI")
m <- bam(EAA_PhenoAge ~ te(Intercept,LinearSlope) + s(family_nb,bs="re",k = 1), 
         family = gaussian, data = individual_coeffs.L, method = "REML")
vis.gam(m, view = c("Intercept", "LinearSlope"), plot.type = "persp",
        theta=-50,phi=-5,
        xlab="Baseline BMI",ylab="Change in BMI")
m <- bam(Adipose ~ te(Intercept,LinearSlope) + s(family_nb,bs="re",k = 1), 
         family = gaussian, data = individual_coeffs.L, method = "REML")
vis.gam(m, view = c("Intercept", "LinearSlope"), plot.type = "persp",
        theta=-50,phi=-5,
        xlab="Baseline BMI",ylab="Change in BMI")

###################################

library(mgcv)
library(mgcViz)
library(gridExtra)
library(ggplot2)

# Create a list to store all the plots.

plots <- list()
outcomes <- c("EAA_PAC", "EAA_HPS","EAA_ProtAge", "EAA_Horvath", "EAA_Hannum", "EAA_PhenoAge", "EAA_Dunedin", "Adipose", "EAA_GrimAge2")
outcome_title <- c("PAC","HPS","ProtAge","Horvath","Hannum","PhenoAge","DunedinPACE","Adipose","GrimAge2")

# Loop over the outcomes and create models and plots for Intercept, Slope, and BMI:
for (i in seq_along(outcomes)) {
  outcome <- outcomes[i]
  title <- outcome_title[i]
  
  
  # Fit GAM models
  m_intercept <- bam(as.formula(paste(outcome, "~ s(Intercept)  + s(family_nb,bs='re',k = 1)")), family = gaussian, data = individual_coeffs.L, method = "REML")
  m_slope <- bam(as.formula(paste(outcome, "~ s(LinearSlope) + Intercept +  s(family_nb,bs='re',k = 1)")), family = gaussian, data = individual_coeffs.L, method = "REML")
  m_bmi <- bam(as.formula(paste(outcome, "~ s(BMI_BS) +  s(family_nb,bs='re',k = 1)")), family = gaussian, data = individual_coeffs.L, method = "REML")
  
  # Create visualization objects
  b_intercept <- getViz(m_intercept)
  b_slope <- getViz(m_slope)
  b_bmi <- getViz(m_bmi)
  
  # Extract edf for labeling
  edf_intercept <- round(summary(m_intercept)$s.table[1, "edf"], 2)
  edf_slope <- round(summary(m_slope)$s.table[1, "edf"], 2)
  edf_bmi <- round(summary(m_bmi)$s.table[1, "edf"], 2)
  
  # Generate and store plots with custom titles
  plots[[paste(outcome, "fitted BMI baseline", sep = "_")]] <- plot(b_intercept, select = 1) +
    labs(y = paste0("Smooth effect of fitted BMI baseline (edf=", edf_intercept, ")"),
         x= "Fitted BMI baseline",
         title = paste0(title," "))
  
  plots[[paste(outcome, "Change in BMI", sep = "_")]] <- plot(b_slope, select = 1) +
    labs(y = paste0("Smooth effect of change in BMI (edf=", edf_slope, ")"),
         x="Change in BMI (units per year)",
         title = paste0(title," "))
  
  plots[[paste(outcome, "BMI at blood sampling", sep = "_")]] <- plot(b_bmi, select = 1) +
    labs(y = paste0("Smooth effect of BMI (edf=", edf_bmi, ")"),
         x="BMI (last measurement)",
         title = paste0(title," "))
}



# Define output directory
setwd("XXX/GAM individual plots")

for(i in 1:length(plots)){
  print(plots[[i]])
  ggsave(paste0("plot",i,".png"),width = 3.8,height=3.8,dpi=500)
}


