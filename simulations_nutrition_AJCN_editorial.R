#R code is modified based on the paper by Tomova et al "Adjustment for energy intake in nutrition research: a causal inference perspective"
#Goal: estimate total effect of sugars, total effect of alcohol, relative effect of sugars vs alcohol 
# using energy parition model(i.e., all-components model), standard multivariate models, and residual model
# Load required packages
library(progress); library(devtools); library(dagitty)

# Update dagitty directly from GitHub to the latest version if required
install_github("jtextor/dagitty/r")

###----------------------------------------------------------------------------------------------
# DATA STRUCTURE

# Two datasets are simulated - without and with proxy confounding by common causes of diet,
# based on DAG and DAGU, respectively. 
###----------------------------------------------------------------------------------------------

# Define data structures for simulation and set up dataframe to store means 

# No confounding
DAG <- dagitty("dag {
                
                NMES  -> GLUC [beta = 0.25]
                CRB   -> GLUC [beta = 0.33]
                FBR   -> GLUC [beta = -0.02]
                UF    -> GLUC [beta = 0.24]
                PRO   -> GLUC [beta = 0.15]
                SF    -> GLUC [beta = 0.175]
                ALC   -> GLUC [beta = 0.09]
                
                U     -> NMES [beta = 0]
                U     -> CRB  [beta = 0]
                U     -> FBR  [beta = 0]
                U     -> SF   [beta = 0]
                U     -> UF   [beta = 0]
                U     -> PRO  [beta = 0]
                U     -> ALC  [beta = 0]

                }") 

# Confounding
DAGU <- dagitty("dag {
                
                NMES  -> GLUC [beta = 0.25]
                CRB   -> GLUC [beta = 0.33]
                FBR   -> GLUC [beta = -0.02]
                UF    -> GLUC [beta = 0.24]
                PRO   -> GLUC [beta = 0.15]
                SF    -> GLUC [beta = 0.175]
                ALC   -> GLUC [beta = 0.09]
                
                U     -> NMES [beta = 0.5]
                U     -> CRB  [beta = 0.25]
                U     -> FBR  [beta = -0.5]
                U     -> SF   [beta = 0.5]
                U     -> UF   [beta = 0.25]
                U     -> PRO  [beta = 0.25]
                U     -> ALC  [beta = 0.5]
                
                }") 

# Set seed to ensure reproducible results
set.seed(96)

# Set up empty data frames to store means later in the simulation
SimulationMeans <- data.frame(
                              model1U=numeric(0), 
                              model1U=numeric(0), 
                              model1U=numeric(0), 
                              model2U=numeric(0),
                              model2U=numeric(0),
                              model2U=numeric(0),
                              model3U=numeric(0),
                              model3U=numeric(0),
                              model3U=numeric(0),
                              model4U=numeric(0),
                              model4U=numeric(0),
                              model5U=numeric(0)
                              )


###----------------------------------------------------------------------------------------------
# DATA SIMULATION

# The simulation involves the following steps:
# 1. Re-scale of all variables based on plausible mean and standard deviation values
# 2. Calculate total and remaining energy intake from the simulated nutrients
# 3. Run the different models that adjust for energy intake: 
#   3.0. The unadjusted model
#   3.1. The energy partition model
#   3.2. The standard model
#   3.3. The nutrient density model, and the multivariable nutrient density model
#   3.4. The residual model
#   3.5. The alternative 'all-components' model
# 4. Save exposure coefficient estimates from each model and store in a vector
###----------------------------------------------------------------------------------------------

# Run 100,000 simulations with 1000 observations of data each

Nsims <- 100000
Nobs  <- 1000



# Create simulation progress bar
pb <- progress_bar$new(total = Nsims, format = ":bar :percent eta: :eta")

# Run simulation
for (i in 1:Nsims) {
  
  SimData  <- simulateSEM(DAG,  N=Nobs)
  SimUData <- simulateSEM(DAGU, N=Nobs)
  
  # Rescale variables based on plausible values, making sure that each target causal effect corresponds to the correct standardised path coefficient before re-scaling:
  
  # 5.0mgdl per 100 calories = 5/25*(125/100) = 0.25 for NMES
  # 3.3mgdl per 100 calories = 3/25*(250/100) = 0.33 for CRB
  # -1.0mgdl per 100 calories = -1/25*(50/100) = -0.02 for FBR
  # 3.0mgdl per 100 calories = 3/25*(200/100) = 0.24 for UF
  # 2.5mgdl per 100 calories = 2.5/25*(150/100) = 0.15 for PRO
  # 3.5mgdl per 100 calories = 3.5/25*(125/100) = 0.175 for SF
  # 4.5mgdl per 100 calories = (4.5/25)*(50/100) = 0.09 for ALC
  
  SimData$NMES   <- SimData$NMES*125+250  
  SimData$CRB    <- SimData$CRB*250+500   
  SimData$FBR    <- SimData$FBR*50+100    
  SimData$UF     <- SimData$UF*200+400    
  SimData$PRO    <- SimData$PRO*150+300   
  SimData$SF     <- SimData$SF*125+275    
  SimData$ALC    <- SimData$ALC*50+175    
  SimData$GLUC   <- SimData$GLUC*25+80
  
  SimUData$NMES  <- SimUData$NMES*125+250
  SimUData$CRB   <- SimUData$CRB*250+500  
  SimUData$FBR   <- SimUData$FBR*50+100  
  SimUData$UF    <- SimUData$UF*200+400  
  SimUData$PRO   <- SimUData$PRO*150+300 
  SimUData$SF    <- SimUData$SF*125+275
  SimUData$ALC   <- SimUData$ALC*50+175
  SimUData$GLUC  <- SimUData$GLUC*25+80
  
  SimUData$TotalEnergy     <- SimUData$NMES + SimUData$CRB + SimUData$FBR + SimUData$UF + SimUData$PRO + SimUData$SF + SimUData$ALC
  

  
### I only consider the case with unmeasurec common cause U
### 1. The energy partition model (same as all components model)
  
  mod1U  <- lm(GLUC ~ NMES + ALC+ CRB+ FBR + UF + PRO + SF , data=SimUData)
  

  mean1U_NMES <- mod1U$coefficients[2]*100
  mean1U_ALC  <- mod1U$coefficients[3]*100
  
  #total effect of NMES
  total_NMES1   <-mod1U$coefficients[2]*100
  
  #total effect of ALC
  total_ALC1    <- mod1U$coefficients[3]*100
  
  # relative effect of NMES vs ALC
  relative1     <-total_NMES1-total_ALC1
  

### 2. Тhe standard model
  mod2U  <- lm(GLUC ~ NMES + TotalEnergy+ CRB+ FBR + UF + PRO + SF , data=SimUData)
  
  mean2U_NMES <- mod2U$coefficients[2]*100
  mean2U_TEI  <- mod2U$coefficients[3]*100
  
  #total effect of NMES
  total_NMES2   <-mean2U_NMES+mean2U_TEI
  
  #total effect of ALC
  total_ALC2     <-mean2U_TEI 
  # relative effect of NMES vs ALC
  relative2      <-mean2U_NMES
  

### 3. The residual model 
  
  
  ResidualU <- lm(NMES ~ TotalEnergy, data=SimUData); summary(ResidualU)
  
  ResidualU1 <- lm(FBR ~ TotalEnergy, data=SimUData); summary(ResidualU1)
  
  ResidualU2 <- lm(UF ~ TotalEnergy, data=SimUData); summary(ResidualU2)
  
  ResidualU3 <- lm(PRO ~ TotalEnergy, data=SimUData); summary(ResidualU3)
  
  ResidualU4 <- lm(SF ~ TotalEnergy, data=SimUData); summary(ResidualU4)
  
  ResidualU5 <- lm(CRB ~ TotalEnergy, data=SimUData); summary(ResidualU5)
  
  ResidualU6 <- lm(ALC ~ TotalEnergy, data=SimUData); summary(ResidualU6)

  
  mod3U <- lm(SimUData$GLUC ~ ResidualU$residuals + SimUData$TotalEnergy+ ResidualU1$residuals+ResidualU2$residuals+ResidualU3$residuals+ResidualU4$residuals+ResidualU5$residuals)
  
  
  alphaU<-ResidualU$coefficients[2]
  alpha1U<-ResidualU1$coefficients[2]
  alpha2U<-ResidualU2$coefficients[2]
  alpha3U<-ResidualU3$coefficients[2]
  alpha4U<-ResidualU4$coefficients[2]
  alpha5U<-ResidualU5$coefficients[2]
  

  aveU  <- mod3U$coefficients[2]*100
  ave1U <- mod3U$coefficients[4]*100
  ave2U <- mod3U$coefficients[5]*100
  ave3U <- mod3U$coefficients[6]*100
  ave4U <- mod3U$coefficients[7]*100
  ave5U <- mod3U$coefficients[8]*100
  aveU_TEI <- mod3U$coefficients[3]*100
  
  #total effect of NMES
  total_NMES3  <- aveU+ aveU_TEI-(aveU*alphaU)-(ave1U*alpha1U)-(ave2U*alpha2U)-(ave3U*alpha3U)-(ave4U*alpha4U)-(ave5U*alpha5U)
  
  #total effect of ALC
  total_ALC3   <-        aveU_TEI-(aveU*alphaU)-(ave1U*alpha1U)-(ave2U*alpha2U)-(ave3U*alpha3U)-(ave4U*alpha4U)-(ave5U*alpha5U)
  
  #relative effect of NMES vs ALC
  relative3<-aveU
  
  
  ### 4. The energy partition model (but fail to specify other dietary variables)
  
  mod4U  <- lm(GLUC ~ NMES  , data=SimUData)
  #total effect of NMES (wrong)
  total_NMES4   <-mod4U$coefficients[2]*100

  mod5U  <- lm(GLUC ~ ALC  , data=SimUData)
  #total effect of ALC (wrong)
  total_ALC4   <-mod5U$coefficients[2]*100

  # relative effect of NMES vs ALC (wrong)
  relative4     <-total_NMES4-total_ALC4
  
  # Save local vectors of means into dataframe:
  SimulationMeans[nrow(SimulationMeans)+1,]  <- c(total_NMES1, total_ALC1, relative1, 
                                                  total_NMES2, total_ALC2, relative2, 
                                                  total_NMES3, total_ALC3, relative3,
                                                  total_NMES4, total_ALC4, relative4 
                                                  )

  rm(aveU, ave1U, ave2U, ave3U, ave4U, ave5U, aveU_TEI, 
     alphaU, alpha1U, alpha2U, alpha3U, alpha4U, alpha5U,
     ResidualU ,  Residual1U,  Residual2U,  Residual3U,  Residual4U,  Residual5U ,    
     total_NMES1, total_ALC1, relative1, 
     total_NMES2, total_ALC2, relative2, 
     total_NMES3, total_ALC3, relative3,
     total_NMES4, total_ALC4, relative4 
     )

  # Display simulation progress
  pb$tick()
  
}


# Create empty dataframes to store point coefficients and simulation intervals for each model
SummaryMeans <- data.frame(model=character(0), 
                             lower=numeric(0), 
                             point=numeric(0), 
                             upper=numeric(0))

# Create a list of model names
modelname <- list("total_NMES1", "total_ALC1", "relative1", 
                  "total_NMES2", "total_ALC2", "relative2", 
                  "total_NMES3", "total_ALC3", "relative3",
                  "total_NMES4", "total_ALC4", "relative4")


# For each model, extract the point estimate and the 95% simulation intervals and store in a vector
for (j in 1:12) {
  
  centiles                             <- c(round(quantile(SimulationMeans[,j], 0.025), digits=2), round(quantile(SimulationMeans[,j], 0.5), digits=2), round(quantile(SimulationMeans[,j], 0.975),digits=2))
  SummaryMeans[nrow(SummaryMeans)+1,]  <- c(modelname[j], unname(as.list(centiles)))
  rm(centiles)
}

# View results
View(SummaryMeans)

# Save results (optional)
write.csv(SummaryMeans, "SummaryMeans_12obs.csv")
#write.csv(SimData, "SimData_osb.csv")
#write.csv(SimUData, "SimUData_obs.csv")

