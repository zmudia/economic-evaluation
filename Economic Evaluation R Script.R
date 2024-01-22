## C1- Economic Evaluation Module LSHTM - R Code
## Author: Manzilina Mudia
## Date: 22 January 2024

# Load necessary package

library(tidyverse)


# Model structure

#![Figure 1. Model Structure](/Users/zmudia/Documents/QMUL Interview/QMUL Interview/Model Structure.png)

# Model characteristics and parameters
start_age <- 40
cycle_length <- 1
time_horizon <- 61

cycles <- c(1:time_horizon)
age <- start_age + cycles

## Willingness to pay threshold based on 27-41% of country's GDP 
## (adopted from [Ochalek, 2018](https://gh.bmj.com/content/3/6/e000964.long))

gdp <- 58000
wtp1 <- 0.27 * gdp
wtp2 <- 0.41 * gdp
wtp3 <- gdp

## Markov model 1 for SVR pathway
state_names_svr <- c("SVR", "Death") # Markov state from SVR to death
states_svr <- length(state_names_svr)


## Markov model 2 for No SVR pathway
state_names_nosvr <- c("Moderate", "Cirrhosis", "Liver Failure", "Death") # Markov state from moderate disease to death
states_nosvr <- length(state_names_nosvr)


# Set up deterministic parameters

## Fixed values
pop <- 1000 # initial population

## Probabilities

### Decision tree probabilities
pSVR_PR <- 0.395 # Probability of having SVR for Pegylated Interferon + Ribavirin
pSVR_LS <- 0.92 # Probabiliy of having SVR for Ledipasvir + Sofosbuvir

### Transition probabilities
tpModCir <- 0.043 # Probability of entering cirrhosis state from moderate disease state
tpCirLiv <- 0.039 # Probability of entering liver failure state from cirrhosis stat
tpDLiv <- 0.13 # Probability of dying from liver failure

### Life tables / death risk
death_risk_matrix <- matrix(c(40, 0.004,
                              50, 0.012,
                              60, 0.027,
                              70, 0.066,
                              80, 0.189), 
                            nrow = 5, ncol = 2, byrow = T)
interval <- findInterval(age, death_risk_matrix[,1])
tpDN <- death_risk_matrix[interval, 2] # Natural death risk


## Cost values
cdrug_PR <- 4245 # cost of common treatment (only at the decision tree)
cdrug_LS <- 40000 # cost of alternative treatment (only at the start of decision tree)

cMod <- 900 # cost of one cycle in moderate disease state
cCir <- 2500 # cost of one cycle in cirrhosis disease statae
cLiv <- 4400 # cost of one cycle in liver failure state

## Utility values
uSVR <- 0.71 # QoL weight for one cycle in SVR state
uMod <- 0.66 # QoL weight for one cycle in moderate disease state
uCir <- 0.55 # QoL weight for one cycle in cirrhosis disease state
uLiv <- 0.45 # QoL weight for one cycle in liver disease state

## Discount 
disc_costs <- 0.035
disc_outcomes <- 0.035
disc_costs_vec <- 1/(1+disc_costs)^cycles
disc_outcomes_vec <- 1/(1+disc_outcomes)^cycles


# Probabilistic Parameters

sims <- 1000 # number of simulations

gen_parameters <- function(){
  
  
  # Decision tree probabilities
  beta_parms <- function(mean, se){
    alpha <- mean * (mean * (1-mean) / (se^2) - 1)
    beta <- (alpha / mean) - alpha
    values <- rbeta(sims, alpha, beta)
    return(values)
  }
  
  pSVR_PR <- beta_parms(mean = 0.395, se = 0.02)
  pSVR_LS <- beta_parms(mean = 0.92, se = 0.03)
  
  
  # Transition probabilities
  tpModCir <- beta_parms(mean = 0.043, se = 0.01)
  tpCirLiv <- beta_parms(mean = 0.039, se = 0.01)
  tpDLiv <- beta_parms(mean = 0.13, se = 0.02)
  
  
  # Utility parameters
  uSVR <- beta_parms(mean = 0.71, se = 0.05)
  uMod <- beta_parms(mean = 0.66, se = 0.054)
  uCir <- beta_parms(mean = 0.55, se = 0.054)
  uLiv <- beta_parms(mean = 0.45, 0.054)
  
  
  # Cost parameters
  gamma_parms <- function(mean, se){
    shape = mean^2 / se^2
    scale = se^2 / mean
    values <- rgamma(sims, shape = shape, scale = scale)
    return(values)
  }
  
  cMod <- gamma_parms(mean = 900, se = 210)
  cCir <- gamma_parms(mean = 2500, se = 1060)
  cLiv <- gamma_parms(mean = 4400, se = 2335)
  
  parms_df <- data.frame(pSVR_LS, pSVR_PR, tpModCir, tpCirLiv, tpDLiv, uMod, uCir, uLiv, cMod, cCir, cLiv)
  return(parms_df)
}

prob_parameters <- gen_parameters()
parm_names <- names(prob_parameters)


# Create model calculation function

run_model <- function(parms){
  
  # Unpack parameters
  for (i in 1:length(parms)) assign(names(parms[i]), as.numeric(parms[i]))
  
  # Seed markov models from decision tree output
  pop_svr_PR <- pop * pSVR_PR
  pop_nosvr_PR <- pop * (1-pSVR_PR)
  pop_svr_LS <- pop * pSVR_LS
  pop_nosvr_LS <- pop * (1-pSVR_LS)
  
  seed_svr_PR <- c(pop_svr_PR, 0)
  seed_svr_LS <- c(pop_svr_LS, 0)
  seed_nosvr_PR <- c(pop_nosvr_PR, 0, 0, 0)
  seed_nosvr_LS <- c(pop_nosvr_LS, 0, 0, 0)
  
  
  # Markov model 1 (SVR pathway)
  
  qaly_values_svr <- c(uSVR, 0)
  
  ## Transition matrix
  transition_matrix_svr <- array(0, dim = c(states_svr, states_svr, time_horizon),
                                 dimnames = list(state_names_svr, state_names_svr))
  
  transition_matrix_svr["SVR", "SVR",] <- 1 - tpDN
  transition_matrix_svr["SVR", "Death",] <- tpDN
  transition_matrix_svr["Death", "Death",] <- 1
  
  ## Markov trace & transitions
  trace_svr_PR <- matrix(0, nrow = time_horizon, ncol = states_svr,
                         dimnames = list(paste("cycle", 1:time_horizon, sep = ""), state_names_svr))
  trace_svr_LS <- trace_svr_PR
  
  trans_svr_PR <- array(0, dim = c(states_svr, states_svr, time_horizon),
                        dimnames = list(state_names_svr, state_names_svr))
  trans_svr_LS <- trans_svr_PR
  
  ## First cycle (using the output from decision tree)
  trace_svr_PR[1,] <- seed_svr_PR %*% transition_matrix_svr[,,1]
  trace_svr_LS[1,] <- seed_svr_LS %*% transition_matrix_svr[,,1]
  
  trans_svr_PR[,,1] <- seed_svr_PR * transition_matrix_svr[,,1] 
  trans_svr_LS[,,1] <- seed_svr_LS * transition_matrix_svr[,,1]
  
  ## Reminder of cycles
  for (i in 2:time_horizon) {
    
    trace_svr_PR[i,] <- trace_svr_PR[i-1,] %*% transition_matrix_svr[,,i]
    trace_svr_LS[i,] <- trace_svr_LS[i-1,] %*% transition_matrix_svr[,,i]
    
    trans_svr_PR[,,i] <- trans_svr_PR[,,i-1] * transition_matrix_svr[,,i]
    trans_svr_LS[,,i] <- trans_svr_LS[,,i-1] * transition_matrix_svr[,,1]
  }
  
  
  # Markov model 2 (No SVR pathway)
  
  state_cost_values_nosvr <- c(cMod, cCir, cLiv, 0)
  qaly_values_nosvr <- c(uMod, uCir, uLiv, 0)
  
  ## Transition matrix
  transition_matrix_nosvr <- array(0, dim = c(states_nosvr, states_nosvr, time_horizon),
                                   dimnames = list(state_names_nosvr, state_names_nosvr))
  
  transition_matrix_nosvr["Moderate", "Moderate",] <- 1 - tpModCir - tpDN
  transition_matrix_nosvr["Moderate", "Cirrhosis",] <- tpModCir
  transition_matrix_nosvr["Moderate", "Death",] <- tpDN
  
  transition_matrix_nosvr["Cirrhosis", "Cirrhosis",] <- 1 - tpCirLiv - tpDN
  transition_matrix_nosvr["Cirrhosis", "Liver Failure",] <- tpCirLiv
  transition_matrix_nosvr["Cirrhosis", "Death",] <- tpDN
  
  transition_matrix_nosvr["Liver Failure", "Liver Failure",] <- 1 - (tpDLiv + tpDN)
  transition_matrix_nosvr["Liver Failure", "Death",] <- tpDLiv + tpDN
  
  transition_matrix_nosvr["Death", "Death",] <- 1
  
  ## Markov trace and transition
  trace_nosvr_PR <- matrix(0, nrow = time_horizon, ncol = states_nosvr,
                           dimnames = list(paste("cycle", 1:time_horizon, sep = ""), state_names_nosvr))
  trace_nosvr_LS <- trace_nosvr_PR
  
  trans_nosvr_PR <- array(0, dim = c(states_nosvr, states_nosvr, time_horizon),
                          dimnames = list(state_names_nosvr, state_names_nosvr))
  trans_nosvr_LS <- trans_nosvr_PR
  
  ## First cycle (output from decision tree)
  trace_nosvr_PR[1,] <- seed_nosvr_PR %*% transition_matrix_nosvr[,,1]
  trace_nosvr_LS[1,] <- seed_nosvr_LS %*% transition_matrix_nosvr[,,1]
  
  trans_nosvr_PR[,,1] <- seed_nosvr_PR * transition_matrix_nosvr[,,1]
  trans_nosvr_LS[,,1] <- seed_nosvr_LS * transition_matrix_nosvr[,,1]
  
  ## Reminder of cycles
  for (i in 2:time_horizon) {
    
    trace_nosvr_PR[i,] <- trace_nosvr_PR[i-1,] %*% transition_matrix_nosvr[,,i]
    trace_nosvr_LS[i,] <- trace_nosvr_LS[i-1,] %*% transition_matrix_nosvr[,,i]
    
    trans_nosvr_PR[,,i] <- trans_nosvr_PR[,,i-1] * transition_matrix_nosvr[,,i]
    trans_nosvr_LS[,,i] <- trans_nosvr_LS[,,i-1] * transition_matrix_nosvr[,,i]
  }
  
  # Check the markov trace sums to 1 across rows
  apply(trace_svr_PR, 1, sum) + apply(trace_nosvr_PR, 1, sum)
  apply(trace_svr_LS, 1, sum) + apply(trace_nosvr_LS, 1, sum)
  
  
  # Life years
  ly_values_svr <- c(1,0)
  ly_values_nosvr <- c(1,1,1,0)
  
  
  ly_vec_svr_PR <- trace_svr_PR %*% ly_values_svr
  ly_vec_nosvr_PR <- trace_nosvr_PR %*% ly_values_nosvr
  ly_vec_total_PR <- ly_vec_svr_PR + ly_vec_nosvr_PR
  
  ly_vec_svr_LS <- trace_svr_LS %*% ly_values_svr
  ly_vec_nosvr_LS <- trace_nosvr_LS %*% ly_values_nosvr
  ly_vec_total_LS <- ly_vec_svr_LS + ly_vec_nosvr_LS
  
  ly_disc_PR <- disc_outcomes_vec %*% ly_vec_total_PR
  ly_disc_LS <- disc_outcomes_vec %*% ly_vec_total_LS
  
  
  # QALYs
  qaly_vec_svr_PR <- trace_svr_PR %*% qaly_values_svr
  qaly_vec_nosvr_PR <- trace_nosvr_PR %*% qaly_values_nosvr
  qaly_vec_total_PR <- qaly_vec_svr_PR + qaly_vec_nosvr_PR
  
  qaly_vec_svr_LS <- trace_svr_LS %*% qaly_values_svr
  qaly_vec_nosvr_LS <- trace_nosvr_LS %*% qaly_values_nosvr
  qaly_vec_total_LS <- qaly_vec_svr_LS + qaly_vec_nosvr_LS
  
  qaly_disc_PR <- disc_outcomes_vec %*% qaly_vec_total_PR
  qaly_disc_LS <- disc_outcomes_vec %*% qaly_vec_total_LS
  
  
  # Costs
  state_cost_vec_nosvr_PR <- trace_nosvr_PR %*% state_cost_values_nosvr
  state_cost_vec_nosvr_LS <- trace_nosvr_LS %*% state_cost_values_nosvr
  
  cost_vec_nosvr_PR <- colSums(state_cost_vec_nosvr_PR) + (pop*cdrug_PR)
  cost_vec_nosvr_LS <- colSums(state_cost_vec_nosvr_LS) + (pop*cdrug_LS)
  
  cost_disc_PR <- disc_costs_vec %*% state_cost_vec_nosvr_PR + (pop*cdrug_PR)
  cost_disc_LS <- disc_costs_vec %*% state_cost_vec_nosvr_LS + (pop*cdrug_LS)
  
  
  # Analysis
  output <- c(inc_cost = cost_disc_LS - cost_disc_PR,
              inc_qaly = qaly_disc_LS - qaly_disc_PR,
              inc_ly = ly_disc_LS - ly_disc_PR) 
  
  return(output)
}


# Passing parameters into the model

# Selecting parameters
det_parms <- c(pSVR_LS, pSVR_PR, tpModCir, tpCirLiv, tpDLiv, uMod, uCir, uLiv, cMod, cCir, cLiv)
names(det_parms) <- parm_names

## Deterministic model
det_results <- run_model(det_parms)
det_results <- round(det_results, 3)

### Plot
det_results_df <- det_results %>% t() %>% as.data.frame

(ICER_plot <- ggplot(det_results_df, aes(x = inc_qaly, y = inc_cost)) +
    geom_point(shape = 21, fill = "red", alpha = 0.6) +
    theme_classic() + 
    geom_abline(aes(intercept = 0, slope = wtp1, colour = "WTP 1: 27% GDP ~ 15660 EGP")) +
    geom_abline(aes(intercept = 0, slope = wtp2, colour = "WTP 2: 41% GDP ~ 23780 EGP")) + 
    labs(x = "QALYs difference", y = "Costs Difference (EGP)", 
         colour = "Willingness to pay threshold (WTP)", caption = "Figure 2. Base case analysis plot on ICER plane") +
    theme(axis.title = element_text(face = "bold"),
          plot.margin = unit(c(0.2, 0.7, 0.2, 0.2), "cm"),
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)),
          plot.caption = element_text(hjust = 0),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.position = c(0.8, 0.3)) +
    guides(colour = guide_legend(override.aes = list(size = 7))) +
    scale_y_continuous(limits = c(0, det_results_df$inc_cost*1.2),
                       labels = scales::label_number(big.mark = ","), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, det_results_df$inc_qaly*1.2),
                       expand = c(0, 0)))


# Probabilistic Sensitivity Analysis

simulation_results <- matrix(NA, nrow = sims, ncol = 3)
colnames(simulation_results) <- c("inc_cost", "inc_ly", "inc_qaly")

for (i in 1:sims) simulation_results[i,] <- run_model(prob_parameters[i, ])
print(colMeans(simulation_results), digits = 3)

## CEAC
lambda <- seq(0, gdp, 100)
lambda_table <- matrix(lambda, ncol = length(lambda), nrow = dim(simulation_results)[1], byrow = T)
inmb_count <- ((simulation_results[,3] * lambda_table) - simulation_results[,1]) > 0
prob_ce <- colMeans(inmb_count)
ceac <- data.frame(lambda, prob_ce)

## Plots
### CE Plane
sim_res_df <- as.data.frame(simulation_results)

(ce_plane <- ggplot(sim_res_df, aes(x = inc_qaly, y = inc_cost)) +
    geom_point(shape = 21, fill = "dodgerblue1", alpha = 0.6) + 
    theme_classic() +
    geom_abline(aes(intercept = 0, slope = wtp1, colour = "WTP 1: 27% GDP ~ 15660 EGP")) +
    geom_abline(aes(intercept = 0, slope = wtp2, colour = "WTP 2: 41% GDP ~ 23780 EGP")) +
    geom_abline(aes(intercept = 0, slope = wtp3, colour = "WTP 3: 1x GDP ~ 58000 EGP")) +
    labs(x = "Incremental QALYs", y = "Incremental costs (EGP)", 
         colour = "Willingness to pay threshold (WTP)", caption = "Figure 3. Cost effectiveness plane") +
    theme(axis.title = element_text(face = "bold"),
          plot.margin = unit(c(0.2, 0.7, 0.2, 0.2), "cm"),
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)),
          plot.caption = element_text(hjust = 0),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.position = c(0.8, 0.2)) +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    scale_y_continuous(limits = c(0, max(sim_res_df$inc_cost)*1.2),
                       labels = scales::label_number(big.mark = ","), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, max(sim_res_df$inc_qaly)*1.2),
                       expand = c(0, 0)))


### CEAC Plot

(ceac_plot <- ggplot(ceac) + 
    geom_line(aes(x = lambda, y = prob_ce)) +
    theme_classic() +
    geom_vline(aes(xintercept = wtp1, colour = "WTP 1: 15660 EGP")) +
    geom_vline(aes(xintercept = wtp2, colour = "WTP 2: 23780 EGP")) +
    geom_vline(aes(xintercept = wtp3, colour = "WTP 3: 58000 EGP")) +
    labs(x = "Willingness to pay threshold (EGP)", y = "Probability of cost-effectiveness", 
         colour = "", caption = "Figure 4. Cost-effectiveness acceptability curve") + 
    theme(axis.title = element_text(face = "bold"),
          plot.margin = unit(c(0.2, 0.7, 0.2, 0.2), "cm"),
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)),
          legend.direction = "horizontal",
          legend.position = c(0.4,1),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          plot.caption = element_text(hjust = 0)) +
    guides(colour = guide_legend(override.aes = list(size = 8))) +
    scale_y_continuous(limits = c(0,1), breaks = c(seq(0, 1, 0.2)), expand = c(0,0)) +
    scale_x_continuous(limits = c(0,wtp3), labels = scales::label_number(big.mark = ","), expand = c(0,0)))
