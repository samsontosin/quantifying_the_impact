## The R codes for Quantifying the impact of Wolbachia releases on dengue infection in Townsville, Australia.

## Required packages needed for the program
library(dplyr)
library(deSolve)
library(ggplot2)
library(car)
library(MASS)
library(fitdistrplus)
library(tidyr)
library(lubridate)
library(readxl)

## Import data file for "Imported monthly dengue case notifications in Townsville"
Imported_cases <- read.csv("dengue_imported2.csv")
Imported_cases
Imported_cases$Date = as.Date(Imported_cases$Date)
Imported_cases

## Import data file for "Locally acquired monthly dengue case notifications in Townsville"
Locally_acquired_cases <- read.csv("dengue_locally22.csv")
Locally_acquired_cases
Locally_acquired_cases$Date = as.Date(Locally_acquired_cases$Date)
Locally_acquired_cases


## Plot the imported dengue cases data with respect to status (Before and After Wolbachia releases)
ggplot2::ggplot(Imported_cases, aes(x=Date, y=Cases, fill=Status)) +
  geom_col() +
  ylab("Dengue cases") +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(size = 0.4, colour = "grey10"),
        text = element_text(size=12,  family="serif"),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.position = "top",
        strip.background =element_rect(fill="royalblue"),
        strip.text = element_text(size = 10, colour = 'white')) + 
  ggtitle("Imported cases")

## Plot the locally acquired dengue cases data with respect to status (Before and After Wolbachia releases)
ggplot2::ggplot(Locally_acquired_cases, aes(x=Date, y=Cases, fill=Status)) +
  geom_col() +
  ylab("Dengue cases") +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(size = 0.4, colour = "grey10"),
        text = element_text(size=12,  family="serif"),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.position = "top",
        strip.background =element_rect(fill="royalblue"),
        strip.text = element_text(size = 10, colour = 'white')) +
  ggtitle("Local cases")

#the initial local cases data ends in 01/01/2017
## Define the dates
dates <- c("02/01/17", "02/01/19", "02/01/19", "01/01/01", "10/01/14", "01/01/80")
wolb_stop_date = as.Date(dates[1], "%m/%d/%y")        ## The Wolbachia release stop date
end_date_im = as.Date(dates[2], "%m/%d/%y")           ## The end date for imported cases in the data
end_date_lc = as.Date(dates[3], "%m/%d/%y")
start_date = as.Date(dates[4], "%m/%d/%y")
post_wolb_date = as.Date(dates[5], "%m/%d/%y")
start_sim_date = as.Date(dates[6], "%m/%d/%y")
cases_period = Locally_acquired_cases %>%
  filter(Date >= start_date, Date <= end_date_lc)

# Determine the number of observations
n_observations = length(Locally_acquired_cases$Date)
n_observations = nrow(Locally_acquired_cases)

# Model parameters
parameters = c(
  activation_rate = 1/5.5, #Gubler 1998 from Ndii 
  recovery_rate = 1/5,        #Ndii et al
  biting_rate_u = 0.3,
  trans_prob_u = 0.2147,
  biting_rate_w = 0.3 * 0.95,
  trans_prob_wh = 0.0052,
  mu = 0.000034,         #0.000034,              # Adeshina et al
  xi_1 = 0.4 * 4 * 12 / 365.25, #rate of importation before Wolbachia per month 0.4
  xi_2 = 1 * 12 / 365.25,  #rate of importation after Wolbachia per month 1
  prob_symp = 1/4,          #Kamtchum-Tatuene et al
  mu_u = 0.043,
  mu_w = 0.068,
  mu_Au = 0.02,
  mu_Aw = 0.02,
  rho_u = 13,                #my paper
  rho_w = 10,               #my paper
  gamma = 0.95,             #walker et al
  phi = 0.95,               #ant et al
  tau_u = 0.11,             #walker and hoffman
  tau_w = 0.11,             #walker et al
  mosq_act_rate = 0.1,      #chowell et al
  mosq_act_wol_rate = 0.1,    #chowel et al
  sigma = 0,
  wolb_rate = 5e3,
  L = 2
)

# Initial conditions
Total_population = 187500
Initial_exposed_im = 0
Initial_exposed_lc = 0
Initial_infected_im_sym = 0
Initial_infected_lc_sym = 0
Initial_cuminc = 0
Initial_recovered = 0
Initial_susceptible = Total_population - Initial_exposed_im - Initial_exposed_lc - Initial_infected_im_sym - Initial_infected_lc_sym - Initial_recovered
Lambda = Initial_susceptible*parameters["mu"]
Initial_Aq_vec_uninfected = 200000
Initial_Sus_vec_uninfected = 175000
Initial_Exp_vec_uninfected = 0
Initial_Inf_vec_uninfected = 0
Initial_Aq_wol_vec_uninfected = 0
Initial_Sus_wol_vec_uninfected = 0
Initial_Exp_wol_vec_uninfected = 0
Initial_Inf_wol_vec_uninfected = 0


# State variables
          #human state variables
state = c(Susceptible = Initial_susceptible,
          Exposed_im = Initial_exposed_im,
          Exposed_lc = Initial_exposed_lc,
          Infected_im_sym = Initial_infected_im_sym,
          Infected_lc_sym = Initial_infected_lc_sym,
          cuminc = Initial_cuminc,
          Recovered = Initial_recovered,
          #mosquito state variables
          Aq_vec = Initial_Aq_vec_uninfected,
          Sus_vec = Initial_Sus_vec_uninfected,
          Exp_vec = Initial_Exp_vec_uninfected,
          Inf_vec = Initial_Inf_vec_uninfected,
          Aq_wol_vec = Initial_Aq_wol_vec_uninfected,
          Sus_wol_vec = Initial_Sus_wol_vec_uninfected,
          Exp_wol_vec = Initial_Exp_wol_vec_uninfected,
          Inf_wol_vec = Initial_Inf_wol_vec_uninfected
)

# Define the dengue importation events
Imp_event1 = data.frame(Date=c(seq.Date(start_sim_date, start_date-1, by = 'month')), Cases=c(rep(0,length(Date))))
Imp_event2 = data.frame(subset(Imported_cases, select = c("Date", "Cases"))) #data.frame(c(times[times>="2001-01-01"]), c(Imported_cases$Cases))
Imp_event2$Cases = Imp_event2$Case * 4
Imported_event = rbind.data.frame(Imp_event1, Imp_event2)
#Imported_event["Cases"][Imported_event["Cases"] != 0] <- 0

# Model function
Dengue_base_model = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    Total_population = Susceptible + Exposed_im + Exposed_lc + Infected_im_sym + Infected_lc_sym + Recovered    # Total human population
    F_u = Sus_vec + Exp_vec + Inf_vec                   #non-wolbachia adult mosquitoes
    F_w = Sus_wol_vec + Exp_wol_vec + Inf_wol_vec       #Wolbachia-infected adult mosquitoes
    FF = F_u + F_w                                      #All adult mosquitoes
    QQ = Aq_vec + Aq_wol_vec                            #All aquatic-staged mosquitoes 
    K = L * Total_population * (cos (2 * pi * (t + 15)/(365.25)) + 1) / 2  #Chowell et al                       #Carrying capacity
    kappa = ifelse((t >  as.numeric(post_wolb_date - start_sim_date)) &&
                                     t <  as.numeric(wolb_stop_date - start_sim_date), wolb_rate, 0)                #Wolbachia mosquitoes importation rate
    
    xi = ifelse(t < as.numeric(start_date - start_sim_date), xi_1, 0)  #importation rate of dengue

    beta1_plus_beta2 = ((biting_rate_u * trans_prob_u * L * Inf_vec) + (biting_rate_w * trans_prob_wh * L * Inf_wol_vec))/Total_population
    
    beta3 =  biting_rate_u * trans_prob_u * (Infected_im_sym + Infected_lc_sym)/Total_population  # with uninfected mosq biting rate
    #beta3 = 0
    beta4 =  biting_rate_w * trans_prob_u * (Infected_im_sym + Infected_lc_sym)/Total_population # with Wolb-mosq biting rate
    #beta4 = 0
    # Calculate the net (instantaneous) change in each state variable
    Susceptible_change = Lambda - beta1_plus_beta2 * Susceptible - mu * Susceptible
    Exposedim_change =  - Exposed_im * activation_rate - mu * Exposed_im + xi
    Exposedlc_change = beta1_plus_beta2 * Susceptible - Exposed_lc * activation_rate - mu * Exposed_lc
    Infectedim_sym_change = Exposed_im * activation_rate - Infected_im_sym * recovery_rate - mu * Infected_im_sym
    Infectedlc_sym_change = Exposed_lc * activation_rate - Infected_lc_sym * recovery_rate - mu * Infected_lc_sym
    Cummulative_inc = Exposed_lc * activation_rate
    Recovered_change = (Infected_im_sym + Infected_lc_sym) * recovery_rate - mu * Recovered
    Aq_vec_change = max(0, ((((rho_u * F_u^2) + (rho_w * ((1-gamma)*F_w^2 + (1-phi)*F_w*F_u)))/FF) * (1-(QQ/K)))) - ((tau_u + mu_Au) * Aq_vec) 
    Sus_vec_change = (tau_u * Aq_vec/2) - (beta3 + mu_u)* Sus_vec
    Exp_vec_change = beta3 * Sus_vec - (mosq_act_rate + mu_u)* Exp_vec
    Inf_vec_change = mosq_act_rate * Exp_vec - mu_u * Inf_vec + sigma * Inf_wol_vec
    Aq_wol_vec_change = (rho_w * (gamma * F_w^2 + phi * F_w*F_u)/FF) * (1-(QQ/K)) - (tau_w + mu_Aw) * Aq_wol_vec + kappa
    Sus_wol_vec_change = (tau_w * Aq_wol_vec/2) - (beta4 + mu_w)* Sus_wol_vec
    Exp_wol_vec_change = beta4 * Sus_wol_vec - (mosq_act_wol_rate + mu_w)* Exp_wol_vec
    Inf_wol_vec_change = mosq_act_wol_rate * Exp_wol_vec - mu_w * Inf_wol_vec - sigma * Inf_wol_vec
    return(list(
      c(
        Susceptible_change,
        Exposedim_change,
        Exposedlc_change,
        Infectedim_sym_change,
        Infectedlc_sym_change,
        Cummulative_inc,
        Recovered_change,
        Aq_vec_change,
        Sus_vec_change,
        Exp_vec_change,
        Inf_vec_change,
        Aq_wol_vec_change,
        Sus_wol_vec_change,
        Exp_wol_vec_change,
        Inf_wol_vec_change
      )
    ))
  })
}

#Change date to numeric variable
Imported_event$t = as.numeric(Imported_event$Date - Imported_event$Date[1])

# Create wrapper that solves the model across time window, segmented by importation times
solve_base_model = function(y, model_func, parms, imports){
  
  full_output = matrix(data = c(NA, y), nrow = 1, ncol = 16, dimnames = list(NULL, c("time", "Susceptible", "Exposed_im", "Exposed_lc", "Infected_im_sym", "Infected_lc_sym", "Recovered", "cuminc", "Aq_vec", "Sus_vec", "Exp_vec", "Inf_vec", "Aq_wol_vec", "Sus_wol_vec", "Exp_wol_vec", "Inf_wol_vec")))
  
  for (i in 1:(nrow(imports)-1)){
    y_initial = full_output[nrow(full_output), -1]
    y_initial["Infected_im_sym"] = y_initial["Infected_im_sym"] + imports$Cases[i]
    output = ode(y_initial, seq(imports$t[i], imports$t[i+1], by=1), model_func, parms)
    full_output = rbind(full_output[-nrow(full_output), ], output)
}
  return(full_output)
}
out = solve_base_model(state, Dengue_base_model, parameters, Imported_event)
out 
  
  # Calculate the prevalence, incidence and cumulative incidence (for comparison with data)
  Wolb_frequency = (out[, "Aq_wol_vec"] + out[, "Sus_wol_vec"]+out[, "Exp_wol_vec"] + out[, "Inf_wol_vec"])/(out[, "Aq_vec"] + out[, "Sus_vec"] +out[, "Exp_vec"] + out[, "Inf_vec"] +out[, "Aq_wol_vec"] + out[, "Sus_wol_vec"]+out[, "Exp_wol_vec"] + out[, "Inf_wol_vec"])
  Prevalence = out[, "Exposed_lc"] + out[, "Infected_lc_sym"]
  Incidence = parameters["prob_symp"] * (out[, "Exposed_lc"]) * parameters["activation_rate"]
  Cumulative_incidence = cumsum(Incidence) + out[1, "Infected_lc_sym"]

  # Append these derived outputs to the main solution
  out = cbind(out, Prevalence, Incidence, Cumulative_incidence, Wolb_frequency)
  colnames(out)[c(17, 18, 19, 20)] = c("Prevalence", "Incidence", "Cumulative_incidence","Wolb_frequency")
  
  out = as.data.frame(out)
  
  Dates = seq.Date(start_sim_date, end_date_lc, by = 'day')
  out$Date = Dates

  out = as.data.frame(out)
  
  # Aggregate the output to monthly output to include the monthly cumulative incidence
  out = out %>% mutate(Month = month(Date), Year = year(Date)) %>%
   mutate(Ini_date = ymd(paste(Year, Month, "01", sep = "-"))) %>%
   group_by(Ini_date) %>%
   filter(Date == max(Date))
  out = as.data.frame(out)
  out
  monthly_incidence = c(out$Cumulative_incidence[1], diff(out$Cumulative_incidence))#  parameters["prob_symp"] * c(out$cuminc[1], diff(out$cuminc))
  
  out = as.data.frame(out)

  # Finally, add in the dates corresponding to the time points in the data
  out$Incidence_monthly = monthly_incidence
  out$Status = c(ifelse(out$time < as.numeric(post_wolb_date - start_sim_date), "Pre-Wolbachia", "Post-Wolbachia"))
  #out$Status = Locally_acquired_cases$Status
  out.init = data.frame(out)
  out.init

#Number of Wolbachia-infected mosquitoes with time
Wol_only = ggplot(out.init, aes(x=Date, y=(Aq_wol_vec + Sus_wol_vec + Exp_wol_vec + Inf_wol_vec))) +
  xlab("Time") +
  ylab(expression(paste("Number of ", italic("Wolbachia-"),"infected mosquitoes")))+
  geom_line() +theme_bw()
Wol_only

#Number of non-Wolbachia mosquitoes with time
Wol_non = ggplot(out.init, aes(x=Date, y=(Aq_vec + Sus_vec + Exp_vec + Inf_vec)))+
  xlab("Time") +
  ylab(expression(paste("Number of non-", italic("Wolbachia"),"mosquitoes")))+
  geom_line() +theme_bw()
Wol_non

#Wolbachia frequency i.e. the ratio of Wolbachia-infected mosquitoes to the total mosquito population with time
Freq_Wol = ggplot(out.init, aes(x=Date, y=(Aq_wol_vec + Sus_wol_vec + Exp_wol_vec + Inf_wol_vec)/(Aq_vec + Sus_vec + Exp_vec + Inf_vec + Aq_wol_vec + Sus_wol_vec + Exp_wol_vec + Inf_wol_vec)))+
  xlab("Time") +
  ylab(expression(paste(italic("Wolbachia "),"frequency")))+
  geom_line() +theme_bw()
Freq_Wol

Incid = ggplot(out.init, aes(x=Date, y=(Incidence_monthly)))+
  xlab("Time") +
  ylab("Dengue incidence")+
  geom_line() +theme_bw()
Incid


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Import data file for "Wolbachia releases dataset in Townsville"
Wolbachia_data <- read_xlsx("Wolbachia_dataset_townsville.xlsx", sheet = 1)

#Convert data into dataframe
Wolbachia_data1 <- data.frame(Wolbachia_data)

class(Wolbachia_data1$Collection.Date)

#Convert the collection date to the "Date" format
Wolbachia_data1$Collection.Date <- as.Date(Wolbachia_data1$Collection.Date,
                                           format = "%Y-%m-%d")

#Order the collection date from the first recorded date of release in all suburbs
Wolbachia_data1 = Wolbachia_data1[order(Wolbachia_data1$Collection.Date),]
Wolbachia_data1

#Compute the number of Wolbachia positive mosquitoes from the percentages given
Wolbachia_data1$num.of.Wol.pos = (Wolbachia_data1$Number.Aedes.aegypti * Wolbachia_data1$perc.wolbachia.positive / 100)

#Convert the Dates to Months
Wolbachia_data1$month <- strftime(Wolbachia_data1$Collection.Date, "%Y-%m")

#Sum up the total aedes mosquitoes captured in each month
Wolbachia_sum_aedes <- aggregate(Wolbachia_data1$Number.Aedes.aegypti, list(Wolbachia_data1$month), sum)

#Sum up the total Wolbachia-infected mosquitoes captured from the aedes mosquitoes in each month 
Wolbachia_sum_wpos <- aggregate(Wolbachia_data1$num.of.Wol.pos, list(Wolbachia_data1$month), sum)

#Create a new dataframe to comprise the number of Wolbachia-infected and the total aedes mosquitoes captured monthly
new_Wolbachia_data <- data.frame(Wolbachia_sum_aedes,Wolbachia_sum_wpos)
new_Wolbachia_data
colnames(new_Wolbachia_data)

#Rename the column names
names(new_Wolbachia_data) <- c('month', 'Num_of_aedes_collected', 'month1', 'Num_of_Wol_pos_aedes')    #month = month1
new_Wolbachia_data

#Compute the monthly percentages of Wolbachia-positive and Wolbachia-negative mosquitoes
new_Wolbachia_data$perc_of_Wol_pos_aedes <- (new_Wolbachia_data$Num_of_Wol_pos_aedes / new_Wolbachia_data$Num_of_aedes_collected * 100)
new_Wolbachia_data$perc_of_Wol_neg_aedes <- (100 - new_Wolbachia_data$perc_of_Wol_pos_aedes)
new_Wolbachia_data

# Finally, add in the dates in the model to corresponding to the time points in the data
out.init1 = out.init[out.init$Date >= "2014-11-01", ]
new_Wolbachia_data$month1 <- out.init1$Date

# Proportion of Wolbachia-infected mosquitoes
new_Wolbachia_data$prop <- new_Wolbachia_data$perc_of_Wol_pos_aedes /100

## Plot of the wolbachia frequency in the model and Wolbachia-positive proportion in the data
library(reshape2)
library(scales)
nnn = ggplot()+ 
  geom_line(data = out.init1, aes(x=Date, y=Wolb_frequency, color = "Model"), size=2) +
  geom_point(data = new_Wolbachia_data, aes(x = month1, y = prop, color = "Data"), size = 2) + 
  xlab("Time") +
 ylab(expression(paste(italic("Wolbachia "),"frequency")))+
  #labs(color = "Legend") +
  scale_colour_manual("Colour", 
                      breaks = c("Model", "Data"),
                      values = c("cyan", "red")) +
  scale_x_date(date_breaks = "1 year", 
               labels=date_format(format = "%Y"),
               limits = as.Date(c('2014-09-01','2019-03-01'))) +
  theme_bw()
nnn

Wol_freq = ggplot(out.init, aes(x=Date, y=Wolb_frequency))+geom_line() +theme_bw()
Wol_freq

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate the negative log-likelihood for our parameter values
out.init2 = out.init[out.init$Date >= "2001-01-01", ]
poisson_negative_log_likelihood_initial = -sum(dpois(cases_period$Cases,
                                                     out.init2$Incidence_monthly, log = TRUE))
poisson_negative_log_likelihood_initial

init.parameters = log(c("trans_prob_u" = 0.22, "trans_prob_wh" = 0.007))

#solve_base_model = function(y, model_func, parms, imports)

nllik = function(t_parameters,
                 state_base = state,
                 times. = Imported_event,
                 func. = Dengue_base_model,
                 parms_base = parameters
){
  trans_prob_u = exp(t_parameters["trans_prob_u"])
  trans_prob_wh = exp(t_parameters["trans_prob_wh"])
  parms_base["trans_prob_u"] = trans_prob_u
  parms_base["trans_prob_wh"] = trans_prob_wh
  print(parms_base)
  print(state_base)
  
  out1 = solve_base_model(state_base,
                         func.,
                         parms_base,
                         times.)
  #out1  = out[out$Date >= "2001-01-01", ]
  # Calculate the prevalence, incidence and cumulative incidence (for comparison with data)
  Wolb_frequency = (out1[, "Aq_wol_vec"] + out1[, "Sus_wol_vec"]+out1[, "Exp_wol_vec"] + out1[, "Inf_wol_vec"])/(out1[, "Aq_vec"] + out1[, "Sus_vec"] +out1[, "Exp_vec"] + out1[, "Inf_vec"] +out1[, "Aq_wol_vec"] + out1[, "Sus_wol_vec"]+out1[, "Exp_wol_vec"] + out1[, "Inf_wol_vec"])
  Prevalence = out1[, "Exposed_lc"] + out1[, "Infected_lc_sym"]
  Incidence = parameters["prob_symp"] * (out1[, "Exposed_lc"]) * parameters["activation_rate"]
  Cumulative_incidence = cumsum(Incidence) + out1[1, "Infected_lc_sym"]
  
  # Append these derived out1puts to the solution
  out1 = cbind(out1, Prevalence, Incidence, Cumulative_incidence, Wolb_frequency)
  colnames(out1)[c(17, 18, 19, 20)] = c("Prevalence", "Incidence", "Cumulative_incidence","Wolb_frequency")
  
  out1 = as.data.frame(out1)
  
  Dates = seq.Date(start_sim_date, end_date_lc, by = 'day')
  out1$Date = Dates
  
  out1 = as.data.frame(out1)
  
  # Aggregate the output to monthly output to include the monthly cumulative incidence
  out1 = out1 %>% mutate(Month = month(Date), Year = year(Date)) %>%
    mutate(Ini_date = ymd(paste(Year, Month, "01", sep = "-"))) %>%
    group_by(Ini_date) %>%
    filter(Date == max(Date))
  out1 = as.data.frame(out1)
  out1
  
  monthly_incidence = c(out1$Cumulative_incidence[1], diff(out1$Cumulative_incidence))# 
  
  out1 = as.data.frame(out1)
  
  out1$Incidence_monthly = monthly_incidence
  out1$Status = c(ifelse(out$time < as.numeric(post_wolb_date - start_sim_date), "Pre-Wolbachia", "Post-Wolbachia"))
  #out1$Status = Locally_acquired_cases$Status
  out2 = out1[out1$Date>="2001-01-01", ]
  out2 = data.frame(out2)
  print(head(out2,5))
    return(-sum(dpois(cases_period$Cases, out2$Incidence_monthly, log=TRUE)))
}

nllik.initial=nllik(init.parameters)
nllik.initial

###### grid search using nll ###########################################################
# trans_prob_u_vec <-  seq(from=0.01, to=0.3, by=0.02)
# trans_prob_wh_vec = seq(from=0.01, to=0.3, by=0.02)
# 
# nll.grid = function(trans_prob_u_vec, trans_prob_wh_vec) {
#   parameter_grid = expand.grid(trans_prob_u = trans_prob_u_vec,
#                                trans_prob_wh = trans_prob_wh_vec
#   )
#   with(parameter_grid,{
#     nll = mapply(function(x, y) {
#       nllik(log(c(
#         "trans_prob_u" = x,
#         "trans_prob_wh" = y
#       )))
#     },
#     trans_prob_u,
#     trans_prob_wh)
# 
#     return(data.frame(trans_prob_u=trans_prob_u,
#                       trans_prob_wh=trans_prob_wh,
#                       NLL = nll))
#   })
# 
# }
# 
# # call the function that will calculate likelihood surface
# nll_grid_results = nll.grid(trans_prob_u_vec,trans_prob_wh_vec)
# 
# 
# # Plot likelihood surface
# nll_plot = nll_grid_results %>%
#   ggplot2::ggplot(aes(
#     x = trans_prob_u,
#     y = trans_prob_wh,
#     fill = NLL,
#     z = NLL
#   )) +
#   geom_tile() +
#   geom_contour(breaks = seq(100, 700, 1000)) +
#   ylab("trans_prob_wh") +
#   xlab("trans_prob_u") +
#   ggtitle("Likelihood surface") +
#   theme_bw()
# 
# # Extract optimum from grid search
# grid_optimum = nll_grid_results[which.min(nll_grid_results$NLL), ]
# grid_optimum
# 
# # Superimpose optimum on likelihood surface plot
# nll_plot <- nll_plot +
#   geom_point(aes(x=grid_optimum$trans_prob_u,
#                  y=grid_optimum$trans_prob_wh),
#              colour="red",
#              size=3)
# # Superimpose optimal point to negative-log-likelihood plot
# # nll_plot <- nll_plot + geom_point(aes(x = optimum_parameters["trans_prob_u"],
# #                                       y = optimum_parameters["trans_prob_wh"]),
# #                                   colour = "yellow",
# #                                   size = 3)
# 
# 
# 
# print(nll_plot)
# #### grid search ends #####################################################################

##Model fitting using the optim function with initial parameters a_u and a_wh

optim_NM <- optim(par = init.parameters,
                  fn = nllik,
                  control = list(maxit=500),
                  method = "Nelder-Mead",
                  hessian = TRUE)
# Inspect solution
optim_NM$par
optim_NM$convergence
# Back-transform parameters
trans_prob_u_optim = exp(optim_NM$par["trans_prob_u"])
trans_prob_wh_optim = exp(optim_NM$par["trans_prob_wh"])
trans_prob_u_optim
trans_prob_wh_optim

optim_NM$value
optim_NM$convergence

optim_NM$hessian
# Covariance matrix
covar_matrix = solve(optim_NM$hessian)
covar_matrix
# We use the covariance matrix to generate standard errors for each parameter
std_errors = sqrt(diag(covar_matrix))
estimates_transformed = data.frame(estimate = optim_NM$par,
                                   std_error = std_errors) %>%
  dplyr::mutate(lower95CI = estimate - 1.96 * std_error,
                upper95CI = estimate + 1.96 * std_error)
estimates_transformed
# Back-transform estimates and confidence intervals
estimates = exp(estimates_transformed[,-2])
estimates
############################################################################################################

# Calculate the optimal solution
optim_params = parameters
optim_params["trans_prob_u"] = trans_prob_u_optim
optim_params["trans_prob_wh"] = trans_prob_wh_optim
#optim_params["wolb_rate"] = 0

optim_state = state

optim_solution = solve_base_model(optim_state,
                                  Dengue_base_model,
                                  optim_params,
                                  Imported_event)
Wolb_frequency = (optim_solution[, "Aq_wol_vec"] + optim_solution[, "Sus_wol_vec"]+optim_solution[, "Exp_wol_vec"] + optim_solution[, "Inf_wol_vec"])/(optim_solution[, "Aq_vec"] + optim_solution[, "Sus_vec"] +optim_solution[, "Exp_vec"] + optim_solution[, "Inf_vec"] +optim_solution[, "Aq_wol_vec"] + optim_solution[, "Sus_wol_vec"]+optim_solution[, "Exp_wol_vec"] + optim_solution[, "Inf_wol_vec"])
Prevalence = optim_solution[, "Exposed_lc"] + optim_solution[, "Infected_lc_sym"]
Incidence = parameters["prob_symp"] * (optim_solution[, "Exposed_lc"]) * parameters["activation_rate"]
Cumulative_incidence = cumsum(Incidence) + optim_solution[1, "Infected_lc_sym"]
#monthly_incidence = optim_solution[, "cuminc"]

# Append these derived optim_solutionputs to the main solution
optim_solution = cbind(optim_solution, Prevalence, Incidence, Cumulative_incidence, Wolb_frequency)
colnames(optim_solution)[c(17, 18, 19, 20)] = c("Prevalence", "Incidence", "Cumulative_incidence","Wolb_frequency")

optim_solution = as.data.frame(optim_solution)

# Finally, add in the dates corresponding to the time points
Dates = seq.Date(start_sim_date, end_date_lc, by = 'day')
optim_solution$Date = Dates

optim_solution = as.data.frame(optim_solution)
#   
# Aggregate the output to monthly output to include the monthly cumulative incidence
optim_solution = optim_solution %>% mutate(Month = month(Date), Year = year(Date)) %>%
  mutate(Ini_date = ymd(paste(Year, Month, "01", sep = "-"))) %>%
  group_by(Ini_date) %>%
  filter(Date == max(Date))
optim_solution = as.data.frame(optim_solution)
optim_solution

monthly_incidence = c(optim_solution$Cumulative_incidence[1], diff(optim_solution$Cumulative_incidence))#  parameters["prob_symp"] * c(optim_solution$cuminc[1], diff(optim_solution$cuminc))

optim_solution = as.data.frame(optim_solution)

optim_solution$Incidence_monthly = monthly_incidence
optim_solution$Status = c(ifelse(optim_solution$time < as.numeric(post_wolb_date - start_sim_date), "Pre-Wolbachia", "Post-Wolbachia"))
optim_solution
optim_solution = optim_solution[optim_solution$Date >= "2001-01-01", ]
out.init$lower50 = qpois(0.25, out.init$Incidence)
out.init$upper50 = qpois(0.75, out.init$Incidence)

out.init$lower95 = qpois(0.025, out.init$Incidence)
out.init$upper95 = qpois(0.975, out.init$Incidence)

optim_solution$Date = as.Date(Locally_acquired_cases$Date)

#optim_solution$R_eff = 0.9 * (optim_solution$Susceptible / Total_population)

optim_solution$lower50 = qpois(0.25, optim_solution$Incidence_monthly)
optim_solution$upper50 = qpois(0.75, optim_solution$Incidence_monthly)

optim_solution$lower95 = qpois(0.025, optim_solution$Incidence_monthly)
optim_solution$upper95 = qpois(0.975, optim_solution$Incidence_monthly)

##############################################################################################################

# Plotting the locally acquired data and model

library(scales)
sammy15<- ggplot(Locally_acquired_cases) +
  geom_col(aes(x=Date, y=Cases, fill=Status), width=50) +
  geom_line(data=optim_solution, aes(x=Date, y=Incidence_monthly, colour="Cumulative monthly incidence"), size=1) +
  geom_ribbon(data=optim_solution, aes(x=Date, ymin=lower50, ymax=upper50), size=1, alpha=0.5) +
  geom_ribbon(data=optim_solution, aes(x=Date, ymin=lower95, ymax=upper95), size=1, alpha=0.2) +
  ylab("Monthly dengue notifications") +
  xlab("Time")+
  scale_color_manual(name = "Model", values = c("red")) +
  scale_x_date(date_breaks = "1 year", 
               labels=date_format("%Y"),
               limits = as.Date(c('2001-01-01','2019-03-01'))) +
  #ggtitle("Dengue local case notifications in Townsville, 2001 - 2017") +
  scale_fill_discrete(
    labels = c(expression(paste("post-wolbachia" = "post-", italic("Wolbachia"))),
               expression(paste("pre-wolbachia"  = "pre-", italic("Wolbachia"))))) +
  annotate(geom = "vline",
           x = Locally_acquired_cases[Locally_acquired_cases$Date=="2014-10-01",1],
           xintercept = Locally_acquired_cases[Locally_acquired_cases$Date=="2014-10-01",1],
           linetype = "dashed") +
  annotate(geom = "text",
           label = expression(paste(italic("Wolbachia "),"introduction")),
           x = Locally_acquired_cases[Locally_acquired_cases$Date=="2014-10-01",1],
           y = c(25),
           angle = 90, 
           vjust = 1) +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(size = 0.4, colour = "grey10"),
        text = element_text(size=12,  family="serif"),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.position = "top",
        strip.background =element_rect(fill="royalblue"),
        strip.text = element_text(size = 10, colour = 'white'))

sammy15

##############################################################################################################

# Computing the corresponding reduction in dengue incidence via Wolbachia intervention 

#dd1 = Locally_acquired_cases[Locally_acquired_cases$Status=="post-wolbachia", ]
#length(dd1[,1])
#dd11 = sum(dd1[,"Cases"])
#dd11
dd2 = optim_solution[optim_solution$Date >= "2014-10-01", ]  # for wolb_rate = 0 (no wolbachia introduction)
length(dd2[,1])
dd22 = sum(dd2[,"Incidence_monthly"])  
dd22

dd3 = optim_solution[optim_solution$Date >= "2014-10-01", ]  # for wolb_rate = 5000 (Wolbachia introduction)
length(dd3[,1])
dd33 = sum(dd3[,"Incidence_monthly"])  
dd33

#dd111 = (dd22 - dd33)/dd22 * 100 # percentage decrease in cases from model values
#dd111

dd111 = (271.6531 - 32.3766)/271.6531
dd111_LCI = (243.955 - 29.6787)/243.955        #95% CI lower 
dd111_UCI = (303.7923 - 35.6059)/303.7923        #95% CI upper
#  The reduction is 89% with CI = 0.88 - 0.91.
##################################################################################################

# Computing and plotting the time-varying reproductive number (R_0).

activation_rate = 1/5.5 #Gubler 1998 from Ndii 
recovery_rate = 1/5        #Ndii et al
biting_rate_u = 0.3
trans_prob_u = 0.2147
biting_rate_w = 0.3 * 0.95
trans_prob_wh = 0.0052
mu = 0.000034         #0.000034,              # Adeshina et al
prob_symp = 1/4          #Kamtchum-Tatuene et al
mu_u = 0.043
mu_w = 0.068
mu_Au = 0.02
mu_Aw = 0.02
rho_u = 13                #my paper
rho_w = 10               #my paper
gamma = 0.95             #walker et al
phi = 0.95               #ant et al
tau_u = 0.11             #walker and hoffman
tau_w = 0.11             #walker et al
mosq_act_rate = 0.1      #chowell et al
mosq_act_wol_rate = 0.1    #chowel et al
sigma = 0
wolb_rate = 5e3
L = 2

#R0 computation and visualization
  alpha1 = (biting_rate_w^2 * trans_prob_wh * mosq_act_wol_rate)/((mu_w + mosq_act_wol_rate)*mu_w)
  alpha2 = (biting_rate_u^2 * trans_prob_u * mosq_act_rate)/((mu_u + mosq_act_rate) * mu)
  alpha = alpha1/alpha2
  R1 = (biting_rate_u^2 * trans_prob_u^2 * mosq_act_rate * activation_rate * optim_solution$Sus_vec / (Total_population)) #numerator uninf mosq
  R11 = ((mu + recovery_rate) * (mu_u + mosq_act_rate) * (mu + activation_rate) *  mu_u)  #denominator uninf mosq
  R2 = (biting_rate_w^2 * trans_prob_wh * trans_prob_u * mosq_act_wol_rate * activation_rate * optim_solution$Sus_wol_vec / (Total_population)) #numerator wolb inf mosq
  R22 = ((mu + recovery_rate) * (mu_w + mosq_act_wol_rate) * (mu + activation_rate) * mu_w)  #denominator wolb inf mosq
  R01 = sqrt(R1/R11)
  R02 = sqrt(R2/R22)
  R0u = R01
  R0 = sqrt(R01^2 + R02^2)
  optim_solution$R_nut = R0
  optim_solution$R_0u = R0u
  optim_solution$wol_prop = optim_solution$Sus_wol_vec/(optim_solution$Sus_vec + optim_solution$Sus_wol_vec)
  optim_solution$R0toRmax = sqrt((1-optim_solution$wol_prop)+ alpha * optim_solution$wol_prop)
  Rvswol_prop = data.frame(optim_solution$wol_prop,optim_solution$R0toRmax)
  
  y.expression <- expression(R(t))
  sammy12 <- ggplot(optim_solution) +
    geom_line(aes(x=Date, y=R_nut), size=1, colour="blue") +
    ylab(y.expression) +
    xlab("Time")+
    scale_x_date(date_breaks = "1 year", 
                 labels=date_format("%Y"),
                 limits = as.Date(c('2001-01-01','2019-03-01'))) +
    #ggtitle("Plot of R0 in the presence of Wolbachia") +
    annotate(geom = "vline",
             x = Locally_acquired_cases[Locally_acquired_cases$Date=="2014-10-01",1],
             xintercept = Locally_acquired_cases[Locally_acquired_cases$Date=="2014-10-01",1],
             linetype = "dashed") +
    annotate(geom = "text",
             label = expression(paste(italic("Wolbachia "),"introduction")),
             x = Locally_acquired_cases[Locally_acquired_cases$Date=="2014-10-01",1],
             y = c(0.4),
             angle = 90,
             vjust = 1) +
    theme(panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid.major = element_line(colour = "grey90"),
          panel.grid.major.x = element_blank(),
          axis.line = element_line(size = 0.4, colour = "grey10"),
          text = element_text(size=12,  family="serif"),
          legend.key = element_rect(fill = "white", colour = "white"),
          legend.position = "top",
          strip.background =element_rect(fill="royalblue"),
          strip.text = element_text(size = 10, colour = 'white'))
  
  sammy12
  
  #####################################################################################################################
  
  ## To account for the changing R(eta)/Rmax with time based on the proportion of Wolbachia infected mosquitoes introduced,
  
  plot(x=optim_solution$wol_prop, y=optim_solution$R0toRmax, ylab = expression(R/R(max)), xlab = expression(paste("Proportion of ", italic("Wolbachia-"),"infected mosquitoes")), col = "red", lwd = 3)
  
  #the mean R0 for pre-wolbachia
  SSS = optim_solution[optim_solution$Status=="Pre-Wolbachia", ]
  mean_of_R_nut_before_wol = mean(SSS$R_nut)
  mean_of_R_nut_before_wol
  
  #the mean R0 for post wolbachia
  TTT = optim_solution[optim_solution$Status=="Post-Wolbachia", ]
  mean_of_R_nut_after_wol = mean(TTT$R_nut)
  mean_of_R_nut_after_wol

############################################################################################################
  
  # Plotting the yearly aggregated dengue case notifications for pre-Wolbachia (from 2001 to 2014) and 
  # post-Wolbachia (from 2014 - 2019) periods for both Imported and locally acquired cases
  
  library(ggsignif)
  library(tidyverse)
  both_cases <- read.csv("dengue_both.csv")
  both_cases
  sam <- ggplot(both_cases, aes(x=relevel(status, "pre-wolbachia"), y=dengue_cases)) +
    geom_boxplot(aes(fill=status)) +
    geom_signif(comparisons = list(c("pre-wolbachia", "post-wolbachia")), 
                map_signif_level=TRUE, tip_length = 0, vjust = 0.2) +
    scale_fill_discrete(
      labels = c(expression(paste("post-wolbachia" = "post-", italic("Wolbachia"), " (2001-2014)")),
                 expression(paste("pre-wolbachia"  = "pre-", italic("Wolbachia"), " (2015-2019)")))) +
    xlab("Status") +
    ylab("dengue cases") +
    theme_bw()
  sam + facet_grid(.~type) + theme(strip.text.x = element_text(size=16, color="black",face="bold"))
  
 