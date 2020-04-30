### Bayesian Covid Tuto script

## Libraries
library(rmarkdown)
library(ggplot2)
library(ggridges)
library(utils)
library(rjags)
library(EnvStats)

## Variables and Data Frames

# File saving path
path <- 'D:/Jonathan/Lescot/Projets/GT_Statistiques/Bayesianism/Bayesian_Statistics/'

# Download data
data <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")
data$dateRep <- as.Date(data$dateRep, format = "%d/%m/%Y")

# German data subset
german_data <- subset(data, data$countriesAndTerritories=="Germany" & data$dateRep >="2020-03-12")
german_data <- german_data[order(german_data$dateRep),]

## Empty df for fatality rate update tracing
rate_df = data.frame(day=integer(), Q25=double(), Q50=double(), Q97=double())

## Bayesian model preparation
# select the data table used to feed the process
data_table <- german_data

# Set initial Beta distribution coeffs correponding to a uniform distribution
coeff_alpha = 1
coeff_beta = 1

## Bayesian updating 
# For loop that iterates on the german-data table and proceed to the modelling steps
# update the beta distribution to set the daily new priors 
# and plot everyday death rate parameter probability distribution
for (row in 1:nrow(data_table)) {
  
  # retrieve the cases and deaths daily quantities
  cases <- data_table[row, "cases"]
  deaths  <- data_table[row, "deaths"]
  
  # DEFINE the model
  covid_model <- "model{
              # Likelihood model for X
              X ~ dbin(p, n)

              # Prior model for p
              p ~ dbeta(a, b)
              }"
  
  # Compile the model
  covid_jags <- jags.model(textConnection(covid_model), 
                           data=list(a=coeff_alpha, 
                                     b=coeff_beta, 
                                     X=deaths, 
                                     n=cases),
                           inits=list(.RNG.name="base::Wichmann-Hill",
                                      .RNG.seed=100))
  
  # SIMULATE the posterior
  covid_sim <- coda.samples(model=covid_jags, 
                            variable.names=c("p"), 
                            n.iter= 1000000)
  
  # Create a df with covid_sim chain
  covid_chain<- data.frame(covid_sim[[1]], iter = 1:1000000)

  
  # Approximate the beta parameters based on the chain quantiles 
  # (rriskDistributions package function)
  coeffs <-  ebeta(covid_chain$p, method = "mle")
  coeff_alpha <- coeffs$parameters[[1]]
  coeff_beta <- coeffs$parameters[[2]]
  

  # Append posterior distribution 2.4%, 50% & 97.5% quantile values to
  # to fatality rate estimation table
  day_number=row
  chain_quant <- quantile(covid_chain$p, c(0.025, 0.5, 0.975))
  Q2 <- chain_quant[[1]]
  Q5 <- chain_quant[[2]]
  Q9 <- chain_quant[[3]]
  quantile_row <- cbind(day_number, Q2, Q5, Q9)
  rate_df <- rbind(rate_df, quantile_row)
  
  # PLOT the simulated posterior
  day <- paste("Jour ", row)
  cas <- paste(cases, "malades")
  morts <- paste(deaths, "morts")
  a <- paste("alpha : ", coeff_alpha)
  b <- paste("beta : ", coeff_beta)
  plot <- ggplot(covid_chain, aes(x = p)) + 
              geom_density(fill="green", alpha=0.4, position = "stack")+
              xlim(0, 0.05) +
              ylim(0, 6000) +
              xlab('Taux de létalité') +
              ylab('Densité') +
              ggtitle(paste0(day, ", ",cas, ", ",morts, ", ",a, ", ",b))
  print(plot)
  file_name <- paste0(path, day, ".png")
  ggsave(filename = file_name)
  }


ggplot(rate_df, aes(x=rate_df$day, y=rate_df$Q5))+
  geom_line(aes(x=rate_df$day, y=rate_df$Q9), color="grey")+
  geom_line(aes(x=rate_df$day, y=rate_df$Q2), color="grey")+
  geom_ribbon(data=rate_df, 
              aes(ymin=rate_df$Q2,ymax=rate_df$Q9), 
                  fill="lightgrey", alpha="0.7")+
  geom_line()+
  ylim(0, 0.036) +
  xlab('Jours') +
  ylab('Taux de létalité') +
  ggtitle('Évolution de la distribution posterior du taux de létalité')
ggsave(filename = paste0(path, 'rate_evolution.png'))
