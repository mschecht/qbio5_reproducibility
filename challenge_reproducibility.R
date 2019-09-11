#----------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------ Functions -----------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------

#Create function for computing dispersion coefficient: var/mean
compute_cd <- function(x) {
  if(any(is.na(x))) {
    message("Input vectors contains NAs. Will ignore.")
    x <- na.omit(x)
  }
  y <- rep((seq_along(x)-1), times = x)
  return(var(y)/mean(y))
}

#Create function for computing mean
#Could be collapsed with the previous function
compute_empirical_lambda <- function(x) {
  if(any(is.na(x))) {
    message("Input vectors contains NAs. Will ignore.")
    x <- na.omit(x)
  }
  y <- rep((seq_along(x)-1), times = x)
  return(mean(y))
}

#----------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------- Analysis -----------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(RMKdiscrete)

#Load tables and rename columns to match
df1 <- read.csv(file = "reproducibility/data/cole_arthropod_data_1946.csv") %>% rename("k_number" = "k_number_of_arthropods")
df2 <- read.csv(file = "reproducibility/data/mitchell_weevil_egg_data_1975.csv") %>% rename("k_number" = "k_number_of_eggs")

#Join the two tables table
df <- full_join(df1, df2, by = "k_number")

#Compute dispersion coefficients for each species
#Sowbugs seem to deviate from Poisson distribution (overdispersed)
dispersion_coefficient <- df %>% select(-1) %>% apply(2, compute_cd)
empirical_l <- df %>% select(-1) %>% apply(2, compute_empirical_lambda)

#Compute number of observations for each species
number_of_observations <- df %>% select(-1) %>% colSums(na.rm = TRUE)

#Estimate expected distribution assuming Poisson distribution
#Adding more columns (k_number, distribution) to prepare for plotting
theoretical_pmf <- empirical_l %>% sapply(function(x) {dpois(x = df$k_number, lambda = x)})
theoretical_poisson_distribution <- theoretical_pmf %>% 
  apply(1, function(x) x*number_of_observations) %>% 
  t %>% as.data.frame %>% 
  mutate("k_number" = df$k_number,
         "distribution" = "Poisson")

#Compute LGP distribution parameters based on hints
#empirical_l is the mean

l1_spiders <- empirical_l[1]
l2_spiders <- 0

l2_sowbugs <- 0.53214
l1_sowbugs <- empirical_l[2]*(1-l2_sowbugs)






data_weevil <- df[,3]/sum(df[,3])
err_function <- function(params) {
  err <- sum((dLGP(df$k_number, theta = params[1], lambda = params[2]) - data_weevil)^2)
}

estimated_params <- optim(par = c(empirical_l[3], 0.1), fn = err_function)$par

plot(df$k_number, df[,3], pch = 16)
lines(dLGP(df$k_number, theta = estimated_params$par[1], lambda = estimated_params$par[2])*number_of_observations[3])







#Prepare data frame for plotting
df_melted <- df %>% reshape2::melt(id = "k_number")
theo_melted <- theoretical_poisson_distribution %>% reshape2::melt(id = c("k_number", "distribution"))

#Plot
ggplot() +
  geom_point(data = df_melted, 
             aes(x = k_number, y = value)) + 
  geom_line(data = theo_melted,
            aes(x = k_number, y = value, 
                group = distribution, col = distribution)) +
  facet_grid(variable~., scales = "free")
