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

#Create a function that estimates LGP parameters using SSE minimization
estimate_LGP_params <- function(data) {
  err_function <- function(params) {
    expected <- dLGP(df$k_number, theta = params[1], lambda = params[2])
    observed <- data
    return(sum((observed-expected)^2))
  }
  return(optim(par = c(1, -0.5), fn = err_function)$par)
}

#----------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------- Analysis -----------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(RMKdiscrete)

#Load tables and rename columns to match
df1 <- read.csv(file = "cole_arthropod_data_1946.csv") %>% rename("k_number" = "k_number_of_arthropods")
df2 <- read.csv(file = "mitchell_weevil_egg_data_1975.csv") %>% rename("k_number" = "k_number_of_eggs")

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
  set_colnames(c("Spiders", "Sowbugs", "Weevil_eggs")) %>% 
  mutate("k_number" = df$k_number,
         "Model" = "Poisson")

#Compute LGP distribution parameters based on hints
#empirical_l is the mean
l1_spiders <- empirical_l[1]
l2_spiders <- 0

l2_sowbugs <- 0.53214
l1_sowbugs <- empirical_l[2]*(1-l2_sowbugs)

#Estimating weevil parameters minimizing SSE
lambdas_weevil <- estimate_LGP_params(na.omit(df[,4]/number_of_observations[3]))
params_table <- data.frame("Spiders" = c(l1_spiders, l2_spiders),
                           "Sowbugs" = c(l1_sowbugs, l2_sowbugs),
                           "Weevil_eggs" = lambdas_weevil)

theoretical_LGP_pmf <- params_table %>% 
  apply(2, function(x) {
    dLGP(df$k_number, theta = x[1], lambda = x[2])
  })

theoretical_LGP_distribution <- theoretical_LGP_pmf %>%
  apply(1, function(x) {
    x*number_of_observations
  })%>% t %>% as.data.frame %>% 
  mutate("k_number" = df$k_number,
         "Model" = "LGP")

theoretical_distributions <- rbind(theoretical_poisson_distribution,
                                   theoretical_LGP_distribution)

#Prepare data frame for plotting
df_melted <- df %>%
  set_colnames(c("k_number", "Spiders", "Sowbugs", "Weevil eggs")) %>% 
  reshape2::melt(id = "k_number")
theo_melted <- theoretical_distributions %>% 
  set_colnames(c("Spiders", "Sowbugs", "Weevil eggs", "k_number", "Model")) %>%
  reshape2::melt(id = c("k_number", "Model"))

#Plot
ggplot() +
  geom_bar(data = df_melted, 
           aes(x = k_number, y = as.numeric(value)),
           fill = "grey80",
           col = "grey40",
           stat = "identity") + 
  geom_line(data = theo_melted,
            aes(x = as.numeric(as.character(k_number)), y = as.numeric(value), 
                group = Model, col = Model,
                linetype = Model)) +
  geom_point(data = theo_melted,
             aes(x = as.numeric(as.character(k_number)), y = as.numeric(value), 
                 group = Model, col = Model,
                 shape = Model)) +
  facet_grid(variable~., scales = "free") +
  ggtitle("Observed and expected distribution") +
  xlab("Counts") +
  ylab("Frequency") +
  theme_bw()

#Compute log-likelihood of each model for all three arthropods
LL_poisson <- sapply(1:3, function(x) {
  rep(theoretical_pmf[1:length(na.omit(df[,x+1])),x], times = na.omit(df[,x+1])) %>% log %>% sum
})

LL_LGP <- sapply(1:3, function(x) {
  rep(theoretical_LGP_pmf[1:length(na.omit(df[,x+1])),x], times = na.omit(df[,x+1])) %>% log %>% sum
})

#Compute LRT statistic, should follow a Chi-Square distribution with 2-1=1 degree of freedom
LRT_statistic <- -2*(LL_poisson-LL_LGP)
pvalue <- 1-pchisq(q = LRT_statistic, df = 1)
#LGP model is better than Poisson for Sowbugs and Weevil eggs (adjusted p < 0.001)