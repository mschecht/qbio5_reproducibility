# Set up environment
#-------------------

# Make list of packages
list.of.packages <- c("tidyverse")

# Subset packages that are not installed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

# Install list of uninstalled packages
if(length(new.packages)) install.packages(new.packages)

# Library all functions
lapply(list.of.packages, library, character.only = TRUE)

# Load data
arthropods <- read_csv("cole_arthropod_data_1946.csv")

weevil <- read_csv("mitchell_weevil_egg_data_1975.csv")

load("pvalueData_PNAS.rda")


# 1) Plot the Poisson distribution with the same mean as the spider counts, along with the data

#Calculate the mean number of total number of observed spiders
total_observation_count_spiders <- sum(arthropods[ ,"C_count_of_boards_with_k_spiders"])
spider_count_mean <- (sum(arthropods$k_number_of_arthropods*arthropods$C_count_of_boards_with_k_spiders))/total_observation_count_spiders

rpois(n = 1000, lambda = spider_count_mean) %>% hist(main = paste("Histogram of" , "Spider Counts"),col="lightblue1", border="blue")

spider_pois <- rpois(n = 1000, lambda = spider_count_mean)



# 2) Plot the Poisson distribution with the same mean as the sowbug counts, along with the data
total_observation_count_sowbugs <- sum(arthropods[ ,"C_count_of_boards_with_k_sowbugs"])
sowbugs_count_mean <- (sum(arthropods$k_number_of_arthropods*arthropods$C_count_of_boards_with_k_sowbugs))/total_observation_count_sowbugs

rpois(n = 1000, lambda = sowbugs_count_mean) %>% hist(main = paste("Histogram of" , "Sowbugs Counts"), col="lightblue1", border="blue")

sowbug_pois <- rpois(n = 1000, lambda = sowbugs_count_mean)



# 3) Plot the Poisson distribution with the same mean as the weevil egg counts, along with the data
total_observation_count_weevil <- sum(weevil[ ,"C_count_of_beans_with_k_eggs"])
weevil_count_mean <- (sum(weevil$k_number_of_eggs*weevil$C_count_of_beans_with_k_eggs))/total_observation_count_weevil

rpois(n = 1000, lambda = weevil_count_mean) %>% hist(main = paste("Histogram of" , "Weevil Eggs Counts"), col="lightblue1", border="blue")

weevil_pois <- rpois(n = 1000, lambda = weevil_count_mean)



# 4) Add a curve to Plot 1) showing the LGP distribution with the parameter hint below for the spider counts

spider_pois %>%
  as_tibble() %>%
  ggplot(aes(x = value)) +
  geom_line(aes(y = ..density.., colour = 'Empirical'), stat = 'density') +
  geom_histogram(aes(y = ..density..), alpha = 0.4)

# 5) Add a curve to Plot 2) showing the LGP distribution with the parameter hint below for the sowbug counts

install.packages("RMKdiscrete")
library(RMKdiscrete)
lambda2 <- 0.53214
lambda1 <- sowbugs_count_mean*(1 - lambda2)
x <- c(0:17)
sowbugsmodel<-dLGP(x, lambda1, lambda2)


# 6) Add a curve to Plot 3) showing the LGP distribution with the parameter hint below for the weevil eg counts

