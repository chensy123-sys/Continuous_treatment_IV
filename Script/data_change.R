library(dplyr);library(tibble);library(tidyr);library(readr)
my_data <- read.table("jtpa_han.tab", header = TRUE, sep = "\t") %>% 
  select(sex, n_hs2, edu, prevearn) %>%
  rename(L = sex, A = edu, Z = n_hs2, Y= prevearn)
head(my_data)
write.csv(my_data,'jpta_han.csv')