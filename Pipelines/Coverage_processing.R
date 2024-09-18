###Coverage
library (stringr)
library (dplyr)
library (tidyr)

avg_cov <- read.table (snakemake@input[[1]], header=T)
 
 # get samples that pass coverage filter
avg_cov <- separate (avg_cov, sample, c ("sample", "replicate"), "_") # split sample column
avg_cov <- spread (avg_cov, replicate, mean)
avg_cov<- rename (avg_cov, rep_2= 3, rep_1= 2)
avg_cov <- mutate (avg_cov, Pass_coverage = ifelse (rep_1 >= 1000 & rep_2 >=1000, "Pass", "Fail")) # filter by coverage 
avg_cov_pass <- filter (avg_cov, Pass_coverage == "Pass")

#add in subject id info 

avg_cov_pass <- mutate (avg_cov_pass, hhsubid = sample)
avg_cov_pass$hhsubid <- avg_cov_pass$hhsubid %>% str_sub ( 1, 6)

write.csv (avg_cov_pass, snakemake@output[[2]], row.names=F, quote=F)
# filter for first sample per individual
avg_cov_pass_first<- avg_cov_pass %>% group_by (hhsubid) %>% filter (sample == min(sample))

write.csv (avg_cov_pass_first,snakemake@output[[1]], row.names=F, quote=F)

