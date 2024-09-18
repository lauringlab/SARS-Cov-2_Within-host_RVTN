library (dplyr)
library (tidyr)
library (gdata)
library (stringr)
library(data.table)

snv_meta <- read.csv ("Results/snv_meta.csv")

## keep snv from individuals with 2 or more samples


snv_count <-  meta %>% group_by (cdc_studyid) %>% distinct (cdc_specid, .keep_all=T) %>% dplyr::count (cdc_studyid) %>% filter (n >1)
snv_2_sample <- filter (snv_meta, cdc_studyid %in% snv_count$cdc_studyid) 
snv_2_sample <- unite (snv_2_sample, hhsubid_pos, c (cdc_studyid, POS), sep="-", remove=F)
snv_2_sample <- unite (snv_2_sample, hhsubid_mutation, c (cdc_studyid, mutation), sep="-", remove=F)


snv_2_sample$date <- as.Date(snv_2_sample$date)

# check for more than 2 alleles
# check for more than 2 alleles 
plus_2_alleles <-snv_2_sample %>%
  distinct (hhsubid_pos, ALT, .keep_all=T) %>% count (hhsubid_pos) %>% filter (n>=2)

snv_2_sample <- filter (snv_2_sample, !(hhsubid_pos %in% plus_2_alleles$hhsubid_pos))
# get data of first sample for each individual

date_first_sample <- snv_2_sample %>% group_by (cdc_studyid) %>% dplyr::slice (which.min (date))%>%
  select (cdc_studyid, date) %>% dplyr::rename (date_first_sample =date)


# separate time of first sample

snv_2_sample <- left_join (snv_2_sample, date_first_sample , by =  "cdc_studyid" ) ## add in date of first sample
#minor_wfabc$SPECDT_1 <- as.Date (minor_wfabc$SPECDT_1, "%Y-%m-%d" ) 
#minor_wfabc$date_first_sample <- as.Date (minor_wfabc$date_first_sample , "%Y-%m-%d" )
snv_2_sample <- mutate (snv_2_sample, days_since_first_sample = date-date_first_sample ) # calculated days since first sample
snv_2_sample$days_since_first_sample <- as.numeric (snv_2_sample$days_since_first_sample)
snv_2_sample <- mutate (snv_2_sample, gen_since_first_sample = 3* days_since_first_sample) # calculated generations since first sample (assume 6hr generation time)


### add in allele count for each mutation - depth * allele frequency
snv_2_sample <- mutate (snv_2_sample , Total_DP  = TOTAL_DP_1+ TOTAL_DP_2, Allele_count = round (avg_freq* Total_DP))

### get frequency for each timepoint (long form)###


#depth df
minor_wfabc_DP <- snv_2_sample %>% ungroup ()  %>% select ( hhsubid_pos, Total_DP, gen_since_first_sample)
minor_wfabc_DP <- minor_wfabc_DP %>% group_by (hhsubid_pos) %>% spread (gen_since_first_sample, Total_DP)

# allele count df
minor_wfabc_allele_count <- snv_2_sample %>% ungroup ()  %>% select ( hhsubid_pos,  Allele_count, gen_since_first_sample)
minor_wfabc_allele_count <- minor_wfabc_allele_count %>% group_by (hhsubid_pos) %>% spread (gen_since_first_sample, Allele_count)


## separate individuals/mutations with different timepoints

#minor_wfabc_allele_count <-select (minor_wfabc_allele_count,-REGION, -POS, -mutation_type, -hhsubid,-PCR_RESULT_1, -onsetdt, -age.new, -sex, -p_seasvx, -Year, -hhid, -reference, -hh_onsetdt, -initial_mutation)
minor_wfabc_ac_no_names <-minor_wfabc_allele_count %>% ungroup ( ) %>%  select (-hhsubid_pos)

#minor_wfabc_DP <- select (minor_wfabc_DP,-REGION, -POS, -mutation_type, -hhsubid,-PCR_RESULT_1, -onsetdt, -age.new, -sex, -p_seasvx, -Year, -hhid, -reference, -hh_onsetdt, -initial_mutation)
minor_wfabc_DP_no_names <- minor_wfabc_DP %>% ungroup ( ) %>%  select (-hhsubid_pos)





# get combinations of time points  
make_combinations <- function(x) {
  
  l <- length(x)
  mylist <- lapply(2:l, function(y) {
    combn(x, y, simplify = FALSE)
  })
  mylist
  unlist(mylist, recursive = FALSE)
}

results <- make_combinations(colnames(minor_wfabc_ac_no_names))

for (i in 1:length(results)) {
  NA_names <- minor_wfabc_allele_count %>% ungroup %>% select ( -results [[i]], -hhsubid_pos)
  NA_names <- list(colnames(NA_names)) %>% unlist (.)
  
  
  if (length (NA_names)==0) {
    ac.df <- minor_wfabc_allele_count %>% filter_at(c(results[[i]]), all_vars(!is.na(.))) %>% 
      select (hhsubid_pos, results[[i]]) 
    dp.df <- minor_wfabc_DP %>% filter_at(c(results[[i]]), all_vars(!is.na(.))) %>% 
      select (hhsubid_pos, results[[i]]) 
  }else {
    ac.df <- minor_wfabc_allele_count %>% filter_at(c(results[[i]]), all_vars(!is.na(.))) %>% 
      filter_at(c(NA_names), all_vars(is.na(.)))%>% 
      select (hhsubid_pos, results[[i]]) %>% ungroup ()
    
    dp.df <- minor_wfabc_DP %>% filter_at(c(results[[i]]), all_vars(!is.na(.))) %>% 
      filter_at(c(NA_names), all_vars(is.na(.)))%>% 
      select (hhsubid_pos, results[[i]])  %>% ungroup ()
  }
  
  if (dim(ac.df)[1] != 0) {
    wfabc_ac.df <- ac.df %>% select(2:ncol(.) ) # get frequency data frame
    depth.df <- dp.df %>% select(2:ncol(.) )  # get coverage data frame
    names <-  ac.df %>% select (1) # get names of individuals and use to make csv file
    generation_col <- colnames (wfabc_ac.df)  # get column names to use in  output files
    generation_name <- paste(generation_col,collapse="_")
    
    
    
    file_name <- paste ("Results/WFABC3/", generation_name, "_names.csv", sep="")
    write.csv(names, file_name, row.names = F, quote=F)
    
    header_info <- c (nrow (wfabc_ac.df), ncol (wfabc_ac.df))# get header info to add to top of csv file
    header_file <- paste ("Results/WFABC3/", generation_name, "_header.txt", sep="")
    write(header_info, header_file)
    
    depth.df <- as.matrix(depth.df)
    wfabc_ac.df <- as.matrix(wfabc_ac.df)
    wfabc_final <- interleave (depth.df, wfabc_ac.df) 
    wfabc_file <- paste ("Results/WFABC3/" , generation_name, ".txt", sep="")
    write.csv (wfabc_final, wfabc_file, row.names=F, quote=F)
  }
}

Ne_ac.df <-data.frame(Col1 = character(),
                      Col2 = character(),
                      Col3 = character ())

Ne_ac.df <- rename (Ne_ac.df, "0" = 1, "3" =2, "6"=3)

Ne_dp.df <-data.frame(Col1 = character(),
                      Col2 = character(),
                      Col3 = character ())
Ne_dp.df <- rename (Ne_dp.df, "0" = 1, "3" =2, "6" =3)


# make input file using allele counts and actual depth 
for (row in 1:nrow(minor_wfabc_ac_no_names)){
  row_wfabc <- minor_wfabc_ac_no_names [row, ]
  row_wfabc <-row_wfabc %>% select_if(~ !any(is.na(.)))
  
  
  if (ncol (row_wfabc) > 2) {
    my_colnames <- as.numeric (colnames(row_wfabc))
    if (my_colnames[2]-my_colnames [1] ==3){
      if (my_colnames[3]-my_colnames [2] ==3){
      row_wfabc <- row_wfabc %>% select ( 1,2, 3) %>% rename ("0" = 1, "3" =2, "6"=3)
      row_wfabc [,1]<- as.character(row_wfabc [,1])
      row_wfabc [,2]<- as.character(row_wfabc [,2])
      row_wfabc [,3]<- as.character(row_wfabc [,3])
      Ne_ac.df <- bind_rows ( Ne_ac.df ,row_wfabc )
    }
  }
  }

for (row in 1:nrow(minor_wfabc_DP_no_names)){
  row_depth <- minor_wfabc_DP_no_names [row, ]
  row_depth <-row_depth %>% select_if(~ !any(is.na(.)))  
  
  if (ncol (row_depth) > 2) {
    my_colnames <- as.numeric (colnames(row_depth))
    if (my_colnames[2]-my_colnames [1] ==3){
      if (my_colnames[3]-my_colnames [2] ==3){
        row_depth <- row_depth %>% select ( 1,2, 3) %>% rename ("0" = 1, "3" =2, "6"=3)
        row_depth [,1]<- as.character(row_depth [,1])
        row_depth [,2]<- as.character(row_depth [,2])
        row_depth [,3]<- as.character(row_depth [,3])
        Ne_dp.df <- bind_rows ( Ne_dp.df ,row_depth )
      }
    }
  }

  
Ne_final_ac <- interleave (Ne_dp.df, Ne_ac.df )
  Ne_file <- paste ("Results/WFABC3/Ne_sample_size_input_3gen.txt", sep="")
  write.csv (Ne_final_ac , Ne_file , row.names=F, quote=F)
  