library (dplyr)

###### Filter SNV and add metadata #######

###Coverage
avg_cov <- read.table("Results/AvgCoverage.all", stringsAsFactors = FALSE, header=T)
avg_cov <- separate (avg_cov, sample, c ("sample", "replicate"), "_") # split sample column
avg_cov <- spread (avg_cov, replicate, mean)
avg_cov<- rename (avg_cov, rep_2= 3, rep_1= 2)
avg_cov <- mutate (avg_cov, Pass_coverage = ifelse (rep_1 >= 1000 & rep_2 >=1000, "Pass", "Fail")) # filter by coverage 

write.csv (avg_cov, "Results/Pass_coverage.csv", quote=F, row.names=F)

avg_cov_pass <- filter (avg_cov, Pass_coverage == "Pass")

snv <- read.table ("Results/all_variants_filtered", header=T)
snv <- filter (snv, sample %in% avg_cov_pass$sample )

# filter out T11075C 

snv <- filter (snv, !(mutation == "T11075C"))

# add metadata 

meta <- read.csv ("Metadata/meta_1000.csv")
snv_meta <- left_join (snv, meta, by == "sample")

# add in amino acid position

snv_meta <- snv_meta %>% mutate (AA_Position = ifelse (GFF_FEATURE == "ORF1a", POS-265, ifelse(
  GFF_FEATURE == "ORF1b", POS-13467, ifelse(GFF_FEATURE == "S", POS-21562, ifelse (
      GFF_FEATURE == "ORF3a", POS-25392, ifelse (GFF_FEATURE == "E", POS-26244, ifelse (
          GFF_FEATURE == "M", POS-26522, ifelse (GFF_FEATURE== "ORF6", POS-27201, ifelse (
            GFF_FEATURE == "ORF7a", POS-27393, ifelse (GFF_FEATURE == "ORF8", POS-27893, ifelse (
              GFF_FEATURE == "N" , POS-28273, "NA"
            ))
          )
            
          )
        )
      ) 
    )
  )
) )
)

snv_meta <- snv_meta %>% mutate (AA_Position = ceiling (as.numeric (AA_Position)/3))

write.csv ("Results/snv_meta.csv", row.names=F, quote=F)