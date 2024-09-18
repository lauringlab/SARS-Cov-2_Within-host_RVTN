library(dplyr)
library (stringr)
library (tibble)

variants_analyze  <- read.table(snakemake@input[[1]],stringsAsFactors=F,comment.char = '#', header=T)
names(variants_analyze)[10:29] <- c('REF_DP_1', 'REF_RV_1', 'REF_QUAL_1', 'ALT_DP_1', 'ALT_RV_1', 'ALT_QUAL_1', 'ALT_FREQ_1', 'TOTAL_DP_1', 'PVAL_1', 'PASS_1','REF_DP_2', 'REF_RV_2', 'REF_QUAL_2', 'ALT_DP_2', 'ALT_RV_2', 'ALT_QUAL_2', 'ALT_FREQ_2', 'TOTAL_DP_2', 'PVAL_2', 'PASS_2')


variants_filtered <- variants_analyze %>% 
  filter(ALT_FREQ_1 < 0.995) %>%
  filter(ALT_FREQ_2 < 0.995) %>%
  mutate(mutation = paste0(REF, POS, ALT)) %>%
  filter(!str_detect(string = ALT, pattern = "\\+") & !str_detect(string = ALT, pattern = "\\-")) %>%
  filter(TOTAL_DP_1 > 400 & PVAL_1 < 0.00001 & ALT_QUAL_1 > 35  & TOTAL_DP_1 >400 & PVAL_2 < 0.00001 & ALT_QUAL_2 > 35 )

# add in mutation type 
# add in average freq
variants_filtered <- variants_filtered %>% mutate(avg_freq = (ALT_FREQ_1 + ALT_FREQ_2)/2) %>%
      mutate (mutation_type = ifelse (REF_AA==ALT_AA, "Syn", "Non"))

#mask sites
mask_sites <- read.table (snakemake@input[[2]], header =T)
variants_filtered <- filter (variants_filtered, !POS %in% mask_sites$POS)

## add sample name
filename <-snakemake@input[[1]]
filename2 <- sub(".merged.tsv", "", filename)

filename_vec <- strsplit(filename2, split = "/")[[1]]

variants_filtered <- dplyr::mutate (variants_filtered, sample = filename_vec[4])

write.table (variants_filtered, snakemake@output[[1]], quote=F, row.names=F)
