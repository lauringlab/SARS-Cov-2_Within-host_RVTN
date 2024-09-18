library (dplyr)
library (tidyr)
library (ggplot2)
library (cowplot)

## preprocessing, stats and figures for number of isnv per sample and isnv frequency

snv_meta <- read.csv("Results/snv_meta.csv")
meta <- read.csv ("Metadata/meta_sample_1000.csv")

# Number of iSNV per sample
no_isnv_samples <- meta %>% filter (! (cdc_specid %in% snv_meta$cdc_specid))%>% 
  select (cdc_specid, age_at_enrollment, vaccinated, clade, Days_post_symp, age_class) %>%
  mutate (n = "0")
no_isnv_samples$vaccinated <- as.factor (no_isnv_samples$vaccinated)   
no_isnv_samples$n <- as.numeric(no_isnv_samples$n)   

isnv_count <- count (snv_meta, cdc_specid, age_at_enrollment, vaccinated, clade, Days_post_symp, age_class)
isnv_count$vaccinated <- as.factor (isnv_count$vaccinated)   
isnv_count <- bind_rows (isnv_count, no_isnv_samples)


## statistics 

#Stats for number of isnv per sample
wilcox.test (n ~ age_class, data=isnv_count)

wilcox.test (n ~ vaccinated, data=isnv_count)

kruskal.test (n ~ clade, data=isnv_count)
dunnTest(n ~ clade, data=isnv_count)

kruskal.test (n ~ Days_post_symp, data=isnv_count)
dunnTest(n ~ Days_post_symp, data=isnv_count)


# Stats for isnv frequency
wilcox.test (avg_freq ~ age_class, data=snv_meta)
wilcox.test (avg_freq ~ vaccinated, data=snv_meta)
wilcox.test (avg_freq ~ mutation_type, data=snv_meta)
kruskal.test (avg_freq ~ clade, data=snv_meta)
kruskal.test (avg_freq ~ Days_post_symp, data=snv_meta)
dunnTest(n ~ Days_post_symp, data=isnv_count)

# Stats for Number of samples per day post symptom onset 
kruskal.test (Days_post_symp ~ clade, data=meta)# not significant
wilcox.test (Days_post_symp ~ vaccinated, data=meta) #not significant
wilcox.test (Days_post_symp ~ age_class, data=meta) # not significant

### Plots

## Figure 2

#isnv freq
snv_freq_plot <- ggplot (snv_meta, aes(avg_freq))+ geom_histogram(binwidth =0.01,  boundary=0)+
  theme_cowplot (12) + xlab ("iSNV frequency") + ylab ("Number of iSNV")+
  theme(axis.title= element_text(size = 11), axis.text = element_text(size =9)) 

#isnv number

snv_count2 <- count (isnv_count, n)
snv_count2 [nrow(snv_count2) + 1,] = c ("0", "135")# number of samples without isnv

snv_count_plot <-ggplot (isnv_count_2, aes (n, nn))+ geom_col()+
  theme_cowplot (12) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n =6 )) +
  xlab ("Number of iSNV per specimen") + ylab ("Number of specimens")+
  theme(axis.title= element_text(size = 11), axis.text = element_text(size =9)) 



#isnv number by covariates
#by clade
snv_count_clade <-ggplot (isnv_count, aes (n))+ geom_density(aes (fill=clade), alpha =0.75)+
  theme_cowplot (12) + xlim (0,12)+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_fill_manual(values=c( "#FDAFA7","#0B2A5C", "#55794A"))+ 
  xlab ("Number of iSNV per sample") + ylab ("Density")+
  theme (legend.position = "none", axis.title = element_text(size=11), axis.text = element_text(size =9))  

#by vaccine status
isnv_count_vaccine <-ggplot (isnv_count, aes (n))+ geom_density(aes (fill=vaccinated), alpha =0.7)+
  theme_cowplot (12) + xlim (0,12)+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_fill_manual(values=c( "#FDAFA7","#0B2A5C"))+ 
  xlab ("Number of iSNV per sample") + ylab ("Density")+
  theme (legend.position = "none", axis.title = element_text(size=11), axis.text = element_text(size =9))  

#by age
isnv_count_age <-ggplot (isnv_count, aes (n))+ geom_density(aes (fill=age_class), alpha =0.7)+
  theme_cowplot (12) + xlim (0,12)+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_fill_manual(values=c( "#FDAFA7","#0B2A5C"))+ 
  xlab ("Number of iSNV per sample") + ylab ("Density")+
  theme (legend.position = "none", axis.title = element_text(size=11), axis.text = element_text(size =9))  

#by day post infection
snv_count_day_plot <-  ggplot (data=subset( isnv_count, !is.na(Days_post_symp)), aes (as.factor (Days_post_symp), n)) + 
  geom_jitter (alpha =0.4, shape =16, size=1) + 
  geom_crossbar(data=time_freq_sum, aes(ymin = n, ymax = n),
                col="red", width = .5, size =0.4)+
  scale_y_continuous(breaks = c (2,4,6,8,10,12))+
  theme_cowplot (12)+
  theme(axis.title= element_text(size = 11), axis.text = element_text(size =9)) +
  xlab ("Days Post Symptom Onset") + ylab ("Number of iSNV")


#combine plots
top_Snv_plot <- plot_grid (snv_freq_plot, snv_count_plot, labels = c ("A", "B"), vjust =0.5)  
density_plot_isnv_count <- plot_grid (isnv_count_vaccine, NULL, isnv_count_age, NULL, isnv_count_clade, labels = c ("C","","D","","E"), 
                                      nrow=1, rel_widths = c(1,.1,1,.1,1), vjust = 0.5)

isnv_count_day_plot <- plot_grid (snv_count_day_plot, vjust=0.5, labels = "F")
plot_grid (NULL,top_Snv_plot, NULL, density_plot_isnv_count, NULL, isnv_count_day_plot, 
           rel_heights = c (.1,1,.1,1,.1,1.3), nrow=6)


##Figure 2S

## isnv frequency by covariates
#Histogram of isnv frequency by mutation type
mutation_type_freq_plot <-  ggplot (data=subset(snv_meta, !is.na(mutation_type)), aes (avg_freq))+ 
  geom_density(aes (group= mutation_type, fill= mutation_type), alpha =0.5)+ 
  theme_cowplot (12) + 
  ylab ("Density") + xlab ("iSNV Frequency") + 
  scale_fill_manual(values=c( "#0B2A5C", "#FDAFA7")) +
  theme(legend.position = "none", axis.title= element_text(size = 11), 
        axis.text = element_text(size =9)) 

#Histogram of isnv frequency by vaccine status 
vaccinated_type_freq_plot <-  ggplot (data=subset(snv_meta, !is.na(vaccinated)), aes (avg_freq))+
  geom_density(aes (group= vaccinated, fill= as.factor (vaccinated)), alpha =0.5)+ 
  theme_cowplot (12) + 
  ylab ("Density") + xlab ("iSNV Frequency")+ 
  scale_fill_manual(values=c( "#0B2A5C", "#FDAFA7")) +
  theme(legend.position = "none", axis.title= element_text(size = 11), 
        axis.text = element_text(size =9)) 

#Histogram of isnv frequency by Age
age_freq_plot <- ggplot (snv_meta, aes (avg_freq))+ 
  geom_density(aes (fill= age_class), alpha =0.7) +
  theme_cowplot (12) + 
  ylab ("Density") + xlab ("iSNV Frequency") +
  scale_fill_manual(values=c( "#0B2A5C", "#FDAFA7")) +
  theme(legend.position = "none", axis.title= element_text(size = 11), 
        axis.text = element_text(size =9)) 

#Histogram of isnv frequency by Clade
clade_type_freq_plot <- ggplot (snv_meta, aes (avg_freq))+ geom_density(fill = "grey")+ facet_grid(cols  = vars (clade)) +
  theme_cowplot (12) + 
  ylab ("Density") + xlab ("iSNV Frequency") +
  theme(axis.title= element_text(size = 11), axis.text = element_text(size =9)) 


#Jitter plot of days post symptom onset and isnv frequency
  #Get average frequency by day
time_freq_sum <- snv_meta %>% group_by (Days_post_symp) %>% summarize (avg_freq = mean(avg_freq))


days_post_sym_freq_plot <- ggplot (snv_meta, aes ( Days_post_symp, avg_freq)) + 
  geom_jitter (alpha =0.4, shape =16, size=1) + 
  geom_crossbar(data=time_freq_sum, aes(ymin = avg_freq, ymax = avg_freq),
                col="red", width = .5, size =0.4) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6))+
  theme_cowplot (12)+
  theme(axis.title= element_text(size = 11), axis.text = element_text(size =9)) +
  xlab ("Days Post Symptom Onset") + ylab ("iSNV Frequency")

# Combine plots to make Figure S2
top_plot <- plot_grid (mutation_type_freq_plot , vaccinated_type_freq_plot, age_freq_plot, labels = c ("A","B", "C"), vjust =0.5, nrow =1 )

plot_grid (NULL, top_plot, NULL,  clade_type_freq_plot, NULL, days_post_sym_freq_plot, nrow = 6, 
           rel_heights = c (0.1, 1, 0.1, 0.8, 0.1, 1), labels = c ("", "", "", "D", "", "E"),vjust = 0.5)


## Plot Number of samples per day by covariate
clade_specimen_per_day_plot <- ggplot (meta, aes ( Days_post_symp))+ 
  geom_bar (aes (fill=clade), position="dodge")+
  theme_cowplot (12) + xlab ("Days Post Symptom Onset") + ylab ("Number of Specimens") +
  scale_fill_manual(values=c( "#FDAFA7","#0B2A5C",  "#55794A"))+ 
  theme (legend.position = "none", axis.title = element_text(size=11), axis.text = element_text(size =9))  


vaccine_specimen_per_day_plot <- ggplot (meta, aes ( Days_post_symp))+ 
  geom_bar (aes (fill= as.factor (vaccinated)), position="dodge") +
  theme_cowplot (12) + xlab ("Days Post Symptom Onset") + ylab ("Number of Specimens")+
  scale_fill_manual(values=c( "#FDAFA7","#0B2A5C"))+ 
  theme (legend.position = "none",  axis.title = element_text(size=11), axis.text = element_text(size =9))  


age_specimen_per_day_plot <- ggplot (meta, aes ( Days_post_symp))+ 
  geom_bar (aes (fill= age_class), position="dodge") +
  theme_cowplot (12) + xlab ("Days Post Symptom Onset") + ylab ("Number of Specimens")+
  scale_fill_manual(values=c( "#FDAFA7","#0B2A5C"))+ 
  theme (legend.position = "none",  axis.title = element_text(size=11), axis.text = element_text(size =9))  

# combine individual plots for samples per day (Figure S3)
bottom <- plot_grid(vaccine_specimen_per_day_plot, age_specimen_per_day_plot , clade_specimen_per_day_plot, labels = c ("A","B","C"),vjust=0.5, nrow=1)
plot_grid (NULL, bottom, rel_heights = c (.1, 1), nrow=2)

