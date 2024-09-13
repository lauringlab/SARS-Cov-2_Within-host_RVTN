library (scico)
library (ggplot2)
library (patchwork)


snv_meta <- read.csv ("Results/snv_meta.csv")
rep_allele_freq_plot <- ggplot (snv_meta, aes (ALT_FREQ_1 , ALT_FREQ_2 )) + geom_point (alpha =.5) + 
  xlab (" Allele Frequency in Replicate 1")+ ylab ("Allele Frequency in Replicate 2") +
  theme_cowplot (12) + theme (legend.position = "none", axis.title= element_text(size = 11), 
                             axis.text = element_text(size =9))


plot_insert <- ggplot (snv_meta, aes (ALT_FREQ_1, ALT_FREQ_2))+ geom_point(size=0.3, alpha =0.5, shape=16)+ 
  theme_cowplot(12) + theme (panel.background = element_rect(fill = "white"), axis.title= element_text(size = 11), 
                             axis.text = element_text(size =7))+
  scale_x_continuous(breaks = c (0, 0.05, 0.1), limits = c(0, 0.1)) + 
  scale_y_continuous(breaks= c (0, 0.05, 0.1),limits = c(0, 0.1)) + 
  xlab ("")+ ylab ("")

snv_rep_plot_inset <- rep_allele_freq_plot + inset_element(plot_insert, 0.52, .00, .99, .47)



meta_sample <- read.csv ("Metadata/meta_sample_1000.csv")
sample_per_person <- count (meta_sample, cdc_studyid)
ggplot ( sample_per_person, aes (as.character(n))) + geom_bar () + 
  xlab ("Number of specimens") + ylab ("Number of individuals") +
  theme_cowplot (12) #+ scale_x_discrete(n.breaks=9)


sample_per_person_plot <- ggplot ( sample_per_person, aes (as.character(n))) + geom_bar () + 
  xlab ("Number of specimens") + ylab ("Number of individuals") +
  theme_cowplot (12)
sample_cum_person <- read.csv ("Results/cum_samples_per_person.csv")



bottom <-plot_grid (sample_per_person_plot, NULL, snv_rep_plot_inset, labels = c ("A","","B"), vjust = 0.5, rel_widths = c(1,0.1,1), nrow=1)
 plot_grid (NULL, bottom, nrow=2, rel_heights = c (0.2,1))
 
 
