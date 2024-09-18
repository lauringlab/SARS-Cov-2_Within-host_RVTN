library (dplyr)
library (ggplot2)
library (cowplot)
library (ggrepel)
library (tidyr)
library (gggenes)

# are there isnv that appear as lineage defining mutations in the next wave
lineage <- read.csv ("References/lineage_defining_mutations.csv")
lineage_pos <- distinct (lineage, AA_Position, .keep_all = T)
AA_lineage_pos<-snv_meta %>% filter (GFF_FEATURE=="S", mutation_type=="Non")
AA_lineage_pos <- inner_join(lineage_pos, AA_lineage_pos, by = c("AA_Position", "clade")) 
write.csv (AA_lineage_pos, "Results/VOC_AA_changes_Spike.csv", row.names=F, quote=F)

AA_lineage_pos <- AA_lineage_pos %>% distinct (cdc_studyid, AA_Position, .keep_all=T)
write.csv (AA_lineage_pos, "Results/VOC_AA_changes_Spike_distinct.csv", row.names=F, quote=F)

AA_lineage_snv <-snv_meta %>% filter (GFF_FEATURE=="S", mutation_type=="Non")
AA_lineage_snv <-AA_lineage_snv %>% inner_join(lineage, by = c("AA_mutation", "clade", "AA_Position")) 


ggplot (snv_meta, aes (Days_post_symp, avg_freq, group = mutation))+ geom_line () + facet_wrap(vars(cdc_studyid))


## selection coefficient results from WFABC
# read in selection coefficients
selection_Coef <- read.csv ("Results/WFABC/selection_coefficients.csv")

# read in other mutational data 

snv_ps <- snv_meta %>% select (POS, REF, ALT, GFF_FEATURE, AA_mutation, mutation_type, cdc_studyid, age_at_enrollment, sex, clade, vaccine_covid_doses, mutation, Days_post_symp, avg_freq) 

snv_ps <-  unite (snv_ps, hhsubid_pos, c (cdc_studyid, POS), sep="-", remove=F) 

snv_ps_sub <- snv_ps %>% distinct (hhsubid_pos, .keep_all=T)


#filter for positive selection coefficients- 95 posterior density cannot cross 0 for any of the 3 sample sizes
selection_pos <-filter (selection_Coef, mean_pos_selection == "TRUE",
                   minus_sd_pos_selection == "TRUE", plus_sd_pos_selection == "TRUE")
selection_pos <-left_join (selection_pos, snv_ps, by = "hhsubid_pos") 
write.csv (selection_pos_sub, "Results/wfabc_pos_selection.csv", row.names=F, quote=F)
write.csv (selection_pos, "Results/wfabc_pos_selection_all_time_points.csv", row.names=F, quote=F)


selection_pos_sub <- selection_pos %>% distinct (hhsubid_pos, .keep_all=T)

selection_pos <-left_join (selection_pos, snv_ps, by = c("hhsubid_pos", "REF","POS","ALT", "cdc_studyid") )

selection_pos_sub <- selection_pos %>% distinct (hhsubid_pos, .keep_all=T) %>% filter (!is.na(GFF_FEATURE))
selection_pos_sub<-selection_pos_sub %>% mutate (non_spike_AA = ifelse (!(GFF_FEATURE=="S"), mutation, NA))


## plot selection coefficients
gene_map <- read.csv ("References/WH1_ORFs.csv")

spike_map <- read.csv ("References/spike_motif_coordinates.csv")
spike_map <- mutate (spike_map, feature = "CDS")

non_spike <- filter (selection_pos_sub, !is.na (non_spike_AA))
lollipop <- ggplot (selection_pos_sub, aes (POS, mean_s))+ geom_point(aes (color=mutation_type))+ 
  geom_segment( aes(x=POS, xend=POS, y=0, yend=mean_s, color = mutation_type)) +
  scale_color_manual(values=c("#481468", "#21A386")) + 
  xlab ("") + ylab ("Selection Coefficient")+ xlab ("Genome Position")+
  xlim (0,30000)+
  theme_classic ()+
  theme (legend.position = "none", axis.title = element_text(size=11), axis.text = element_text(size =9)) +
  #, text = element_text(size = 10) , axis.text.x = element_text(size =9), axis.text.y = element_text(size =9)) +
  geom_label_repel(data=non_spike, aes (label = AA_mutation), min.segment.length = 0, size=2)

lollipop_spike <- ggplot (data= subset (selection_pos_sub, GFF_FEATURE=="S"), aes (POS, mean_s))+ geom_point(aes (color=mutation_type))+ 
  geom_segment( aes(x=POS, xend=POS, y=0, yend=mean_s, color = mutation_type)) +
  scale_color_manual(values=c("#481468", "#21A386")) + xlim (21000, 25500)+
  ylab ("Selection Coefficient")+ xlab ("Genome Position")+
  theme_classic ()+theme (legend.position = "none", , axis.title = element_text(size=11), axis.text = element_text(size =9)) + 
  #text = element_text(size = 10), 
  # axis.text.x = element_text(size =9), axis.text.y = element_text(size =9)) +
  geom_label_repel( aes (label = AA_mutation),min.segment.length = 0, size=2)

genome_plot_inset <- ggplot(gene_map, aes(xmin = start, xmax = end, y = feature))+
  geom_gene_arrow(fill="gray46")+
  geom_gene_label(aes(label= GFF_FEATURE), color="grey95")+
  xlim (0,30000)+
  theme_genes () + theme (legend.position = "none", axis.text.x = element_blank(), 
                          axis.ticks.x = element_blank(), axis.title.y= element_blank(), axis.text.y=element_blank(), plot.background = element_rect(fill='transparent', color=NA))


spike_plot_inset <- ggplot(spike_map, aes(xmin = start, xmax = end, y = feature))+
  geom_gene_arrow(fill="gray46")+
  geom_gene_label(aes(label= Motif), color="grey95")+
  xlim (21000, 25500)+
  theme_genes () + theme (legend.position = "none", axis.text.x = element_blank(), 
                          axis.ticks.x = element_blank(), axis.title.y= element_blank(), axis.text.y=element_blank(), plot.background = element_rect(fill='transparent', color=NA))

genome_pos_sel <- plot_grid(genome_plot_inset, lollipop, ncol = 1, rel_heights = c(.25, 1)) ## add gene/region coordinates
spike_pos_sel <- plot_grid(spike_plot_inset, lollipop_spike, ncol = 1, rel_heights = c(.25, 1))

top <- plot_grid(genome_pos_sel, spike_pos_sel, labels = c('A', 'B'), rel_widths = c(1,.8))


## allele trajectory of selected isnv 
selection_pos_cds <- filter (selection_pos, !is.na (GFF_FEATURE))


trajectory_plot <- ggplot (selection_pos_cds, aes (Days_post_symp, avg_freq))+ 
  geom_line (aes (color = mutation_type, group=hhsubid_pos))+
  xlab ("Days post symptom onset (DPS)") + ylab ("iSNV frequency")+
scale_color_manual(values=c("#481468", "#21A386"))+
  theme_bw()+ theme (legend.position = "none", axis.title = element_text(size=11), axis.text = element_text(size =9))+
  facet_wrap (vars (cdc_studyid), nrow = 2) 
  
linkage <- filter (selection_pos_cds, cdc_studyid == "103403")

linkage_plot <-ggplot (linkage, aes (Days_post_symp, avg_freq))+ 
  geom_line (aes (color = mutation_type))+
  xlab ("DPS") + ylab ("iSNV frequency")+
  scale_color_manual(values=c("#481468", "#21A386"))+
  theme_bw()+ theme (legend.position = "none", axis.title = element_text(size=11), axis.text = element_text(size =9))+
  facet_wrap (vars (mutation), nrow = 3) 

bottom <- plot_grid(trajectory_plot, linkage_plot, rel_widths =c (1, 0.3), labels = c ("C", "D"))

## combine selection coefficient plot and allele trajectory plot
plot_grid (NULL, top, NULL, bottom, rel_heights = c (.1, 1, .1,1), nrow = 4)



## plot population level selection coefficient from "Fitness effects of mutations to SARS-CoV-2 proteins"

pop_sel_coef <- read.csv("Results/selection_coefficient_bloom.csv")
pop_sel_coef <- pop_sel_coef %>% mutate (spike = ifelse (Gene=="S", "yes", "no"))

ggplot (pop_sel_coef, aes (selection_coefficient, E50))+ geom_point (aes (color=mutation_type, shape = spike), size =2)+
  xlab ("Within-host selection coefficent") + ylab ("population level selection coefficient") + 
  scale_color_manual(values=c("#481468", "#21A386")) +theme_cowplot (12)
