library (dplyr)
library(tidyr)
library(broom)
library (ggplot2)
library(grid)
library(gridExtra)
library (cowplot)
library (FSA)

## get average divergence rate for each individual using the linear regression method. 
## includes stats and plots


##WholeGenome
meta <- read.csv ("Metadata/meta_1000.csv")
DataWithinDiv <- read.csv ("Results/SampleDivPerSiteGenomeAllSamples.csv")# read in divergence rate per specimen

# add in studyid for each specimen
DataWithinRates <- meta %>% select (cdc_studyid, cdc_specid) %>% right_join(DataWithinDiv, by = "cdc_specid")


DataWithinRates <- DataWithinRates %>% 
  filter (!is.na(DPI)) %>% # filter out asymptomatic cases
  group_by(cdc_studyid,  mutation_type) %>% 
  do(tidy(lm(DivPerSite ~ DPI, .))) %>% #perform linear regression for each indivdual
  filter(term=="DPI") %>%
  dplyr::rename(MeanDivPerSitePerDay=estimate,
                SEDivPerSitePerDay=std.error)

# add in meta data on individual (clade, age, vaccination)
meta_ind <- distinct (meta, cdc_studyid, .keep_all=T) %>%
  select (cdc_studyid,  age_at_enrollment, clade, vaccinated) %>% 
  mutate (age_class= ifelse (age_at_enrollment<18, "Child", "Adult"))

DataWithinRates <- left_join (DataWithinRates, meta_ind, by="cdc_studyid")

# add in number of samples for each individual
sample_count <- count (meta, cdc_studyid) %>% rename ("Num_samples" = "n") 
DataWithinRates <- left_join (DataWithinRates, sample_count, by="cdc_studyid")

## Complete above analysis for each gene, divergence rates for each sample  need to be completed first

DataWithinDivGene<- snv_meta  %>% 
  filter (!is.na (GFF_FEATURE)) %>%
  filter (cdc_specid %in% meta_qpcr$cdc_specid)%>%
  group_by(cdc_specid, mutation_type, GFF_FEATURE) %>%
  summarise(Div=sum(avg_freq))


# Identify samples that are not represented in this dataframe
# because they have zero variants.
# Add them to the dataframe in the form of a dummy value.

sample_no_ISNV_Gene <- meta %>% 
  filter(!(cdc_specid %in% DataWithinDivGene$cdc_specid))%>%
  dplyr::select(cdc_specid) %>%
  mutate( mutation_type="Non",GFF_FEATURE = "S", Div=0)

DataWithinDivGene <- rbind (DataWithinDivGene, sample_no_ISNV_Gene )

# Fill in zero values for samples, genes, and mutation classes
# that have zero reported variants using the complete function.
DataWithinDivGene <- DataWithinDivGene %>% ungroup() %>%
  complete(cdc_specid, mutation_type, GFF_FEATURE,  fill=list(Div=0))

# Add metadata about the length of each coding sequence.
DataWithinDivGene<- left_join (DataWithinDivGene, CodingSequenceLengths, by = "GFF_FEATURE" )



# Add metadata about the proportion of sites of each mutation type.
DataWithinDivGene <- left_join(DataWithinDivGene, AvailableSites, 
                               by=c("mutation_type"))


# Normalize sample divergence based on the number of available sites.
DataWithinDivGene$Length <- as.numeric (DataWithinDivGene$Length )

DataWithinDivGene<- DataWithinDivGene%>%
  mutate(DivPerSite=Div/(Length*PercentSites))

# Normalize sample divergence based on days post infection
meta_sub <- meta %>% ungroup () %>% select ( cdc_specid, Days_post_symp, age_at_enrollment, vaccinated, clade, clade_who, covqpcr_load)
DataWithinDivGene <- left_join(DataWithinDivGene, meta_sub, by = "cdc_specid" ) %>% 
  mutate (DPI = Days_post_symp +2 )

DataWithinDivGene <- DataWithinDivGene %>% mutate (DivPerDay =  DivPerSite/DPI)



DataWithinRatesGene <- meta %>% select (cdc_studyid, cdc_specid) %>% right_join(DataWithinDivGene, by = "cdc_specid")

DataWithinRatesGene <- DataWithinRatesGene %>% 
  filter (!is.na(DPI)) %>% 
  group_by(cdc_studyid,  mutation_type, GFF_FEATURE) %>%
  do(tidy(lm(DivPerSite ~ DPI, .))) %>%
  filter(term=="DPI") %>%
  dplyr::rename(MeanDivPerSitePerDay=estimate,
                SEDivPerSitePerDay=std.error)

## plots 

ggplot (DataWithinRates, aes (x= as.factor (Num_samples),y= MeanDivPerSitePerDay, group=mutation_type)) + 
  geom_point(position=position_jitterdodge(dodge.width=.95), aes (color=mutation_type ), show.legend=F, alpha = 0.6) +
  scale_color_manual(values=c("#77427F", "#21A386"))+
  theme_cowplot ()+ 
  theme(axis.title= element_text(size = 11), axis.text = element_text(size =9)) +
  xlab ("Days post infection") + ylab ("Divergence rate")



vax_rate_plot <- ggplot (DataWithinRates, aes (as.factor (vaccinated), MeanDivPerSitePerDay))+
   geom_boxplot2 (aes(fill= mutation_type), width.errorbar = 0.35, show.legend = F)+
  theme_cowplot (12) + #scale_color_manual(values=c("#77427F", "#21A386")) +
  scale_fill_manual(values=c("#77427F", "#21A386")) +
  xlab ("Vaccinated")+ ylab ("Div/Site/Day")+
  theme(axis.title= element_text(size = 11), axis.text = element_text(size =9)) 
  

age_rate_plot <- ggplot (DataWithinRates, aes (age_class, MeanDivPerSitePerDay))+
  geom_boxplot2 (aes(fill= mutation_type), width.errorbar = 0.35, show.legend = F)+
  theme_cowplot (12) + #scale_color_manual(values=c("#77427F", "#21A386")) +
  scale_fill_manual(values=c("#77427F", "#21A386")) +
  xlab ("Age")+ ylab ("Div/Site/Day")+
  theme(axis.title= element_text(size = 11), axis.text = element_text(size =9)) 

clade_rate_plot <- ggplot (DataWithinRates, aes (clade, MeanDivPerSitePerDay))+
  geom_boxplot2 (aes(fill= mutation_type), width.errorbar = 0.35, show.legend = F)+
  theme_cowplot (12) + #scale_color_manual(values=c("#77427F", "#21A386")) +
  scale_fill_manual(values=c("#77427F", "#21A386")) +
  xlab ("Clade")+ ylab ("Div/Site/Day")+
  theme(axis.title= element_text(size = 11), axis.text = element_text(size =9)) 


sample_rate_plot <- ggplot (DataWithinRates, aes (as.factor (Num_samples), MeanDivPerSitePerDay))+
  geom_boxplot2 (aes(fill= mutation_type), width.errorbar = 0.35, show.legend = F)+
  theme_cowplot (12) + #scale_color_manual(values=c("#77427F", "#21A386")) +
  scale_fill_manual(values=c("#77427F", "#21A386")) +
  xlab ("Number of specimens")+ ylab ("Div/Site/Day")+
  theme(axis.title= element_text(size = 11), axis.text = element_text(size =9)) 


DataWithinRatesGene$GFF_FEATURE <- ordered(DataWithinRatesGene$GFF_FEATURE, levels = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "N"))

gene_rate_plot <- ggplot (DataWithinRatesGene, aes (GFF_FEATURE, MeanDivPerSitePerDay, group=mutation_type))+
  #geom_boxplot2 (aes(fill= mutation_type), width.errorbar = 0.35, show.legend = F)+
  geom_point(position=position_jitterdodge(dodge.width=.95), aes (color=mutation_type ), show.legend=F, alpha = 0.6) +
  theme_cowplot (12) + scale_color_manual(values=c("#77427F", "#21A386")) +
  xlab ("Gene") + ylab ("Div/Day/Site" ) +
  theme(axis.title= element_text(size = 11), axis.text = element_text(size =9)) 

top_plot <- plot_grid (vax_rate_plot,age_rate_plot,clade_rate_plot, nrow=1, rel_widths = c(1,1,1.2), labels= c("A", "B", "C"), vjust=0.5 )

plot_grid (NULL,top_plot, NULL, gene_rate_plot, NULL, sample_rate_plot,NULL, nrow=7, 
           rel_heights = c(.1, 1, .1, 1, .1, 1.3, .1), labels = c ("","","","D","","E", ""), vjust =0.5)
sample_rate_plot 


## stats

DataWithinRatesNon <- filter (DataWithinRates, mutation_type == "Non")
DataWithinRatesSyn <- filter (DataWithinRates, mutation_type == "Syn")

wilcox.test(MeanDivPerSitePerDay ~ vaccinated, data=DataWithinRates) 
wilcox.test(MeanDivPerSitePerDay ~ vaccinated, data=DataWithinRatesNon ) 
wilcox.test(MeanDivPerSitePerDay ~ vaccinated, data=DataWithinRatesSyn) 


wilcox.test(MeanDivPerSitePerDay~ age_class, data=DataWithinRates) 
wilcox.test(MeanDivPerSitePerDay ~ age_class, data=DataWithinRatesNon ) 
wilcox.test(MeanDivPerSitePerDay ~ age_class, data=DataWithinRatesSyn) 

kruskal.test(MeanDivPerSitePerDay ~ clade, data = DataWithinRates)
kruskal.test(MeanDivPerSitePerDay ~ clade, data = DataWithinRatesNon)
kruskal.test(MeanDivPerSitePerDay ~ clade, data = DataWithinRatesSyn)

kruskal.test(MeanDivPerSitePerDay ~ Num_samples, data = DataWithinRates)
kruskal.test(MeanDivPerSitePerDay ~ Num_samples, data = DataWithinRatesNon)
kruskal.test(MeanDivPerSitePerDay ~ Num_samples, data = DataWithinRatesSyn)

DataWithinRatesGeneNon <- filter (DataWithinRatesGene, mutation_type == "Non")
DataWithinRatesGeneSyn <- filter (DataWithinRatesGene, mutation_type == "Syn")


kruskal.test(MeanDivPerSitePerDay ~ GFF_FEATURE, data = DataWithinRatesGene)
kruskal.test(MeanDivPerSitePerDay ~ GFF_FEATURE, data = DataWithinRatesGeneNon)
kruskal.test(MeanDivPerSitePerDay ~ GFF_FEATURE, data = DataWithinRatesGeneSyn)
