library (dplyr)
library(tidyr)
library(broom)
library (ggplot2)
library(grid)
library(gridExtra)
library (cowplot)
library (FSA)

## calculation of divergence rate using the point method, statistics and figures. For highest ct value sample and all samples, by gene and whole genome

meta <- read.csv ("Metadata/meta_1000.csv")

# Import information on the number of available sites of each mutation type.
AvailableSites <- read.csv("References/AvailableSitesByType.csv",stringsAsFactors = FALSE)

# Import information on coding sequence lengths.
CodingSequenceLengths <- read.csv("References/CodingSequenceLengths.csv", stringsAsFactors = FALSE)

## divergence rate for all samples - whole genome
Genome_length <- CodingSequenceLengths %>% summarise (Length = sum (Length))


## get average isnv frequency for wach specimen
DataWithinDiv<- snv_meta  %>% 
  filter (!is.na (GFF_FEATURE)) %>%
  group_by(cdc_specid, mutation_type) %>%
  summarise(Div=sum(avg_freq))


# Identify samples that are not represented in this dataframe
# because they have zero variants.
# Add them to the dataframe in the form of a dummy value.

sample_no_ISNV <- meta %>% ungroup () %>% 
  filter(!(cdc_specid %in% DataWithinDiv$cdc_specid))%>%
  dplyr::select(cdc_specid) %>%
  mutate( mutation_type="Non",Div=0)

DataWithinDiv <- rbind (DataWithinDiv, sample_no_ISNV )

# Fill in zero values for  mutation classes that have zero reported variants using the complete function.
DataWithinDiv <- DataWithinDiv %>% ungroup() %>%
  complete(cdc_specid, mutation_type, fill=list(Div=0))

# Add metadata about the length of each coding sequence.
DataWithinDiv<- mutate (DataWithinDiv, Length = "29006")
                        

# Add metadata about the proportion of sites of each mutation type.
DataWithinDiv <- left_join(DataWithinDiv, AvailableSites, 
                           by=c("mutation_type"))


# Normalize sample divergence based on the number of available sites.
DataWithinDiv$Length <- as.numeric (DataWithinDiv$Length )

DataWithinDiv <- DataWithinDiv%>%
  mutate(DivPerSite=Div/(Length*PercentSites))

# Normalize sample divergence based on days post infection
meta_sub <- meta %>%  select ( cdc_specid, Days_post_symp)
DataWithinDiv <- left_join(DataWithinDiv, meta_sub, by = "cdc_specid" )  %>%
  mutate (DPI = Days_post_symp +2 )

DataWithinDiv <- DataWithinDiv %>% mutate (DivPerDay =  DivPerSite/DPI)

write.csv(DataWithinDiv, "Results/SampleDivPerSiteGenomeAllSamples.csv",
          quote=FALSE, row.names=FALSE)





# use only highest viral load sample per individual
meta$covqpcr_load <- as.numeric (meta$covqpcr_load)
meta_qpcr <- meta %>% group_by (cdc_studyid) %>% 
  slice_max (covqpcr_load) %>% # chose sample with highest viral load per person
  slice_min (Days_post_symp) # chose earliest timepoint if multiple samples have identical viral loads

DataWithinDivqPCR<- snv_meta  %>% 
  filter (!is.na (GFF_FEATURE)) %>%
  filter (cdc_specid %in% meta_qpcr$cdc_specid)%>%
  group_by(cdc_specid, mutation_type) %>%
  summarise(Div=sum(avg_freq))


# Identify samples that are not represented in this dataframe
# because they have zero variants.
# Add them to the dataframe in the form of a dummy value.

sample_no_ISNV_qPCR <- meta_qpcr %>% ungroup () %>% 
  filter(!(cdc_specid %in% DataWithinDivqPCR$cdc_specid))%>%
  dplyr::select(cdc_specid) %>%
  mutate( mutation_type="Non",Div=0)

DataWithinDivqPCR <- rbind (DataWithinDivqPCR, sample_no_ISNV_qPCR )

# Fill in zero values for samples, genes, and mutation classes
# that have zero reported variants using the complete function.
DataWithinDivqPCR <- DataWithinDivqPCR %>% ungroup() %>%
  complete(cdc_specid, mutation_type, fill=list(Div=0))

# Add metadata about the length of each coding sequence.
DataWithinDivqPCR<- mutate (DataWithinDivqPCR, Length = "29006")


# Add metadata about the proportion of sites of each mutation type.
DataWithinDivqPCR <- left_join(DataWithinDivqPCR, AvailableSites, 
                           by=c("mutation_type"))


# Normalize sample divergence based on the number of available sites.
DataWithinDivqPCR$Length <- as.numeric (DataWithinDivqPCR$Length )

DataWithinDivqPCR <- DataWithinDivqPCR%>%
  mutate(DivPerSite=Div/(Length*PercentSites))

# Normalize sample divergence based on days post infection
meta_sub_qpcr <- meta_qpcr %>% ungroup () %>% select ( cdc_specid, Days_post_symp, age_at_enrollment, vaccinated, clade, clade_who, covqpcr_load)
DataWithinDivqPCR <- left_join(DataWithinDivqPCR, meta_sub_qpcr, by = "cdc_specid" ) %>% 
  mutate (DPI = Days_post_symp +2 )

DataWithinDivqPCR <- DataWithinDivqPCR %>% mutate (DivPerDay =  DivPerSite/DPI)

write.csv(DataWithinDivqPCR, "Results/SampleDivPerSiteGenomeHighestViralLoad.csv",
          quote=FALSE, row.names=FALSE)

##without high divergence samples

high_divergence_samples <- DataWithinDivqPCR %>% filter (DivPerDay > 1e-5)
DataWithinDivqPCRh<- filter ( DataWithinDivqPCR , (!cdc_specid %in% high_divergence_samples$cdc_specid))

### divergence by gene

meta_qpcr <- meta %>% group_by (cdc_studyid) %>% 
  slice_max (covqpcr_load) %>% # chose sample with highest viral load per person
  slice_min (Days_post_symp) # chose earliest timepoint if multiple samples have identical viral loads

DataWithinDivGene<- snv_meta  %>% 
  filter (!is.na (GFF_FEATURE)) %>%
  filter (cdc_specid %in% meta_qpcr$cdc_specid)%>%
  group_by(cdc_specid, mutation_type, GFF_FEATURE) %>%
  summarise(Div=sum(avg_freq))


# Identify samples that are not represented in this dataframe
# because they have zero variants.
# Add them to the dataframe in the form of a dummy value.

sample_no_ISNV_Gene <- meta_qpcr %>% ungroup () %>% 
  filter(!(cdc_specid %in% DataWithinDivGene$cdc_specid))%>%
  dplyr::select(cdc_specid) %>%
  mutate( mutation_type="Non",GFF_FEATURE = "S", Div=0)

DataWithinDivGene <- rbind (DataWithinDivGene, sample_no_ISNV_Gene )

# Fill in zero values for samples, genes, and mutation classes
# that have zero reported variants using the complete function.
DataWithinDivGene <- DataWithinDivGene %>% ungroup() %>%
  complete(cdc_specid, mutation_type, GFF_FEATURE,  fill=list(Div=0))

# Add metadata about the length of each coding sequence.
CodingSequenceLengths <- dplyr::rename (CodingSequenceLengths, GFF_FEATURE = Attr) %>% select (-Seqname)
DataWithinDivGene<- left_join (DataWithinDivGene, CodingSequenceLengths, by = "GFF_FEATURE" )



# Add metadata about the proportion of sites of each mutation type.
DataWithinDivGene <- left_join(DataWithinDivGene, AvailableSites, 
                               by=c("mutation_type"))


# Normalize sample divergence based on the number of available sites.
DataWithinDivGene$Length <- as.numeric (DataWithinDivGene$Length )

DataWithinDivGene<- DataWithinDivGene%>%
  mutate(DivPerSite=Div/(Length*PercentSites))

# Normalize sample divergence based on days post infection
meta_sub_qpcr <- meta_qpcr %>% ungroup () %>% select ( cdc_specid, Days_post_symp, age_at_enrollment, vaccinated, clade, clade_who, covqpcr_load)
DataWithinDivGene <- left_join(DataWithinDivGene, meta_sub_qpcr, by = "cdc_specid" ) %>% 
  mutate (DPI = Days_post_symp +2 )

DataWithinDivGene <- DataWithinDivGene %>% mutate (DivPerDay =  DivPerSite/DPI)

write.csv(DataWithinDivGene, "Results/SampleDivPerSiteGeneHighestViralLoad.csv",
          quote=FALSE, row.names=FALSE)



#genome wide plot using all samples
mutation.labs <- c("Nonsynonymous", "Synonymous")
names(mutation.labs) <- c("Non", "Syn")

ggplot (DataWithinDiv, aes (x= DPI,y= DivPerDay)) + 
  geom_point (alpha =0.6, shape = 16, aes (color=mutation_type), show.legend = FALSE)+
  scale_color_manual(values=c("#77427F", "#21A386"))+
  theme_bw ()+ 
  theme(axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab ("Days post infection") + ylab ("Divergence rate") +
  facet_wrap (vars(mutation_type), labeller = labeller(mutation_type = mutation.labs))


#genome wide plot using highest viral load sample plotted by various samples
## non vs syn 
full_sum <- DataWithinDivqPCR %>% filter (!is.na (DivPerDay)) %>%
  filter (!is.na (vaccinated)) %>%
  group_by(mutation_type)%>% 
  summarise (mean = mean (DivPerDay), n= n(), sd =sd (DivPerDay), SE = sd/sqrt(n))

ggplot (DataWithinDivqPCR, aes (x=mutation_type, y= DivPerDay, group=mutation_type))+ 
  geom_jitter (alpha =0.6, aes (color=mutation_type),show.legend = FALSE) +
  geom_crossbar(data=full_sum, aes(y= mean, ymin = mean, ymax = mean), linewidth=.3, width = .5)+ 
  #geom_jitter (alpha =0.6, aes (color=mutation_type),show.legend = FALSE) +
  theme_cowplot (12) + scale_color_manual(values=c("#77427F", "#21A386")) +
  xlab ("Mutation type") + ylab ("Div/Day/Site" ) 
  
ggplot (DataWithinDivqPCR, aes (x=mutation_type, y= DivPerDay, group=mutation_type))+ 
  geom_crossbar(data=full_sum, aes(y= mean, ymin = mean, ymax = mean), size=.4, width = .5)+ 
  geom_jitter (alpha =0.6, aes (color=mutation_type),show.legend = FALSE) +
  theme_cowplot (12) + scale_color_manual(values=c("#77427F", "#21A386")) +
  xlab ("Mutation type") + ylab ("Div/Day/Site" ) +ylim (0,2.5e-6)

## vaccine

vaccine_sum <- DataWithinDivqPCR %>% filter (!is.na (DivPerDay)) %>%
  filter (!is.na (vaccinated)) %>%
  group_by(vaccinated, mutation_type)%>% 
  summarise (mean = mean (DivPerDay), n= n(), sd =sd (DivPerDay), SE = sd/sqrt(n))

ggplot (DataWithinDivqPCR, aes (x=as.factor (vaccinated), y= DivPerDay, group=mutation_type))+ 
  geom_point(position=position_jitterdodge(dodge.width=.95), aes (color=mutation_type ), show.legend=F, alpha = 0.6) +
  geom_crossbar(data=vaccine_sum, aes(y= mean, ymin = mean, ymax = mean), size=.4, width = .5, position = position_dodge(width =.95)) + 
# geom_point(position=position_jitterdodge(dodge.width=.95), aes (color=mutation_type ), show.legend=F, alpha = 0.6) +
  theme_cowplot (12) + scale_color_manual(values=c("#77427F", "#21A386")) +
  xlab ("Vaccination") + ylab ("Div/Day/Site" ) 



vvaccine_mean_plot <- ggplot (vaccine_sum, aes (vaccinated, mean))+ 
  geom_pointrange (aes (ymin = mean-SE, ymax = mean+SE, color= mutation_type), 
                   position = position_dodge(width =0.4), show.legend = FALSE) +
  scale_color_manual(values=c("#77427F", "#21A386"))+
  theme_cowplot  (12 )+
  xlab ("Vaccination") + ylab ("")

# subtype
subtype_sum <- DataWithinDivqPCR %>% filter (!is.na (DivPerDay)) %>%
  group_by(clade, mutation_type)%>% 
  summarise (mean = mean (DivPerDay), n= n(), sd =sd (DivPerDay), SE = sd/sqrt(n))

subtype_mean_plot <- ggplot (subtype_sum, aes (clade, mean) )+ 
  geom_pointrange (aes (ymin = mean-SE, ymax = mean+SE, color= mutation_type), 
                   position = position_dodge(width =0.4), show.legend = FALSE) +
  scale_color_manual(values=c("#77427F", "#21A386"))+
  scale_x_discrete(labels=c("0"="No", "1"= "Yes"))+
  theme_cowplot  (12 )+
  xlab ("Subtype") + ylab ("")

ggplot (DataWithinDivqPCR, aes (x=as.factor (clade), y= DivPerDay, group=mutation_type))+ 
  geom_point(position=position_jitterdodge(dodge.width=.95), aes (color=mutation_type ), show.legend=F, alpha = 0.6) +
  geom_crossbar(data=subtype_sum, aes(y= mean, ymin = mean, ymax = mean), size=.4, width = .5, position = position_dodge(width =.95)) + 
  #geom_point(position=position_jitterdodge(dodge.width=.95), aes (color=mutation_type ), show.legend=F, alpha = 0.6) +
  theme_cowplot (12) + scale_color_manual(values=c("#77427F", "#21A386")) +
  xlab ("Clade") + ylab ("Div/Day/Site" ) 



# age


ggplot (DataWithinDivqPCR, aes (x= age_at_enrollment,y= DivPerDay)) + 
  geom_point (alpha =0.6, shape = 16, aes (color=mutation_type), show.legend = FALSE)+
  scale_color_manual(values=c("#77427F", "#21A386"))+
  theme_bw ()+ 
  theme(axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab ("Age") + ylab ("Divergence rate") +
  facet_wrap (vars(mutation_type), labeller = labeller(mutation_type = mutation.labs))

### gene

segment_rate_plot <-  ggplot (DataWithinDivGene, aes (GFF_FEATURE, DivPerDay, fill = mutation_type)) + 
  geom_boxplot ( aes (GFF_FEATURE), show.legend = FALSE) +
  scale_fill_manual(values=c("#77427F", "#21A386"))+
  scale_x_discrete(labels=c("NA_"="NA"))+
  theme_cowplot  (12 )+
  xlab ("Age") + ylab ("")


DataWithinDivGene <- filter ( DataWithinDivGene, !is.na(DivPerDay))
segment_sum <-  DataWithinDivGene %>% filter (!is.na (DivPerDay)) %>% 
  group_by(GFF_FEATURE, mutation_type)%>% 
  summarise (mean = mean (DivPerDay), n= n(), sd =sd (DivPerDay), SE = sd/sqrt(n))


segment_sum$GFF_FEATURE <- ordered(segment_sum$GFF_FEATURE, levels = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "N"))

segment_rate_mean_plot <-  ggplot (  segment_sum, aes (GFF_FEATURE, mean, fill = mutation_type)) + 
  geom_pointrange (aes (ymin = mean-SE, ymax = mean+SE, color= mutation_type), position = position_dodge(width =0.4), show.legend = FALSE) +
  scale_color_manual(values=c("#77427F", "#21A386"))+
  scale_x_discrete(labels=c("NA_"="NA"))+
  theme_cowplot  (12 )+
  xlab ("Segment") + ylab ("")


ggplot (DataWithinDivGene, aes (x=GFF_FEATURE, y= DivPerDay, group=mutation_type))+ 
  geom_point(position=position_jitterdodge(dodge.width=.95), aes (color=mutation_type ), show.legend=F, alpha = 0.6) +
  geom_crossbar(data=segment_sum, aes(y= mean, ymin = mean, ymax = mean), size=.4, width = .5, position = position_dodge(width =.95)) + 
  #geom_point(position=position_jitterdodge(dodge.width=.95), aes (color=mutation_type ), show.legend=F, alpha = 0.6) +
  theme_cowplot (12) + scale_color_manual(values=c("#77427F", "#21A386")) +
  xlab ("Gene") + ylab ("Div/Day/Site" ) + ylim (0,1e-5)


### stats 
DataWithinDivNon <- filter (DataWithinDiv, mutation_type == "Non")
DataWithinDivSyn <- filter (DataWithinDiv, mutation_type == "Syn")

kruskal.test(DivPerDay ~ DPI, data = DataWithinDiv)
kruskal.test(DivPerDay ~ DPI, data = DataWithinDivNon)
kruskal.test(DivPerDay ~ DPI, data = DataWithinDivSyn)



DataWithinDivqPCRNon <- filter (DataWithinDivqPCR, mutation_type == "Non")
DataWithinDivqPCRSyn <- filter (DataWithinDivqPCR, mutation_type == "Syn")

wilcox.test(DivPerDay ~ vaccinated, data=DataWithinDivqPCR) 
wilcox.test(DivPerDay ~ vaccinated, data=DataWithinDivqPCRNon ) 
wilcox.test(DivPerDay ~ vaccinated, data=DataWithinDivqPCRSyn) 


DataWithinDivqPCR<- DataWithinDivqPCR %>% mutate (age_class =  ifelse (age_at_enrollment >=18 , "Adult", "Child"))
DataWithinDivqPCRNon<- DataWithinDivqPCRNon %>% mutate (age_class =  ifelse (age_at_enrollment >=18 , "Adult", "Child"))
DataWithinDivqPCRSyn<- DataWithinDivqPCRSyn %>% mutate (age_class = ifelse (age_at_enrollment >=18 , "Adult", "Child"))


wilcox.test(DivPerDay ~ age_class, data=DataWithinDivqPCR) 
wilcox.test(DivPerDay ~ age_class, data=DataWithinDivqPCRNon ) 
wilcox.test(DivPerDay ~ age_class, data=DataWithinDivqPCRSyn) 

kruskal.test(DivPerDay ~ clade, data = DataWithinDivqPCR)
kruskal.test(DivPerDay ~ clade, data = DataWithinDivqPCRNon)
kruskal.test(DivPerDay ~ clade, data = DataWithinDivqPCRSyn)

dunnTest(DivPerDay ~ clade, data = DataWithinDivqPCR)
dunnTest(DivPerDay ~ clade, data = DataWithinDivqPCRNon)
dunnTest(DivPerDay ~ clade, data = DataWithinDivqPCRSyn)





DataWithinDivGeneNon <- filter (DataWithinDivGene, mutation_type == "Non")
DataWithinDivGeneSyn<- filter (DataWithinDivGene, mutation_type == "Syn")


kruskal.test(DivPerDay ~ GFF_FEATURE, data = DataWithinDivGene)
kruskal.test(DivPerDay ~ GFF_FEATURE, data = DataWithinDivGeneNon)
kruskal.test(DivPerDay ~ GFF_FEATURE, data = DataWithinDivGeneSyn)

dunnTest(DivPerDay ~ GFF_FEATURE, data = DataWithinDivGene)
dunnTest(DivPerDay ~ GFF_FEATURE, data = DataWithinDivGeneNon)
dunnTest(DivPerDay ~ GFF_FEATURE, data = DataWithinDivGeneSyn)
