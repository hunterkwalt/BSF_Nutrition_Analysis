### Written by Hunter K. Walt
### Final Version: 14 April 2025

#set working directory
wd <- "/mnt/md0/bsf_protein_study/github_repo/BSF_Nutrition_Analysis"
setwd(wd) 

#libraries
library(phyloseq)
library(vegan)
library(RColorBrewer)
library(ggplot2)
library(stats)
library(rstatix)
library(ggpubr)
library(microbiome)
library(PerformanceAnalytics)
library(tidyverse)
library(RVAideMemoire)
library(data.table)
library(clusterProfiler)

### necessary files to reproduce analysis ###
kb <- "community_files/kraken_biom_genus.json" #kraken biom file - this is generated from the output of kracken2/bracken
genes <- "community_files/normalized_counts.csv" # Full set of normalized counts from BSF
deg <- "community_files/DE_genes.txt" #List of differentially expressed BSF Genes
mcts <- "community_files/microbial_normalized_cts.csv" # Full set of normalized counts from microbes
mdegs <- "community_files/microbe_DE.csv" #List of differentially expressed microbial Genes
bsf_degann <- "community_files/BSF_gene_annotations.tsv" #annotations of BSF DEGs from NCBI
mic_degann <- "community_files/microbial_annotations.tsv" #annotations of Microbial sig correlated DEGs from eggnog

#import biom file from kraken-biom
data <-import_biom(kb, parseFunction=parse_taxonomy_default)
data@tax_table@.Data <- substring(data@tax_table@.Data, 4)
#rename columns to species level
colnames(data@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#get rid of some contaminant/misassigned reads (Only for genus and species)
data <- subset_taxa(data, (!Genus %in% "Homo"))
data <- subset_taxa(data, (!Genus %in% "Callorhinchus"))

#make sample data for phyloseq object
col <- names(data.frame(data@otu_table))
features <- factor(c(rep("Frass", 9), rep("Larvae", 9), "Control", 
                     rep("Frass", 9), rep("Larvae", 9), "Control",
                     rep("Frass", 9), rep("Larvae", 9), "Control",
                     rep("Frass", 9), rep("Larvae", 9), "Control"), levels = c("Larvae", "Frass", "Control"))
features2 <- factor(c(rep("Gainesville", 19), rep("5P:1C", 19),rep("1P:5C", 19), rep("1P:1C", 19)), 
                    levels = c("Gainesville", "1P:1C", "1P:5C", "5P:1C"))

#change this based on what taxonomic level you imported from Bracken. e.g., if you imported genus level data,
##change "_._report_bracken_families" to "_._report_bracken_genuses" 
feature3 <- gsub("_._report_bracken_genuses","", colnames(data@otu_table@.Data))
features3 <- gsub("_report_bracken_genuses","", feature3) #change this one to the correct taxonomy, too.

#create sample data for phyloseq object
samdf <- data.frame(col, features, features2, features3)
colnames(samdf)<-c("Data","Sample","Diet", "Enclosure")
row.names(samdf) <- samdf[,1]
samdf <- samdf[,2:4]
data@sam_data<-sample_data(samdf)

data <- subset_samples(data, Diet != "Gainesville")

data@sam_data$Diet <- factor(data@sam_data$Diet, levels = c("1P:5C", "1P:1C", "5P:1C"))

#alpha diveristy
data.alpha <- estimate_richness(data)
shannon <- data.alpha[,"Shannon"]
shannon_div <- data.frame(row.names = rownames(data.alpha), shannon, features[which(features2 != "Gainesville")], 
                          factor(features2[which(features2 != "Gainesville")], 
                                 levels = c( "1P:5C","1P:1C", "5P:1C")), 
                          features3[which(features2 != "Gainesville")])
colnames(shannon_div) <- c("Shannon", "Sample","Diet", "Enclosure")

#drop control values: no statistical power
shannon_no_ctl <- shannon_div[!(shannon_div$Sample=="Control"),]

#plot boxes
b <- ggplot(shannon_no_ctl, aes(x=Sample,y=Shannon, fill=Diet)) + geom_boxplot() 
alpha_div <- b + ggtitle("Alpha Diversity") + ylab("Shannon Diversity Index")

#look at distribution
hist(shannon_no_ctl$Shannon)
shannon_frass <- subset(shannon_no_ctl, Sample == "Frass")
shannon_larvae <- subset(shannon_no_ctl, Sample == "Larvae")
#write.csv(shannon_larvae, quote = F, file="/mnt/md0/bsf_protein_study/larvae_alpha_div.csv")


l_pairwise <- shannon_larvae %>% t_test(Shannon ~ Diet)
f_pairwise <- shannon_frass %>% t_test(Shannon ~ Diet)

l_pairwise <- l_pairwise %>% add_xy_position(x = "Diet")
f_pairwise <- f_pairwise %>% add_xy_position(x = "Diet")



#plot larvae and frass separately
b_l <- ggboxplot(shannon_larvae, x="Diet",y="Shannon", fill="Diet", palette = c("#0D7760", "goldenrod", "cornflowerblue"), font.label = list(size = 13)) 
alpha_div_l <- b_l + ggtitle("Alpha Diversity: Larvae") + ylab("Shannon Diversity") + 
  theme(plot.title =  element_text(size = 13), axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),  legend.position = "none") +
  stat_pvalue_manual(l_pairwise, label = "p.adj.signif", tip.length = 0.01)


b_f <- ggboxplot(shannon_frass, x="Diet",y="Shannon", fill="Diet", palette= c("#0D7760", "goldenrod", "cornflowerblue"), font.label = list(size = 13)) 
alpha_div_f <- b_f + ggtitle("Alpha Diversity: Frass") + ylab("Shannon Diversity") + 
  theme(plot.title = element_text(size = 13), axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), legend.position = "right", legend.text=element_text(size=12)) + stat_pvalue_manual(f_pairwise, label = "p.adj.signif", tip.length = 0.01)

#check normality
shapiro.test(shannon_no_ctl$Shannon) 
shapiro.test(shannon_frass$Shannon)
shapiro.test(shannon_larvae$Shannon)

#two way anova for frass and larvae separately
shan.aov_l <- aov(Shannon ~ Diet, data = shannon_larvae)
summary(shan.aov_l)

shan.aov_f <- aov(Shannon ~ Diet, data = shannon_frass)
summary(shan.aov_f)

#tukey post hoc to test for pairwise differences
TukeyHSD(shan.aov_l)
TukeyHSD(shan.aov_f)


#subset larvae and frass samples. 
larvae <- subset_samples(data, Sample == "Larvae")
frass <- subset_samples(data, Sample == "Frass")
control <- subset_samples(data, Sample == "Control")
#lpf <- subset_samples(frass, Diet == "Low Protein")

#make dataset for relative abundance
ps.prop_larvae <- transform_sample_counts(larvae, function(otu) otu/sum(otu))
ps.prop_frass <- transform_sample_counts(frass, function(otu) otu/sum(otu))
ps.prop_ctl <- transform_sample_counts(control, function(otu) otu/sum(otu))


### make PCoA of larvae and frass individually.
ord.nmds.bray <- ordinate(physeq = ps.prop_larvae, method="PCoA", distance="bray")
plot_ordination(ps.prop_larvae, ord.nmds.bray, color="Diet", title="Bray PCoA") #+ geom_point(aes(size = 1)) + scale_size_continuous(guide = "none")
pcoa_larvae <- plot_ordination(ps.prop_larvae, ord.nmds.bray, title="Beta Diversity: Larvae", color = "Diet",) + geom_point(size = 2) + scale_size_continuous(guide = "none")+
  xlab("PCo1 (52%)")+ ylab("PCo2 (24.1%)") + theme(text = element_text(size = 12), legend.position = "none")   + scale_color_manual(values = c( "goldenrod", "#0D7760","cornflowerblue")) + stat_ellipse(aes(color=Diet)) 

ord.nmds.bray <- ordinate(physeq = ps.prop_frass, method="PCoA", distance="bray")
plot_ordination(ps.prop_frass, ord.nmds.bray, color="Diet", title="Bray PCoA") #+ geom_point(aes(size = 1)) + scale_size_continuous(guide = "none")
pcoa_frass <- plot_ordination(ps.prop_frass, ord.nmds.bray, title="Beta Diversity: Frass", color = "Diet",) + geom_point(size = 2) + scale_size_continuous(guide = "none")+
  xlab("PCo1 (85.2%)")+ ylab("PCo2 (6.4%)") + theme(text = element_text(size = 12), legend.text=element_text(size=12))  + 
  scale_color_manual(values = c( "goldenrod", "#0D7760","cornflowerblue")) + stat_ellipse(aes(color=Diet)) 
#pcoa_frass



#how many taxa make up more than 1% mean?
ldf <- psmelt(ps.prop_larvae)
gtone <- ldf %>% 
  group_by(Genus) %>% 
  summarize(Mean = mean(Abundance *100)) %>%
  arrange(-Mean) 


fdf <- psmelt(ps.prop_frass)
fgtone <- fdf %>% 
  group_by(Genus) %>% 
  summarize(Mean = mean(Abundance *100)) %>%
  arrange(-Mean) 

cdf <- psmelt(ps.prop_ctl)
cgtone <- cdf %>% 
  group_by(Genus) %>% 
  summarize(Mean = mean(Abundance *100)) %>%
  arrange(-Mean) 


lgt_sub <- subset(gtone, Mean >= 1)
fgt_sub <- subset(fgtone, Mean >= 1)
cgt_sub <- subset(cgtone, Mean >= 1)

# Make Taxonomic Abundance Charts
col_tax <- c("Acetobacter" = "#88CCEE", "Arsenophonus" = "#CC6677", "Clavulicium" ="#0072B2", "Endomicrobium" = "#DDCC77","Enterobacter" = "#117733", "Enterococcus" = "firebrick2", "Haemophilus" = "#332288", 
             "Klebsiella" = "#AA4499", "Lactobacillus" = "#44AA99", "Magnusiomyces" = "goldenrod" , "Myroides" = "#6699CC", "Pichia" = "#661100", "Serratia" = "#999933" ,
             "Rhodotorula" = "#D55E00", "Staphylococcus" = "#888888", "Streptococcus" = "plum1", "Other (<1%)" = "grey", "Neocatenulostroma" = "palegreen2" , "Yueomyces" = "grey25")

larv_fact <- c("Acetobacter", "Arsenophonus" , "Clavulicium" ,"Enterobacter", 
               "Klebsiella", "Lactobacillus", "Magnusiomyces", "Pichia",
               "Rhodotorula", "Staphylococcus", "Streptococcus", "Other (<1%)")

for_gg <- psmelt(ps.prop_larvae)
for_gg$Genus[which(!for_gg$Genus %in% lgt_sub$Genus)] <- "Other (<1%)"
#change factor levels
for_gg$Genus <- factor(for_gg$Genus, levels = larv_fact)
#species
#for_gg$Sample <- gsub("_report_bracken_species","", for_gg$Sample)

#genus
for_gg$Sample <- gsub("_report_bracken_genuses","", for_gg$Sample)

#family
for_gg$Sample <- gsub("_report_bracken_families","", for_gg$Sample)

#make pct level bar chart
ncols <- 12
mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(ncols) 
pct_larvae <- ggplot(data = for_gg, 
                     mapping = aes_string(x = "Sample",y = "Abundance"))  + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(face = "italic", size = 11),
        text = element_text(size = 12)) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="fill") +
  facet_wrap(~Diet, scales="free_x") + scale_fill_manual(values = col_tax) + scale_color_manual(values = col_tax)+ ggtitle("Relative Abundance: Larvae")

pct_larvae


###### Make relative abundance plot for frass ########
for_gg <- psmelt(ps.prop_frass)
for_gg$Genus[which(!for_gg$Genus %in% fgt_sub$Genus)] <- "Other (<1%)"

frass_fact <- c(fgt_sub$Genus[fgt_sub$Genus %>% order()], "Other (<1%)")
#change factor levels
for_gg$Genus <- factor(for_gg$Genus, levels = frass_fact)

#species
#for_gg$Sample <- gsub("_report_bracken_species","", for_gg$Sample)

#genus
for_gg$Sample <- gsub("_report_bracken_genuses","", for_gg$Sample)

#family
for_gg$Sample <- gsub("_report_bracken_families","", for_gg$Sample)
#species
#for_gg$Sample <- gsub("_report_bracken_species","", for_gg$Sample)

#genus
for_gg$Sample <- gsub("_report_bracken_genuses","", for_gg$Sample)

#family
for_gg$Sample <- gsub("_report_bracken_families","", for_gg$Sample)

#make pct level bar chart
ncols <- 30
mycolors <- colorRampPalette(brewer.pal(30,"Paired"))(ncols) 
pct_frass <- ggplot(data = for_gg, 
                    mapping = aes_string(x = "Sample",y = "Abundance"))  + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.text = element_text(face = "italic", , size = 11),
        text = element_text(size = 12) , legend.position = "right") + 
  ggtitle("Relative Abundance: Frass") +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="fill") +
  facet_wrap(~Diet, scales="free_x") + scale_fill_manual(values = col_tax) + scale_color_manual(values = col_tax)

pct_frass



################ CONTROL REL ABUNDANCE (For supplemental Figure) ##################
#get top ten taxa fpr larvae
top20 <- names(sort(taxa_sums(control), decreasing=TRUE))[1:15]
ps.top20 <- transform_sample_counts(control, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
b = plot_bar(ps.top20, x="Sample", fill="Family") #+ facet_wrap(~Diet, scales="free_x")
taxa_bar_genus <- b+geom_bar(aes(fill=Family), stat="identity", position = "stack") +
  scale_fill_brewer(palette = "Set1") 

col_tax <- c("Acetobacter" = "#88CCEE", "Arsenophonus" = "#CC6677", "Clavulicium" ="#0072B2", "Endomicrobium" = "#DDCC77","Enterobacter" = "#117733", "Haemophilus" = "#332288", 
             "Klebsiella" = "#AA4499", "Lactobacillus" = "#44AA99", "Magnusiomyces" = "goldenrod" , "Myroides" = "#6699CC", "Pichia" = "#661100", "Serratia" = "#999933" ,
             "Rhodotorula" = "#D55E00", "Staphylococcus" = "#888888", "Streptococcus" = "plum1")


for_gg <- psmelt(ps.top20)


#make pct level bar chart
ncols <- 15
mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(ncols) 
pct_control <- ggplot(data = for_gg, 
                      mapping = aes_string(x = "Sample",y = "Abundance"))  + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
                                                                                   legend.text = element_text(face = "italic", size = 11),text = element_text(size = 12))+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="fill") +
  facet_wrap(~Diet, scales="free_x") + scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors)+ ggtitle("Relative Abundance: No Larvae Control")

pct_control
#######################################################

#make full figure
diversity_fig <- alpha_div_l+ alpha_div_f + pcoa_larvae + pcoa_frass + pct_larvae + pct_frass + patchwork::plot_layout(nrow = 3, ncol = 2, guides = "keep") + patchwork::plot_annotation(tag_levels = "A")
diversity_fig

#OUTPUT PLOT
# width = 10
# height = (9/16) * width
# svg(filename = "/mnt/md0/bsf_protein_study/microbiome_figures/diversity_figure.svg", width = 9, height = 9.5)
# diversity_fig
# dev.off()


#Output control plot
# png("/mnt/md0/bsf_protein_study/microbiome_figures/control_community.png", width = 7, height = 5, units = "in", res = 300 )
# pct_control
# dev.off()
pct_control

########## BETA DIVERSITY STATISTICAL ANALYSIS #########
#vegan analysis
factors_to_test <- c("Diet") #replace with the factors in your data you want to test
ps.prop_larvae <- transform_sample_counts(larvae, function(otu) otu/sum(otu))
ps.prop_frass <- transform_sample_counts(frass, function(otu) otu/sum(otu))
#1: extract the distance
braydist_larvae <- phyloseq::distance(ps.prop_larvae, "bray")
braydist_frass <- phyloseq::distance(ps.prop_frass, "bray")

#2: make the data frame
sampledf_larvae <- data.frame(sample_data(ps.prop_larvae))
ps.prop_larvae@sam_data<-sample_data(sampledf_larvae)

sampledf_frass <- data.frame(sample_data(ps.prop_frass))
ps.prop_frass@sam_data<-sample_data(sampledf_frass)


#look for homogeneity within groups. Insig. means that the groups are homogenous.
beta_larvae <- betadisper(braydist_larvae, sampledf_larvae$Diet)
permutest(beta_larvae)

beta_frass <- betadisper(braydist_frass, sampledf_frass$Diet)
permutest(beta_frass)

plot(beta_larvae)
plot(beta_frass)


## PERMANOVA
sig_larvae <- adonis2(braydist_larvae ~ Diet, data = sampledf_larvae, permutations = 100000)
sig_frass <- adonis2(braydist_frass ~ Diet, data = sampledf_frass, permutations = 100000)

pairwise.perm.manova(braydist_larvae, sampledf_larvae$Diet, p.method = "bonferroni")
pairwise.perm.manova(braydist_frass, sampledf_frass$Diet, p.method = "bonferroni")

######## CORRELATIONS ######

################ Correlate microbes and genes ###########
cts <- read.table(genes, header = T, sep = ",")
degs <- read.table(deg, header = F)

# ADD DEG GROUPS TO EACH DEG
degs$group <-NA
degs$group[1:165] <- "A"
degs$group[166:199] <- "D"
degs$group[200:351] <- "B"
degs$group[352:399] <- "C"
names(degs) <- c("GeneID","Group")

#get only differentiall expressed genes
colnames(cts)[1] <- "gene"
deg_sub <- cts[which(cts$gene %in% degs$GeneID),]
#deg_sub <- deg_sub[,-c(2:10)]
deg_long <- deg_sub %>% pivot_longer(!gene, values_to = "Normalized_Counts")
deg_wide <- deg_long %>% pivot_wider(names_from = gene, values_from = "Normalized_Counts")
names(deg_wide)[1] <- "Sample"


#Prepare data for correlations: Make matrices using diversity (shannon) and the normalized DEG counts.
row.names(shannon_larvae) <- gsub("_report_bracken_genuses", "", row.names(shannon_larvae))
x <- as.matrix(shannon_larvae[,1])
y<- as.matrix(deg_wide[,2:ncol(deg_wide)])

# Make correlation table
correlation.table <- associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1, p.adj.method = "bonferroni")
correlation.table$Correlation[which(correlation.table$p.adj > 0.05)] <- 0
#correlation.table <- subset(correlation.table, p.adj <= 0.05)
knitr::kable(head(correlation.table))
#write.csv(correlation.table, file = "/mnt/md0/bsf_protein_study/no_gainesville/correlation_table.csv", quote = F)


########### CORRELATE MICROBIAL DEGS TO BSF DEGS ############
mgc <- read.csv(mcts)
mde <- read.csv(mdegs)

#get only differentiall expressed genes
colnames(mgc)[1] <- "gene"
mdeg_sub <- mgc[which(mgc$gene %in% mde$genes),]
#deg_sub <- deg_sub[,-c(2:10)]
mdeg_long <- mdeg_sub %>% pivot_longer(!gene, values_to = "Normalized_Counts")
mdeg_wide <- mdeg_long %>% pivot_wider(names_from = gene, values_from = "Normalized_Counts")
names(mdeg_wide)[1] <- "Sample"
mdeg_wide <- subset(mdeg_wide, select =-c(hpl7_b_DN4673_c0_g1_i1))

#prepare data for correlations: make matrcies from DEG counts
x <- as.matrix(mdeg_wide[,2:ncol(mdeg_wide)])

y<- as.matrix(deg_wide[,2:ncol(deg_wide)])

# make correlation table
correlation.table <- associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1, p.adj.method = "bonferroni")
correlation.table$Correlation[which(correlation.table$p.adj > 0.05)] <- 0
#correlation.table <- subset(correlation.table, p.adj <= 0.05)
knitr::kable(head(correlation.table))
#write.csv(correlation.table, file = "/mnt/md0/bsf_protein_study/no_gainesville/correlation_table_DEG.csv", quote = F)


###### Read in alternative names for correlated transcripts
bsfdeg <- read.delim(bsf_degann, header = T, sep = "\t")
micdeg <- read.delim(mic_degann, header = F, sep = "\t")
correlation.table$Bac_description <- NA
correlation.table$BSF_description <- NA


for(i in 1:nrow(correlation.table)){
  correlation.table$Bac_description[i] <- micdeg$V9[which(correlation.table$X1[i] == micdeg$V1)] 
}

for(i in 1:nrow(correlation.table)){
  correlation.table$BSF_description[i] <- bsfdeg$Description[which(substr(correlation.table$X2[i], 4,nchar(toString(correlation.table$X2[i]))) == bsfdeg$NCBI.GeneID)] 
}

#make with line break
sttest <- unique(paste(correlation.table$X2, correlation.table$BSF_description, sep = "\n"))
#change to named list
names(sttest) <- gsub("\n.*", "",sttest)

#make plot
cor_bac_bsf = correlation.table %>%
  ggplot(aes(Bac_description, X2, fill = Correlation)) +
  geom_hline(yintercept = 1:50, color = "lightgray", alpha = 0.7, linetype = "dotted") +
  geom_vline(xintercept = 1:50, color = "lightgray", alpha = 0.7, linetype = "dotted") +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "steelblue1", high = "firebrick1", mid = "white", 
                       na.value = "grey", midpoint = 0, limit = c(-1,1), 
                       space = "Lab", name = NULL) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  geom_text(aes(Bac_description, X2, label = round(Correlation, 2)), color = "black") +
  #labs(x = NULL, y = NULL, title = "Distance (Filtering)") +
  theme_bw() +
  #geom_vline(xintercept = 8.5, color = "blue", linetype = "dashed") +
  #geom_vline(xintercept = 17.5, color = "blue", linetype = "dashed") +
  #geom_vline(xintercept = 25.5, color = "blue", linetype = "dashed") +
  #geom_hline(yintercept = 8.5, color = "blue", linetype = "dashed") +
  #geom_hline(yintercept = 17.5, color = "blue", linetype = "dashed") +
  #geom_hline(yintercept = 25.5, color = "blue", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9.5, hjust = 1, 
                                   face = "plain"),
        axis.text.y = element_text(size = 9.5, face = "plain"),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        #axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), 
        legend.position = "none") +
  ylab("BSF Gene")+
  xlab("Microbial Transcript\n(EggNOG Orthogroup)") +
  scale_y_discrete(labels = sttest) +
  coord_flip() + ggtitle("Spearman Correlations: Microbial DEGs vs. BSF DEGs")
cor_bac_bsf

# png(filename = "/mnt/md0/bsf_protein_study/microbiome_figures/correlation_microbeDEvsBSFde.png", width = 10, height = 5, res = 600, units = "in")
# cor_bac_fung
# dev.off()
