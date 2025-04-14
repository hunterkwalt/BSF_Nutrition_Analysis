### Written by Hunter K. Walt 
### Final Version: 14 April 2025

library(tximportData)
library(tximport)
library(DESeq2)
library(SummarizedExperiment)
library(knitr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library("pheatmap")
library("RColorBrewer")
library(ggrepel)
library(BiocParallel)
library(data.table)
library(grid)
library(gridExtra)
library(dplyr)
library(DEGreport)
library(clusterProfiler)
library(RCurl)
library(topGO)

#This script takes a few files as input. Please change these to where the files exist in your directory.
samps<-"/mnt/md0/bsf_protein_study/assemblies/txomes/individ_cdhit/cdhit_quantified/abundance/larvae_samps_no_gainesville.txt" #sample names
abun <- "/mnt/md0/bsf_protein_study/assemblies/txomes/larvae_txome/cdhit_90/final_quants/abundance" #directory with Kallisto abundance
egnog <- "/mnt/md0/bsf_protein_study/assemblies/txomes/larvae_txome/cdhit_90/annotation/best_hit_annotations_final.tsv" #Emapper annotations
gtog <-"/mnt/md0/bsf_protein_study/assemblies/txomes/larvae_txome/cdhit_90/annotation/gene2go.tsv" # gene ID to go term
g1_anns <- "/mnt/md0/bsf_protein_study/microbiome_figures/tables/annotations/group_1_larvae_annotations.tsv" #file with manually curated functional terms for cluster 1 genes
g2_anns <- "/mnt/md0/bsf_protein_study/microbiome_figures/tables/annotations/group_2_larvae_annotations.tsv" #file with manually curated functional terms for cluster 2 genes
g3_anns <- "/mnt/md0/bsf_protein_study/microbiome_figures/tables/annotations/group_3_larvae_annotations.tsv" #file with manually curated functional terms for cluster 3 genes

#make sample files 
l_samples= read.table(samps, header = F)
colnames(l_samples) <- "sample"

#path to abundance files
files_larvae <- file.path(abun, paste0(l_samples$sample, ".tsv"))


#name the files by sample
names(files_larvae) <- paste0(l_samples$sample)

#make files
txi.larvae.tsv <- tximport(files_larvae, type = "kallisto", txIn = T, txOut = T)
head(txi.larvae.tsv$counts)

#add metadata as factors
features <- factor(c(rep("5P:1C", 9),rep("1P:5C", 9), rep("1P:1C", 9)), levels = c("1P:5C","1P:1C","5P:1C") )

#make metadata
sampleTable_l <- data.frame(Diet = features)
rownames(sampleTable_l) <- colnames(txi.larvae.tsv$counts)

#make DEseq object
dds_l <- DESeqDataSetFromTximport(txi.larvae.tsv, sampleTable_l, ~Diet)
levels(dds_l$Diet)
dds_l <- DESeq(dds_l, parallel = T, BPPARAM = MulticoreParam(48))


######## LOAD IN ANNOTATION DATA FOR ENRICHMENT ANALYSIS ########

#read in data
#eggnog_output file
eggnog_data_larvae <- read.table(egnog, header = F, sep = "\t", quote = "")

#get columns 1 (query) and 9 (KO terms)
kegg_data_larvae <- eggnog_data_larvae[c(1,12)]
colnames(kegg_data_larvae) <- c("X.query","KEGG_ko")

# clean up by removing the "ko:" in front of every KO term
kegg_data_larvae$KEGG_ko <- gsub( "ko:", "", as.character(kegg_data_larvae$KEGG_ko))

# expand, since some genes/proteins will have multiple assigned KO terms
kegg_larvae <- data.table(kegg_data_larvae)
kegg_larvae <- kegg_larvae[, list(KEGG.KO = unlist(strsplit(KEGG_ko , ","))), by = X.query]

# select the needed columns
kegg_larvae_final <- kegg_larvae[,c(2,1)]


######## go #############
#read in data and functions
geneID2GO_l <- readMappings(file = gtog, IDsep = ",")
geneNames <- names(geneID2GO_l)

getGeneUniverse <- function(deseqRes){
  vgs = as.numeric(deseqRes$padj)
  names(vgs) = rownames(deseqRes)
  vgs
}

sel_pval <- function(allScore){ return(allScore < 0.05)}

######## KEGG PATHWAY ANNOTATIONS #########
#get columns 1 (query) and 9 (KO terms)
kegg_data_larvae_map <- eggnog_data_larvae[c(1,13)]
colnames(kegg_data_larvae_map) <- c("X.query","KEGG_Pathway")

# clean up by removing the "ko:" in front of every KO term
#kegg_data_larvae$KEGG_ko <- gsub( "ko:", "", as.character(kegg_data_larvae$KEGG_ko))

# expand, since some genes/proteins will have multiple assigned KO terms
kegg_larvae_map <- data.table(kegg_data_larvae_map)
kegg_larvae_map <- kegg_larvae_map[, list(KEGG_Pathway = unlist(strsplit(KEGG_Pathway , ","))), by = X.query]
kegg_larvae_map <- kegg_larvae_map[grep("map", kegg_larvae_map$KEGG_Pathway),]
# select the needed columns
kegg_larvae_final_map <- kegg_larvae_map[,c(2,1)]
#kegg_larvae_final_map$X.query <- gsub('.{2}$',"",kegg_larvae_final_map$X.query)


############## GO FUNCTIONS ##############
make_GO <- function(x, geneNames){
  #genes <- read.table(x, header = F)
  #genes <- genes$V1 
  #get gene counts
  #get the gene counts
  gl <- factor(as.integer(geneNames %in% x))
  names(gl) <- geneNames
  return(gl)
  #head(gl)
}


prep_go <- function(goEnrichment){
  goEnrichment$pval <- as.numeric(goEnrichment$pval)
  goEnrichment <- goEnrichment[goEnrichment$pval<0.05,]
  goEnrichment <- goEnrichment[,c("GO.ID","Term","pval")]
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
  return(goEnrichment)
}

plotgo <- function(goEnrichment, title){
  ggplot(goEnrichment[1:10,], aes(x=Term, y=-log10(as.numeric(pval)), fill = (as.numeric(pval)))) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    #xlab("Biological Process") +
    #ylab("Enrichment") +
    ggtitle(title) +
    #scale_y_continuous(breaks = round(seq(0, max(-log10(as.numeric(goEnrichment$pval))), by = 2), 1)) +
    scale_fill_gradient(low="red",high="blue", name = "p-value") +
    #theme_cowplot() +
    theme_bw(base_size=12) +
    theme(
      #legend.position='none',
      legend.background=element_rect(),
      plot.title=element_text(angle=0, size=11, face="bold", vjust=1),
      axis.text.x=element_text(angle=0, size=10, hjust=1.10),
      axis.text.y=element_text(angle=0, size=10),
      axis.title=element_text(size=14, face="bold"),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.key=element_blank(),     #removes the border
      legend.key.size=unit(0.25, "cm"),      #Sets overall area/size of the legend
      legend.text=element_text(size=10),  #Text size
      title=element_text(size=12)) +
    #guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip()
}

############## CLUSTER ANALYSIS ###############
# Perform likelihood ratio test for differential expression
dds <- DESeq(dds_l, test="LRT", reduced = ~ 1)
full_res <- results(dds)

#microbial counts
mcts <- counts(dds, normalized=TRUE)
#write.csv(mcts, file = "/mnt/md0/bsf_protein_study/no_gainesville/microbial_normalized_cts.csv")

# Quality Control (QC) - PCA
vsdata <- vst(dds, blind=FALSE)
p <- plotPCA(vsdata, intgroup = "Diet")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("Diet")])
colnames(df) <- "Diet"
row.names(df) <- colnames(vsdata)
pheatmap(assay(vsdata)[select,], cluster_rows=F, show_rownames=F,
         cluster_cols=T, annotation_col=df)


# Quality Control (QC) - Dispersion estimate
par(mfrow=c(1,1))
plotDispEsts(dds)

# Quality Control (QC) - MA plot
#plotMA(full_res, main="LRT")

#filter out significantly differentially expressed genes 
sigs <- na.omit(full_res)
sigs <- sigs[sigs$padj < 0.05, ]
sigs <- as.data.frame(sigs)
#sigs <- sigs[order(sigs$padj, decreasing = F), ]
sigs <- subset(sigs, baseMean > 5)

rlg <- rlog(dds, blind=FALSE)
rlog_mat <- assay(rlg)

selected_genes <- rownames(sigs)
cluster_rlog <- rlog_mat[selected_genes, ]

#make clusters
clusters_larvae <- degPatterns(cluster_rlog, metadata = sampleTable_l, time="Diet", col=NULL)
cluster_groups_larvae <- clusters_larvae$df
g1 <- subset(cluster_groups_larvae, cluster == 1)
g2 <- subset(cluster_groups_larvae, cluster == 2)
g3<- subset(cluster_groups_larvae, cluster == 3)
#g4 <- subset(cluster_groups, cluster == 4)
write.csv(cluster_groups_larvae, file = "/mnt/md0/bsf_protein_study/no_gainesville/microbe_DE.csv", row.names = TRUE)

#make factor levels correct
#clusters_larvae[["normalized"]]$Diet <- factor(clusters_larvae[["normalized"]]$Diet, levels = c("1P:5C", "1P:1C", "5P:1C"))

#make new names for plot facet labels
cluster_labels <- c("Cluster 1: 439 Transcripts", "Cluster 2: 22 Transcripts", "Cluster 3: 37 Transcripts")
names(cluster_labels) <- c("1", "2", "3")

# Use the data yourself for custom figures
larvae_patterns <- ggplot(clusters_larvae[["normalized"]],
                          aes(Diet, value)) +
  geom_boxplot(outliers = F, color = "#0D7760", fill = "green", alpha = 0.25) + ggtitle("Patterns of DE") + ylab("Z-score of Gene Abundance") +
  theme(plot.title = element_text(face = "bold"), axis.title = element_text(size = 12), 
        #axis.text.x = element_text(angle = 75, hjust = 1,1, vjust = 1, size = 12), 
        strip.text = element_text(size = 11)) +
  geom_jitter(color = "darkgray", alpha = 0.5, size = 1) +
  # change the method to make it smoother
  #geom_smooth(aes(group=Diet), method = "lm") + 
  facet_wrap(cluster~., labeller = labeller(cluster = cluster_labels), dir = "v") +geom_line(aes(group = genes), linewidth = 0.25, color = "gray", alpha = 0.25) 


######## MAKE PLOT OF TAXA DE ##########
#annotate lrt table
lrt_annotations <- eggnog_data_larvae[which(eggnog_data_larvae$V1 %in% cluster_groups_larvae$genes),] 

lrt_annotations$Group <- NA

lrt_annotations$Group[which(lrt_annotations$V1 %in% g1$genes)] <-  "Group1"
lrt_annotations$Group[which(lrt_annotations$V1 %in% g2$genes)] <-  "Group2"
lrt_annotations$Group[which(lrt_annotations$V1 %in% g3$genes)] <-  "Group3"

#write.csv(lrt_annotations, file = "/mnt/md0/bsf_protein_study/microbiome_figures/lrt_functional_table.csv")

lrt_tax <- lrt_annotations %>% ggplot( aes(x= fct_infreq(V6))) + geom_bar() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1, size = 14), text = element_text(size = 14)) +
  labs(x = "Most Specific Annotation", y = "Count") + ggtitle("")


lrt_tax_df <- data.frame(table(lrt_annotations$V6)[which(table(lrt_annotations$V6) > 1)])
lrt_tax_df$Var1 <- gsub(".*\\|", "", lrt_tax_df$Var1)
colnames(lrt_tax_df) <- c("Maximum Annotation", "Count")
lrt_tax_df<- lrt_tax_df[2:nrow(lrt_tax_df), ]
lrt_tax_df$Domain <- "Bacteria"
lrt_tax_df$Domain[which(lrt_tax_df$`Maximum Annotation` == "Fungi" |lrt_tax_df$`Maximum Annotation` == "Eukaryota" | lrt_tax_df$`Maximum Annotation` == "Opisthokonta" )] = "Eukaryota"

lrt_taxa_plot <- ggplot(lrt_tax_df, aes(x=reorder(`Maximum Annotation`, -Count), y = Count, fill = Domain)) + geom_col() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1, size = 12), text = element_text(size = 16)) +
  labs(x = "Most Specific Annotation", y = "DE Transcript Count") + ggtitle("DE Transcripts by Max Taxonomic Annotation") + scale_fill_manual(values = c("darkblue", "#0D7760"))

width = 10.2
height = (3/4)*width
#ggsave(lrt_taxa_plot,filename = "/mnt/md0/bsf_protein_study/microbiome_figures/lrt_taxa_barplot.pdf", width = width, height = height)


####### GROUP 1 GO ##########

g1_l <- make_GO(g1$genes,geneNames)

GOdata <- new("topGOdata", ontology= "BP", allGenes = g1_l,
              annot = annFUN.gene2GO, gene2GO = geneID2GO_l, nodeSize = 5)

allGO = usedGO(object = GOdata)
GOresult = runTest(GOdata, statistic = "fisher", algorithm = "classic")

allRes <- GenTable(GOdata, pval = GOresult, topNodes = length(allGO))

g1_go =  subset(allRes, as.numeric(pval) <= 0.05)

#write.table(g1_go, file = "/mnt/md0/bsf_protein_study/microbiome_figures/tables/group1_larvae_go_BP.tsv", sep = "\t", quote = F, row.names = T)

g1_bp <- prep_go(g1_go)
g1_bp_plot <- plotgo(g1_bp, "GO Enrichment: BP")

g1_mf <- prep_go(g1_go)
g1_mf_plot <- plotgo(g1_mf, "GO Enrichment: MF")

g1_cc <- prep_go(g1_go)
g1_cc_plot <- plotgo(g1_cc, "GO Enrichment: CC")

######## group 2 ###########
g2_l <- make_GO(g2$genes,geneNames)

GOdata <- new("topGOdata", ontology= "BP", allGenes = g2_l,
              annot = annFUN.gene2GO, gene2GO = geneID2GO_l, nodeSize = 5)

allGO = usedGO(object = GOdata)
GOresult = runTest(GOdata, statistic = "fisher", algorithm = "classic")

allRes <- GenTable(GOdata, pval = GOresult, topNodes = length(allGO))

g2_go =  subset(allRes, as.numeric(pval) <= 0.05)

####### Group 3 ###########

g3_l <- make_GO(g3$genes,geneNames)

GOdata <- new("topGOdata", ontology= "BP", allGenes = g3_l,
              annot = annFUN.gene2GO, gene2GO = geneID2GO_l, nodeSize = 5)

allGO = usedGO(object = GOdata)
GOresult = runTest(GOdata, statistic = "fisher", algorithm = "classic")

allRes <- GenTable(GOdata, pval = GOresult, topNodes = length(allGO))

g3_go =  subset(allRes, as.numeric(pval) <= 0.05)

#write.table(g3_go, file = "/mnt/md0/bsf_protein_study/microbiome_figures/tables/group3_larvae_go_BP.tsv", sep = "\t", quote = F, row.names = T)

g3_bp <- prep_go(g3_go)
g3_bp_plot <- plotgo(g3_bp, "GO Enrichment: BP")

g3_mf <- prep_go(g3_go)
g3_mf_plot <- plotgo(g3_mf, "GO Enrichment: MF")

g3_cc <- prep_go(g3_go)
g3_cc_plot <- plotgo(g3_cc, "GO Enrichment: CC")


########## KEGG Pathway #######

#get enrichment

#no enrichment in cluster1
# g1_pathway <- enricher(g1$genes, TERM2GENE=kegg_larvae_final_map, 
#                        pvalueCutoff = 0.05, pAdjustMethod = "BH", 
#                        qvalueCutoff = 0.05, minGSSize = 5)

g2_pathway <- enricher(g2$genes, TERM2GENE=kegg_larvae_final_map, 
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                       qvalueCutoff = 0.05, minGSSize = 5)


g3_pathway <- enricher(g3$genes, TERM2GENE=kegg_larvae_final_map, 
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                       qvalueCutoff = 0.05, minGSSize = 5)




#match KEGG IDs to common names
name_list <- c()
for (i in 1:length(row.names(subset(g2_pathway@result, p.adjust <= 0.05)))){
  #a <- enr_res_lpf@result[i,1]
  u <- "http://togows.dbcls.jp/entry/kegg-pathway/"
  r <- g2_pathway@result[i,1]
  l <- "/name"
  f <- paste0(u,r,l) 
  print(f)
  name <- getURL(f)
  name_list[i] <- name}
name_list <- gsub('.{1}$',"",name_list)
g2_pathway@result$Description[1:1:length(row.names(subset(g2_pathway@result, p.adjust <= 0.05)))] <- name_list


name_list <- c()
for (i in 1:length(row.names(subset(g3_pathway@result, p.adjust <= 0.05)))){
  #a <- enr_res_lpf@result[i,1]
  u <- "http://togows.dbcls.jp/entry/kegg-pathway/"
  r <- g3_pathway@result[i,1]
  l <- "/name"
  f <- paste0(u,r,l) 
  print(f)
  name <- getURL(f)
  name_list[i] <- name}
name_list <- gsub('.{1}$',"",name_list)
g3_pathway@result$Description[1:1:length(row.names(subset(g3_pathway@result, p.adjust <= 0.05)))] <- name_list

#barplot(g1_pathway)
barplot(g2_pathway)
barplot(g3_pathway)

######### KO enrichment #########

#Get KO enrichment
g1_ko <- enricher(g1$genes, TERM2GENE=kegg_larvae_final, 
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, minGSSize = 5)

g2_ko <- enricher(g2$genes, TERM2GENE=kegg_larvae_final, 
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, minGSSize = 5)


g3_ko <- enricher(g3$genes, TERM2GENE=kegg_larvae_final, 
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, minGSSize = 5)


# match KO terms to common terms
name_list <- c()
for (i in 1:length(row.names(subset(g1_ko@result, p.adjust <= 0.05)))){
  #a <- enr_res_lpf@result[i,1]
  u <- "http://togows.dbcls.jp/entry/kegg-pathway/"
  r <- g1_ko@result[i,1]
  l <- "/name"
  f <- paste0(u,r,l) 
  print(f)
  name <- getURL(f)
  name_list[i] <- name}
name_list <- gsub('.{1}$',"",name_list)
g1_ko@result$Description[1:1:length(row.names(subset(g1_ko@result, p.adjust <= 0.05)))] <- name_list



name_list <- c()
for (i in 1:length(row.names(subset(g2_ko@result, p.adjust <= 0.05)))){
  #a <- enr_res_lpf@result[i,1]
  u <- "http://togows.dbcls.jp/entry/kegg-pathway/"
  r <- g2_pathway@result[i,1]
  l <- "/name"
  f <- paste0(u,r,l) 
  print(f)
  name <- getURL(f)
  name_list[i] <- name}
name_list <- gsub('.{1}$',"",name_list)
g2_pathway@result$Description[1:1:length(row.names(subset(g2_ko@result, p.adjust <= 0.05)))] <- name_list


name_list <- c()
for (i in 1:length(row.names(subset(g3_ko@result, p.adjust <= 0.05)))){
  #a <- enr_res_lpf@result[i,1]
  u <- "http://togows.dbcls.jp/entry/kegg-pathway/"
  r <- g3_ko@result[i,1]
  l <- "/name"
  f <- paste0(u,r,l) 
  print(f)
  name <- getURL(f)
  name_list[i] <- name}
name_list <- gsub('.{1}$',"",name_list)
g3_ko@result$Description[1:1:length(row.names(subset(g3_ko@result, p.adjust <= 0.05)))] <- name_list

# barplot(g1_ko)
# barplot(g2_ko)
# barplot(g3_ko)

### Annotate DE genes
g1_annotations <- subset(lrt_annotations, Group == "Group1")
g1_annotations$Pathway <- ""
g1_annotations$KO <- ""
g1_annotations$BP <- ""


g2_annotations <- subset(lrt_annotations, Group == "Group2")
g2_annotations$Pathway <- ""
g2_annotations$KO <- ""
g2_annotations$BP <- ""

g3_annotations <- subset(lrt_annotations, Group == "Group3")
g3_annotations$Pathway <- ""
g3_annotations$KO <- ""
g3_annotations$BP <- ""

#loops to populate annotations for enrichment

### GO ###
for(x in 1:nrow(g1_go)){
  g1_annotations$BP[grep(g1_go$GO.ID[x], g1_annotations$V10)] = paste0(g1_annotations$BP[grep(g1_go$GO.ID[x], g1_annotations$V10)], g1_go[x,2], sep=",")
}

#no GO enrichment for g2
# for(x in 1:nrow(g2_go)){
#   g2_annotations$BP[grep(g2_go$GO.ID[x], g2_annotations$V10)] = paste0(g2_annotations$BP[grep(g2_go$GO.ID[x], g2_annotations$V10)], g2_go[x,2], sep=",")
# }

for(x in 1:nrow(g3_go)){
  g3_annotations$BP[grep(g3_go$GO.ID[x], g3_annotations$V10)] = paste0(g3_annotations$BP[grep(g3_go$GO.ID[x], g3_annotations$V10)], g3_go[x,2], sep=",")
}


### KO ###
for(x in 1:nrow(subset(g1_ko@result, p.adjust <= 0.05))){
  g1_annotations$KO[grep(g1_ko@result$ID[x], g1_annotations$V12, ignore.case = T)] = paste0(g1_annotations$KO[grep(g1_ko@result$ID[x], g1_annotations$V12, ignore.case = T)], g1_ko@result[x,2], sep=",")
}


for(x in 1:nrow(subset(g2_ko@result, p.adjust <= 0.05))){
  g2_annotations$KO[grep(g2_ko@result$ID[x], g2_annotations$V12, ignore.case = T)] = paste0(g2_annotations$KO[grep(g2_ko@result$ID[x], g2_annotations$V12, ignore.case = T)], g2_ko@result[x,2], sep=",")
}

for(x in 1:nrow(subset(g3_ko@result, p.adjust <= 0.05))){
  g3_annotations$KO[grep(g3_ko@result$ID[x], g3_annotations$V12, ignore.case = T)] = paste0(g3_annotations$KO[grep(g3_ko@result$ID[x], g3_annotations$V12, ignore.case = T)], g3_ko@result[x,2], sep=",")
}

### Pathway ###
# for(x in 1:nrow(subset(g1_pathway@result, p.adjust <= 0.05))){
#   g1_annotations$KO[grep(g1_pathway@result$ID[x], g1_annotations$V13, ignore.case = T)] = paste0(g1_annotations$KO[grep(g1_pathway@result$ID[x], g1_annotations$V13, ignore.case = T)], g1_ko@result[x,2], sep=",")
# }


for(x in 1:nrow(subset(g2_pathway@result, p.adjust <= 0.05))){
  g2_annotations$Pathway[grep(g2_pathway@result$ID[x], g2_annotations$V13, ignore.case = T)] = paste0(g2_annotations$Pathway[grep(g2_pathway@result$ID[x], g2_annotations$V13, ignore.case = T)], g2_pathway@result[x,2], sep=",")
}

for(x in 1:nrow(subset(g3_pathway@result, p.adjust <= 0.05))){
  g3_annotations$Pathway[grep(g3_pathway@result$ID[x], g3_annotations$V13, ignore.case = T)] = paste0(g3_annotations$Pathway[grep(g3_pathway@result$ID[x], g3_annotations$V13, ignore.case = T)], g3_pathway@result[x,2], sep=",")
}

#subset important
g1_sub <- g1_annotations[,c("V1","V5", "V6","V8","V21","V12","KO","V13","Pathway","V10","BP")]
g2_sub <- g2_annotations[,c("V1","V5", "V6","V8","V21","V12","KO","V13","Pathway","V10","BP")]
g3_sub <- g3_annotations[,c("V1","V5", "V6","V8","V21","V12","KO","V13","Pathway","V10","BP")]

#label_columns
colnames(g1_sub) <- c("genes", "Tax_Annotation", "Max_Annotation", "Description","PFAM", "KO_All","KO_sig","Pathway_all","Pathway_sig","GO_BP","GO_sig")
colnames(g2_sub) <- c("genes", "Tax_Annotation", "Max_Annotation", "Description","PFAM", "KO_All","KO_sig","Pathway_all","Pathway_sig","GO_BP","GO_sig")
colnames(g3_sub) <- c("genes", "Tax_Annotation", "Max_Annotation", "Description","PFAM", "KO_All","KO_sig","Pathway_all","Pathway_sig","GO_BP","GO_sig")

#write out tables for figure term labeling
#write.table(g1_sub, file="/mnt/md0/bsf_protein_study/microbiome_figures/tables/group_1_larvae_annotations.tsv", sep = "\t", quote = F, row.names = F)
#write.table(g2_sub, file="/mnt/md0/bsf_protein_study/microbiome_figures/tables/group_2_larvae_annotations.tsv", sep = "\t", quote = F, row.names = F)
#write.table(g3_sub, file="/mnt/md0/bsf_protein_study/microbiome_figures/tables/group_3_larvae_annotations.tsv", sep = "\t", quote = F, row.names = F)

#read in tables after labeling figure terms

####### THESE WERE READ IN AFTER MANUAL CURATION OF SUMMARY TERMS
g1_fig <- read.table(file = g1_anns, sep = "\t",header = T, quote = "")
g2_fig <-read.table(file = g2_anns, sep = "\t", header = T, quote = "")
g3_fig <-read.table(file = g3_anns,, sep = "\t", header = T, quote = "")

#add values for uninformative annotations and ambiguous annotations
g1_fig$Figure_term[g1_fig$Figure_term == ""] <- "Other"
g2_fig$Figure_term[g2_fig$Figure_term == ""] <- "Other"
g3_fig$Figure_term[g3_fig$Figure_term == ""] <- "Other"

g1_fig$Figure_term[g1_fig$Description == "-" & g1_fig$PFAM == "-"  ] <- "No Annotation"
g2_fig$Figure_term[g2_fig$Description == "-" & g2_fig$PFAM == "-" ] <- "No Annotation"
g3_fig$Figure_term[g3_fig$Description == "-" & g3_fig$PFAM == "-" ] <- "No Annotation"

#colors for figure
col_func <- c("Bacterial Outer Membrane Protein" = "#88CCEE", "Carbohydrate Metabolism" = "#CC6677", "Carbohydrate Transport" ="#0072B2", "Cellular Stress" = "#DDCC77","Cellular Respiration" = "#117733", "Nucleotide Metabolism" = "#332288", 
              "Other" = "darkgrey", "Protein Metabolism" = "#44AA99", "Translation" = "goldenrod" , "Transcription" = "plum1", "No Annotation" = "gray40", "Serratia" = "#999933" ,
              "Rhodotorula" = "#D55E00", "Staphylococcus" = "#888888", "Streptococcus" = "#6699CC")

### MAKE PIE CHARTS ##

library(ggrepel)
g1_fig_counts <- data.frame(table(g1_fig$Figure_term))
colnames(g1_fig_counts) <- c("Function", "Count")
g1_fig_counts$Function <- forcats::fct_relevel(g1_fig_counts$Function, "Other", after = Inf)
g1_fig_counts$Function <- forcats::fct_relevel(g1_fig_counts$Function, "No Annotation", after = Inf)
pie1 <- g1_fig_counts %>% ggplot( aes(x="", y=Count, fill = Function)) + geom_bar(stat="identity", width=1, color="white") + 
  coord_polar("y", start=0)+theme_void() + scale_fill_manual(values = col_func) + 
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 14)) + 
  geom_text(aes(x=1.75, label = Count), position = position_stack(vjust = 0.5)) # the x axis label had to be labeled 

g2_fig_counts <- data.frame(table(g2_fig$Figure_term))
colnames(g2_fig_counts) <- c("Function", "Count")
g2_fig_counts$Function <- forcats::fct_relevel(g2_fig_counts$Function, "Other", after = Inf)
g2_fig_counts$Function <- forcats::fct_relevel(g2_fig_counts$Function, "No Annotation", after = Inf)
pie2 <- g2_fig_counts %>% ggplot( aes(x="", y=Count, fill = Function)) + geom_bar(stat="identity", width=1, color="white") + 
  coord_polar("y", start=0)+theme_void() + scale_fill_manual(values = col_func) + 
  theme(legend.position = "none") + geom_text(aes(x=1.6, label = Count), position = position_stack(vjust = 0.5))

g3_fig_counts <- data.frame(table(g3_fig$Figure_term))
colnames(g3_fig_counts) <- c("Function", "Count")
g3_fig_counts$Function <- forcats::fct_relevel(g3_fig_counts$Function, "Other", after = Inf)
g3_fig_counts$Function <- forcats::fct_relevel(g3_fig_counts$Function, "No Annotation", after = Inf)
pie3 <- g3_fig_counts %>% ggplot( aes(x="", y=Count, fill = Function)) + geom_bar(stat="identity", width=1, color="white") + 
  coord_polar("y", start=0)+theme_void() + scale_fill_manual(values = col_func)+ 
  theme(legend.position = "none") + geom_text(aes(x=1.6, label = Count), position = position_stack(vjust = 0.5))


#rename first column of annotation table
colnames(g1_fig)[colnames(g1_fig) == 'Transcript'] <- "genes"
colnames(g2_fig)[colnames(g2_fig) == 'Transcript'] <- "genes"
colnames(g3_fig)[colnames(g3_fig) == 'Transcript'] <- "genes"

#join data together
cluster_normalized <- clusters_larvae[["normalized"]] 

cluster_2_normalized <- cluster_normalized %>%
  left_join(g1_fig %>% dplyr::select(genes, Figure_term), by = "genes") %>%
  mutate(colored = Figure_term) %>%
  dplyr::select(-Figure_term)



# Function to generate plot for a given cluster and sheet number
generate_plot <- function(cluster_num, cluster_ann, my_colors, show_x_axis = TRUE) {
  figure_curation <- cluster_ann
  
  cluster_normalized <- cluster_normalized %>%  # Ensure this dataframe is available in your environment
    filter(cluster == cluster_num)
  
  cluster_2_normal<- cluster_normalized %>%
    left_join(cluster_ann %>% dplyr::select(genes, Figure_term), by = "genes") %>%
    mutate(colored = Figure_term) %>%
    dplyr::select(-Figure_term)
  
  # Adjust diet treatment labels
  cluster_2_normal$Diet <- factor(cluster_2_normal$Diet, levels = c("1P:5C", "1P:1C", "5P:1C"),
                                  labels = c("1P:5C", "1P:1C", "5P:1C"))
  
  p <- ggplot(cluster_2_normal, aes(Diet, value)) + 
    geom_line(aes(group = genes, color = colored, y = value, x = Diet), alpha = 0.2) +  # Add transparent lines connecting same genes
    #geom_jitter(data = subset(cluster_2_normal, colored == "other"), aes(color = colored), size =.75,
    geom_jitter(data = subset(cluster_2_normal, colored == "Other" | colored == "No Annotation"), aes(color = colored), alpha = 0.4, size = 1.5) +
    geom_jitter(data = subset(cluster_2_normal, colored != "Other" & colored != "No Annotation"), aes(color = colored), alpha = 1, size = 1.5) + scale_color_manual(values=col_func) +
    #            position = position_jitter()) +  # Add grey jittered points behind
    #geom_jitter(data = subset(cluster_2_normal, colored != "other"), aes(color = colored), size =.75,
    #            position = position_jitter()) +  # Add colorful jittered points in front
    scale_color_manual(values = my_colors) +  # Manual color scale without labels
    geom_boxplot(fill = "transparent", color = "black", alpha = 0) + 
    labs(title = "", x = "Diet", y = "", caption = "") +  # Customize titles
    theme_minimal() +                                  # Apply a minimal theme
    theme(
      legend.position = "none",                    # Remove legend
      plot.margin = margin(0.5, 0.5, 0.5, 0.5),
      axis.title.x = element_text(size = 12.5),      # Adjust x-axis title size
      axis.text.x = element_text(size = 12)        # Adjust x-axis tick label size
    ) + 
    guides(color = guide_legend(override.aes = list(alpha = 1)))  # Ensure full opacity in legend
  
  return(p)
}


plot1 <- generate_plot(1, g1_fig, col_func) + 
  #ggtitle("Group A: 439 Transcripts")+ 
  theme(axis.title.x = element_blank())
plot2 <- generate_plot(2, g2_fig, col_func) + 
  #ylab("Z-Score of Gene Abundance") + 
  #ggtitle("Group B: 22 Transcripts") + 
  theme(axis.title.x = element_blank())
plot3 <- generate_plot(3, g3_fig, col_func) +
  #ggtitle("Group C: 37 Transcripts")
  theme(axis.title.x = element_blank())

#primary plot
plot1 + pie1 + plot2 + pie2 +plot3 + pie3 + patchwork::plot_layout(ncol = 2, nrow = 3, guides = "collect") + patchwork::plot_annotation(tag_levels = "A", title = "Microbial Patterns of DE")

# Adjust titles
yleft <- textGrob("Z-Score of Gene Abundance", rot = 90, gp = gpar(fontsize = 15, fontfamily = "Arial"))
top1 <- textGrob("Group A - 439 Transcripts", rot = 0, gp = gpar(fontsize = 18, fontfamily = "Arial"), x = unit(0.70, "npc"), just = "right")
top2 <- textGrob("Group B - 22 Transcripts", rot = 0, gp = gpar(fontsize = 18, fontfamily = "Arial"), x = unit(0.70, "npc"), just = "right")
top3 <- textGrob("Group C: 37 Transcripts", rot = 0, gp = gpar(fontsize = 18, fontfamily = "Arial"), x = unit(0.7, "npc"), just = "right")
#top4 <- textGrob("Cluster 4 - 35 Transcripts", rot = 0, gp = gpar(fontsize = 20, fontfamily = "Helvetica"), x = unit(0.7, "npc"), just = "right")

# Arrange plots and pie charts
pie1_nl <- pie1 + theme(legend.position = "none") #make piechart without legend
comb1 <- grid.arrange(plot1, pie1_nl, ncol = 2, nrow = 1, widths = c(1, 1), top = top1)
comb2 <- grid.arrange(plot2, pie2, ncol = 2, nrow = 1, widths = c(1, 1), top = top2)
comb3 <- grid.arrange(plot3, pie3, ncol = 2, nrow = 1, widths = c(1, 1), top = top3)


#get common legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#extract legend. This legend as all necessary labels for the pie charts, this may not always be true.
my_legend <- g_legend(pie1)

#final plot
grid.arrange(ggarrange(comb1,comb2,comb3, labels = c("A","B","C"), ncol = 1, nrow = 3, font.label = list(face = "plain")), left = yleft, right = my_legend)

#output plots
# svg('/mnt/md0/bsf_protein_study/microbiome_figures/final_figures/updated_names_microbial_DE.svg', height = 8, width = 9)
# grid.arrange(ggarrange(comb1,comb2,comb3, labels = c("A","B","C"), ncol = 1, nrow = 3, font.label = list(face = "plain")), left = yleft, right = my_legend)
# dev.off()
# 
# png('/mnt/md0/bsf_protein_study/microbiome_figures/final_figures/updated_names_microbial_DE.png', height = 8, width = 9, units = "in", res = 300)
# grid.arrange(ggarrange(comb1,comb2,comb3, labels = c("A","B","C"), ncol = 1, nrow = 3, font.label = list(face = "plain")), left = yleft, right = my_legend)
# dev.off()

