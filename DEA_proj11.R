

# ----------------------------- #
#        LOAD LIBRARIES        #
# ----------------------------- #
library(DESeq2)
library(tidyverse)
library(readxl)
library(topGO)
library(clusterProfiler)
library(enrichplot)
library(data.table)
library(lintr)
library(lattice)
library(ggplot2)

# ----------------------------- #
#     LOAD COUNT DATASET       #
# ----------------------------- #
count_data <- read_excel("/Users/michaelsmith/Desktop/Project_11/grad_students/counts_data.xlsx", sheet = 1)
count_data <- as.data.frame(count_data)
rownames(count_data) <- count_data[, 1]
count_data <- count_data[, -1]

# Filter out lowly expressed genes
keep <- rowSums(count_data) >= 10
count_data <- count_data[keep, ]

# ----------------------------- #
#     LOAD SAMPLE METADATA     #
# ----------------------------- #
colData <- read_excel("/Users/michaelsmith/Desktop/Project_11/grad_students/meta_data.xlsx")
colData <- as.data.frame(colData)
rownames(colData) <- colData[, 1]
colData$Diet <- factor(colData$Diet, levels = c("carb_biased", "balanced", "protein_biased"))

# Ensure column names in count_data match rownames in colData
stopifnot(all(colnames(count_data) %in% rownames(colData)))
stopifnot(all(colnames(count_data) == rownames(colData)))

# ----------------------------- #
#   DESEQ2 DIFFERENTIAL TEST   #
# ----------------------------- #
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = colData, design = ~ Diet)
dds$Diet <- relevel(dds$Diet, ref = "balanced")

# Normalize counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, "/Users/michaelsmith/Desktop/Project_11/no_gainesville/normalized_counts.csv")

# Filter low count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Regularized log transformation
rlg <- rlog(dds, blind = FALSE)
rlog_mat <- assay(rlg)

# Likelihood Ratio Test
dds <- DESeq(dds, test = "LRT", reduced = ~ 1)
res <- results(dds)
write.csv(as.data.frame(res), "/Users/michaelsmith/Desktop/Project_11/no_gainesville/larva_all_LRT_genes.csv")

# ----------------------------- #
#       QUALITY CONTROL        #
# ----------------------------- #
vsdata <- vst(dds, blind = FALSE)
plotPCA(vsdata, intgroup = "Diet")
plotDispEsts(dds)
plotMA(res, main = "LRT")

# ----------------------------- #
#      SIGNIFICANT GENES       #
# ----------------------------- #
sigs <- na.omit(res)
sigs <- subset(sigs, padj < 0.05 & baseMean > 5)
write.csv(sigs, "/Users/michaelsmith/Desktop/Project_11/no_gainesville/larva_sig_LRT_genes.csv")

# ----------------------------- #
#      CLUSTERING GENES        #
# ----------------------------- #
top_sig <- sigs %>% arrange(padj) %>% head(1000)
selected_genes <- rownames(top_sig)
cluster_rlog <- rlog_mat[selected_genes, ]
write.csv(cluster_rlog, "/Users/michaelsmith/Desktop/Project_11/no_gainesville/larva_sig_normalized_counts.csv")

clusters <- degPatterns(cluster_rlog, metadata = colData, time = "Diet", col = NULL)
cluster_groups <- clusters$df
write.csv(cluster_groups, "/Users/michaelsmith/Desktop/Project_11/no_gainesville/larava_clusters.csv")

# ----------------------------- #
#       GO ENRICHMENT          #
# ----------------------------- #
gene_to_go_mapping_pre <- read_excel("/Users/michaelsmith/Desktop/Project_11/benchtop/GO_universe.xlsx")
gene_to_go_mapping <- lapply(gene_to_go_mapping_pre$GO, function(x) unlist(strsplit(x, ",")))
names(gene_to_go_mapping) <- gene_to_go_mapping_pre$gene

go_data <- read_excel("/Users/michaelsmith/Desktop/Project_11/benchtop/GO_universe.xlsx")
go_data$Group <- as.numeric(go_data$Group)
go_sorted <- go_data %>% select(gene, Group) %>% arrange(Group)
topgo <- go_sorted$Group
names(topgo) <- go_sorted$gene

for (go_category in c("BP")) {
  my_go_data <- new("topGOdata",
                    description = paste("GOtest", go_category, sep = "_"),
                    ontology = go_category,
                    geneSel = function(x) x == 1,
                    allGenes = topgo,
                    gene2GO = gene_to_go_mapping,
                    annot = annFUN.gene2GO,
                    nodeSize = 5
  )
  
  result_weight_fisher <- runTest(my_go_data, algorithm = "weight01", statistic = "fisher")
  result_weight_output <- GenTable(my_go_data,
                                   weight_fisher = result_weight_fisher,
                                   orderBy = "weight_fisher",
                                   topNodes = length(score(result_weight_fisher))
  )
  
  result_weight_output <- result_weight_output[result_weight_output$weight_fisher < 0.05, ]
  result_weight_output$category <- go_category
  result_weight_output$Term <- gsub(" ", "_", result_weight_output$Term)
  
  write.csv(result_weight_output, "/Users/michaelsmith/Desktop/Project_11/benchtop/GO_bp.csv", row.names = FALSE)
}

# ----------------------------- #
#     KEGG ENRICHMENT (KO)     #
# ----------------------------- #
ko_lines <- readLines("/Users/michaelsmith/Downloads/ko00001.keg")
ko_df <- data.frame(Lines = ko_lines)
ko_filtered <- ko_df[grep("^D", ko_df$Lines), ]
ko_terms <- data.frame(
  term = substr(ko_filtered$Lines, 8, 13),
  name = sub(".*;\\s", "", ko_filtered$Lines)
)

eggnog_data <- read_excel("/Users/michaelsmith/Desktop/Project_11/final_anal/universe_eggnog.xlsx", sheet = 1)
colnames(eggnog_data) <- as.character(eggnog_data[2, ])
eggnog_data <- eggnog_data[-(1:2), ]
eggnog_data <- eggnog_data[1:(nrow(eggnog_data) - 3), ]

kegg_data <- eggnog_data[c(1, 12)]  # KO column
kegg_data$KEGG_ko <- gsub("ko:", "", as.character(kegg_data$KEGG_ko))
kegg_data <- data.table(kegg_data)
kegg_data <- kegg_data[, list(KEGG_ko = unlist(strsplit(KEGG_ko, ","))), by = query]
kegg_final <- kegg_data[, c(2, 1)]

# Enrichment for KO
protein_IDs <- read_excel("/Users/michaelsmith/Desktop/Project_11/benchtop/GO_universe.xlsx")$DEGs
protein_IDs <- as.character(protein_IDs)

enr_res <- enricher(protein_IDs, TERM2GENE = kegg_final, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10)
enr_res@result$Description <- ko_terms$name[match(enr_res@result$ID, ko_terms$term)]

write.table(enr_res, "/Users/michaelsmith/Desktop/Project_11/benchtop/korth.csv")
jpeg("/Users/michaelsmith/Desktop/Project_11/benchtop/korth.jpg", width = 9, height = 9, units = "in", res = 1200)
dotplot(enr_res, showCategory = 30)
dev.off()

# ----------------------------- #
#    KEGG ENRICHMENT (PATH)    #
# ----------------------------- #
kpath_lines <- readLines("/Users/michaelsmith/Downloads/ko00001.keg")
kpath_df <- data.frame(Lines = kpath_lines)
kpath_filtered <- kpath_df[grep("^C", kpath_df$Lines), ]
kpath_terms <- data.frame(
  term = substr(kpath_filtered$Lines, 6, 10),
  name = substr(kpath_filtered$Lines, 12, nchar(kpath_filtered$Lines) - 14)
)
kpath_terms$term <- paste0("map", kpath_terms$term)

kp_data <- eggnog_data[c(1, 13)]  # KEGG Pathway column
kp_data <- data.table(kp_data)
kp_data <- kp_data[, list(KEGG_Pathway = unlist(strsplit(KEGG_Pathway, ","))), by = query]
kp_data <- kp_data[!grepl("^k", KEGG_Pathway)]
kp_final <- kp_data[, c(2, 1)]

enr_res <- enricher(protein_IDs, TERM2GENE = kp_final, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10)
enr_res@result$Description <- kpath_terms$name[match(enr_res@result$ID, kpath_terms$term)]

write.table(enr_res, "/Users/michaelsmith/Desktop/Project_11/benchtop/kpath.csv")
jpeg("/Users/michaelsmith/Desktop/Project_11/benchtop/kpath.jpg", width = 9, height = 9, units = "in", res = 1200)
dotplot(enr_res, showCategory = 30)
dev.off()

