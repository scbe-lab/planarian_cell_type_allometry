# SIZES DGE ANALYSIS

# Setup
dir <- "/mnt/sda/alberto/colabos/sizes/"
setwd(dir)

source("code/R_code/r_functions/sourcefolder.R")
sourceFolder("code/R_code/r_functions/",recursive = TRUE)

library(Matrix)
library(topGO)
library(readr)
library(ggplot2)
library(ggrepel)
library(ggvenn)
library(openxlsx)

# Load necessary data
sizes_Idents <- read.delim2("data/sizes_Idents.csv", sep = "," , header = TRUE)
sizes_genes <- read.delim2("data/sizes_genes.csv", sep = "," , header = TRUE)
sizes_leiden_col <- read.delim2("/mnt/sda/alberto/projects/smed_cisreg/data/std_ref_cellAnnotation/leiden_3_colors.csv", sep = ",")
sizes_Idents$color <- translate_ids(x=sizes_Idents$leiden_3,dict = sizes_leiden_col)

sizes_ctypes <- unique(sizes_Idents[,c("leiden_3","leiden_3_names","broad_names","color")])
sizes_ctypes <- merge(
  sizes_ctypes,
  as.data.frame(table(sizes_Idents$leiden_3)),
  by.x = 1,
  by.y = 1,
  all.x = TRUE
)

rosetta <- read.delim2("~/projects/smed_rink_gene_annot/20230915_Smed_Rink_Simplified_Annotation_Table.tsv",header = TRUE,sep="\t")

# Load the matrix
sizes_X <- readMM(file = "data/matrix/matrix.mtx.gz")
colnames(sizes_X) <-  read.table("data/matrix/barcodes.tsv")[,1]
rownames(sizes_X) <-  read.table("data/matrix/features.tsv")[,1]
sizes_X <- sizes_X[rownames(sizes_X) %in% sizes_genes$X,colnames(sizes_X) %in% sizes_Idents$X]

# Create table of counts per condition
sizes_psbulk_cond_rep <- pseudobulk_cond_rep(
  x = sizes_X,
  identities = sizes_Idents$leiden_3_names,
  conditions = sizes_Idents$Size,
  replicates = sizes_Idents$Library
)

# Clean sampletable
sizes_sampletable <- clean_sampletable(sizes_psbulk_cond_rep$sampletable)
sizes_matrix <- 
  sizes_psbulk_cond_rep$matrix[
    rownames(sizes_psbulk_cond_rep$matrix) %in% rosetta$gene[rosetta$gene_type == "hconf"],
    colnames(sizes_psbulk_cond_rep$matrix) %in% sizes_sampletable$sample
  ]

# SINGLE CELL ANALYSIS OVER ALL CELL TYPES
sizes_DGE_all <- list()
for(i in unique(sizes_sampletable$ctype)){
  print(paste0("starting with cell type ",i))
  sizes_DGE_all[[i]] <-
    deseq_pseudobulk(
      count_matrix = sizes_matrix,
      samples_info = sizes_sampletable[,-2],
      celltype = i,
      filter_by = "pvalue", p_threshold = 0.05,
      contrast_info = c("condition","L","S"),
      plot_results = FALSE, min_passing_samples = 2,
      min_counts_per_sample = 1,
      keep_dubious = FALSE
    )
  print(paste0("done cell type ",i))
}

# METAPLOT DATA FRAME
sizes_diffreg <- data.frame(
  num_diff = sapply(
    sizes_DGE_all,
    function(x){
      a = x$diffgenes
      if(is.na(a[1])){
        b = 0
      } else{
        b = length(a)
      }
      return(b)
    }
  )
)

sizes_diffreg <- 
  merge(
    sizes_diffreg,
    sizes_ctypes,
    by.x = 0,
    by.y = "leiden_3_names",
    all.x = TRUE
  )
colnames(sizes_diffreg)[1] <- "ctype"
sizes_diffreg <- sizes_diffreg[-c(grep("unannotated",sizes_diffreg$ctype)),]

## PANEL A: META PLOT NUM DIFF GENES AND CLUSTER SIZE
p1 <- ggplot(sizes_diffreg, aes(x = log(Freq), y = num_diff, label = ctype, color = color)) +
  geom_point() +
  labs(
    title = "Sizes cluster size vs # DiffReg",
    x = "cluster size (log n cells)",
    y = "no. diffreg genes"
  ) +
  scale_color_identity()+
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  geom_text_repel(cex=3) + labs(title = "Cluster size vs # DiffReg")

p1

# PANEL B: venn diagram basal & epidermal
sizes_DGE_basal_epid <-
  list(
    `basal cells` = sizes_DGE_all$`basal cells`$diffgenes,
    epidermis = sizes_DGE_all$epidermis$diffgenes
  )

p2 <- ggvenn(
  sizes_DGE_basal_epid,
  fill_color = c(sizes_ctypes$color[sizes_ctypes$leiden_3_names == "basal cells"], sizes_ctypes$color[sizes_ctypes$leiden_3_names == "epidermis"]),
  fill_alpha = 0.5,
  stroke_color = "grey80",
  stroke_alpha = 0.1,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 4,
  text_color = "black",
  text_size = 4,
  label_sep = ",",
  count_column = NULL,
  show_outside = c("auto", "none", "always"),
  auto_scale = FALSE
)
p2

basal_diffgenes <- data.frame(
  gene = sizes_DGE_all$`basal cells`$diffgenes,
  type = ifelse(
    sizes_DGE_all$`basal cells`$diffgenes %in%
      rownames(
        sizes_DGE_all$`basal cells`$res[
          sizes_DGE_all$`basal cells`$res$log2FoldChange > 0,
        ]
      ),
    "up","down"
  )
)

epid_diffgenes <- data.frame(
  gene = sizes_DGE_all$epidermis$diffgenes,
  type = ifelse(
    sizes_DGE_all$epidermis$diffgenes %in%
      rownames(
        sizes_DGE_all$epidermis$res[
          sizes_DGE_all$epidermis$res$log2FoldChange > 0,
        ]
      ),
    "up","down"
  )
)

# PANEL C: VOLCANO PLOTS
sizes_DGE_basal <-
  deseq_pseudobulk(
    count_matrix = sizes_matrix,
    samples_info = sizes_sampletable[,-2],
    celltype = "basal cells",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","L","S"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
p3 <- sizes_DGE_basal$res %>%
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "red")) +
  guides(color = "none") + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "basal cells DEGs")
p3

sizes_DGE_epid <-
  deseq_pseudobulk(
    count_matrix = sizes_matrix,
    samples_info = sizes_sampletable[,-2],
    celltype = "epidermis",
    filter_by = "pvalue", p_threshold = 0.05,
    contrast_info = c("condition","L","S"),
    plot_results = TRUE, min_passing_samples = 2,
    min_counts_per_sample = 1,
    keep_dubious = FALSE
  )
p4 <- sizes_DGE_epid$res %>% 
  ggplot(aes(x = log2FoldChange, y = -log(pvalue), color = factor(pvalue < 0.05))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "red")) +
  guides(color = "none") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75))+
  ggtitle(label = "epidermis DEGs")
p4

# PANEL D: GENE ONTOLOGIES
smed_id_GO <- readMappings("~/projects/smed_cisreg/outputs/gene_annotation/smed_GOs.tsv")

list_diffregs_forGO <-
  list(
    basal_down = sizes_DGE_basal$diffgenes[
      sizes_DGE_basal$diffgenes %in% 
        rownames(sizes_DGE_basal$res[sizes_DGE_basal$res$log2FoldChange<0,])
    ],
    basal_up =  sizes_DGE_basal$diffgenes[
      sizes_DGE_basal$diffgenes %in% 
        rownames(sizes_DGE_basal$res[sizes_DGE_basal$res$log2FoldChange>0,])
    ],
    epid_down =  sizes_DGE_epid$diffgenes[
      sizes_DGE_epid$diffgenes %in% 
        rownames(sizes_DGE_epid$res[sizes_DGE_epid$res$log2FoldChange<0,])
    ],
    epid_up = sizes_DGE_epid$diffgenes[
      sizes_DGE_epid$diffgenes %in% 
        rownames(sizes_DGE_epid$res[sizes_DGE_epid$res$log2FoldChange>0,])
    ]
  )

sizes_diffreg_GOs <- getGOs(
  genelist = list_diffregs_forGO,
  gene_universe= rownames(sizes_matrix),
  gene2GO = smed_id_GO
)

sizes_diffreg_GOs$GOtable$basal_down$brief_term <- abrevi(sizes_diffreg_GOs$GOtable$basal_down$Term)
p5 <- 
  ggplot(
    sizes_diffreg_GOs$GOtable$basal_down[1:14,], 
    aes(
      x = reorder(brief_term, -log10(classicFisher)),
      y = Significant/Expected,
      fill = -log10(classicFisher)
    )
  )+
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors = sequential_hcl(7,"PinkYl")[1:6], name = "-log10(pvalue)") +
  coord_flip() +
  theme_classic() +
  labs(x = "GO term", y = "Significant/Expected ratio",
       title = "Basal Downreg GOs") + scale_y_reverse()+
  scale_x_discrete(position = "top")+
  geom_text(aes(label=brief_term, y = max(Significant/Expected)/2))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
        )

sizes_diffreg_GOs$GOtable$basal_up$brief_term <- abrevi(sizes_diffreg_GOs$GOtable$basal_up$Term)
p6 <- 
  ggplot(
    sizes_diffreg_GOs$GOtable$basal_up[1:15,], 
    aes(
      x = reorder(brief_term, -log10(classicFisher)),
      y = Significant/Expected,
      fill = -log10(classicFisher)
    )
  )+
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors = sequential_hcl(7,"PinkYl")[1:6], name = "-log10(pvalue)") +
  coord_flip() +
  theme_classic() +
  labs(x = "GO term", y = "Significant/Expected ratio",
       title = "Basal Upreg GOs")+
  geom_text(aes(label=brief_term, y = max(Significant/Expected)/2))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

sizes_diffreg_GOs$GOtable$epid_down$brief_term <- abrevi(sizes_diffreg_GOs$GOtable$epid_down$Term)
p7 <- 
  ggplot(
    sizes_diffreg_GOs$GOtable$epid_down[1:14,], 
    aes(
      x = reorder(brief_term, -log10(classicFisher)),
      y = Significant/Expected,
      fill = -log10(classicFisher)
    )
  )+
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors = sequential_hcl(7,"PinkYl")[1:6], name = "-log10(pvalue)") +
  coord_flip() +
  theme_classic() +
  labs(x = "GO term", y = "Significant/Expected ratio",
       title = "Epidermis Downreg GOs") + scale_y_reverse()+
  scale_x_discrete(position = "top")+
  geom_text(aes(label=brief_term, y = max(Significant/Expected)/2))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

sizes_diffreg_GOs$GOtable$epid_up$brief_term <- abrevi(sizes_diffreg_GOs$GOtable$epid_up$Term)
p8 <- 
  ggplot(
    sizes_diffreg_GOs$GOtable$epid_up[1:15,], 
    aes(
      x = reorder(brief_term, -log10(classicFisher)),
      y = Significant/Expected,
      fill = -log10(classicFisher)
    )
  )+
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors =  sequential_hcl(7,"PinkYl")[1:6], name = "-log10(pvalue)") +
  coord_flip() +
  theme_classic() +
  labs(x = "GO term", y = "Significant/Expected ratio",
       title = "Epidermis Upreg GOs")+
  geom_text(aes(label=brief_term, y = max(Significant/Expected)/2))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

# PDFs of FIGURES
pdf(
  file = "graphics/6A.pdf",
  width = 7,
  height = 6
)
print(p1)
dev.off()

pdf(
  file = "graphics/6B.pdf",
  width = 5,
  height = 5
)
print(p2)
dev.off()

pdf(
  file = "graphics/6C1.pdf",
  width = 4,
  height = 7
)
print(p3)
dev.off()

pdf(
  file = "graphics/6C2.pdf",
  width = 4,
  height = 7
)
print(p4)
dev.off()

pdf(
  file = "graphics/6D1.pdf",
  width = 5,
  height = 4
)
print(p5)
dev.off()

pdf(
  file = "graphics/6D2.pdf",
  width = 5,
  height = 4
)
print(p6)
dev.off()

pdf(
  file = "graphics/6D3.pdf",
  width = 5,
  height = 4
)
print(p7)
dev.off()

pdf(
  file = "graphics/6D4.pdf",
  width = 5,
  height = 4
)
print(p8)
dev.off()

# Save lists of diffregs
write.table(
  basal_diffgenes,
  file = "outputs/basal_diffreg.tsv",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)
write.table(
  epid_diffgenes,
  file = "outputs/epid_diffreg.tsv",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)

wb <- createWorkbook()
for(i in names(sizes_DGE_all)){
  r <- sizes_DGE_all[[i]]$res
  r <- r[order(r$pvalue),]
  r$DEG <- ifelse(r$pvalue < 0.05, TRUE, FALSE)
  r$up_down <- ifelse(
    r$DEG == TRUE & r$log2FoldChange > 0,
    "up-regulated",
    ifelse(
      r$DEG == TRUE & r$log2FoldChange < 0,
      "down-regulated",
      "none"
    )
  )
  sheetname <- substr(i,1,30)
  addWorksheet(wb, sheetName = sheetname)
  writeData(wb,sheetname,r,rowNames = TRUE)
}
saveWorkbook(wb,file = "outputs/all_DGE_analyses.xlsx")

# Save objects
saveRDS(
  sizes_DGE_all,
  "outputs/sizes_DGE_all_celltypes.RDS"
)

saveRDS(
  sizes_psbulk_cond_rep,
  "outputs/sizes_DGE_pseudobulk_output.RDS"
)

saveRDS(
  sizes_sampletable,
  "outputs/sizes_DGE_pseudobulk_cond_rep_sampletable.RDS"
)

saveRDS(
  sizes_matrix,
  "outputs/sizes_DGE_pseudobulk_cond_rep_matrix.RDS"
)
