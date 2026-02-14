# ==============================================================================
# Part 1: Seurat & CellChat Analysis
# ==============================================================================

library(Seurat)
library(dplyr)
library(CellChat)
library(ggplot2)
library(patchwork)
library(future)

# 1. Load and Filter Data ------------------------------------------------------

# PLACEHOLDER: Update path to the Tabula Sapiens RDS file
seurat_object <- readRDS("path/to/tabular_sapiens.rds")

# Define tissues of interest
tissues_of_interest <- c("Blood", "Bone_Marrow", "Lymph_Node", "Spleen", "Thymus")

# Filter Seurat object
seurat_filtered <- subset(seurat_object, subset = organ_tissue %in% tissues_of_interest)
seurat_filtered <- NormalizeData(seurat_filtered)

# QC / Summary Stats
total_cells <- ncol(seurat_filtered)
cat(sprintf("Total cells (filtered): %d\n\n", total_cells))

# Counts (and percentages) per tissue
per_tissue <- seurat_filtered@meta.data %>%
  count(organ_tissue, name = "n") %>%
  arrange(desc(n)) %>%
  mutate(pct = round(100 * n / sum(n), 2))

print(per_tissue)

# Check for missing tissues
missing <- setdiff(tissues_of_interest, per_tissue$organ_tissue)
if (length(missing) > 0) {
  cat("\nNo cells found for:", paste(missing, collapse = ", "), "\n")
}

# 2. Initialize CellChat -------------------------------------------------------

data.input <- GetAssayData(seurat_filtered, slot = "data", assay = "RNA")
meta <- seurat_filtered@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_ontology_class")

# 3. Custom Database Curation --------------------------------------------------

# Manually input external association data (e.g., GWAS or Proteomics data)
df <- vroom("path/to/association_data.txt")

# Load CellChat Database
CellChatDB <- CellChatDB.human
CellChatDB_gene_info <- CellChatDB$geneInfo

# Merge external data with CellChat DB info
merged_df <- merge(df, CellChatDB_gene_info, by.x = "UniProt_ID", by.y = "EntryID.uniprot", all.x = TRUE)

# Filter interactions based on symbols in the external data
CellChatDB_filtered <- CellChatDB$interaction %>%
  filter(ligand.symbol %in% merged_df$Symbol | receptor.symbol %in% merged_df$Symbol)

# Manually correct/reverse direction for JAG1_CD46 interaction
if("CD46_JAG1" %in% rownames(CellChatDB_filtered)){
  CellChatDB_filtered["JAG1_CD46", ] <- CellChatDB_filtered["CD46_JAG1", ]
  CellChatDB_filtered["JAG1_CD46", "interaction_name"] <- "JAG1_CD46"
  CellChatDB_filtered["JAG1_CD46", "ligand"] <- "JAG1"
  CellChatDB_filtered["JAG1_CD46", "receptor"] <- "CD46"
  CellChatDB_filtered["JAG1_CD46", "interaction_name_2"] <- "JAG1_CD46"
  CellChatDB_filtered["JAG1_CD46", "ligand.symbol"] <- "JAG1"
  CellChatDB_filtered <- CellChatDB_filtered[-which(rownames(CellChatDB_filtered) == "CD46_JAG1"), ]
}

# Update CellChat object with filtered DB
CellChatDB.use <- CellChatDB
CellChatDB.use$interaction <- CellChatDB_filtered
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

# 4. Run CellChat Analysis -----------------------------------------------------

future::plan("multisession", workers = 2) # Enable parallel processing

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probability
cellchat <- computeCommunProb(cellchat, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 50)

# Extract network dataframe
df.net <- subsetCommunication(cellchat)

# Aggregate network
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Adjust P-values
p_values <- cellchat@netP$prob
adjusted_pvalues <- p.adjust(p_values, method = "BH")
cellchat@netP$prob <- adjusted_pvalues

# 5. CellChat Visualization ----------------------------------------------------

# Basic Bubble plot
netVisual_bubble(cellchat, remove.isolate = FALSE, font.size = 5)

# Bubble plot for specific pathways
netVisual_bubble(cellchat, signaling = c("BAFF","APRIL","CD40"), remove.isolate = FALSE)

# Custom Visualization of Top Interactions
df_top5 <- df.net %>%
  group_by(interaction_name) %>%
  top_n(5, wt = prob) %>%
  ungroup() %>%
  mutate(source_target = paste(source, "->", target))

# Plot
ggplot(df_top5, aes(x = interaction_name_2, y = source_target, fill = prob)) +
  geom_point(size = 3.75, shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(x = "Interaction Name", y = "Source -> Target", fill = "Prob") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_line(color = "grey", size = 0.2)
  )


# ==============================================================================
# Part 2: NicheNet Analysis
# ==============================================================================

library(nichenetr)
library(vroom)
library(tidyverse)

# 1. Load NicheNet Networks ----------------------------------------------------

organism <- "human"

if(organism == "human"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
} else {
  # Mouse options kept for flexibility
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
}

lr_network <- lr_network %>% distinct(from, to)

# 2. Define Receiver and Ligands -----------------------------------------------

# Define receiver cell type
seurat_filtered <- SetIdent(seurat_filtered, value = "cell_ontology_class")
receiver <- c("plasma cell")

# Get expressed genes in receiver
expressed_genes_receiver <- get_expressed_genes(receiver, seurat_filtered, pct = 0.05)
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

# Define ligands of interest (APRIL/BAFF)
potential_ligands <- c("TNFSF13","TNFSF13B")

# 3. Load Disease Signature (Target Genes) -------------------------------------

# PLACEHOLDER: Update path to your differential expression file
mm_scrna_seq <- vroom("path/to/mm_scrnaseq_deg.txt")

# Filter for genes of interest (e.g., MGUS stage)
mm_scrna_seq_filtered <- mm_scrna_seq %>%
  group_by(disease_stage) %>%
  arrange(`pval_adj (FDR)`) %>%
  ungroup()

mgusreceiver_genes <- mm_scrna_seq_filtered %>%
  filter(disease_stage == "MGUS") %>%
  pull(gene) %>%
  unique()

# Filter for genes present in ligand-target matrix
mgusreceiver_genes <- mgusreceiver_genes %>% .[. %in% rownames(ligand_target_matrix)]

# 4. Predict Ligand Activity ---------------------------------------------------

ligand_activities <- predict_ligand_activities(geneset = mgusreceiver_genes,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>%
  arrange(-aupr_corrected) %>%
  mutate(rank = rank(desc(aupr_corrected)))

# 5. Visualization -------------------------------------------------------------

# Get top weighted links
active_ligand_target_links_df <- potential_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = mgusreceiver_genes,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

# Prepare visualization
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33)

order_ligands <- intersect(potential_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])

# Generate Heatmap
make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke", high = "purple")
