#################################################################################
# 
# Figure S7 
# 
# Software: scran::cyclone v1.20.1 
# Outcome: Assign each single cell to a discrete cell cycle stage (G1, S, G2.M).
#################################################################################


library(scran)

#------ Figure S7B: Cell cycle stage assignment

# import mouse cycle cycle marker genes
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
package="scran"))

# import SCE object
sce = readRDS(file = "./data/sce.rds")
    
# get EMSEMBL gene ID
library(biomaRt)
ensembl_mus <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
biomart_res <- getLDS(attributes="external_gene_name", 
                      filters = "external_gene_name", 
                      values = rownames(sce), mart=ensembl_mus, 
                      attributesL="ensembl_gene_id", martL=ensembl_mus)

    
# Using Ensembl IDs to match up with the annotation in 'mm.pairs'.
rowData(sce) = data.frame(Symbol = rownames(sce), 
        Ensembl = biomart_res$Gene.stable.ID[match(rownames(sce), biomart_res$Gene.name)])
assignments <- cyclone(sce, 
                       mm.pairs, gene.names=rowData(sce)$Ensembl)
    
# Figure S7B
tab =  as.data.frame(table(assignments$phases, clust_assign_k_18))
ggplot(tab, aes(x = clust_assign_k_18, y = Freq, fill = Var1)) + 
  geom_bar(position="fill", stat="identity") +
  labs(title = "Cyclone predicted Cell cycle phase\n", x = "\nTopics", y = "Frequency\n") +
  theme_classic()
    