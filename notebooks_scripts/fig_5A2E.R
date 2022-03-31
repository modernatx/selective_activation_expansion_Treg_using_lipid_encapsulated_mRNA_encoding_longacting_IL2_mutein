#################################################################################
# 
# Figure 5 A-E
#
#################################################################################

#-- Prepare for UMAP visualization

# Import filtered expression data
sce = readRDS(file = "./data/sce.rds")
    

library("Rtsne")
library("uwot")
library("ggplot2")
library(umap)

# pip install geosketch
geosketch <- import('geosketch')


# import fastTopics and UMAP results
pdata <- colData(sce)
fit_k = readRDS("./output/fastTopics/fit_k_18.rds")
fit.umap = readRDS(paste0("./output/umap/umap_topics_18.rds"))
L_sub = fit_k$L[which(pdata$macro != "Discard"),]
pdata_sub = pdata[which(pdata$macro != "Discard"),]
custom.config = umap.defaults

custom.config$random_state = 222
sketch.indices <- geosketch$gs(L_sub, as.integer(nrow(L_sub)*0.2), one_indexed = TRUE)
sketch.indices <- as.numeric(sketch.indices)

fit.umap = umap(L_sub[sketch.indices,], config = custom.config)

fit.densmap <- densmap(L_sub[sketch.indices,], dens_frac = 0.5, dens_lambda = 0.5)

df = data.frame(fit.umap$layout,
                LIB0000 = pdata_sub[sketch.indices,]$LIB0000,
               condition_day = factor(pdata_sub[sketch.indices,]$condition_day, 
                                      levels =c("utx_0", "M_2", "M_4", "WT_2", "WT_4"),
                                      labels =c("utx_0", "M_2", "M_4", "WT_2", "WT_4")),
               condition_day_rep = pdata_sub[sketch.indices,]$condition_day_rep,
               clust = pdata_sub[sketch.indices,]$clust_assign_k_18,
               clust_update = pdata_sub[sketch.indices,]$clust_assign_k_18_update,
               macro = pdata_sub[sketch.indices,]$macro,
               macro_label = pdata_sub[sketch.indices,]$macro_label)

    df$LIB0000 = as.factor(df$LIB0000)
    df$clust = as.factor(df$clust)
    colnames(df)[1:2] = c("umap.1", "umap.2")


#--- Figure 5A
#- UMAP with the big families

df$macro_label = factor(df$macro_label, 
                        levels = c("Tregs", "Activated T cells", "Other T cells", "g/d NKT", "Discard"))
pp = ggplot(df, aes(x = umap.1, y = umap.2, col = macro_label)) +
    geom_point(size = .1, alpha = .9) + #facet_wrap(~condition_day, ncol = 2) +
    ylab("Dimension 2") + xlab("Dimension 1") +
    ggtitle("UMAP") + theme_classic() + 
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

#cairo_ps(filename = "./figs/fig_5A_umap_all.eps",
#         width = 6, height = 5, pointsize = 12,
#         fallback_resolution = 300)
#pp
#dev.off()


#--- Figure 5B
df = pdata[which(pdata$macro_label != "Discard"),]

df$macro_label = droplevels(df$macro_label)
tab = table(df$condition_day, df$macro_label)
tab_plot = data.frame(freq = as.vector(tab))
tab_plot$condition_day = rep(rownames(tab), length(unique(df$macro_label)))
tab_plot$macro = factor(rep(colnames(tab), each = nrow(tab)))
tab_plot$condition_day = factor(tab_plot$condition_day,
                               levels = c("WT_4", "WT_2", "utx_0","M_2","M_4"))
tab_plot$macro = droplevels(tab_plot$macro)
tab_plot$macro = factor(tab_plot$macro, 
                        levels = c("Activated T cells", "Other T cells", "Treg", "g/d NKT"))

library(scales)
show_col(hue_pal()(4))

pp = ggplot(tab_plot, aes(fill=macro, y=freq, x=condition_day)) + 
    geom_bar(position="fill", stat="identity") +
    ggtitle("Macro cluster proportions") + 
    xlab("") + ylab("Proportion per condition/day") +
    scale_fill_manual(values = hue_pal()(4)[c(2,3,1,4)]) +
    theme_classic()

#cairo_ps(filename = "./figs/fig_5B_umap_all_freq.eps",
#         width = 6, height = 5, pointsize = 12,
#         fallback_resolution = 300)
#pp
#dev.off()




#--- Figure 5C
diff_count_res_k_18 = readRDS("./output/fastTopics/diff_count_res_k_18.rds")

markers_diffexp = c("Egr2", "Cxcr5", "Bcl6", "Atp5c1", "Cox4i1", "Ndufv3", "Gpr83", 
                    "Ikzf2", "Lrrc32", "Nusap1", "Stmn1", "Ccna2", "Pimreg", "Top2a", 
                    "Gzmb", "Icos", "Il1rl1")
markers_diffexp = str_to_sentence(markers_diffexp)
markers_diffexp = markers_diffexp[order(markers_diffexp)]

dmat = diff_count_res_k_18$beta[rownames(diff_count_res_k_18$beta) %in% markers_diffexp,]

cls = c(c("5","6","16","17","18","15","1"))
micro = c("aTreg", "mTreg", "OxTreg", "OxTreg", "ST2 Treg.1", "ST2 Treg.2", "Tfr")
dmat = dmat[,paste0("k",cls)]

library(pheatmap)
library(viridis)


# List with colors for each annotation.
cols = RColorBrewer::brewer.pal(n = 8, name = 'Dark2')
mat_colors <- list(group = c(cols[2], cols[1], cols[3], cols[3], cols[6], cols[4], cols[5]))
names(mat_colors$group) <- micro

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(dmat, n = 11)

library(sparklyr)
dmat_fac = dmat %>% as.data.frame() %>%
    mutate(k5_c = as.numeric(cut(k5, breaks = mat_breaks, include.lowest = TRUE)),
           k6_c = as.numeric(cut(k6, breaks = mat_breaks, include.lowest = TRUE)),
           k16_c = as.numeric(cut(k16, breaks = mat_breaks, include.lowest = TRUE)),
           k17_c = as.numeric(cut(k17, breaks = mat_breaks, include.lowest = TRUE)),
           k18_c = as.numeric(cut(k18, breaks = mat_breaks, include.lowest = TRUE)),
           k15_c = as.numeric(cut(k15, breaks = mat_breaks, include.lowest = TRUE)),
           k1_c = as.numeric(cut(k1, breaks = mat_breaks, include.lowest = TRUE))) %>%
    select(c('k5_c','k6_c','k16_c','k17_c','k18_c','k15_c','k1_c')) %>%
    as.matrix()
    
mat_col <- data.frame(group = micro)
rownames(mat_col) <- colnames(dmat_fac)
    

pp = pheatmap(dmat_fac,
              cluster_rows = TRUE, cluster_cols = FALSE,
              color = rev(plasma(10)),
              annotation_col    = mat_col,
              annotation_colors = mat_colors,
              fontsize=10, fontsize_rows = 10, fontsize_cols = 10)    

#cairo_ps(filename = "./fig_5C_Treg_heatmap_selected_degenes.eps",
#         width = 6, height = 7, pointsize = 12,
#         fallback_resolution = 150)
#pp
#dev.off()



#--- Figure 5D
df = readRDS(file = "./output/umap/df_geosub_20pnt.rds")
cols = RColorBrewer::brewer.pal(n = 8, name = 'Dark2')
df$condition_day = factor(df$condition_day, 
                          levels = c("utx_0", "M_2", "M_4", "WT_2", "WT_4"),
                             labels = c("UTX", "M_2", "M_4", "WT_2", "WT_4"))
pp = df %>% ggplot(aes(x = umap.1, y = umap.2, col = clust_update)) +
        geom_point(size = .5, alpha = .5) + #facet_wrap(~condition_day, ncol = 3) +
        ylab("Dimension 2") + xlab("Dimension 1") +
        theme_classic() +
        scale_color_manual(values = c("13" = scales::alpha("gray85", .5), 
                                      "3" = scales::alpha("gray85", .5), 
                                      "12_14" = scales::alpha("gray85", .5), 
                                      "11" = scales::alpha("gray85", .5), 
                                      "7_9" = scales::alpha("gray85", .5), 
                                      "10" = scales::alpha("gray85", .5),
                                      "4" = scales::alpha("gray85", .5),
                                      "6" = cols[1], "5" = cols[2], "16_17" = cols[3], 
                                      "15" = cols[4], "1" = cols[5], "18" = cols[6])) +
         guides(color = "none")
tmp_pp = df %>% filter(macro =="Macro3") %>%
    ggplot(aes(x = umap.1, y = umap.2, col = clust_update)) +
        geom_point(size = .5, alpha = .5) + theme_classic() + 
        scale_color_manual(name = "Tregs", 
                           values = c("6" = cols[1], "5" = cols[2], "16_17" = cols[3], 
                                      "15" = cols[4], "1" = cols[5], "18" = cols[6]))

#cairo_ps(filename = "./fig_5D_umap_Treg_all.eps",
#         width = 6, height = 5, pointsize = 12,
#         fallback_resolution = 300)
#tmp_pp
#dev.off()

