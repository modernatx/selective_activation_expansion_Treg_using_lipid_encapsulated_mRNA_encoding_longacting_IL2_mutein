#--- fastTopics v0.5-25
                       
library(fastTopics)

sce <- readRDS(file = "./data/sce.rds")
counts <- counts(sce)
pdata <- colData(sce)

# input to fastTopics is cell by gene
counts = t(counts)

for (i in seq(5,18)) {
    message("\n", "fitting k=", i)
    message("fit_poisson_nmf...")
    fit0 <- fit_poisson_nmf(tmp,k=i,init_method = "topicscore",numiter = 50,
                           method = "scd",verbose = FALSE,
                           control = list(extrapolate = TRUE))
    message("poisson2multinom...")
    fit_k <- poisson2multinom(fit0)
    saveRDS(fit_k, file = paste0("./output/fastTopics/fit_k_",i,".rds"))
    
    message("\n", "umap")    
    fit_k = readRDS(paste0("./output/fastTopics/fit_k_",i,".rds"))
    clust_assign = apply(fit_k$L, 1, which.max)

    fit.umap = umap(fit_k$L)
    saveRDS(fit.umap, paste0("./output/umap/umap_topics_",i,".rds"))
}    
    
#--- plot out UMAP
for (i in seq(5,18)) {
    fit_k = readRDS(paste0("./output/fastTopics/fit_k_",i,".rds"))
    clust_assign = apply(fit_k$L, 1, which.max)
    fit.umap = readRDS(paste0("./output/umap/umap_topics_",i,".rds"))
    
    L = fit_k$L
    df = data.frame(fit.umap$layout,
               L,
               LIB0000 = pdata_combined$LIB0000,
               condition_day = factor(pdata_combined$condition_day, 
                                      levels =c("UTX_0", "M_2", "M_4", "WT_2", "WT_4"),
                                      labels =c("UTX_0", "M_2", "M_4", "WT_2", "WT_4"),),
               condition_day_rep = pdata_combined$condition_day_rep,
               clust = clust_assign)
    df$LIB0000 = as.factor(df$LIB0000)
    df$clust = as.factor(df$clust)
    colnames(df)[1:2] = c("umap.1", "umap.2")

    png(file.path(paste0("./figs/umap_topics_k_",i,".png")), 
        width = 1000, heigh=950)
    pp = ggplot(df, aes(x = umap.1, y = umap.2, col = clust)) +
            geom_point(size = .6) + #facet_wrap(~condition_day, ncol = 2) +
            ylab("Dimension 2") + xlab("Dimension 1") +
            ggtitle(paste0("k=",i))
    print(pp)
    dev.off()
    
    png(file.path(paste0("./figs/umap_topics_k_",i,"_batch.png")), 
        width = 1200, heigh=1200)
    pp = ggplot(df, aes(x = umap.1, y = umap.2, col = condition_day)) +
            geom_point(size = .6) + facet_wrap(~condition_day_rep, ncol = 4) +
            ylab("Dimension 2") + xlab("Dimension 1") +
            ggtitle(paste0("k=",i))
    print(pp)
    dev.off()
    
    png(file.path(paste0("./figs/umap_topics_k_",i,"_batch_cluster_tx.png")), 
        width = 900, heigh=1200)
    pp = ggplot(subset(df, condition_day != "utx_0"), 
                aes(x = umap.1, y = umap.2, col = clust)) +
            geom_point(size = .6) + facet_wrap(~condition_day_rep, nrow = 4) +
            ylab("Dimension 2") + xlab("Dimension 1") +
            ggtitle(paste0("k=",i))
    print(pp)
    dev.off()

    png(file.path(paste0("./figs/umap_topics_k_",i,"_batch_cluster_utx.png")), 
        width = 700, heigh=600)
    pp = ggplot(subset(df, condition_day == "utx_0"), aes(x = umap.1, y = umap.2, col = clust)) +
            geom_point(size = .6) + facet_wrap(~condition_day_rep, nrow = 2) +
            ylab("Dimension 2") + xlab("Dimension 1") +
            ggtitle(paste0("k=",i))
    print(pp)
    dev.off()
    
                    
    png(file.path(paste0("./figs/umap_topics_k_",i,
                      "_batch_cluster_combined.png")), 
        width = 1300, heigh=900)
    pp = ggplot(df, aes(x = umap.1, y = umap.2, col = clust)) +
            geom_point(size = .6) + facet_wrap(~condition_day, nrow = 2) +
            ylab("Dimension 2") + xlab("Dimension 1") +
            ggtitle(paste0("k=",i))
    print(pp)
    dev.off()
}
    


             
             

                              
