plot_gsea_emap <- function(deg_df, category_num = 50, p_val_adj_cutoff = 0.01) { # nolint
    library(dplyr)
    library(tidyverse)
    library(org.Hs.eg.db, quietly = TRUE)
    library(AnnotationDbi)
    library(ggplot2)
    create_genelist <- function(df) {
        genelist <- df$avg_log2FC
        names(genelist) <- df$entrez
        return(genelist)
    }
    if (typeof(deg_df) == data.frame) {
        deg_df %<>%
            tibble::rownames_to_column("gene")
    }
    deg_df %>%
        mutate(
            entrez = mapIds(
                org.Hs.eg.db, .$gene, "ENTREZID",
                "SYMBOL"
            )
        ) %>%
        dplyr::filter(
            (!isNA(entrez)) & (p_val_adj < 0.01)
        ) %>%
        dplyr::select(entrez, avg_log2FC) %>%
        arrange(desc(avg_log2FC)) %>%
        create_genelist() %>%
        clusterProfiler::gseGO(
            OrgDb = org.Hs.eg.db,
            ont = "ALL", verbose = TRUE
        ) %>%
        enrichplot::pairwise_termsim() %>%
        enrichplot::emapplot(
            showCategory = 120,
            color = "enrichmentScore"
        ) +
        theme(
            legend.justification = "bottom"
        ) +
        scale_fill_gradient2(
            low = "#0000b2",
            mid = "white",
            high = "#e20000"
        ) %>%
        return()
}
