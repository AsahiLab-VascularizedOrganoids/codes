estimate_trajectory <- function(seurat, root_cell_ident, output_path, stretch = 0) { # nolint
    message("@ Start executing estimate_trajectory >>>")
    library(ggplot2, quietly = TRUE)
    library(Seurat, quietly = TRUE)
    library(SeuratDisk, quietly = TRUE)
    library(dplyr, quietly = TRUE)
    library(slingshot, quietly = TRUE)
    library(SingleCellExperiment, quietly = TRUE)
    if (is.null(seurat@reductions$umap)) {
        stop(paste(
            "[ERROR] Could not find dimention reduction: UMAP",
            "Please reduct dimention with UMAP"
        ))
    }
    if (is.null(seurat@meta.data$seurat_clusters)) {
        stop(paste(
            "[ERROR] Could not find seurat cluster.",
            "Please execute clustering analysis."
        ))
    }
    message(
        "Start estimating trajetory"
    )
    if (typeof(root_cell_ident) == "double") {
        message(
            paste(
                "[Info] Use seurat clusters:",
                root_cell_ident
            )
        )
    } else if (typeof(root_cell_ident) == "character") {
        message(
            paste(
                "[Info] Use cell type:",
                root_cell_ident
            )
        )
        root_cell_ident <-
            seurat@meta.data %>%
            dplyr::filter(
                celltype %in% root_cell_ident
            ) %>%
            dplyr::distinct(seurat_clusters) %>%
            pull()
    }
    message(
        "[Info] Estimating trajectory pathway with slingshot..."
    )
    sds <-
        slingshot(
            Embeddings(seurat, "umap"),
            clusterLabels = seurat$seurat_clusters,
            start.clus = root_cell_ident,
            stretch = 0
        )
    message(
        "<<< [Info] Finished."
    )
    trajectory_info <-
        list(
            root_cell_ident = root_cell_ident,
            pseudotime = TrajectoryUtils::averagePseudotime(sds),
            lineages = sds@metadata$lineages
        )
    if (!stringr::str_ends(output_path, ".trajc")) {
        output_path <-
            paste0(
                output_path,
                ".trajc"
            )
    }
    trajectory_info %>%
        saveRDS(output_path)
    message(
        paste(
            "<<< [Info] Saved to",
            output_path
        )
    )
    return(seurat)
}
