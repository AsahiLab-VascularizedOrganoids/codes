plot_trajectory <- function(seurat,
                            trajectory_file_path,
                            features, out_lambda = NULL,
                            combine = TRUE, ncol = NULL,
                            nrow = NULL,
                            arrow_style = arrow(
                                length = unit(0.30, "cm"), type = "closed"
                            ),
                            curve_width = 0.5,
                            curve_color = "#002d70",
                            curvature = 0.2,
                            legend_size = 5,
                            viridis_option = "D",
                            plot_cluster = FALSE,
                            cluster_pre_label = "") {
    library(Seurat)
    library(SeuratDisk)
    library(SingleCellPipeline)
    library(slingshot)
    library(tidyverse)
    library(TrajectoryUtils)
    library(dplyr)
    library(ggplot2)
    trajectory_info <-
        readRDS(
            trajectory_file_path
        )
    plts <-
        features %>%
        as.list() %>%
        lapply(
            function(feature) {
                if (feature == "celltype") {
                    plt <-
                        seurat %>%
                        DimPlot(
                            raster = FALSE,
                            group.by = "seurat_clusters"
                        )
                    data <-
                        plt$data %>%
                        tibble::rownames_to_column("cell") %>%
                        left_join(
                            seurat@meta.data %>%
                                tibble::rownames_to_column("cell") %>%
                                select(cell, celltype, color.celltype),
                            by = "cell"
                        )
                } else {
                    plt <-
                        seurat %>%
                        set_default_assay("SCT") %>%
                        FeaturePlot(
                            features = feature,
                            raster = FALSE
                        )
                    dimplt <-
                        seurat %>%
                        DimPlot(
                            raster = FALSE,
                            group.by = "seurat_clusters"
                        )
                    data <-
                        plt$data %>%
                        tibble::rownames_to_column("cell") %>%
                        left_join(
                            dimplt$data %>%
                                tibble::rownames_to_column("cell") %>%
                                select(cell, seurat_clusters),
                            by = "cell"
                        )
                }
                data <-
                    data %>%
                    mutate(
                        coordinate = select(
                            .,
                            starts_with("UMAP")
                        ) %>%
                            rowMeans(),
                        seurat_clusters = as.double(seurat_clusters) - 1
                    ) %>%
                    left_join(
                        trajectory_info$pseudotime %>%
                            as.data.frame() %>%
                            tibble::rownames_to_column("cell") %>%
                            rename_all(
                                vars(
                                    c("cell", "pseudotime")
                                )
                            ),
                        by = "cell"
                    )
                path_data <-
                    trajectory_info$lineages %>%
                    lapply(function(lineage) {
                        tibble::tibble(
                            from = lineage[1:length(lineage) - 1],
                            to = lineage[2:length(lineage)]
                        ) %>%
                            return()
                    }) %>%
                    bind_rows() %>%
                    distinct(from, to) %>%
                    mutate(
                        from = as.double(from),
                        to = as.double(to)
                    ) %>%
                    left_join(
                        data %>%
                            select(seurat_clusters, pseudotime, coordinate) %>%
                            group_by(seurat_clusters) %>%
                            summarize_all(mean),
                        by = c(from = "seurat_clusters")
                    ) %>%
                    dplyr::rename(
                        from_x = pseudotime,
                        from_y = coordinate
                    ) %>%
                    left_join(
                        data %>%
                            select(seurat_clusters, pseudotime, coordinate) %>%
                            group_by(seurat_clusters) %>%
                            summarize_all(mean),
                        by = c(to = "seurat_clusters")
                    ) %>%
                    dplyr::rename(
                        to_x = pseudotime,
                        to_y = coordinate
                    )
                result <-
                    data %>%
                    ggplot(
                        aes(
                            x = pseudotime,
                            y = coordinate,
                            color = !!sym(feature)
                        )
                    ) +
                    geom_point(
                        size = 0.4
                    ) +
                    geom_curve(
                        data = path_data,
                        aes(
                            x = from_x,
                            y = from_y, xend = to_x, yend = to_y
                        ),
                        colour = curve_color,
                        size = curve_width,
                        curvature = curvature,
                        arrow = arrow_style
                    ) +
                    theme_minimal() +
                    theme(
                        legend.justification = "bottom",
                        plot.background =
                            element_rect(
                                fill = "white", colour = "white"
                            ),
                        panel.background =
                            element_rect(
                                colour = "white"
                            )
                    ) +
                    ggtitle(feature)
                if (feature == "celltype") {
                    result <-
                        result +
                        guides(
                            color = guide_legend(
                                override.aes = list(size = legend_size)
                            )
                        ) +
                        scale_color_manual(
                            values =
                                data %>%
                                    arrange(celltype) %>%
                                    distinct(color.celltype) %>%
                                    pull()
                        )
                } else {
                    result <-
                        result +
                        viridis::scale_color_viridis(
                            option = viridis_option
                        )
                }
                if (plot_cluster) {
                    annotationdata <-
                        data %>%
                        select(
                            seurat_clusters, pseudotime,
                            coordinate
                        ) %>%
                        group_by(seurat_clusters) %>%
                        summarize_all(mean)

                    for (seurat_cluster in annotationdata$seurat_clusters) {
                        target_data <-
                            annotationdata %>%
                            filter(seurat_clusters == seurat_cluster)
                        result <-
                            result +
                            annotate(
                                "text",
                                x = target_data$pseudotime,
                                y = target_data$coordinate,
                                label = paste(cluster_pre_label, seurat_cluster)
                            )
                    }
                }
                return(result)
            }
        )
    if (combine) {
        plts <-
            plts %>%
            patchwork::wrap_plots(
                ncol = ncol,
                nrow = nrow
            )
    }
    if (is.null(out_lambda)) {
        return(plts)
    } else {
        out_lambda(plts)
        return(seurat)
    }
}
