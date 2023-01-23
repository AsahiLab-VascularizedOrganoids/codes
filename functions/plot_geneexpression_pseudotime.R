plot_geneexpression_pseudotime <- function(assigned_seurat, features, assay = "SCT", pseudotime_slingshot = NULL, target_celltype = c(), target_samples = c(), scale.free = TRUE, facet.wrap = TRUE, layout.nrows = NULL, layout.ncols = NULL, out_lambda = NULL) { # nolint
    library(Seurat)
    library(tidyverse)
    library(dplyr)
    library(ggplot2)
    if (length(target_samples) == 0) {
        target_samples <-
            assigned_seurat@meta.data %>%
            dplyr::distinct(sample) %>%
            pull()
    }
    if (length(target_celltype) == 0) {
        target_celltype <-
            assigned_seurat@meta.data %>%
            dplyr::distinct(celltype) %>%
            pull()
    }
    if (scale.free) {
        scale.param <- "free"
    } else {
        scale.param <- "fixed"
    }
    if (!is.null(pseudotime_slingshot)) {
        pseudo_data <- pseudotime_slingshot
    } else {
        pseudo_data <-
            assigned_seurat@reductions$umap@misc$trajectory$pseudotime
        if (is.null(pseudo_data)) {
            stop("Cound not find pseudotime data. Please give pseudotime_slingshot or run estimate_trajectory")
        }
    }
    plt <-
        features %>%
        as.list() %>%
        lapply(
            function(target_feature) {
                plt <-
                    assigned_seurat %>%
                    set_default_assay(assay = assay) %>%
                    FeaturePlot(
                        features = target_feature,
                        raster = FALSE
                    )
                plt$data %>%
                    rename(
                        expression = dplyr::contains(target_feature)
                    ) %>%
                    tibble::rownames_to_column("cell") %>%
                    left_join(
                        pseudo_data,
                        by = "cell"
                    ) %>%
                    dplyr::left_join(
                        assigned_seurat@meta.data %>%
                            tibble::rownames_to_column("cell") %>%
                            select(cell, sample, celltype, color.sample),
                        by = "cell"
                    ) %>%
                    mutate(
                        expression = as.double(expression),
                        feature = target_feature
                    ) %>%
                    filter(
                        (sample %in% target_samples) &
                            (celltype %in% target_celltype)
                    ) %>%
                    select(
                        cell, celltype, sample, color.sample,
                        pseudotime, feature, expression
                    )
            }
        ) %>%
        bind_rows() %>%
        ggplot(
            aes(
                x = pseudotime,
                y = expression,
                color = sample,
                label = sample
            )
        ) +
        stat_smooth(
            method = "gam", na.rm = TRUE
        ) +
        facet_wrap(
            scales = "free",
            . ~ feature
        ) +
        scale_color_manual(
            values =
                assigned_seurat@meta.data %>%
                    arrange(sample) %>%
                    distinct(color.sample) %>%
                    pull()
        ) +
        scale_y_continuous(limits = c(0, NA)) +
        theme(
            legend.justification = "bottom",
            strip.text = element_text(
                face = "italic",
                colour = "black"
            ),
            panel.background = element_rect(
                fill = "white",
                color = "black"
            ),
            axis.ticks = element_line(
                colour = "black"
            ),
            panel.grid = element_line(
                colour = "#dbdbdb"
            ),
            strip.background = element_rect(
                fill = "white",
                colour = "black"
            )
        ) +
        xlim(0, NA) +
        ylim(0, NA) +
        xlab("Pseudotime") +
        ylab("Expression")
    if (facet.wrap) {
        plt <- plt + facet_wrap(. ~ feature,
            scales = scale.param,
            nrow = layout.nrows, ncol = layout.ncols
        )
    } else {
        plt <- plt + facet_grid(. ~ feature, scales = scale.param)
    }

    if (is.null(out_lambda)) {
        return(plt)
    }
    out_lambda(plt)
    return(assigned_seurat)
}
