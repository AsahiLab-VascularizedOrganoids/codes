plot_expression_cds <- function(cds, features, color_low = "gray", color_high = "blue") { # nolint
    pacman::p_load(
        "monocle3", "Seurat", "tidyverse",
        "SeuratDisk", "scran", "igraph",
        "SeuratWrappers"
    )
    data <-
        SummarizedExperiment::assays(cds[row.names(subset(
            SummarizedExperiment::rowData(cds),
            rownames(cds) %in% c("SOX2", "TOP2A")
        )), ])$counts %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("cell") %>%
        tidyr::pivot_longer(
            cols = c(-cell), # nolint
            names_to = "gene",
            values_to = "expression"
        ) %>%
        dplyr::left_join(
            coordinate %>% # nolint
                tibble::rownames_to_column("cell"),
            by = "cell"
        )
    data %>%
        dplyr::group_split(gene) %>% # nolint
        lapply(
            function(df) {
                df %>%
                    ggplot2::ggplot(
                        ggplot2::aes(
                            x = UMAP_1, # nolint
                            y = UMAP_2, # nolint
                            color = expression
                        )
                    ) +
                    ggplot2::geom_point() +
                    ggplot2::scale_color_gradient2(
                        low = color_low,
                        high = color_high
                    ) +
                    ggplot2::theme_minimal() +
                    ggplot2::theme(
                        plot.background = ggplot2::element_rect(
                            fill = "white",
                            colour = "white"
                        ),
                        legend.justification = "bottom"
                    ) +
                    ggplot2::labs(
                        title = df %>%
                            dplyr::distinct(gene) %>% # nolint
                            dplyr::pull()
                    )
            }
        ) %>%
        patchwork::wrap_plots() %>%
        return()
}
