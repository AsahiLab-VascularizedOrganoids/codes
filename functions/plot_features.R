plot_features <- function(seurat, assay = SeuratObject::DefaultAssay(seurat), manual_markers = NULL, out_lambda = NULL, cols = NULL) { # nolint
    library(tidyverse)
    library(Seurat)
    if (is.null(manual_markers)) {
        manual_markers <-
            seurat@meta.data %>%
            tibble::rownames_to_column("cell") %>%
            dplyr::distinct(celltype, marker)
    }
    if (is.null(cols)) {
        cols <- c(
            "#e9e9e9",
            "#00109c"
        )
    }
    plt <-
        manual_markers %>%
        mutate(
            marker = .$marker %>%
                str_split(",")
        ) %>%
        unnest(marker) %>%
        filter(
            !is.na(marker)
        ) %>%
        arrange(celltype) %>%
        mutate(
            marker = factor(
                marker,
                levels = as.data.frame(marker) %>%
                    distinct(marker) %>%
                    pull()
            )
        ) %>%
        group_split(marker) %>%
        lapply(
            lambda(
                data:
                seurat %>%
                    set_default_assay(assay) %>%
                    FeaturePlot(
                        features = data %>%
                            mutate(
                                marker = as.character(marker)
                            ) %>%
                            distinct(marker) %>%
                            pull(),
                        cols = cols,
                        order = TRUE
                    ) +
                    ggtitle(
                        data$marker,
                        subtitle =
                            data %>%
                                select(celltype) %>%
                                pull() %>%
                                str_c(collapse = ", ")
                    ) +
                    theme(
                        plot.title =
                            element_text(
                                hjust = 0, vjust = 0,
                                face = "italic"
                            ),
                        legend.justification = "bottom"
                    )
            )
        ) %>%
        patchwork::wrap_plots()
    if (is.null(out_lambda)) {
        return(plt)
    } else {
        out_lambda(plt)
    }
}
