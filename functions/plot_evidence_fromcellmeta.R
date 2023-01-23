plot_evidence_fromcellmeta <- function(assigned, out_lambda = NULL) {
    library(ggplot2, quietly = TRUE)
    library(Seurat, quietly = TRUE)
    library(SeuratDisk, quietly = TRUE)
    library(dplyr, quietly = TRUE)
    library(tidyr, quietly = TRUE)
    library(stringr, quietly = TRUE)
    library(tidyverse, quietly = TRUE)
    if (!is_assigned(assigned)) {
        stop(
            paste(
                "No cell type has been assigned_seurat.",
                "Run 'assign_celltype' function first."
            )
        )
        return(NULL)
    }
    plt <-
        assigned@meta.data %>%
        tibble::rownames_to_column("cell") %>%
        distinct(celltype, marker, color.celltype) %>%
        mutate(
            marker = .$marker %>%
                str_split(",")
        ) %>%
        unnest(marker) %>%
        arrange(celltype) %>%
        mutate(
            marker = factor(
                marker,
                levels = as.data.frame(marker) %>%
                    distinct(marker) %>%
                    pull()
            )
        ) %>%
        dplyr::filter(
            !is.na(marker)
        ) %>%
        group_split(marker) %>%
        lapply(
            lambda(
                data:
                assigned %>%
                    set_default_assay("SCT") %>%
                    FeaturePlot(
                        features = data %>%
                            mutate(
                                marker = as.character(marker)
                            ) %>%
                            distinct(marker) %>%
                            pull(),
                        cols = c(
                            "#e9e9e9",
                            data %>%
                                distinct(color.celltype) %>%
                                pull() %>%
                                first()
                        ),
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
       return(assigned)
    }
}
