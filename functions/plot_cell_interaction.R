plot_cell_interaction <- function(means_output, cell_meta, target_gene_pattern, title, output, width = 15, height = 15, units = "cm") {
    library(circlize)
    library(tidyverse)
    data <-
        read_tsv(
            means_output,
            show_col_types = FALSE
        ) %>%
        filter(
            str_detect(interacting_pair, target_gene_pattern)
        ) %>%
        select(
            dplyr::contains("|")
        ) %>%
        stack() %>%
        mutate(
            ind = as.character(
                ind
            )
        ) %>%
        tidyr::separate(
            ind,
            c("celltype1", "celltype2"),
            sep = "\\|"
        ) %>%
        select(celltype1, celltype2, values)
    metadata <-
        read_tsv(
            cell_meta,
            show_col_types = FALSE
        ) %>%
        filter(
            celltype %in% (
                data %>%
                    select(celltype1) %>%
                    pull()
            )
        ) %>%
        select(celltype, color)
    circos.clear()
    circos.par(start.degree = 90)
    png(output, width = width, height = height, units = units, res = 600)
    par(
        oma = c(0, 0, 1, 0)
    )
    data %>%
        chordDiagram(
            grid.col = structure(
                metadata$color,
                names = metadata$celltype
            ),
            scale = TRUE,
            link.arr.type = "big.arrow"
        )
    title(title)
    dev.off()
    circos.clear()
}