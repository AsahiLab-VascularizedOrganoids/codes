estimate_cell_interaction <- function(seurat, assay, workdir) { # nolint
    library(Seurat)
    library(SeuratDisk)
    library(SingleCellPipeline)
    library(biomaRt)
    library(org.Hs.eg.db, quietly = TRUE)
    library(tidyverse)
    maptoens <- function(df) {
        org.Hs.eg.db %>%
            AnnotationDbi::mapIds(
                keys = df %>%
                    dplyr::select("Gene") %>%
                    dplyr::pull(),
                keytype = "SYMBOL",
                column = "ENSEMBL"
            ) %>%
            tibble::as_tibble() %>%
            pull() %>%
            return()
    }
    setwd(workdir)
    if (!dir.exists("./resources")) {
        dir.create("./resources")
    }
    if (!dir.exists("./resources/meta")) {
        dir.create("./resources/meta")
    }
    if (!dir.exists("./resources/counts")) {
        dir.create("./resources/counts")
    }
    for (target.sample in seurat@meta.data %>%
        dplyr::distinct(sample) %>%
        dplyr::pull()) {
        message(paste(
            "reading",
            target.sample,
            sep = " "
        ))
        subsetted <-
            seurat %>%
            subset(
                subset = sample == target.sample
            )
        Seurat::GetAssay(assay)@counts %>%
            as.data.frame() %>%
            tibble::rownames_to_column("Gene") %>%
            dplyr::mutate(
                Gene = maptoens(.)
            ) %>%
            filter(
                !is.na(Gene)
            ) %>%
            readr::write_tsv(
                paste(
                    "./resources/counts",
                    paste(
                        target.sample,
                        "tsv",
                        sep = "."
                    ),
                    sep = "/"
                )
            )
        subsetted@meta.data %>%
            tibble::rownames_to_column("Cell") %>%
            dplyr::rename(
                cell_type = celltype
            ) %>%
            dplyr::select(Cell, cell_type) %>%
            write_tsv(
                paste(
                    "./resources/meta",
                    paste(
                        target.sample,
                        "tsv",
                        sep = "."
                    ),
                    sep = "/"
                )
            )
    }
    message(
        paste(
            ">> Please run codes on bash.",
            paste("cd", workdir),
            "cellphonedb method analysis\n
            ./resources/meta/~~sample~~.tsv\n
            ./resources/counts/~~sample~~.tsv\n
            "
        )
    )
}