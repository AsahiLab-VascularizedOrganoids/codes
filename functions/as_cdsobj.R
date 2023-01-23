as_cdsobj <- function(seurat) {
    pacman::p_load(
        monocle3, Seurat, tidyverse, # nolint
        SeuratDisk, scran, igraph, # nolint
        SeuratWrappers # nolint
    )
}
