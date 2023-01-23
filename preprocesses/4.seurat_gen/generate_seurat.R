pacman::p_load(
    Seurat, SeuratDisk, SingleCellPipeline,
    tidyverse
)

sources <-
    list(
        Fetal = "~/resources/Fetal/counted/filtered_counts.tsv",
        hETV2_hCO = "~/resources/hETV2_hCO/counted/filtered_counts.tsv",
        hETV2_vhCO = "~/resources/hETV2_vhCO/counted/filtered_counts.tsv",
        HUVEC_hCO_d65 = "~/resources/HUVEC_hCO_d65/counted/filtered_counts.tsv",
        HUVEC_vhCO_d65 = "~/resources/HUVEC_vhCO_d65/counted/filtered_counts.tsv",
        HUVEC_hCO_d100 = "~/resources/HUVEC_hCO_d100/counted/filtered_counts.tsv",
        HUVEC_vhCO_d100 = "~/resources/HUVEC_vhCO_d100/counted/filtered_counts.tsv"
    )

insert_samplemeta <- function(seurat, samplename) {
    seurat@meta.data$sample <- samplename
    seurat[["percent.mt"]] <-
        seurat %>%
        seurat::PercentageFeatureSet(
            pattern = "^MT-"
        )
    return(seurat)
}

integrateall <- function(objectlist) {
    library(dplyr)
    features <-
        objectlist %>%
        Seurat::SelectIntegrationFeatures(
            nfeatures = 3000
        )
    objectlist %>%
        Seurat::PrepSCTIntegration(
            anchor.features = features
        ) %>%
        Seurat::FindIntegrationAnchors(
            normalization.method = "SCT",
            anchor.features = features
        ) %>%
        Seurat::IntegrateData(normalization.method = "SCT") %>%
        return()
}

names(sources) %>%
    as.list() %>%
    lapply(
        function(target_sample) {
            seurat <-
                sources[[target_sample]] %>%
                read.table(
                    sep = "\t", row.names = 1, header = TRUE
                ) %>%
                CreateSeuratObject(
                    project = "vhCOs",
                    assay = "RNA",
                    names.delim = "_",
                    meta.data = NULL,
                    min.cells = 3,
                    min.features = 200
                ) %>%
                insert_samplemeta(sample) %>%
                SCTransform() %>%
                subset(
                    (percent.mt < 10.0) &
                        (nFeature_RNA > 1000)
                )
            SaveH5Seurat(
                seurat,
                paste0(
                    "~/resources/",
                    target_sample
                )
            )
            return(seurat)
        }
    ) %>%
    integrateall() %>%
    RunPCA() %>%
    RunUMAP(reduction = "pca", dims = 1:20) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 2.0) %>%
    SaveH5Seurat(
        "~/resources/__integrated"
    )
