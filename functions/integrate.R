# Integrate multiple sample data into seuratobject.
# [WARNING!] This process requires plenty of time and memory.
# @
# > project <- "vascular"
# > basepath <-
# >      "/home/yuyasato/work2/_organoid_scRNA/raw/matrix/__filtered/__fetal"
# > samples <- list(
# >     fetalPCW14 = list(
# >         path = paste(basepath, "fetalPCW14_0rsydy7.seurat", sep = "/")
# >     ),
# >     fetalPCW15 = list(
# >         path = paste(basepath, "fetalPCW15_0rsydy7.seurat", sep = "/")
# >     ),
# >     fetalPCW16 = list(
# >         path = paste(basepath, "fetalPCW16_GSE162170.seurat", sep = "/")
# >     ),
# >     fetalPCW19 = list(
# >         path = paste(basepath, "fetalPCW19_0rsydy7.seurat", sep = "/")
# >     ),
# >     fetalPCW20 = list(
# >         path = paste(basepath, "fetalPCW20_GSE162170.seurat", sep = "/")
# >     ),
# >     fetalPCW21 = list(
# >         path = paste(basepath, "fetalPCW21_GSE162170.seurat", sep = "/")
# >     ),
# >     fetalPCW24 = list(
# >         path = paste(basepath, "fetalPCW24_GSE162170.seurat", sep = "/")
# >     )
# > )
# > resolution <- 2.0
# > markers <- "/home/yuyasato/work2/_organoid_scRNA/ref/__markers_refined.seurat"


integrate <-
    function(samples, project_name, outputdir, marker_path, max_dim = 20, resolution = 1.5, force_run = FALSE) { # nolint
        library(Seurat)
        library(SeuratDisk)
        library(SeuratData)
        library(tidyverse)
        as_fullpath <- function(path) {
            return(paste(outputdir, project_name, path, sep = "/"))
        }
        log_output <- function(content) {
            if (!dir.exists(as_fullpath("logs"))) {
                dir.create(as_fullpath("logs"))
            }
            if (!file.exists(as_fullpath("logs/__Integration.log"))) {
                cat(
                    paste(
                        paste(
                            "[",
                            Sys.time(),
                            "]",
                            sep = ""
                        ),
                        "Integration process started"
                    ),
                    sep = "\n",
                    file = paste(
                        as_fullpath("logs/__Integration.log")
                    ),
                    append = FALSE
                )
            }
            cat(
                paste(
                    paste(
                        "[",
                        Sys.time(),
                        "]",
                        sep = ""
                    ),
                    content
                ),
                sep = "\n",
                file = paste(
                    as_fullpath("logs/__Integration.log")
                ),
                append = TRUE
            )
        }

        if (!file.exists(marker_path)) {
            stop("Cound not find marker path. Aborting.")
        }
        output <- list(
            seurat = as_fullpath("seurat"),
            markers = as_fullpath("dist/non_biased_markers.tsv"),
            cluster_plt = as_fullpath("dist/cluster.png"),
            cluster_each_plt = as_fullpath("dist/cluster_each.png"),
            merged_plt = as_fullpath("dist/merged.png"),
            feature_plt = as_fullpath("dist/feature.png"),
            feature_each_plt = as_fullpath("dist/cluster_feature.png"),
            feature_dot_plt = as_fullpath("dist/dotplot.png")
        )
        # Prepare
        seurats <- list()
        gc()
        distdir <- paste(outputdir, project_name, sep = "/")
        if (dir.exists(distdir)) {
            if (force_run) {
                unlink(distdir, recursive = TRUE)
            } else {
                message_text <- paste(
                    "The disignated directory",
                    paste(outputdir, project_name, sep = "/"),
                    "is already exists. Aborting."
                )
                log_output(paste(
                    "[ERROR]:",
                    message_text
                ))
                stop(message_text)
            }
        }
        dir.create(paste(outputdir, project_name, sep = "/"))
        dir.create(as_fullpath("dist"))
        log_output("@ Each process start")
        for (sample in names(samples)) {
            if (!file.exists(samples[[sample]])) {
                msg <-
                    paste(
                        "[ERROR] Given path:",
                        samples[[sample]],
                        "does not exists."
                    )
                log_output(msg)
                stop(msg)
            }
            if (!(stringr::str_ends(samples[[sample]], ".seurat"))) {
                msg <- "[ERROR] Data should be seurat file."
                log_output(msg)
                stop(msg)
            }
        }
        for (sample in names(samples)) {
            log_output(paste(
                ">start loading ",
                sample
            ))
            if (stringr::str_ends(samples[[sample]], ".seurat")) {
                seurats[[sample]] <-
                    SeuratDisk::LoadH5Seurat(
                        samples[[sample]]
                    )
                log_output(">> Finished loading Seurat")
            }
            log_output(paste(
                ">finished loading ",
                sample
            ))
        }
        log_output("@ Each process end")

        log_output("Selecting Integration Features...")
        features <- SelectIntegrationFeatures(
            object.list = seurats,
            verbose = TRUE
        )
        log_output("Prepare Integration...")
        seurats <-
            PrepSCTIntegration(
                object.list = seurats,
                anchor.features = features,
                verbose = TRUE
            )

        log_output("Finding Integration anchors...")
        anchors <- Seurat::FindIntegrationAnchors(
            object.list = seurats,
            anchor.features = features,
            normalization.method = "SCT",
            verbose = TRUE
        )
        rm(seurats)
        rm(features)
        log_output("Integrating... new.assay.name=integrated")
        merged <- Seurat::IntegrateData(
            anchorset = anchors,
            normalization.method = "SCT",
            new.assay.name = "integrated",
            verbose = TRUE
        )
        rm(anchors)

        SeuratObject::DefaultAssay(merged) <- "integrated"
        log_output(
            paste(
                "Running PCA: Calc MAX Dims:",
                max_dim
            )
        )
        merged <- Seurat::RunPCA(
            merged,
            npcs = max_dim,
            verbose = TRUE
        )
        log_output(
            paste(
                "Running UMAP: Calc Compressing Dims:",
                max_dim
            )
        )
        merged <- Seurat::RunUMAP(
            merged,
            reduction = "pca",
            dims = 1:max_dim,
            verbose = TRUE
        )
        log_output(paste0(
            "Finding k-neighbors@dim=1:", resolution
        ))
        merged <- Seurat::FindNeighbors(
            merged,
            verbose = TRUE
        )
        log_output(paste(
            "Finding clusters. Resolution:",
            resolution
        ))
        merged <- Seurat::FindClusters(
            merged,
            resolution = resolution,
            verbose = TRUE
        )
        log_output("< Finished.")
        log_output("Starting Dim plot")


        # Clustering plot
        cluster_plt <- Seurat::DimPlot(
            merged,
            reduction = "umap",
            label = TRUE,
            repel = TRUE,
            raster = FALSE
        ) + NoLegend()
        ggplot2::ggsave(output$cluster_plt, cluster_plt,
            width = 20, height = 20, units = "cm"
        )
        rm(cluster_plt)
        ###################
        # Merged plot
        merged_plt <-
            Seurat::DimPlot(
                merged,
                reduction = "umap",
                group.by = "sample",
                raster = FALSE
            ) +
            ggplot2::ggtitle("") +
            ggplot2::theme(
                legend.justification = "bottom"
            )
        ggplot2::ggsave(output$merged_plt, merged_plt,
            width = 22, height = 20, units = "cm"
        )
        rm(merged_plt)


        ggplot2::ggsave(
            output$cluster_each_plt,
            Seurat::DimPlot(
                merged,
                reduction = "umap", split.by = "sample"
            ),
            width = 100, height = 20, units = "cm"
        )
        log_output("Starting Feature plot")

        ggplot2::ggsave(
            output$feature_plt,
            VlnPlot(
                merged,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, group.by = "sample",
                pt.size = 0
            ),
            width = 40, height = 40, units = "cm"
        )

        ggplot2::ggsave(
            output$feature_each_plt,
            VlnPlot(
                merged,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, group.by = "seurat_clusters",
                pt.size = 0
            ),
            width = 80, height = 40, units = "cm"
        )
        log_output("Starting Dotplot")
        marker <-
            readr::read_tsv(marker_path) %>%
            dplyr::select(gene) %>%
            dplyr::pull()
        color <- "white"
        dp <- Seurat::DotPlot(
            merged,
            features = marker,
            scale.min = 0,
            cols = c("#ffffac", "#0034ad"),
            col.min = 0,
            assay = "SCT"
        ) +
            theme(
                rect = element_rect(colour = color, fill = color),
                text = element_text(color = "#000000"),
                plot.background = element_rect(colour = color, fill = color),
                panel.background = element_rect(colour = color, fill = color),
                panel.border =
                    element_rect(fill = NA, colour = "#1b1b1b", size = 1),
                panel.grid.major = element_line(colour = "grey60"),
                panel.grid.minor = element_line(colour = "grey30"),
                legend.key = element_rect(colour = color, fill = color),
                axis.text = element_text(colour = "#000000"),
                axis.text.x.bottom = element_text(
                    size = 18,
                    hjust = .5, vjust = .5
                ),
                axis.text.y.left = element_text(
                    size = 18
                ),
                axis.title.y.left = element_text(
                    size = 20,
                    margin = unit(c(0, 1, 0, 1), "lines")
                ),
                axis.title.x.bottom = element_text(
                    size = 20,
                    margin = unit(c(1, 0, 1, 0), "lines")
                )
            ) +
            xlab("Gene Features") +
            ylab("Cluster") + coord_flip()

        ggplot2::ggsave(
            output$feature_dot_plt,
            dp,
            width = 90, height = 50, units = "cm", dpi = 200
        )

        log_output("Saving Seurat object")
        SeuratDisk::SaveH5Seurat(
            object = merged,
            filename = output$seurat,
            overwrite = TRUE
        )

        log_output("Finding All markers")
        Seurat::FindAllMarkers(
            merged,
            verbose = TRUE
        ) %>%
            readr::write_tsv(
                file = output$markers
            )
        log_output("!!Completed All Processes")
    }
