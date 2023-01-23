plot_umap <- function(seurat, group.by, split.by = NULL, disable.color = "gray", arrow_color = "#060056", point.size = .2, out_lambda = NULL, trajectory_file_path = NULL, raster = TRUE) { # nolint
  pacman::p_load(Seurat, SeuratDisk, tidyverse)
  coordinates <-
    Embeddings(seurat, "umap") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    left_join(
      seurat@meta.data %>%
        dplyr::select(
          !!sym(group.by),
          !!sym(paste0("color.", group.by))
        ) %>%
        tibble::rownames_to_column("cell"),
      by = "cell"
    ) %>%
    mutate(
      label = !!sym(group.by)
    )
  colors <-
    coordinates %>%
    arrange(!!sym(group.by)) %>%
    dplyr::distinct(!!sym(paste0("color.", group.by))) %>%
    pull() %>%
    as.character()
  if (!is.null(split.by)) {
    targets <-
      seurat@meta.data %>%
      distinct(!!sym(split.by)) %>%
      pull()
    levels <-
      targets %>%
      levels()
    if (is.null(levels)) {
      levels <- targets
    }
    coordinates <-
      coordinates %>%
      dplyr::left_join(
        seurat@meta.data %>%
          dplyr::select(
            !!sym(split.by)
          ) %>%
          tibble::rownames_to_column("cell"),
        by = "cell"
      )
    coordinates <-
      coordinates %>%
      dplyr::distinct(!!sym(split.by)) %>%
      pull() %>%
      as.list() %>%
      lapply(
        function(target) {
          coordinates %>%
            mutate(
              split_key = factor(
                target,
                levels = levels
              ),
              label = if_else(
                !!sym(split.by) == target,
                true = as.character(!!sym(group.by)),
                false = "disabled"
              ) %>%
                factor(
                  levels = c(
                    "disabled",
                    coordinates %>%
                      distinct(!!sym(group.by)) %>%
                      pull() %>%
                      levels()
                  )
                ),
              keyarrange = !!sym(split.by) %>%
                factor(
                  levels = c(
                    target,
                    coordinates %>%
                      filter(
                        !!sym(split.by) != target
                      ) %>%
                      distinct(!!sym(split.by)) %>%
                      pull()
                  )
                )
            )
        }
      ) %>%
      bind_rows()
    colors <- c(
      disable.color,
      colors
    )
  }
  p <-
    coordinates %>%
    arrange(label) %>%
    ggplot(
      aes(
        x = UMAP_1,
        y = UMAP_2,
        color = label
      )
    ) +
    geom_point(size = point.size) +
    scale_color_manual(
      values = colors
    ) +
    theme_minimal() +
    theme(
      legend.justification = "bottom",
      plot.background = element_rect(
        fill = "white",
        colour = "white"
      )
    ) +
    guides(
      colour = guide_legend(
        title = group.by,
        override.aes = list(size = 5)
      )
    )
  if (!is.null(split.by)) {
    p <-
      p +
      facet_wrap(. ~ split_key)
  }
  if (!is.null(trajectory_file_path)) {
    lineage_coordinates <-
      Embeddings(seurat@reductions$umap) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell") %>%
      left_join(
        seurat@meta.data %>%
          tibble::rownames_to_column("cell") %>%
          select(cell, seurat_clusters),
        by = "cell"
      ) %>%
      tibble::column_to_rownames("cell") %>%
      group_by(seurat_clusters) %>%
      dplyr::summarise(
        UMAP_x = mean(UMAP_1),
        UMAP_y = mean(UMAP_2)
      ) %>%
      dplyr::mutate(
        seurat_clusters = as.character(seurat_clusters)
      )
    lineages <-
      readRDS(trajectory_file_path)$lineages %>%
      lapply(
        function(lineage) {
          lineage_df <- tibble()
          cnt <- 1
          for (lin in as.integer(lineage)) {
            if (cnt < length(lineage)) {
              lineage_df <-
                lineage_df %>%
                rbind(
                  list(
                    "from" = lineage[cnt],
                    "to" = lineage[cnt + 1]
                  ) %>%
                    as_tibble()
                )
            }
            cnt <- cnt + 1
          }
          return(lineage_df)
        }
      ) %>%
      bind_rows() %>%
      distinct(from, to) %>%
      left_join(
        lineage_coordinates %>%
          dplyr::rename(
            UMAP_from_x = UMAP_x,
            UMAP_from_y = UMAP_y
          ),
        by = c("from" = "seurat_clusters")
      ) %>%
      left_join(
        lineage_coordinates %>%
          dplyr::rename(
            UMAP_to_x = UMAP_x,
            UMAP_to_y = UMAP_y
          ),
        by = c("to" = "seurat_clusters")
      )
    arrow_style <- ggplot2::arrow(
      angle = 20,
      ends = "last",
      type = "closed",
      length = grid::unit(5, "pt")
    )
    p <-
      p +
      ggplot2::geom_segment(
        data = lineages,
        aes(
          x = UMAP_from_x, y = UMAP_from_y,
          xend = UMAP_to_x, yend = UMAP_to_y
        ),
        colour = "#161616",
        arrow = arrow_style,
        size = 0.8
      )
  }

  if (is.null(out_lambda)) {
    return(p)
  }
  out_lambda(p)
  return(seurat)
}
