ggsave_stream <- function(plot, filename, width, height, units = "cm") {
    ggplot2::ggsave(
        filename,
        plot = plot,
        device = NULL,
        path = NULL,
        scale = 1,
        width = width,
        height = height,
        units = units,
        dpi = 300,
        limitsize = TRUE,
        bg = NULL,
    )
}