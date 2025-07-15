# R/plots.R
#' Plot Species Abundance and Biomass Dynamics for a Hypernode
#'
#' @param hypernode_dir Character. Path to the top-level hypernode directory.
#' @param final_time    Numeric. Maximum time to plot (default: Inf).
#' @param output_dir    Character. Where to save plots (if NULL, does not save).
#' @return A ggplot object of species abundance over time.
#' @import dplyr tidyr ggplot2 yaml RColorBrewer
#' @export
plot_species_dynamics <- function(hypernode_dir, final_time = Inf, output_dir = NULL) {
  # 1) Read config YAML (not used here but could supply labels/colors)
  cfg_file <- list.files(file.path(hypernode_dir, "config"),
                         pattern = "\\.ya?ml$", full.names = TRUE)[1]
  cfg      <- yaml::read_yaml(cfg_file)

  # 2) Locate trace file
  analysis_dir <- file.path(hypernode_dir,
                            paste0(basename(hypernode_dir), "_analysis"))
  trace_dir  <- if (dir.exists(analysis_dir)) analysis_dir else hypernode_dir
  trace_file <- list.files(trace_dir,
                           pattern = "\\.trace$",
                           recursive = TRUE,
                           full.names = TRUE)[1]

  # 3) Read and melt
  tr <- utils::read.table(trace_file, header = TRUE) %>%
    tidyr::gather(key = "Places", value = "Marking", -Time)

  # 4) Filter only abundance columns (n_<species>)
  df <- tr %>%
    dplyr::filter(grepl("^n_", Places), Time <= final_time) %>%
    dplyr::mutate(
      Species = sub("^n_", "", Places)
    )

  # 5) Plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Time, y = Marking,
                                        color = Species, group = Species)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(
      title = "Species Abundance Dynamics",
      x     = "Time",
      y     = "Cell Count",
      color = "Species"
    ) +
    ggplot2::theme_minimal(base_size = 14)

  # 6) Save if requested
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(file.path(output_dir, "species_dynamics.png"),
                    plot = p, width = 7, height = 4)
    ggplot2::ggsave(file.path(output_dir, "species_dynamics.pdf"),
                    plot = p, width = 7, height = 4)
  }

  return(p)
}



#' Plot Extracellular Metabolite Concentrations Over Time
#'
#' @param hypernode_dir Character. Path to the hypernode directory.
#' @param output_dir    Character. Where to save plots (if NULL, does not save).
#' @return A ggplot object of metabolite concentration dynamics.
#' @import dplyr tidyr ggplot2 yaml RColorBrewer
#' @export
plot_metabolite_dynamics <- function(hypernode_dir, output_dir = NULL) {
  # 1) Read config YAML
  cfg_file <- list.files(file.path(hypernode_dir, "config"), pattern = "\\.ya?ml$", full.names = TRUE)[1]
  cfg      <- yaml::read_yaml(cfg_file)
  mets     <- cfg$boundary_metabolites

  # Colors
  cols <- RColorBrewer::brewer.pal(max(3, length(mets)), "Dark2")[1:length(mets)]
  names(cols) <- mets

  # 2) Locate trace
  analysis_dir <- file.path(hypernode_dir, paste0(basename(hypernode_dir), "_analysis"))
  trace_dir    <- if (dir.exists(analysis_dir)) analysis_dir else hypernode_dir
  trace_file   <- list.files(trace_dir, pattern = "\\.trace$", recursive = TRUE, full.names = TRUE)[1]

  # 3) Read and tidy
  tr <- utils::read.table(trace_file, header = TRUE) %>%
    tidyr::gather(key = "Places", value = "Marking", -Time)

  # 4) Filter metabolites
  df <- tr %>%
    dplyr::filter(Places %in% mets) %>%
    dplyr::rename(Concentration = Marking)

  # 5) Plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Time, y = Concentration, color = Places, group = Places)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(
      title = "Extracellular Metabolite Dynamics",
      x = "Time", y = "Concentration", color = "Metabolite"
    ) +
    ggplot2::theme_minimal(base_size = 14)

  # 6) Save
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(file.path(output_dir, "metabolite_dynamics.png"), p, width = 7, height = 4)
    ggplot2::ggsave(file.path(output_dir, "metabolite_dynamics.pdf"), p, width = 7, height = 4)
  }

  return(p)
}



#' Plot Species Biomass Dynamics for a Hypernode
#'
#' @param hypernode_dir Character. Path to the top-level hypernode directory.
#' @param final_time    Numeric. Maximum time to plot (default: Inf).
#' @param output_dir    Character. Where to save plots (if NULL, does not save).
#' @return A ggplot object of species biomass over time.
#' @import dplyr tidyr ggplot2 yaml
#' @export
plot_biomass_dynamics <- function(hypernode_dir, final_time = Inf, output_dir = NULL) {
  # 1) Locate the .trace file (in _analysis if present)
  analysis_dir <- file.path(hypernode_dir, paste0(basename(hypernode_dir), "_analysis"))
  trace_dir    <- if (dir.exists(analysis_dir)) analysis_dir else hypernode_dir
  trace_file   <- list.files(trace_dir, pattern="\\.trace$", recursive=TRUE, full.names=TRUE)[1]

  # 2) Read & pivot
  tr <- utils::read.table(trace_file, header = TRUE) %>%
    tidyr::gather(key="Places", value="Marking", -Time)

  # 3) Keep only biomass columns (those beginning with "biomass_e_")
  df <- tr %>%
    dplyr::filter(grepl("^biomass_e_", Places), Time <= final_time) %>%
    dplyr::mutate(
      Species = sub("^biomass_e_", "", Places)
    )

  # 4) Plot (same colour mapping by Species)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Time, y = Marking,
                                        color = Species, group = Species)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(
      title = "Time-dependent Average Biomass",
      x     = "Time (h)",
      y     = "Biomass (pg)",
      color = "Species"
    ) +
    ggplot2::theme_minimal(base_size = 14)

  # 5) Save if requested
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(file.path(output_dir, "biomass_dynamics.png"), plot = p,
                    width = 7, height = 4)
    ggplot2::ggsave(file.path(output_dir, "biomass_dynamics.pdf"), plot = p,
                    width = 7, height = 4)
  }

  return(p)
}


