#' Build and analyse a hyper‑node automatically
#'
#' A single entry‑point that recreates the full *minimal‑doublet* workflow
#' for any community. The user only needs to supply:
#'
#' * a **PNPRO template** (e.g. `blank.PNPRO`)
#' * a **YAML configuration** in the format already used by the demos
#' * optional **custom *.mat** files placed in a folder you point to
#'
#' The function looks for `.mat` models in four places **in order**:
#' 1. If `model_name` in the YAML **already ends in `.mat`** and exists → use as‑is.
#' 2. `mat_dir` (default "models/") – user‑supplied models for this project.
#' 3. `system.file("MATmodels", …, package = "epimodFBAfunctions")` –
#'    the built‑ins shipped with the package.
#' 4. Error if not found.
#'
#' All generated artefacts live in `base_dir/hypernodes/<hypernode_name>/`.
#'
#' @param hypernode_name  Character. Logical name for the run (no ext).
#' @param config_yaml     Path to *config_<name>.yaml*.
#' @param pnpro_template  PNPRO skeleton.  File name alone is looked up
#'                        inside *inst/petri_nets_library/* of the package.
#' @param mat_dir         Directory where the user dropped custom *.mat
#'                        models (default "models").
#' @param base_dir        Working directory (default = `getwd()`).
#' @param overwrite       Logical. Remove/replace existing outputs? (FALSE)
#' @param debug           Passed to `epimod::model.analysis()`.
#' @return Invisibly a list of created sub‑paths.
#' @export
build_hypernode <- function(hypernode_name,
                            config_yaml,
                            pnpro_template = "blank.PNPRO",
                            mat_dir        = "models",
                            base_dir       = getwd(),
                            overwrite      = FALSE,
                            debug          = FALSE) {

  # ------------------------------------------------------------
  # 0)   Resolve templates shipped in inst/ ---------------------
  # ------------------------------------------------------------
  if (!fs::is_absolute_path(pnpro_template)) {
    pnpro_template <- system.file("petri_nets_library", pnpro_template,
                                  package = "epimodFBAfunctions")
  }
  if (pnpro_template == "")
    stop("PNPRO template not found: ", pnpro_template)

  tmpl_cpp <- system.file("templates", "general_functions_template.cpp",
                          package = "epimodFBAfunctions")
  tmpl_R   <- system.file("templates", "functions_hypernode_template.R",
                          package = "epimodFBAfunctions")

  # ------------------------------------------------------------
  # 1)   Create directory tree ---------------------------------
  # ------------------------------------------------------------
  hyper_root <- fs::path(base_dir, "hypernodes", hypernode_name)
  if (fs::dir_exists(hyper_root) && !overwrite)
    stop("Folder already exists. Set overwrite = TRUE to rebuild: ", hyper_root)
  fs::dir_create(hyper_root, recurse = TRUE)

  paths <- list(
    config = fs::path(hyper_root, "config"),
    src    = fs::path(hyper_root, "src"),
    output = fs::path(hyper_root, "output"),
    biounit= fs::path(hyper_root, "biounits"),
    gen    = fs::path(hyper_root, "gen")
  )
  purrr::walk(paths, fs::dir_create, recurse = TRUE)

  # copy YAML next to the run for provenance
  fs::file_copy(config_yaml, paths$config, overwrite = TRUE)
  cfg_path <- fs::path(paths$config, fs::path_file(config_yaml))
  cfg      <- yaml::read_yaml(cfg_path)

  # ------------------------------------------------------------
  # 2)   Resolve *.mat files -----------------------------------
  # ------------------------------------------------------------
  find_mat <- function(name) {
    # a) explicit path?
    if (grepl("\\.mat$", name, ignore.case = TRUE) && fs::file_exists(name))
      return(fs::path_abs(name))
    # b) user project folder
    user_path <- fs::path(mat_dir, paste0(name, ".mat"))
    if (fs::file_exists(user_path))
      return(fs::path_abs(user_path))
    # c) package built‑ins
    pkg_path <- system.file("MATmodels", paste0(name, ".mat"),
                            package = "epimodFBAfunctions")
    if (pkg_path != "") return(pkg_path)
    stop("MAT model not found (searched mat_dir & package): ", name)
  }

  model_names <- vapply(cfg$cellular_units, `[[`, character(1), "model_name")
  mat_paths   <- vapply(model_names, find_mat, character(1))

  # collect param vectors
  biomass_params <- lapply(cfg$cellular_units, `[[`, "biomass")
  pop_params     <- lapply(cfg$cellular_units, `[[`, "population")
  init_counts    <- as.numeric(vapply(cfg$cellular_units, `[[`, character(1), "initial_count"))

  # ------------------------------------------------------------
  # 3)   Build biounit models + write population parameters -----
  # ------------------------------------------------------------
  message("▶ building biounit models …")
  biounit_models <- epimodFBAfunctions::make_biounit_models(
    mat_paths, biomass_params, pop_params, init_counts)

  epimodFBAfunctions::write_population_params(
    biounit_models,
    fs::path(paths$config, "population_parameters.csv"))

  # ------------------------------------------------------------
  # 4)   Process models, project boundaries, validate PN --------
  # ------------------------------------------------------------
  message("▶ processing models …")
  purrr::walk(biounit_models, function(m)
    epimodFBAfunctions::process_model(m, hypernode_name))

  message("▶ projecting boundary reactions …")
  epimodFBAfunctions::project_boundary_reactions(
    biounit_models      = biounit_models,
    boundary_metabolites= cfg$boundary_metabolites,
    out_dir             = paths$output,
    hypernode_name      = hypernode_name)

  epimodFBAfunctions::validate_pnpro(
    pnpro2validate    = pnpro_template,
    hypernode_root    = hyper_root,
    biounit_models    = biounit_models,
    boundary_metabolites = cfg$boundary_metabolites,
    out_dir           = paths$output,
    hypernode_name    = hypernode_name)

  epimodFBAfunctions::generate_pnpro(
    arc_df   = readr::read_csv(fs::path(paths$output, "repaired_arcs.csv"), show_col_types = FALSE),
    pnpro_out= fs::path(paths$biounit, "..", "..", "petri_nets_library", paste0(hypernode_name, ".PNPRO")))

  # ------------------------------------------------------------
  # 5)   Exchange bounds (default + diet) ----------------------
  # ------------------------------------------------------------
  message("▶ computing exchange bounds …")
  epimodFBAfunctions::run_full_ex_bounds(
    hypernode_name        = hypernode_name,
    biounit_models        = biounit_models,
    projected_base_ub     = cfg$fba_upper_bound,
    projected_base_lb     = cfg$fba_lower_bound,
    not_projected_base_lb = cfg$background_met * cfg$volume,
    not_projected_base_ub = 1000)

  if (!is.null(cfg$exchange_bounds))
    private_adjust_bounds(cfg, biounit_models, hyper_root, cfg$volume)

  # ------------------------------------------------------------
  # 6)   Generate source from templates ------------------------
  # ------------------------------------------------------------
  epimodFBAfunctions::generate_cpp_from_arcs(
    arcs_csv     = fs::path(paths$output, "repaired_arcs.csv"),
    cpp_template = tmpl_cpp,
    output_cpp   = fs::path(paths$src, paste0("general_functions_", hypernode_name, ".cpp")),
    volume       = cfg$volume,
    cell_density = cfg$cell_density)

  epimodFBAfunctions::generate_R_from_pnpro(
    pnpro_file = fs::path(paths$biounit, "..", "..", "petri_nets_library", paste0(hypernode_name, ".PNPRO")),
    r_template = tmpl_R,
    biounit_models = biounit_models,
    output_r   = fs::path(paths$src, paste0("functions_", hypernode_name, ".R")))

  # ------------------------------------------------------------
  # 7)   epimod generation & analysis --------------------------
  # ------------------------------------------------------------
  fba_files <- vapply(biounit_models, function(m)
    fs::path(hyper_root, "biounits", m$FBAmodel, m$txt_file), character(1))

  epimod::model.generation(
    net_fname        = fs::path(paths$biounit, "..", "..", "petri_nets_library", paste0(hypernode_name, ".PNPRO")),
    transitions_fname= fs::path(paths$src, paste0("general_functions_", hypernode_name, ".cpp")),
    fba_fname        = fba_files)

  fs::dir_create(paths$gen)
  purrr::walk(c(".solver", ".def", ".fbainfo", ".net", ".PlaceTransition"), function(ext) {
    src <- fs::path(base_dir, paste0(hypernode_name, ext))
    if (fs::file_exists(src)) fs::file_move(src, paths$gen)
  })

  epimod::model.analysis(
    solver_fname     = fs::path(paths$gen, paste0(hypernode_name, ".solver")),
    parameters_fname = fs::path(paths$config, "initial_data.csv"),
    functions_fname  = fs::path(paths$src, paste0("functions_", hypernode_name, ".R")),
    debug            = debug,
    i_time = 0, f_time = 10, s_time = 1, atol = 1e-6, rtol = 1e-6,
    fba_fname        = fba_files,
    user_files       = c(
      fs::path(paths$config, "population_parameters.csv"),
      fs::path(paths$gen, paste0(hypernode_name, ".fbainfo")),
      fs::path(paths$output, "ub_bounds_projected.csv"),
      fs::path(paths$output, "ub_bounds_not_projected.csv")))

  invisible(paths)
}

# ---------------------------------------------------------------------
# internal helper: diet‑based bound adjustment -------------------------
# ---------------------------------------------------------------------
private_adjust_bounds <- function(cfg, biounit_models, hypernode_root, volume) {
  not_proj_csv <- fs::path(hypernode_root, "output", "ub_bounds_not_projected.csv")
  proj_csv     <- fs::path(hypernode_root, "output", "ub_bounds_projected.csv")

  not_proj_df <- readr::read_csv(not_proj_csv, col_names = c("reaction", "FBAmodel", "upper_bound"), show_col_types = FALSE)
  proj_df     <- readr::read_csv(proj_csv,     col_names = c("reaction", "FBAmodel", "upper_bound"), show_col_types = FALSE)

  init_counts <- vapply(biounit_models, `[[`, numeric(1), "initial_count")
  names(init_counts) <- vapply(biounit_models, `[[`, character(1), "FBAmodel")

  rb <- cfg$exchange_bounds
  rb$reaction_r      <- paste0(rb$reaction, "_r")
  rb$background_met  <- as.numeric(rb$value)

  updated <- purrr::map_dfr(seq_len(nrow(rb)), function(i) {
    rxn <- rb$reaction_r[i]
    subset_df <- dplyr::filter(not_proj_df, .data$reaction == rxn)
    if (nrow(subset_df) == 0) return(NULL)
    orgs <- unique(subset_df$FBAmodel)
    tot  <- sum(init_counts[orgs])
    tibble::tibble(
      reaction   = rxn,
      FBAmodel   = orgs,
      upper_bound= abs((rb$value[i] * volume) / tot)
    )
  })

  not_proj_df <- dplyr::left_join(not_proj_df, updated,
                                  by = c("reaction", "FBAmodel"),
                                  suffix = c("", "_diet")) %>%
    dplyr::mutate(upper_bound = dplyr::coalesce(.data$upper_bound_diet, .data$upper_bound)) %>%
    dplyr::select(-.data$upper_bound_diet)

  # apply to projected bounds
  proj_df <- dplyr::mutate(proj_df,
    upper_bound = dplyr::if_else(stringr::str_detect(.data$reaction, "_r$"),
                                 rb$background_met[match(stringr::str_remove(.data$reaction, "_r$"), rb$reaction)],
                                 .data$upper_bound))

  readr::write_csv(not_proj_df, not_proj_csv)
  readr::write_csv(proj_df,     proj_csv)
}

