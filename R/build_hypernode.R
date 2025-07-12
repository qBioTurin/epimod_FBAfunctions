#' Build and analyse a hyper‑node automatically
#'
#' One function = the whole workflow.  It reproduces the original
#' *minimal‑doublet* demo but works for any community size and lets the
#' user mix built‑in and custom *.mat* models.
#'
#' Required user files *only*:
#'   • a **YAML** configuration (`config_yaml`)
#'   • a PNPRO **template** (`pnpro_template`)
#'   • any custom *.mat* files (placed in `mat_dir`)
#'
#' Everything else (templates, built‑in MAT models, lookup tables…) is
#' retrieved with `system.file()`.
#'
#' Paths resolved:
#'   * PNPRO template
#'       1. absolute / relative path as given
#'       2. `inst/extdata/<file>` shipped with the package
#'   * MAT models
#'       1. explicit path ending in `.mat`
#'       2. `<mat_dir>/<name>.mat`
#'       3. `inst/MATmodels/<name>.mat`
#'
#' All artefacts are written to:
#'   `<base_dir>/hypernodes/<hypernode_name>/`
#' and the generated PNPRO is saved to
#'   `<base_dir>/petri_nets_library/<hypernode_name>.PNPRO`.
#'
#' @param hypernode_name Character.  Label for the run (no extension).
#' @param config_yaml    Path to the YAML file.
#' @param pnpro_template Path to a PNPRO skeleton.  If just a filename it
#'                       is searched in the package’s `extdata/` folder.
#' @param mat_dir        Folder containing user *.mat* files (default "models").
#' @param base_dir       Root working directory (default = `getwd()`).
#' @param overwrite      Logical.  Delete any previous run with the same
#'                       name before rebuilding (FALSE).
#' @param debug          Passed to `epimod::model.analysis()`.
#' @return (Invisibly) a list of created sub‑paths.
#' @export
build_hypernode <- function(hypernode_name,
                            config_yaml,
                        		boundary_conditions_file,
                        		initial_data,
                            mat_dir        = "models",
                            base_dir       = getwd(),
                            overwrite      = FALSE,
                            debug          = FALSE) {
  # helper -------------------------------------------------------------
  abs_path <- function(x) fs::path_abs(fs::path_expand(x))

  # 0)  Resolve PNPRO template ----------------------------------------
  pnpro_template <- system.file("extdata", "blank.PNPRO",
                                package = "epimodFBAfunctions")
  
  if (pnpro_template == "")
    stop("PNPRO template not found: ", pnpro_template)
    
  pnpro_template <- abs_path(pnpro_template)

  tmpl_cpp <- system.file("templates", "general_functions_template.cpp",
                          package = "epimodFBAfunctions")
  tmpl_R   <- system.file("templates", "functions_hypernode_template.R",
                          package = "epimodFBAfunctions")

  # 1) Directory scaffold ---------------------------------------------
  hyper_root <- fs::path(base_dir, "hypernodes", hypernode_name)
  if (fs::dir_exists(hyper_root) && !overwrite)
    stop("Run already exists, set overwrite = TRUE: ", hyper_root)
	if (fs::dir_exists(hyper_root))
		fs::dir_delete(hyper_root)
  fs::dir_create(hyper_root)

  paths <- list(
    config = fs::path(hyper_root, "config"),
    src    = fs::path(hyper_root, "src"),
    output = fs::path(hyper_root, "output"),
    biounit= fs::path(hyper_root, "biounits"),
    gen    = fs::path(hyper_root, "gen")
  )
  purrr::walk(paths, fs::dir_create, recurse = TRUE)

  # local petri_nets_library (mirrors old script)
  petri_lib <- fs::path(base_dir, "petri_nets_library")
  fs::dir_create(petri_lib)

  # 2) Load YAML (+ optional JSON) ------------------------------------
  fs::file_copy(initial_data, paths$config, overwrite = TRUE)
  fs::file_copy(config_yaml, paths$config, overwrite = TRUE)
  cfg_path <- fs::path(paths$config, fs::path_file(config_yaml))
  cfg      <- yaml::read_yaml(cfg_path)

  # ----- boundary_conditions_file argument --------------------------
  if (missing(boundary_conditions_file) || is.null(boundary_conditions_file))
    stop("`boundary_conditions_file` (JSON) must be supplied explicitly.")

  if (!fs::file_exists(boundary_conditions_file))
    stop("boundary_conditions_file not found: ", boundary_conditions_file)

  # copy next to YAML (re-run reproducibility)
  fs::file_copy(boundary_conditions_file, paths$config, overwrite = TRUE)
  bc_path <- fs::path(paths$config, fs::path_file(boundary_conditions_file))

  # merge JSON content into cfg, JSON values take precedence
  bc_json <- jsonlite::fromJSON(bc_path, simplifyVector = TRUE)
  cfg     <- utils::modifyList(cfg, bc_json)

  # 3) Resolve *.mat ---------------------------------------------------
  find_mat <- function(name) {
    if (grepl("\\.mat$", name, ignore.case = TRUE) && fs::file_exists(name))
      return(abs_path(name))
    user_path <- fs::path(mat_dir, paste0(name, ".mat"))
    if (fs::file_exists(user_path))
      return(abs_path(user_path))
    pkg_path <- system.file("MATmodels", paste0(name, ".mat"),
                            package = "epimodFBAfunctions")
    if (pkg_path != "") return(pkg_path)
    stop("MAT model not found: ", name)
  }

  model_names <- vapply(cfg$cellular_units, `[[`, character(1), "model_name")
  mat_paths   <- vapply(model_names, find_mat, character(1))
  biomass_params <- lapply(cfg$cellular_units, `[[`, "biomass")
  pop_params     <- lapply(cfg$cellular_units, `[[`, "population")
  init_counts    <- as.numeric(vapply(cfg$cellular_units, `[[`, character(1), "initial_count"))

  # 4) Build & process models -----------------------------------------
  message("▶ building biounit models …")
  biounit_models <- epimodFBAfunctions::make_biounit_models(mat_paths, biomass_params, pop_params, init_counts)
  epimodFBAfunctions::write_population_params(biounit_models, fs::path(paths$config, "population_parameters.csv"))
  
	purrr::walk(
		biounit_models,
		epimodFBAfunctions::process_model,
		hypernode_name = hypernode_name,
		mat_dir        = mat_dir,   
		base_dir       = base_dir
	)

  # 5) Boundary projection + PN repair --------------------------------
	epimodFBAfunctions::project_boundary_reactions(
		biounit_models       = biounit_models,
		boundary_metabolites = cfg$boundary_metabolites,
		out_dir              = paths$output,
		hypernode_name       = hypernode_name,
		base_dir             = base_dir   
	)

  epimodFBAfunctions::validate_pnpro(pnpro_template, hyper_root, biounit_models, cfg$boundary_metabolites, paths$output, hypernode_name)

  repaired_pn <- fs::path(petri_lib, paste0(hypernode_name, ".PNPRO"))
  epimodFBAfunctions::generate_pnpro(readr::read_csv(fs::path(paths$output, "repaired_arcs.csv"), show_col_types = FALSE), repaired_pn)

  # 6) Exchange bounds -------------------------------------------------
  epimodFBAfunctions::run_full_ex_bounds(hypernode_name, biounit_models, cfg$fba_upper_bound, cfg$fba_lower_bound, cfg$background_met * cfg$volume, 1000)
  if (!is.null(cfg$exchange_bounds))
    private_adjust_bounds(cfg, biounit_models, hyper_root, cfg$volume)

	# 7) Emit C++ / R helpers -------------------------------------------

	## ── resolve template files **inside the installed package** ─────────
	##    (mirrors how we located blank.PNPRO above)
	tmpl_cpp <- system.file("templates", "general_functions_template.cpp",
		                      package = "epimodFBAfunctions")
	tmpl_R   <- system.file("templates", "functions_hypernode_template.R",
		                      package = "epimodFBAfunctions")

	if (tmpl_cpp == "" || tmpl_R == "")
		stop("Template files not found inside package ‘epimodFBAfunctions’")

	## ── generate helper files ───────────────────────────────────────────
	epimodFBAfunctions::generate_cpp_from_arcs(
		fs::path(paths$output, "repaired_arcs.csv"),
		tmpl_cpp,
		fs::path(paths$src, paste0("general_functions_", hypernode_name, ".cpp")),
		cfg$volume,
		cfg$cell_density
	)

	epimodFBAfunctions::generate_R_from_pnpro(
		repaired_pn,
		tmpl_R,
		biounit_models,
		fs::path(paths$src, paste0("functions_", hypernode_name, ".R"))
	)


  # 8) epimod model generation + analysis -----------------------------
  fba_files <- vapply(biounit_models, function(m)
    fs::path(hyper_root, "biounits", m$FBAmodel, m$txt_file), character(1))

  epimod::model.generation(net_fname = repaired_pn,
                           transitions_fname = fs::path(paths$src, paste0("general_functions_", hypernode_name, ".cpp")),
                           fba_fname = fba_files)

  fs::dir_create(paths$gen)
  purrr::walk(c(".solver", ".def", ".fbainfo", ".net", ".PlaceTransition"), function(ext) {
    src <- fs::path(base_dir, paste0(hypernode_name, ext))
    if (fs::file_exists(src)) fs::file_move(src, paths$gen)
  })
  

  epimod::model.analysis(solver_fname     = fs::path(paths$gen, paste0(hypernode_name, ".solver")),
                         parameters_fname = fs::path(paths$config, fs::path_file("initial_data.csv")),
                         functions_fname  = fs::path(paths$src, paste0("functions_", hypernode_name, ".R")),
                         debug = debug, i_time = 0, f_time = 10, s_time = 1, atol = 1e-6, rtol = 1e-6,
                         fba_fname = fba_files,
                         user_files = c(fs::path(paths$config, "population_parameters.csv"),
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

