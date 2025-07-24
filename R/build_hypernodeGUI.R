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
build_hypernodeGUI <- function(hypernode_name,
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
  petri_lib <- fs::path(hyper_root, "petri_net")
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

  model_names <- vapply(
    cfg$cellular_units,
    function(u) as.character(u[["model_name"]]),
    character(1)
  )
  mat_paths     <- vapply(model_names, find_mat, character(1))
  biomass_params<- lapply(cfg$cellular_units, `[[`, "biomass")
  pop_params    <- lapply(cfg$cellular_units, `[[`, "population")
  # — coerce to numeric
  init_counts <- vapply(
    cfg$cellular_units,
    function(u) as.numeric(u[["initial_count"]]),
    numeric(1)
  )


  # 4) Build & process models -----------------------------------------
  message("▶ building biounit models …")
  biounit_models <- epimodFBAfunctions::make_biounit_models(mat_paths, biomass_params, pop_params, init_counts)
  epimodFBAfunctions::write_population_params(biounit_models, fs::path(paths$config, "population_parameters.csv"))
  
  message("▶ copying prebuilt biounits …")
  fs::dir_create(paths$biounit, recurse = TRUE)  # now exists
  prebuilt_dirs <- fs::dir_ls(mat_dir, type = "directory", recurse = FALSE)
  for (src in prebuilt_dirs) {
    model_name <- fs::path_file(src)
    dest       <- fs::path(paths$biounit, model_name)
    message("   • ", model_name)
    fs::dir_copy(src, dest, overwrite = TRUE)
  }

  
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

	# 6) Exchange bounds  ----------------------------------------------------------
	epimodFBAfunctions::run_full_ex_bounds(
		hypernode_name             = hypernode_name,
		biounit_models             = biounit_models,
		projected_lower_bound      = cfg$projected_lower_bound,
		projected_upper_bound      = cfg$projected_upper_bound,
		not_projected_lower_bound  = cfg$not_projected_lower_bound,
		not_projected_upper_bound  = cfg$not_projected_upper_bound,
		volume                     = cfg$volume,
		base_dir                   = base_dir
	)


	#────────────────────────────────────────────────────────────────────────────
	#  DIET-ADJUSTMENT  — aggiorna i bound con i valori di cfg$exchange_bounds
	#────────────────────────────────────────────────────────────────────────────
	if (!is.null(cfg$exchange_bounds)) {

		# 1) path ai due nuovi CSV
		proj_csv   <- fs::path(hyper_root, "output", "ub_bounds_projected.csv")
		nproj_csv  <- fs::path(hyper_root, "output", "non_projected_bounds.csv")

		# 2) template vuoto
		empty_df <- tibble::tibble(
		  reaction        = character(),
		  FBAmodel        = character(),
		  background_conc = double(),
		  system_volume   = double()
		)

		# 3) leggi (o inizializza) i CSV
		proj_df <- if (fs::file_exists(proj_csv))   readr::read_csv(
		              proj_csv, show_col_types = FALSE,
		              col_types = readr::cols(
		                reaction        = readr::col_character(),
		                FBAmodel        = readr::col_character(),
		                background_conc = readr::col_double(),
		                system_volume   = readr::col_double()
		              )
		           ) else empty_df

		nproj_df <- if (fs::file_exists(nproj_csv)) readr::read_csv(
		              nproj_csv, show_col_types = FALSE,
		              col_types = readr::cols(
		                reaction        = readr::col_character(),
		                FBAmodel        = readr::col_character(),
		                background_conc = readr::col_double(),
		                system_volume   = readr::col_double()
		              )
		            ) else empty_df

		# 4) Se entrambi vuoti non c’è nulla da fare
		if (nrow(proj_df) == 0 && nrow(nproj_df) == 0) return(invisible())

		# 5) initial_count per ripartire i consumi tra i modelli
		init_counts <- vapply(biounit_models, `[[`, numeric(1), "initial_count")
		names(init_counts) <- vapply(biounit_models, function(x) model_dir(x$FBAmodel), character(1))

		# 6) loop sugli exchange_bounds del JSON
		exch <- cfg$exchange_bounds
		if (!all(c("reaction", "value") %in% names(exch)))
		  stop("exchange_bounds must contain fields 'reaction' and 'value'.")

		for (i in seq_len(nrow(exch))) {

		  base_rxn <- exch$reaction[i]
		  val      <- as.numeric(exch$value[i])

		  # segno del valore → direzione
		  if (val < 0) {
		    target_rxn <- paste0(base_rxn, "_r")   # reverse
		  } else if (val > 0) {
		    target_rxn <- paste0(base_rxn, "_f")   # forward
		  } else {
		    next                                   # 0 => nessuna modifica
		  }

		  ## ---------- PROJECTED -------------------------------------------------
		  proj_df$background_conc[proj_df$reaction == target_rxn] <- abs(val)

		  ## ---------- NON-PROJECTED --------------------------------------------
		  rows <- which(nproj_df$reaction == target_rxn)
		  if (length(rows)) {
		    # modelli che partecipano a questa reazione
		    orgs <- unique(nproj_df$FBAmodel[rows])
		    tot  <- sum(init_counts[orgs])
		    # quota per organismo ∝ initial_count
		    for (org in orgs) {
		      share <- abs(val) * (init_counts[org] / tot)
		      nproj_df$background_conc[rows & nproj_df$FBAmodel == org] <- share
		    }
		  }
		}

		# 7) scrivi i CSV aggiornati
		readr::write_csv(proj_df,  proj_csv)
		readr::write_csv(nproj_df, nproj_csv)

		cat(
		  "✔ Diet-adjusted bounds written to:\n",
		  "  - ", proj_csv,  "\n",
		  "  - ", nproj_csv, "\n", sep = ""
		)
	}

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
  
  invisible(paths)
}

