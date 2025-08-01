#' General wrapper for epimod::model.analysis (GUI) with console debug
#'
#' @export
model_analysis_GUI <- function(
  paths,
  hypernode_name,
  debug_solver = FALSE,
  i_time       = 0,
  f_time       = 10,
  s_time       = 1,
  atol         = 1e-6,
  rtol         = 1e-6,
  fba_fname,
  user_files   = character(),
  volume       = getwd()
) {
  message("[DEBUG] Starting model_analysis_GUI for hypernode: ", hypernode_name)

  # ───────────────────────────────────
  # 0) Attach & validate
  # ───────────────────────────────────
  if (!"package:epimod" %in% search()) {
    message("[DEBUG] Loading epimod package")
    library(epimod)
  }
  message("[DEBUG] Validating inputs")
  stopifnot(
    is.character(paths), length(paths) >= 4,
    all(c("gen","config","src","output") %in% names(paths)),
    is.character(hypernode_name), nzchar(hypernode_name),
    is.logical(debug_solver), length(debug_solver) == 1,
    is.numeric(i_time), is.numeric(f_time), is.numeric(s_time),
    is.numeric(atol), is.numeric(rtol),
    is.character(fba_fname), length(fba_fname) >= 1,
    all(file.exists(fba_fname)),
    is.character(user_files)
  )

  # ───────────────────────────────────
  # 1) Build core paths
  # ───────────────────────────────────
  solver_fname     <- fs::path(paths["gen"], paste0(hypernode_name, ".solver"))
  parameters_fname <- fs::path(paths["config"], "initial_data.csv")
  orig_fun_fname   <- fs::path(paths["src"],   paste0("functions_", hypernode_name, ".R"))
  message("[DEBUG] solver file: ", solver_fname)
  message("[DEBUG] parameters file: ", parameters_fname)
  message("[DEBUG] original functions.R: ", orig_fun_fname)

  # ───────────────────────────────────
  # 2) Ensure files exist
  # ───────────────────────────────────
  needed <- c(solver_fname, parameters_fname, orig_fun_fname)
  missing <- needed[!file.exists(needed)]
  if (length(missing)) {
    stop("Missing required analysis files: ", paste(basename(missing), collapse = ", "))
  }
  message("[DEBUG] All core files exist")

  # ───────────────────────────────────
  # 3) Read GUI snapshot YAML
  # ───────────────────────────────────
  gui_yaml_path <- fs::path(paths["config"], paste0(hypernode_name, "_gui.yaml"))
  message("[DEBUG] Looking for GUI YAML at: ", gui_yaml_path)
  if (!file.exists(gui_yaml_path)) {
    stop("Missing GUI‐snapshot YAML: ", gui_yaml_path)
  }
  gui_yml <- yaml::read_yaml(gui_yaml_path)
  message("[DEBUG] Read GUI YAML successfully")

  # ───────────────────────────────────
  # 3a) mu_max CSV  (high precision)
  # ───────────────────────────────────
  mu_defs <- gui_yml$cellular_units %||% list()
  if (length(mu_defs) > 0) {
    mu_df <- data.frame(
      Model  = vapply(mu_defs, function(u) u$model_name, character(1)),
      mu_max = vapply(mu_defs, function(u) u$mu_max, numeric(1)),
      stringsAsFactors = FALSE
    )
    mu_csv <- fs::path(paths["config"], "mu_max_values_gui.csv")
    message("[DEBUG] Writing mu_max CSV: ", mu_csv)
		write.csv(mu_df, mu_csv, row.names = FALSE, quote = FALSE)   # ← tolto digits

    user_files <- c(user_files, mu_csv)
  }

  # ───────────────────────────────────
  # 3b) population_parameters.csv
  # ───────────────────────────────────

	pop_csv <- fs::path(paths["config"], "population_parameters.csv")

	if (!file.exists(pop_csv)) {             # ← aggiungi questo if
		pop_mat <- do.call(rbind, lapply(pop_defs, function(u)
		  c(u$population$starv, u$population$dup, u$population$death)))

		message("[DEBUG] Writing population parameters CSV: ", pop_csv)
		write.table(pop_mat, pop_csv, sep = ",", row.names = FALSE,
		            col.names = FALSE, quote = FALSE)
		user_files <- c(user_files, pop_csv)
	} else {
		message("[DEBUG] Keeping existing population_parameters.csv")
	}



  # ───────────────────────────────────
  # 5) Prepare GUI-patched R stub
  # ───────────────────────────────────
  gui_fun_fname <- fs::path(paths["src"], paste0("functions_", hypernode_name, "_gui.R"))
  file.copy(orig_fun_fname, gui_fun_fname, overwrite = TRUE)
  fun_lines <- readLines(gui_fun_fname)

  # locate y_ini / yini.names
  ini_idx   <- grep("^\\s*y_ini\\s*<-\\s*c\\(", fun_lines)
  names_idx <- grep("^\\s*yini\\.names\\s*<-\\s*c\\(", fun_lines)
  if (!length(ini_idx) || !length(names_idx))
    stop("Cannot find y_ini or yini.names definitions in stub")

  raw_names <- strsplit(sub("^.*c\\((.*)\\).*", "\\1", fun_lines[names_idx]), ",")[[1]]
  name_vec  <- gsub("['\"]", "", trimws(raw_names))
  val_vec   <- trimws(strsplit(sub("^.*c\\((.*)\\).*", "\\1", fun_lines[ini_idx[1]]), ",")[[1]])

  # ------------------------------------------------------------------
  #  Fallback: ricostruisco cellular_units se mancante
  # ------------------------------------------------------------------
  if (length(gui_yml$cellular_units %||% list()) == 0 &&
      !is.null(gui_yml$models)) {

    biomass_tags <- grep("^biomass_e_", name_vec, value = TRUE)
    suff         <- sub("^biomass_e_", "", biomass_tags)
    if (length(suff) < length(gui_yml$models))
      suff <- c(suff, sprintf("u%02d", seq_len(length(gui_yml$models))))[seq_len(length(gui_yml$models))]

    gui_yml$cellular_units <- mapply(function(m, lab) {
      dat <- gui_yml$models[[m]]
      list(
        model_name      = m,
        label           = lab,
        mu_max          = dat$params$mu_max %||% 1,
        biomass         = list(
          max  = dat$params$bioMax,
          mean = dat$params$bioMean,
          min  = dat$params$bioMin
        ),
        population      = list(
          starv = dat$params$starv,
          dup   = dat$params$dup,
          death = dat$params$death
        ),
        initial_biomass = dat$initial_biomass,
        initial_count   = dat$population
      )
    },
    m   = names(gui_yml$models),
    lab = suff,
    SIMPLIFY = FALSE)
  }

  # ------------------------------------------------------------------
  # 5a) boundary metabolites
  # ------------------------------------------------------------------
  bc_list <- gui_yml$simulation$boundary_concentrations %||% list()
  for (met in names(bc_list)) {
    pos <- which(name_vec == met)
    if (length(pos)) val_vec[pos] <- as.character(bc_list[[met]])
  }

  # ------------------------------------------------------------------
  # 5b) initial biomass
  # ------------------------------------------------------------------
  for (unit in gui_yml$cellular_units %||% list()) {
    if (is.null(unit$label) || unit$label == "")           # label fallback
      unit$label <- substr(gsub("[^a-z]", "", tolower(unit$model_name)), 1, 4)

    target <- paste0("biomass_e_", unit$label)
    pos    <- which(name_vec == target)
    if (length(pos))
      val_vec[pos] <- as.character(unit$initial_biomass %||% unit$biomass$mean)
  }

  # ------------------------------------------------------------------
  # 5c) initial population
  # ------------------------------------------------------------------
  for (unit in gui_yml$cellular_units %||% list()) {
    target <- paste0("n_", unit$label)
    pos    <- which(name_vec == target)
    if (length(pos))
      val_vec[pos] <- as.character(unit$initial_count)
  }

  # write patched stub
  fun_lines[ini_idx[1]] <- sub(
    "^\\s*(y_ini\\s*<-\\s*c\\().*(\\).*)$",
    paste0("\\1", paste(val_vec, collapse = ", "), "\\2"),
    fun_lines[ini_idx[1]]
  )
  writeLines(fun_lines, gui_fun_fname)


  # ───────────────────────────────────
  # 7) Call epimod::model.analysis
  # ───────────────────────────────────
  message("[DEBUG] Launching epimod::model.analysis …")
  results <- epimod::model.analysis(
    solver_fname     = solver_fname,
    parameters_fname = parameters_fname,
    functions_fname  = gui_fun_fname,
    debug            = TRUE,
    i_time           = i_time,
    f_time           = f_time,
    s_time           = s_time,
    atol             = atol,
    rtol             = rtol,
    fba_fname        = fba_fname,
    user_files       = user_files,
    volume           = volume
  )
  message("[DEBUG] epimod::model.analysis() returned")

  invisible(results)
}

