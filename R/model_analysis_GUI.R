#' General Wrapper for epimod::model.analysis (GUI) with console debug
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

  ## 0) Attach & validate
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

  ## 1) Build core paths
  solver_fname     <- fs::path(paths["gen"], paste0(hypernode_name, ".solver"))
  parameters_fname <- fs::path(paths["config"], "initial_data.csv")
  orig_fun_fname   <- fs::path(paths["src"],   paste0("functions_", hypernode_name, ".R"))
  message("[DEBUG] solver file: ", solver_fname)
  message("[DEBUG] parameters file: ", parameters_fname)
  message("[DEBUG] original functions.R: ", orig_fun_fname)

  ## 2) Ensure files exist
  needed <- c(solver_fname, parameters_fname, orig_fun_fname)
  missing <- needed[!file.exists(needed)]
  if (length(missing)) {
    stop("Missing required analysis files: ", paste(basename(missing), collapse = ", "))
  }
  message("[DEBUG] All core files exist")

  ## 3) Read the GUI snapshot YAML
  gui_yaml_path <- fs::path(paths["config"], paste0(hypernode_name, "_gui.yaml"))
  message("[DEBUG] Looking for GUI YAML at: ", gui_yaml_path)
  if (!file.exists(gui_yaml_path)) {
    stop("Missing GUI‐snapshot YAML: ", gui_yaml_path)
  }
  gui_yml <- yaml::read_yaml(gui_yaml_path)
  message("[DEBUG] Read GUI YAML successfully")

  ## 3a) Write mu_max_values_gui.csv from GUI YAML
  mu_defs <- gui_yml$cellular_units %||% list()
  if (length(mu_defs) > 0) {
    mu_df <- data.frame(
      Model  = vapply(mu_defs, function(u) u$model_name, character(1)),
      mu_max = vapply(mu_defs, function(u) u$mu_max, numeric(1)),
      stringsAsFactors = FALSE
    )
    mu_csv <- fs::path(paths["config"], "mu_max_values_gui.csv")
    message("[DEBUG] Writing mu_max CSV: ", mu_csv)
    write.csv(mu_df, mu_csv, row.names = FALSE, quote = FALSE)
    user_files <- c(user_files, mu_csv)
  }

  ## 3b) Write population_parameters.csv from GUI YAML (no header, no model names)
  pop_defs <- gui_yml$cellular_units %||% list()
  if (length(pop_defs) > 0) {
    pop_mat <- do.call(rbind, lapply(pop_defs, function(u) {
      c(u$population$starv, u$population$dup, u$population$death)
    }))
    pop_csv <- fs::path(paths["config"], "population_parameters.csv")
    message("[DEBUG] Writing population parameters CSV: ", pop_csv)
    write.table(pop_mat, pop_csv,
                sep = ",", row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    # replace any original population_parameters.csv in user_files
    user_files <- c(setdiff(user_files, fs::path(paths["config"], "population_parameters.csv")), pop_csv)
  }

  ## 4) Fix reverse‐bounds CSV volume if needed
  rev_csv <- fs::path(paths["output"], "non_projected_reverse_background_met_gui.csv")
  if (file.exists(rev_csv)) {
    df_r    <- read.csv(rev_csv, stringsAsFactors = FALSE)
    gui_vol <- gui_yml$simulation$system_parameters$volume
    if (!identical(df_r$volume[1], gui_vol)) {
      df_r$volume <- gui_vol
      write.csv(df_r, rev_csv, row.names = FALSE, quote = FALSE)
    }
  }

    ## 5) Prepare GUI-patched R stub with full y_ini patching
  gui_fun_fname <- fs::path(paths["src"], paste0("functions_", hypernode_name, "_gui.R"))
  message("[DEBUG] Copying functions.R to stub: ", gui_fun_fname)
  file.copy(orig_fun_fname, gui_fun_fname, overwrite = TRUE)

  fun_lines <- readLines(gui_fun_fname)
  message("[DEBUG] Read ", length(fun_lines), " lines from stub")

  # locate y_ini and yini.names
  ini_idx   <- grep("^\\s*y_ini\\s*<-\\s*c\\(", fun_lines)
  names_idx <- grep("^\\s*yini\\.names\\s*<-\\s*c\\(", fun_lines)
  if (!length(ini_idx) || !length(names_idx)) {
    stop("Cannot find y_ini or yini.names definitions in stub")
  }

  # parse & clean the names vector
  raw_names <- strsplit(sub("^.*c\\((.*)\\).*", "\\1", fun_lines[names_idx]), ",")[[1]]
  name_vec  <- gsub("['\"]", "", trimws(raw_names))
  message("[DEBUG] Cleaned yini.names: ", paste(name_vec, collapse = ", "))

  # parse the y_ini values
  val_vec   <- trimws(strsplit(sub("^.*c\\((.*)\\).*", "\\1", fun_lines[ini_idx[1]]), ",")[[1]])
  message("[DEBUG] Original y_ini values: ", paste(val_vec, collapse = ", "))

  # 5a) boundary metabolites go in the last positions
  bm_defs    <- gui_yml$boundary_metabolites %||% character()
  new_bounds <- as.character(unname(gui_yml$simulation$boundary_concentrations))
  if (length(bm_defs)) {
    nvals <- length(val_vec)
    val_vec[(nvals - length(new_bounds) + 1):nvals] <- new_bounds
    message("[DEBUG] After boundary patch: ", paste(val_vec, collapse = ", "))
  }

  # 5b) patch initial biomass (from each cellular_units$initial_biomass)
  for (unit in gui_yml$cellular_units %||% list()) {
    lbl    <- unit$label
    target <- paste0("biomass_e_", lbl)
    pos    <- which(name_vec == target)
    if (length(pos)) {
      cnt        <- as.character(unit$initial_biomass %||% unit$biomass$mean)
      val_vec[pos] <- cnt
      message("[DEBUG] Patching ", target, " at position ", pos, " to ", cnt)
    } else {
      message("[DEBUG] biomass entry ", target, " not found, skipping")
    }
  }

  # 5c) patch initial population (as before)
  for (unit in gui_yml$cellular_units %||% list()) {
    target <- paste0("n_", unit$label)
    pos    <- which(name_vec == target)
    if (length(pos)) {
      cnt        <- as.character(unit$initial_count)
      val_vec[pos] <- cnt
      message("[DEBUG] Patching ", target, " at position ", pos, " to ", cnt)
    } else {
      message("[DEBUG] n_ entry ", target, " not found, skipping")
    }
  }

  # reconstruct & write the new y_ini line
  old_y_line <- fun_lines[ini_idx[1]]
  new_y_line <- sub(
    "^\\s*(y_ini\\s*<-\\s*c\\().*(\\).*)$",
    paste0("\\1", paste(val_vec, collapse = ", "), "\\2"),
    old_y_line
  )
  fun_lines[ini_idx[1]] <- new_y_line
  message("[DEBUG] OLD y_ini line: ", old_y_line)
  message("[DEBUG] NEW y_ini line: ", new_y_line)

  writeLines(fun_lines, gui_fun_fname)
  message("[DEBUG] Wrote patched R stub to: ", gui_fun_fname)


  ## 6) Prepare GUI‐patched C++ stub
  cpp_files <- list.files(paths["src"], "\\.cpp$", full.names = TRUE)
  if (length(cpp_files)) {
    orig_cpp <- cpp_files[[1]]
    gui_cpp  <- fs::path(
      paths["src"],
      paste0(tools::file_path_sans_ext(basename(orig_cpp)), "_gui.cpp")
    )
    file.copy(orig_cpp, gui_cpp, overwrite = TRUE)
    cpp_lines <- readLines(gui_cpp)
    vol_val <- gui_yml$simulation$system_parameters$volume
    den_val <- gui_yml$simulation$system_parameters$cell_density
    cpp_lines <- sub("double V = .*?;", sprintf("double V = %g;", vol_val), cpp_lines)
    cpp_lines <- sub("long long int delta = .*?;", sprintf("long long int delta = %g;", den_val), cpp_lines)
    writeLines(cpp_lines, gui_cpp)
  }

  ## 7) Call epimod::model.analysis with GUI stubs
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

