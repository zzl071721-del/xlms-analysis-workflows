suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

# ============================================================
# Generic pLink batch summary from filtered spectra CSV files, split by run ID.
#
# For mixed-run report folders like:
#   reports_10001_10002_10003
#   reports_10011_10012_10013
#
# Output:
# 1) Per report_dir / link_type / run_id:
#    <report_dir>/plink_stats_out_batch_per_run_spectra_only/<link_type>/<run_id>/stats.csv
#    <report_dir>/plink_stats_out_batch_per_run_spectra_only/<link_type>/<run_id>/csm_rows.csv
#    <report_dir>/plink_stats_out_batch_per_run_spectra_only/<link_type>/<run_id>/peptide_pair_rows.csv
#    <report_dir>/plink_stats_out_batch_per_run_spectra_only/<link_type>/<run_id>/protein_inter_rows.csv
#    <report_dir>/plink_stats_out_batch_per_run_spectra_only/<link_type>/<run_id>/protein_intra_rows.csv
#
# 2) Per link_type summary:
#    <report_dir>/plink_stats_out_batch_per_run_spectra_only/<link_type>/all_runs_stats.csv
#
# 3) Global summary:
#    <base_dir>/batch_plink_spectra_summary_per_run.csv
#
# New:
# - Skip already processed link_type if all_runs_stats.csv exists
# - Set overwrite = TRUE to force re-run
# ============================================================

# -------------------------
# USER SETTINGS
# -------------------------
# Root directory containing pLink report folders. Change this before running,
# or define the PLINK_RESULTS_DIR environment variable.
base_dir <- Sys.getenv("PLINK_RESULTS_DIR", unset = "path/to/plink/results")
link_types <- c("cross-linked", "mono-linked", "loop-linked", "regular")


report_pattern <- "^reports_(\\d+_){1,}\\d+$"

# Set TRUE to overwrite previously generated summaries.
overwrite <- FALSE

# -------------------------
# helpers
# -------------------------
pair_key2 <- function(x, y, sep = "||") {
  paste0(pmin(x, y), sep, pmax(x, y))
}

get_acc <- function(x) {
  str_match(x, "^[a-z]{2}\\|([^|]+)\\|")[, 2]
}

read_any_csv <- function(path) {
  readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
}

pick_one_file <- function(files) {
  files[which.max(file.info(files)$mtime)]
}

extract_run_id_from_string <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[is.na(x) | x == ""] <- NA_character_
  
  if (all(is.na(x))) return(x)
  
  s <- basename(x)
  
  rid1 <- str_match(s, "_(\\d+)\\.[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+\\.(?i:dta)$")[, 2]
  rid2 <- str_match(s, "_(\\d+)\\.(?:[0-9]+\\.){2,6}[0-9]+\\.(?i:dta)$")[, 2]
  rid  <- dplyr::coalesce(rid1, rid2)
  
  rid3 <- str_match(s, "(\\d+)(?i)\\.(raw|d|mgf|mzml|mzxml|wiff|wiff2|tsf|tdf)$")[, 1]
  rid3 <- str_match(rid3, "(\\d+)")[, 2]
  rid  <- dplyr::coalesce(rid, rid3)
  
  rid <- ifelse(
    is.na(rid),
    str_replace(s, "(?i)\\.(raw|d|mgf|mzml|mzxml|wiff|wiff2|tsf|tdf|dta)$", ""),
    rid
  )
  
  rid
}

parse_first_protein_pair <- function(proteins_str) {
  proteins_str <- str_trim(proteins_str)
  
  if (is.na(proteins_str) || proteins_str == "") {
    return(tibble(
      first_pair = NA_character_,
      acc_a = NA_character_,
      acc_b = NA_character_
    ))
  }
  
  first_pair <- str_split_fixed(proteins_str, "/", 2)[, 1]
  first_pair <- str_trim(first_pair)
  
  if (!str_detect(first_pair, "-")) {
    return(tibble(
      first_pair = first_pair,
      acc_a = NA_character_,
      acc_b = NA_character_
    ))
  }
  
  a <- str_split_fixed(first_pair, "-", 2)[, 1] |> str_trim()
  b <- str_split_fixed(first_pair, "-", 2)[, 2] |> str_trim()
  
  tibble(
    first_pair = first_pair,
    acc_a = get_acc(a),
    acc_b = get_acc(b)
  )
}

prep_spectra <- function(spectra_raw, report_dir, link_type, spectra_path) {
  df <- spectra_raw
  
  names(df) <- names(df) %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("[\\(\\)]", "_") %>%
    str_replace_all("__+", "_")
  
  pick <- function(cands) {
    hit <- intersect(names(df), cands)
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  
  c_title    <- pick(c("Title", "title"))
  c_score    <- pick(c("Score", "score"))
  c_evalue   <- pick(c("Evalue", "evalue", "EValue", "E_Value", "E-value", "E_value"))
  c_pep      <- pick(c("Peptide", "peptide"))
  c_linker   <- pick(c("Linker", "linker"))
  c_mods     <- pick(c("Modifications", "modifications"))
  c_proteins <- pick(c("Proteins", "Protein", "proteins", "protein"))
  
  need <- c(c_title, c_score, c_evalue, c_pep, c_linker, c_mods)
  if (any(is.na(need))) {
    stop(paste0(
      "Missing required columns in spectra file: ", spectra_path,
      "\nNeed Title / Score / Evalue / Peptide / Linker / Modifications"
    ))
  }
  
  df %>%
    mutate(
      report_dir   = report_dir,
      link_type    = link_type,
      spectra_path = spectra_path,
      Title = .data[[c_title]],
      Score = suppressWarnings(as.numeric(.data[[c_score]])),
      Evalue = suppressWarnings(as.numeric(.data[[c_evalue]])),
      Peptide = .data[[c_pep]],
      Linker = .data[[c_linker]],
      Modifications = dplyr::coalesce(.data[[c_mods]], "null"),
      Proteins = if (!is.na(c_proteins)) .data[[c_proteins]] else NA_character_,
      run_id = extract_run_id_from_string(.data[[c_title]]),
      csm_key = paste0(spectra_path, "||", Title),
      pep_a_raw = str_split_fixed(Peptide, "-", 2)[, 1],
      pep_b_raw = {
        tmp <- str_split_fixed(Peptide, "-", 2)
        ifelse(tmp[, 2] == "", tmp[, 1], tmp[, 2])
      },
      pep_pair_key = paste0(
        pair_key2(pep_a_raw, pep_b_raw),
        "||", Modifications,
        "||", Linker
      )
    )
}

read_one_linktype <- function(base_dir, report_dir, link_type) {
  rd <- file.path(base_dir, report_dir)
  
  sp <- list.files(rd, full.names = TRUE, recursive = FALSE) |>
    keep(~ str_detect(
      basename(.x),
      paste0("^result_.*\\.filtered_", link_type, "_spectra\\.csv$")
    ))
  
  spectra_raw <- tibble()
  spectra_path <- NA_character_
  
  if (length(sp) > 0) {
    spectra_path <- pick_one_file(sp)
    spectra_raw <- read_any_csv(spectra_path)
  }
  
  list(
    spectra_raw = spectra_raw,
    spectra_path = spectra_path
  )
}

process_one_linktype <- function(base_dir, report_dir, link_type, overwrite = FALSE) {
  rd <- file.path(base_dir, report_dir)
  out_root <- file.path(rd, "plink_stats_out_batch_per_run_spectra_only", link_type)
  dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
  
  done_file <- file.path(out_root, "all_runs_stats.csv")
  
  # -------------------------
  # skip if already done
  # -------------------------
  if (!overwrite && file.exists(done_file)) {
    message("    skip existing: ", report_dir, " / ", link_type)
    old_stats <- read_any_csv(done_file) %>%
      mutate(
        csm = suppressWarnings(as.numeric(csm)),
        peptide_pairs = suppressWarnings(as.numeric(peptide_pairs)),
        protein_inter = suppressWarnings(as.numeric(protein_inter)),
        protein_intra = suppressWarnings(as.numeric(protein_intra))
      )
    return(old_stats)
  }
  
  x <- read_one_linktype(base_dir, report_dir, link_type)
  
  if (nrow(x$spectra_raw) == 0) {
    stats_empty <- tibble(
      report_dir = report_dir,
      run_id = NA_character_,
      link_type = link_type,
      csm = 0,
      peptide_pairs = 0,
      protein_inter = 0,
      protein_intra = 0
    )
    write_csv(stats_empty, done_file)
    return(stats_empty)
  }
  
  spectra <- prep_spectra(x$spectra_raw, report_dir, link_type, x$spectra_path)
  
  run_ids <- spectra %>%
    filter(!is.na(run_id), run_id != "") %>%
    distinct(run_id) %>%
    pull(run_id)
  
  if (length(run_ids) == 0) {
    run_ids <- "ALL"
  }
  
  stats_runs <- purrr::map_dfr(run_ids, function(rid) {
    out_dir_run <- file.path(out_root, rid)
    dir.create(out_dir_run, showWarnings = FALSE, recursive = TRUE)
    
    sp_run <- if (rid == "ALL") {
      spectra
    } else {
      spectra %>% filter(run_id == rid)
    }
    
    csm_rows <- sp_run %>%
      arrange(Evalue, desc(Score)) %>%
      distinct(csm_key, .keep_all = TRUE)
    
    pep_rows <- sp_run %>%
      arrange(Evalue, desc(Score)) %>%
      distinct(pep_pair_key, .keep_all = TRUE)
    
    protein_inter_rows <- tibble()
    protein_intra_rows <- tibble()
    
    if (link_type == "cross-linked") {
      prot_pairs <- csm_rows %>%
        transmute(csm_key, Proteins) %>%
        bind_cols(purrr::map_dfr(.$Proteins, parse_first_protein_pair)) %>%
        filter(!is.na(acc_a), !is.na(acc_b)) %>%
        mutate(
          protein_pair_key = pair_key2(acc_a, acc_b),
          is_inter = acc_a != acc_b,
          is_intra = acc_a == acc_b
        ) %>%
        group_by(protein_pair_key) %>%
        summarise(
          acc_a = dplyr::first(acc_a),
          acc_b = dplyr::first(acc_b),
          first_pair = dplyr::first(first_pair),
          n_csm_support = n_distinct(csm_key),
          .groups = "drop"
        ) %>%
        mutate(
          is_inter = acc_a != acc_b,
          is_intra = acc_a == acc_b
        )
      
      protein_inter_rows <- prot_pairs %>% filter(is_inter)
      protein_intra_rows <- prot_pairs %>% filter(is_intra)
    }
    
    stats <- tibble(
      report_dir = report_dir,
      run_id = rid,
      link_type = link_type,
      csm = nrow(csm_rows),
      peptide_pairs = nrow(pep_rows),
      protein_inter = if (link_type == "cross-linked") nrow(protein_inter_rows) else 0,
      protein_intra = if (link_type == "cross-linked") nrow(protein_intra_rows) else 0
    )
    
    write_csv(csm_rows,           file.path(out_dir_run, "csm_rows.csv"))
    write_csv(pep_rows,           file.path(out_dir_run, "peptide_pair_rows.csv"))
    write_csv(protein_inter_rows, file.path(out_dir_run, "protein_inter_rows.csv"))
    write_csv(protein_intra_rows, file.path(out_dir_run, "protein_intra_rows.csv"))
    write_csv(stats,              file.path(out_dir_run, "stats.csv"))
    
    stats
  })
  
  write_csv(stats_runs, done_file)
  stats_runs
}

process_one_report <- function(report_dir, base_dir, link_types, overwrite = FALSE) {
  message("Processing report: ", report_dir)
  
  purrr::map_dfr(link_types, function(lt) {
    message("  link_type: ", lt)
    process_one_linktype(base_dir, report_dir, lt, overwrite = overwrite)
  })
}

# -------------------------
# RUN BATCH
# -------------------------
all_report_dirs <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
target_report_dirs <- all_report_dirs[str_detect(all_report_dirs, report_pattern)]

if (length(target_report_dirs) == 0) {
  stop("No report folders matched: ", report_pattern)
}

all_stats <- purrr::map_dfr(target_report_dirs, function(rd) {
  tryCatch(
    process_one_report(rd, base_dir, link_types, overwrite = overwrite),
    error = function(e) {
      message("Failed: ", rd, " ; reason: ", e$message)
      tibble(
        report_dir = rd,
        run_id = NA_character_,
        link_type = NA_character_,
        csm = NA_real_,
        peptide_pairs = NA_real_,
        protein_inter = NA_real_,
        protein_intra = NA_real_
      )
    }
  )
})

all_stats_clean <- all_stats %>%
  filter(!is.na(link_type)) %>%
  arrange(report_dir, run_id, factor(link_type, levels = link_types))

out_all <- file.path(base_dir, "batch_plink_spectra_summary_per_run.csv")
write_csv(all_stats_clean, out_all)

message("Done.")
message("Global summary written to: ", out_all)
