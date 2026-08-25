# ============================================================
#  PELSA batch limma report — DIA-NN pr_matrix
#
#  Batch input:
#    report.pr_matrix_batch1.tsv
#    report.pr_matrix_batch2.tsv
#    report.pr_matrix_batch3.tsv
#
#  For each file:
#    detect each compound condition automatically
#    run compound vs DMSO
#    n = 4 treatment, n = 4 DMSO
# ============================================================

# -----------------------------
# 0) Install packages if missing
# -----------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma")
}

for (pkg in c("ggplot2", "ggrepel", "patchwork", "dplyr", "reshape2", "scales", "stringr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(dplyr)
  library(reshape2)
  library(scales)
  library(stringr)
})

# ============================================================
# 1) Paths and global settings
# ============================================================

IN_DIR <- Sys.getenv(
  "PELSA_INPUT_DIR",
  unset = "path/to/dia_nn_reports"
)
OUT_DIR <- Sys.getenv(
  "PELSA_OUTPUT_DIR",
  unset = file.path(IN_DIR, "batch_results")
)

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

INPUT_FILES <- list.files(
  IN_DIR,
  pattern = "^report\\.pr_matrix.*\\.tsv$",
  full.names = TRUE
)

if (length(INPUT_FILES) == 0) {
  stop("No report.pr_matrix*.tsv files found in: ", IN_DIR)
}

DMSO_PATTERN <- "DMSO"

ADJ_P_CUTOFF <- 0.05
LOG2_FC_CUTOFF <- 1.6

# For n = 4 treatment and n = 4 DMSO, complete case is acceptable.
# TRUE means peptide must be quantified in all 8 samples for that comparison.
USE_COMPLETE_CASE <- TRUE


# Optional gene highlight
# Example: "TARGET_GENE"
# Set to "" to disable
HIGHLIGHT_GENE <- ""
HIGHLIGHT_P_THRESH <- 0.05

# ============================================================
# 2) Colors
# ============================================================

ORANGE <- "#F97316"
RED    <- "#E74C3C"
BLUE   <- "#5B9BD5"
GREY   <- "#CCCCCC"

# ============================================================
# 3) Helper functions
# ============================================================

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

extract_condition_from_raw <- function(x) {
  x <- basename(x)
  x <- sub("\\.d$", "", x)
  
  # Supported examples:
  #   study_COMPOUND_A_R1.d
  #   study_Comp1_R2.d
  #   study_CP1_R3.d
  #   study_DMSO1_R4.d
  condition <- str_match(
    x,
    regex(
      "(?:^|_)(COMPOUND_[A-Z0-9]+|COMPOUND[0-9]+|COMP[0-9]+|CP[0-9]+|DMSO[0-9]*)(?:_|\\.|$)",
      ignore_case = TRUE
    )
  )[, 2]
  
  condition
}

extract_replicate_from_raw <- function(x) {
  x <- basename(x)
  x <- sub("\\.d$", "", x)
  
  rep <- str_extract(x, "R[0-9]+")
  rep[is.na(rep)] <- paste0("R", seq_len(sum(is.na(rep))))
  
  rep
}

short_sample_name <- function(x) {
  condition <- extract_condition_from_raw(x)
  replicate <- extract_replicate_from_raw(x)
  
  out <- paste0(condition, "_", replicate)
  out <- make.unique(out, sep = "_")
  out
}

detect_intensity_columns <- function(df) {
  int_cols <- grep("\\.d$", colnames(df), value = TRUE)
  
  if (length(int_cols) == 0) {
    int_cols <- grep("PELSA|Slot|\\.raw|\\.d", colnames(df), value = TRUE)
  }
  
  if (length(int_cols) == 0) {
    stop("No intensity columns detected.")
  }
  
  int_cols
}

detect_compounds_in_file <- function(df, dmso_pattern = DMSO_PATTERN) {
  int_cols <- detect_intensity_columns(df)
  
  conditions <- extract_condition_from_raw(int_cols)
  conditions <- unique(conditions)
  conditions <- conditions[!is.na(conditions) & conditions != ""]
  
  dmso_conditions <- conditions[grepl(dmso_pattern, conditions, ignore.case = TRUE)]
  treat_conditions <- setdiff(conditions, dmso_conditions)
  
  if (length(dmso_conditions) == 0) {
    stop("No DMSO condition detected. Conditions found: ", paste(conditions, collapse = ", "))
  }
  
  if (length(treat_conditions) == 0) {
    stop("No treatment condition detected. Conditions found: ", paste(conditions, collapse = ", "))
  }
  
  treat_conditions
}

# ============================================================
# 4) Data loading and normalization
# ============================================================

load_normalize <- function(path,
                           compound_condition,
                           dmso_pattern = DMSO_PATTERN,
                           exclude_pat = NULL,
                           use_complete_case = TRUE,
                           min_valid_treat = 3,
                           min_valid_dmso = 3) {
  
  message("Reading: ", path)
  message("Compound: ", compound_condition)
  
  df <- read.delim(
    path,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_cols <- c(
    "Protein.Group",
    "Protein.Ids",
    "Protein.Names",
    "Genes",
    "Modified.Sequence",
    "Precursor.Charge",
    "Precursor.Id"
  )
  
  missing_required <- setdiff(required_cols, colnames(df))
  if (length(missing_required) > 0) {
    stop("Missing required columns: ", paste(missing_required, collapse = ", "))
  }
  
  df$sequence_charge <- paste0(df$`Modified.Sequence`, "_", df$`Precursor.Charge`)
  
  # Remove duplicated precursor entries
  df <- df[!duplicated(df$sequence_charge), ]
  
  int_cols <- detect_intensity_columns(df)
  
  if (!is.null(exclude_pat)) {
    int_cols <- int_cols[!grepl(exclude_pat, int_cols, perl = TRUE, ignore.case = TRUE)]
  }
  
  conditions <- extract_condition_from_raw(int_cols)
  
  t_cols <- int_cols[conditions == compound_condition]
  d_cols <- int_cols[grepl(dmso_pattern, conditions, ignore.case = TRUE)]
  
  message("Treatment columns: ", length(t_cols))
  message("DMSO columns:      ", length(d_cols))
  
  if (length(t_cols) != 4) {
    warning("Treatment column number is not 4 for ", compound_condition, ": n = ", length(t_cols))
  }
  
  if (length(d_cols) != 4) {
    warning("DMSO column number is not 4: n = ", length(d_cols))
  }
  
  if (length(t_cols) < 2) {
    stop("Too few treatment columns detected for: ", compound_condition)
  }
  
  if (length(d_cols) < 2) {
    stop("Too few DMSO columns detected.")
  }
  
  all_cols <- c(t_cols, d_cols)
  
  for (col in all_cols) {
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
    df[[col]][!is.na(df[[col]]) & df[[col]] < 1] <- NA
  }
  
  if (use_complete_case) {
    keep_idx <- complete.cases(df[, all_cols])
  } else {
    n_treat_valid <- rowSums(!is.na(df[, t_cols, drop = FALSE]))
    n_dmso_valid  <- rowSums(!is.na(df[, d_cols, drop = FALSE]))
    
    keep_idx <- n_treat_valid >= min_valid_treat &
      n_dmso_valid >= min_valid_dmso
  }
  
  pc <- df[keep_idx, ]
  
  message("Peptides after filtering: ", nrow(pc))
  
  if (nrow(pc) == 0) {
    stop("No peptides left after filtering. Try USE_COMPLETE_CASE <- FALSE.")
  }
  
  mat <- as.matrix(pc[, all_cols])
  mat <- log2(mat)
  
  # Median normalization inside this one comparison: compound n=4 + DMSO n=4
  meds <- apply(mat, 2, median, na.rm = TRUE)
  global_med <- median(meds, na.rm = TRUE)
  nrm <- sweep(mat, 2, meds - global_med)
  
  colnames(nrm) <- all_cols
  
  list(
    df = df,
    pc = pc,
    nrm = nrm,
    all_cols = all_cols,
    t_cols = t_cols,
    d_cols = d_cols,
    sample_names = short_sample_name(all_cols)
  )
}

# ============================================================
# 5) limma
# ============================================================

run_limma <- function(dat) {
  
  group <- factor(
    c(
      rep("treat", length(dat$t_cols)),
      rep("dmso",  length(dat$d_cols))
    ),
    levels = c("dmso", "treat")
  )
  
  design <- model.matrix(~0 + group)
  colnames(design) <- c("dmso", "treat")
  
  fit <- lmFit(dat$nrm, design)
  cont <- makeContrasts(treat - dmso, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont))
  
  res <- topTable(
    fit2,
    coef = 1,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  )
  
  extra_cols <- intersect(
    c(
      "Protein.Group",
      "Protein.Ids",
      "Protein.Names",
      "Genes",
      "First.Protein.Description",
      "Proteotypic",
      "Stripped.Sequence",
      "Modified.Sequence",
      "Precursor.Charge",
      "Precursor.Id",
      "Precursor.Mz",
      "sequence_charge"
    ),
    colnames(dat$pc)
  )
  
  out <- cbind(
    dat$pc[, extra_cols, drop = FALSE],
    logFC     = res$logFC,
    P.Value   = res$P.Value,
    adj.P.Val = res$adj.P.Val,
    AveExpr   = res$AveExpr
  )
  
  out$neg_log10_p <- -log10(pmax(out$P.Value, 1e-300))
  
  out
}

# ============================================================
# 6) Plot functions
# ============================================================

plot_volcano <- function(rf, report_title, highlight_df = NULL) {
  
  sig <- !is.na(rf$adj.P.Val) &
    !is.na(rf$logFC) &
    rf$adj.P.Val < ADJ_P_CUTOFF &
    abs(rf$logFC) >= LOG2_FC_CUTOFF
  
  p <- ggplot() +
    geom_point(
      data = rf[!sig, ],
      aes(logFC, neg_log10_p),
      color = GREY,
      alpha = 0.35,
      size = 0.7
    ) +
    geom_point(
      data = rf[sig, ],
      aes(logFC, neg_log10_p),
      color = RED,
      alpha = 0.85,
      size = 2
    )
  
  if (any(sig)) {
    sig_df <- rf[sig, ]
    sig_df$charge <- sig_df$Precursor.Charge
    
    sig_df$label <- ifelse(
      !is.na(sig_df$Genes) & sig_df$Genes != "",
      paste0(sig_df$Genes, "\n", sig_df$Modified.Sequence, " (+", sig_df$charge, ")"),
      paste0(sig_df$Modified.Sequence, " (+", sig_df$charge, ")")
    )
    
    p <- p +
      geom_text_repel(
        data = sig_df,
        aes(logFC, neg_log10_p, label = label),
        color = "#7B0000",
        size = 3.3,
        fontface = "bold",
        lineheight = 0.85,
        max.overlaps = 30,
        min.segment.length = 0.1,
        segment.color = "grey70",
        segment.size = 0.3
      )
  }
  
  if (!is.null(highlight_df) && nrow(highlight_df) > 0) {
    
    highlight_df$charge <- highlight_df$Precursor.Charge
    
    highlight_df$label <- paste0(
      highlight_df$Genes,
      "\n",
      highlight_df$Modified.Sequence,
      " (+",
      highlight_df$charge,
      ")"
    )
    
    p <- p +
      geom_point(
        data = highlight_df,
        aes(logFC, neg_log10_p),
        fill = ORANGE,
        color = "#7c2d00",
        shape = 21,
        size = 2.8,
        stroke = 0.6
      ) +
      geom_text_repel(
        data = highlight_df,
        aes(logFC, neg_log10_p, label = label),
        color = "#7c2d00",
        size = 3.2,
        fontface = "bold",
        lineheight = 0.85,
        max.overlaps = 40,
        min.segment.length = 0.1,
        segment.color = "#7c2d00",
        segment.size = 0.3
      )
  }
  
  p +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4) +
    geom_vline(
      xintercept = c(-LOG2_FC_CUTOFF, LOG2_FC_CUTOFF),
      linetype = "dashed",
      color = "grey65",
      linewidth = 0.4
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey65",
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = sprintf(
        "adj.P < %.2f; |log2FC| >= %.1f; n = %d",
        ADJ_P_CUTOFF,
        LOG2_FC_CUTOFF,
        sum(sig)
      ),
      hjust = -0.05,
      vjust = 1.5,
      size = 5.2,
      color = "grey30"
    ) +
    labs(
      x = expression(log[2] ~ "Fold Change"),
      y = expression(-log[10] ~ "P-value"),
      title = report_title
    ) +
    theme_classic(base_size = 17) +
    theme(
      plot.title = element_text(face = "bold", size = 19),
      aspect.ratio = 1,
      axis.text = element_text(size = 15)
    )
}

plot_id_counts <- function(dat) {
  
  df <- dat$df
  all_cols <- dat$all_cols
  pc <- dat$pc
  
  cnts <- sapply(all_cols, function(c) sum(!is.na(df[[c]])))
  cnts <- c(cnts, `All_quant` = nrow(pc))
  
  labels <- c(dat$sample_names, "All_quant")
  labels <- make.unique(labels, sep = "_")
  
  cols <- c(rep(BLUE, length(all_cols)), "#E36C09")
  
  bdf <- data.frame(
    label = factor(labels, levels = labels),
    count = as.numeric(cnts),
    col = cols
  )
  
  ggplot(bdf, aes(label, count, fill = col)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.4) +
    geom_text(
      aes(label = sprintf("%.1fk", count / 1000)),
      vjust = -0.3,
      size = 4.0
    ) +
    scale_fill_identity() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      x = NULL,
      y = "# peptide precursors",
      title = "IDs per run"
    ) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", size = 17),
      axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 13)
    )
}

plot_cv_violin <- function(dat) {
  
  nrm <- dat$nrm
  t_cols <- dat$t_cols
  d_cols <- dat$d_cols
  all_cols <- dat$all_cols
  
  cv_group <- function(cols) {
    lin <- 2^nrm[, match(cols, all_cols), drop = FALSE]
    
    cv <- apply(
      lin,
      1,
      function(x) {
        sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
      }
    ) * 100
    
    cv[is.finite(cv) & cv > 0 & cv < 200]
  }
  
  vdf <- rbind(
    data.frame(group = "Control",   cv = cv_group(d_cols)),
    data.frame(group = "Treatment", cv = cv_group(t_cols))
  )
  
  vdf$group <- factor(vdf$group, levels = c("Control", "Treatment"))
  
  med_df <- vdf %>%
    group_by(group) %>%
    summarise(med = median(cv, na.rm = TRUE), .groups = "drop")
  
  fill_vals <- c("Control" = BLUE, "Treatment" = RED)
  
  ggplot(vdf, aes(group, cv, fill = group)) +
    geom_violin(alpha = 0.75, width = 0.6, color = NA) +
    geom_boxplot(
      width = 0.15,
      outlier.shape = NA,
      color = "grey30",
      linewidth = 0.4
    ) +
    geom_text(
      data = med_df,
      aes(group, med + 3, label = sprintf("%.1f%%", med)),
      size = 5.2,
      fontface = "bold"
    ) +
    scale_fill_manual(values = fill_vals) +
    scale_y_continuous(limits = c(0, 80)) +
    labs(
      x = NULL,
      y = "CV (%)",
      title = "Peptide CV"
    ) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", size = 17),
      legend.position = "none",
      axis.text = element_text(size = 14)
    )
}

plot_correlation <- function(dat) {
  
  labels <- dat$sample_names
  labels <- make.unique(labels, sep = "_")
  
  corr <- cor(dat$nrm, use = "pairwise.complete.obs")
  dimnames(corr) <- list(labels, labels)
  
  mdf <- melt(corr)
  mdf$label <- sprintf("%.3f", mdf$value)
  
  ggplot(mdf, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    geom_text(
      aes(label = label),
      size = ifelse(ncol(corr) <= 8, 4.6, 3.5)
    ) +
    scale_fill_gradientn(
      colors = c("#d73027", "#fee08b", "#1a9850"),
      limits = c(0.9, 1.0),
      name = "Pearson r",
      oob = scales::squish
    ) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    coord_fixed() +
    labs(
      x = NULL,
      y = NULL,
      title = "Sample correlation"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", size = 17),
      axis.text = element_text(size = 11),
      legend.key.height = unit(0.5, "cm"),
      legend.text = element_text(size = 11)
    )
}

plot_sample_boxplot <- function(dat) {
  
  m <- dat$nrm
  labels <- dat$sample_names
  labels <- make.unique(labels, sep = "_")
  colnames(m) <- labels
  
  mdf <- melt(m)
  colnames(mdf) <- c("Peptide", "Sample", "log2Intensity")
  
  ggplot(mdf, aes(Sample, log2Intensity)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3) +
    labs(
      x = NULL,
      y = "log2 intensity after median normalization",
      title = "Intensity distribution"
    ) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", size = 17),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 13)
    )
}

# ============================================================
# 7) One comparison
# ============================================================

run_one_comparison <- function(input_file, compound_condition) {
  
  file_tag <- tools::file_path_sans_ext(basename(input_file))
  file_tag <- sub("^report\\.pr_matrix_", "", file_tag)
  
  comparison_label <- paste0(compound_condition, "_vs_DMSO")
  
  message("\n============================================================")
  message("File: ", basename(input_file))
  message("Running comparison: ", comparison_label)
  message("============================================================\n")
  
  comp_out_dir <- file.path(
    OUT_DIR,
    safe_filename(file_tag),
    safe_filename(comparison_label)
  )
  
  if (!dir.exists(comp_out_dir)) {
    dir.create(comp_out_dir, recursive = TRUE)
  }
  
  report_title <- paste0(compound_condition, " vs DMSO PELSA")
  output_pdf <- paste0(safe_filename(comparison_label), "_PELSA_report.pdf")
  
  dat <- load_normalize(
    path = input_file,
    compound_condition = compound_condition,
    dmso_pattern = DMSO_PATTERN,
    exclude_pat = EXCLUDE_PATTERN,
    use_complete_case = USE_COMPLETE_CASE,
    min_valid_treat = MIN_VALID_TREAT,
    min_valid_dmso = MIN_VALID_DMSO
  )
  
  meta <- data.frame(
    Sample = dat$sample_names,
    Raw.File = basename(dat$all_cols),
    Group = c(
      rep("TREAT", length(dat$t_cols)),
      rep("DMSO",  length(dat$d_cols))
    ),
    stringsAsFactors = FALSE
  )
  
  message("\nSample metadata:")
  print(meta, row.names = FALSE)
  
  message("\nTreatment samples: ", length(dat$t_cols))
  message("DMSO samples:      ", length(dat$d_cols))
  message("Peptides used:     ", nrow(dat$pc), "\n")
  
  rf <- run_limma(dat)
  
  highlight_df <- NULL
  
  if (!is.null(HIGHLIGHT_GENE) && HIGHLIGHT_GENE != "") {
    highlight_df <- rf[
      grepl(HIGHLIGHT_GENE, rf$Genes, ignore.case = TRUE) &
        rf$P.Value < HIGHLIGHT_P_THRESH,
    ]
  }
  
  n_hits <- sum(
    rf$adj.P.Val < ADJ_P_CUTOFF & abs(rf$logFC) >= LOG2_FC_CUTOFF,
    na.rm = TRUE
  )
  message(
    "Hits at adj.P < ", ADJ_P_CUTOFF,
    " and |log2FC| >= ", LOG2_FC_CUTOFF,
    ": ", n_hits
  )
  
  p_volcano     <- plot_volcano(rf, report_title, highlight_df = highlight_df)
  p_ids         <- plot_id_counts(dat)
  p_cv          <- plot_cv_violin(dat)
  p_correlation <- plot_correlation(dat)
  p_boxplot     <- plot_sample_boxplot(dat)
  
  page <- (
    p_volcano |
      ((p_ids / p_cv) | (p_correlation / p_boxplot))
  ) +
    plot_layout(widths = c(1.05, 1.3)) +
    plot_annotation(
      title = paste("PELSA QC & Analysis —", report_title),
      subtitle = sprintf("Treatment = %s  |  Control = DMSO", compound_condition),
      theme = theme(
        plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5, color = "grey40")
      )
    )
  
  # -----------------------------
  # Save figures
  # -----------------------------
  
  out_main <- file.path(comp_out_dir, output_pdf)
  
  ggsave(
    out_main,
    plot = page,
    width = 24,
    height = 13,
    units = "in",
    device = "pdf"
  )
  
  out_volcano <- file.path(comp_out_dir, sub("\\.pdf$", "_volcano.pdf", output_pdf))
  out_ids     <- file.path(comp_out_dir, sub("\\.pdf$", "_ids.pdf", output_pdf))
  out_cv      <- file.path(comp_out_dir, sub("\\.pdf$", "_cv.pdf", output_pdf))
  out_corr    <- file.path(comp_out_dir, sub("\\.pdf$", "_correlation.pdf", output_pdf))
  out_box     <- file.path(comp_out_dir, sub("\\.pdf$", "_boxplot.pdf", output_pdf))
  
  ggsave(out_volcano, plot = p_volcano,     width = 10, height = 10, units = "in", device = "pdf")
  ggsave(out_ids,     plot = p_ids,         width = 10, height = 6,  units = "in", device = "pdf")
  ggsave(out_cv,      plot = p_cv,          width = 6,  height = 6,  units = "in", device = "pdf")
  ggsave(out_corr,    plot = p_correlation, width = 8,  height = 7,  units = "in", device = "pdf")
  ggsave(out_box,     plot = p_boxplot,     width = 10, height = 6,  units = "in", device = "pdf")
  
  message("Saved report:      ", out_main)
  message("Saved volcano:     ", out_volcano)
  message("Saved IDs:         ", out_ids)
  message("Saved CV:          ", out_cv)
  message("Saved correlation: ", out_corr)
  message("Saved boxplot:     ", out_box)
  
  # -----------------------------
  # Export tables
  # -----------------------------
  
  keep <- intersect(
    c(
      "Protein.Group",
      "Protein.Ids",
      "Protein.Names",
      "Genes",
      "First.Protein.Description",
      "Proteotypic",
      "Stripped.Sequence",
      "Modified.Sequence",
      "Precursor.Charge",
      "Precursor.Id",
      "Precursor.Mz",
      "sequence_charge",
      "logFC",
      "P.Value",
      "adj.P.Val",
      "AveExpr",
      "neg_log10_p"
    ),
    colnames(rf)
  )
  
  rf_ordered <- rf[order(rf$adj.P.Val, rf$P.Value), ]
  
  all_out <- rf_ordered[, keep, drop = FALSE]
  all_out$Hit <- ifelse(
    !is.na(all_out$adj.P.Val) &
      !is.na(all_out$logFC) &
      all_out$adj.P.Val < ADJ_P_CUTOFF &
      abs(all_out$logFC) >= LOG2_FC_CUTOFF,
    "hit",
    "no hit"
  )
  all_out$Comparison <- comparison_label
  all_out$Input.File <- basename(input_file)
  all_out <- all_out[, c("Input.File", "Comparison", "Hit", setdiff(colnames(all_out), c("Input.File", "Comparison", "Hit")))]
  
  hits_out <- all_out[all_out$Hit == "hit", ]
  
  protein_summary <- all_out %>%
    group_by(Input.File, Comparison, Protein.Group, Genes, Protein.Names) %>%
    summarise(
      n_peptides = n(),
      n_hit_peptides = sum(Hit == "hit", na.rm = TRUE),
      best_adj.P.Val = min(adj.P.Val, na.rm = TRUE),
      best_P.Value = min(P.Value, na.rm = TRUE),
      median_logFC = median(logFC, na.rm = TRUE),
      mean_logFC = mean(logFC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(best_adj.P.Val, best_P.Value)
  
  all_file <- file.path(comp_out_dir, sub("\\.pdf$", "_all_peptides.tsv", output_pdf))
  hit_file <- file.path(comp_out_dir, sub("\\.pdf$", "_hits.tsv", output_pdf))
  protein_file <- file.path(comp_out_dir, sub("\\.pdf$", "_protein_summary.tsv", output_pdf))
  meta_file <- file.path(comp_out_dir, sub("\\.pdf$", "_sample_metadata.tsv", output_pdf))
  
  write.table(all_out, all_file, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(hits_out, hit_file, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(protein_summary, protein_file, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(meta, meta_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  message("All-peptides table saved: ", all_file, "  n = ", nrow(all_out))
  message("Hit table saved:          ", hit_file, "  n = ", nrow(hits_out))
  message("Protein summary saved:    ", protein_file, "  n = ", nrow(protein_summary))
  message("Metadata saved:           ", meta_file)
  
  list(
    input_file = basename(input_file),
    comparison = comparison_label,
    compound = compound_condition,
    all_peptides = all_out,
    hits = hits_out,
    protein_summary = protein_summary,
    meta = meta,
    out_dir = comp_out_dir,
    n_peptides = nrow(all_out),
    n_hits = nrow(hits_out)
  )
}

# ============================================================
# 8) Batch run all files
# ============================================================

message("\nInput files:")
print(basename(INPUT_FILES))

all_results <- list()

for (input_file in INPUT_FILES) {
  
  message("\n############################################################")
  message("Processing file: ", basename(input_file))
  message("############################################################\n")
  
  tmp_df <- read.delim(
    input_file,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    nrows = 5
  )
  
  compounds <- detect_compounds_in_file(tmp_df, dmso_pattern = DMSO_PATTERN)
  
  message("Detected compounds: ", paste(compounds, collapse = ", "))
  
  for (compound in compounds) {
    result_name <- paste0(
      tools::file_path_sans_ext(basename(input_file)),
      "__",
      compound,
      "_vs_DMSO"
    )
    
    all_results[[result_name]] <- run_one_comparison(
      input_file = input_file,
      compound_condition = compound
    )
  }
}

# ============================================================
# 9) Combined output across all compounds
# ============================================================

combined_all <- bind_rows(lapply(all_results, function(x) x$all_peptides))
combined_hits <- bind_rows(lapply(all_results, function(x) x$hits))
combined_protein_summary <- bind_rows(lapply(all_results, function(x) x$protein_summary))

batch_summary <- bind_rows(lapply(all_results, function(x) {
  data.frame(
    Input.File = x$input_file,
    Comparison = x$comparison,
    Compound = x$compound,
    N.Peptides = x$n_peptides,
    N.Hits = x$n_hits,
    Out.Dir = x$out_dir,
    stringsAsFactors = FALSE
  )
}))

combined_all_file <- file.path(OUT_DIR, "pelsa_batch_all_peptides_combined.tsv")
combined_hits_file <- file.path(OUT_DIR, "pelsa_batch_hits_combined.tsv")
combined_protein_file <- file.path(OUT_DIR, "pelsa_batch_protein_summary_combined.tsv")
batch_summary_file <- file.path(OUT_DIR, "pelsa_batch_run_summary.tsv")

write.table(combined_all, combined_all_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(combined_hits, combined_hits_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(combined_protein_summary, combined_protein_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(batch_summary, batch_summary_file, sep = "\t", row.names = FALSE, quote = FALSE)

message("\n============================================================")
message("Batch done.")
message("Combined all peptides:     ", combined_all_file)
message("Combined hits:             ", combined_hits_file)
message("Combined protein summary:  ", combined_protein_file)
message("Batch run summary:         ", batch_summary_file)
message("============================================================\n")

print(batch_summary)
