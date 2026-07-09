#!/usr/bin/env Rscript
# ============================================================================
# mtDNA genotype concordance: Array-vs-Sequence and Imputed-vs-Sequence
# using the Matthews Correlation Coefficient (MCC)
#
# Pure base R (no CRAN packages required, e.g. no vcfR).
# ============================================================================

# ---------------------------------------------------------------------------
# 0. USER SETTINGS -- edit these paths / names for your data
# ---------------------------------------------------------------------------
array_vcf_path    <- "ArraySamples_rename.vcf.gz"      # or "array.vcf.gz"
imputed_vcf_path  <- "ArraySamples_rename.vcf.gz"    # or "imputed.vcf.gz"
sequence_vcf_path <- "ArraySubsetSamples.recode.vcf.gz"   # or "sequence.vcf.gz"  (truth set)

# mtDNA contig name in your VCFs. Leave NULL to auto-detect common names
# ("MT", "chrM", "chrMT", "M"). Set explicitly if auto-detect fails,
# e.g. mito_chr <- "chrM"
mito_chr <- NULL

out_dir <- "mtdna_concordance_results"

# ---------------------------------------------------------------------------
# 1. VCF READER (base R only)
# ---------------------------------------------------------------------------
read_vcf <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con)
  
  header_idx <- which(startsWith(lines, "#CHROM"))
  if (length(header_idx) == 0) stop("No #CHROM header line found in: ", path)
  
  col_names <- strsplit(sub("^#", "", lines[header_idx]), "\t")[[1]]
  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[nzchar(trimws(data_lines))]
  if (length(data_lines) == 0) stop("No variant records found in: ", path)
  
  split_fields <- strsplit(data_lines, "\t")
  n_cols <- length(col_names)
  lens <- lengths(split_fields)
  if (any(lens != n_cols)) {
    bad <- which(lens != n_cols)
    stop("Malformed VCF line(s) in ", path, " at record(s): ",
         paste(head(bad, 5), collapse = ", "))
  }
  
  mat <- matrix(unlist(split_fields), ncol = n_cols, byrow = TRUE)
  colnames(mat) <- col_names
  
  fixed_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  fix <- as.data.frame(mat[, fixed_cols, drop = FALSE], stringsAsFactors = FALSE)
  fix$POS <- as.integer(fix$POS)
  
  sample_names <- setdiff(col_names, c(fixed_cols, "FORMAT"))
  if (length(sample_names) == 0) stop("No sample columns found in: ", path)
  
  format_col <- mat[, "FORMAT"]
  geno_block <- mat[, sample_names, drop = FALSE]
  
  format_split <- strsplit(format_col, ":")
  gt_field_idx <- vapply(format_split, function(f) {
    idx <- which(f == "GT")
    if (length(idx) == 0) NA_integer_ else idx[1]
  }, integer(1))
  
  gt_only <- matrix(NA_character_, nrow = nrow(geno_block), ncol = ncol(geno_block),
                    dimnames = dimnames(geno_block))
  for (i in seq_len(nrow(geno_block))) {
    idx <- gt_field_idx[i]
    if (is.na(idx)) next
    vals <- strsplit(geno_block[i, ], ":")
    gt_only[i, ] <- vapply(vals, function(v) if (length(v) >= idx) v[idx] else NA_character_,
                           character(1))
  }
  
  list(fix = fix, gt = gt_only, samples = sample_names)
}

# ---------------------------------------------------------------------------
# 2. Restrict to mtDNA contig
# ---------------------------------------------------------------------------
filter_mito <- function(vcf_obj, mito_chr = NULL) {
  chrom <- vcf_obj$fix$CHROM
  if (is.null(mito_chr)) {
    patterns <- c("^MT$", "^chrM$", "^chrMT$", "^M$")
    keep <- Reduce(`|`, lapply(patterns, function(p) grepl(p, chrom, ignore.case = TRUE)))
    if (!any(keep)) {
      stop("Could not auto-detect a mitochondrial contig among: ",
           paste(unique(chrom), collapse = ", "),
           ". Set 'mito_chr' explicitly.")
    }
  } else {
    keep <- chrom == mito_chr
    if (!any(keep)) stop("No records found with CHROM == '", mito_chr, "'")
  }
  vcf_obj$fix <- vcf_obj$fix[keep, , drop = FALSE]
  vcf_obj$gt  <- vcf_obj$gt[keep, , drop = FALSE]
  vcf_obj
}

# ---------------------------------------------------------------------------
# 3. Standardize genotypes to a single haploid allele call
#    mtDNA is uniparental/haploid, so "0/0","0|0","0" -> "0"; "1/1" -> "1"
#    Mixed calls (e.g. "0/1", heteroplasmic) and missing calls -> NA
# ---------------------------------------------------------------------------
standardize_gt <- function(gt_mat) {
  gt_vec <- as.vector(gt_mat)
  alleles_list <- strsplit(gt_vec, "[/|]")
  call <- vapply(alleles_list, function(a) {
    if (length(a) == 0 || all(is.na(a))) return(NA_character_)
    a <- a[!is.na(a) & a != "."]
    if (length(a) == 0) return(NA_character_)
    ua <- unique(a)
    if (length(ua) == 1) ua else NA_character_   # ambiguous/heteroplasmic -> NA
  }, character(1))
  matrix(call, nrow = nrow(gt_mat), ncol = ncol(gt_mat), dimnames = dimnames(gt_mat))
}

# ---------------------------------------------------------------------------
# 4. Build a variant_id x sample call matrix for a VCF path
# ---------------------------------------------------------------------------
load_mito_calls <- function(path, mito_chr = NULL) {
  vcf <- read_vcf(path)
  vcf <- filter_mito(vcf, mito_chr)
  calls <- standardize_gt(vcf$gt)
  variant_id <- paste(vcf$fix$CHROM, vcf$fix$POS, vcf$fix$REF, vcf$fix$ALT, sep = "_")
  if (any(duplicated(variant_id))) {
    warning(sum(duplicated(variant_id)), " duplicated variant IDs in ", path,
            " -- keeping first occurrence of each.")
    keep <- !duplicated(variant_id)
    calls <- calls[keep, , drop = FALSE]
    variant_id <- variant_id[keep]
  }
  rownames(calls) <- variant_id
  list(calls = calls, samples = vcf$samples, variant_id = variant_id)
}

# ---------------------------------------------------------------------------
# 5. Matthews Correlation Coefficient (generalized multiclass form, Gorodkin 2004)
#    Reduces to the standard 2x2 MCC formula when the confusion matrix is binary.
# ---------------------------------------------------------------------------
mcc_from_confusion <- function(conf_matrix) {
  conf_matrix <- as.matrix(conf_matrix)
  storage.mode(conf_matrix) <- "double"  # avoid integer overflow on large cohorts (c_*s, s^2, etc.)
  s <- sum(conf_matrix)
  if (s == 0) return(NA_real_)
  c_ <- sum(diag(conf_matrix))
  t_k <- rowSums(conf_matrix)
  p_k <- colSums(conf_matrix)
  numerator <- c_ * s - sum(t_k * p_k)
  denominator <- sqrt((s^2 - sum(p_k^2)) * (s^2 - sum(t_k^2)))
  if (denominator == 0) return(NA_real_)
  numerator / denominator
}

mcc_from_calls <- function(truth, predicted) {
  ok <- !is.na(truth) & !is.na(predicted)
  truth <- truth[ok]; predicted <- predicted[ok]
  if (length(truth) < 2) return(list(mcc = NA_real_, n = length(truth), conf = NULL))
  lv <- sort(unique(c(truth, predicted)))
  tt <- factor(truth, levels = lv)
  pp <- factor(predicted, levels = lv)
  conf <- table(Truth = tt, Predicted = pp)
  list(mcc = mcc_from_confusion(conf), n = length(truth), conf = conf)
}

# ---------------------------------------------------------------------------
# 6. Compare one platform (array or imputed) against sequence (truth)
# ---------------------------------------------------------------------------
compare_platforms <- function(seq_obj, other_obj, label) {
  common_variants <- intersect(rownames(seq_obj$calls), rownames(other_obj$calls))
  common_samples  <- intersect(seq_obj$samples, other_obj$samples)
  
  if (length(common_variants) == 0 || length(common_samples) == 0) {
    warning("[", label, "] No overlapping variants/samples with the sequence VCF.")
    return(list(label = label, overall_mcc = NA_real_, n_pairs = 0,
                n_variants = 0, n_samples = 0, confusion = NULL,
                per_variant = data.frame()))
  }
  
  seq_sub   <- seq_obj$calls[common_variants, common_samples, drop = FALSE]
  other_sub <- other_obj$calls[common_variants, common_samples, drop = FALSE]
  
  overall <- mcc_from_calls(as.vector(seq_sub), as.vector(other_sub))
  
  per_variant <- do.call(rbind, lapply(common_variants, function(v) {
    res <- mcc_from_calls(seq_sub[v, ], other_sub[v, ])
    n_concordant <- if (res$n > 0) sum(seq_sub[v, ] == other_sub[v, ], na.rm = TRUE) else NA
    data.frame(variant_id = v, n_samples_compared = res$n,
               n_concordant = n_concordant,
               concordance_rate = ifelse(res$n > 0, n_concordant / res$n, NA),
               mcc = res$mcc, stringsAsFactors = FALSE)
  }))
  
  list(label = label,
       overall_mcc = overall$mcc,
       n_pairs = overall$n,
       n_variants = length(common_variants),
       n_samples = length(common_samples),
       confusion = overall$conf,
       per_variant = per_variant)
}

# ---------------------------------------------------------------------------
# 7. RUN
# ---------------------------------------------------------------------------
run_mtdna_concordance <- function(array_vcf_path, imputed_vcf_path, sequence_vcf_path,
                                  mito_chr = NULL, out_dir = "mtdna_concordance_results") {
  
  message("Loading sequence (truth) VCF ...")
  seq_obj <- load_mito_calls(sequence_vcf_path, mito_chr)
  message("Loading array VCF ...")
  array_obj <- load_mito_calls(array_vcf_path, mito_chr)
  message("Loading imputed VCF ...")
  imputed_obj <- load_mito_calls(imputed_vcf_path, mito_chr)
  
  message("Comparing array vs sequence ...")
  res_array <- compare_platforms(seq_obj, array_obj, "Array_vs_Sequence")
  message("Comparing imputed vs sequence ...")
  res_imputed <- compare_platforms(seq_obj, imputed_obj, "Imputed_vs_Sequence")
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  summary_tbl <- data.frame(
    comparison        = c(res_array$label, res_imputed$label),
    n_variants        = c(res_array$n_variants, res_imputed$n_variants),
    n_samples         = c(res_array$n_samples, res_imputed$n_samples),
    n_genotype_pairs  = c(res_array$n_pairs, res_imputed$n_pairs),
    overall_mcc       = c(res_array$overall_mcc, res_imputed$overall_mcc),
    stringsAsFactors  = FALSE
  )
  
  write.csv(summary_tbl, file.path(out_dir, "overall_mcc_summary.csv"), row.names = FALSE)
  write.csv(res_array$per_variant, file.path(out_dir, "per_variant_mcc_array_vs_sequence.csv"),
            row.names = FALSE)
  write.csv(res_imputed$per_variant, file.path(out_dir, "per_variant_mcc_imputed_vs_sequence.csv"),
            row.names = FALSE)
  
  png(file.path(out_dir, "overall_mcc_comparison.png"), width = 800, height = 600, res = 120)
  bp <- barplot(summary_tbl$overall_mcc, names.arg = summary_tbl$comparison,
                ylim = c(-1, 1), col = c("#4C72B0", "#DD8452"),
                ylab = "Matthews Correlation Coefficient",
                main = "mtDNA genotype concordance vs sequence data")
  abline(h = 0, lty = 2, col = "grey40")
  is_na_mcc <- is.na(summary_tbl$overall_mcc)
  label_vals <- ifelse(is_na_mcc, "NA", round(summary_tbl$overall_mcc, 3))
  label_y    <- ifelse(is_na_mcc, 0, summary_tbl$overall_mcc)
  label_pos  <- ifelse(is_na_mcc | summary_tbl$overall_mcc >= 0, 3, 1)
  text(bp, label_y, labels = label_vals, pos = label_pos)
  dev.off()
  
  cat("\n================ mtDNA genotype concordance (MCC) ================\n")
  print(summary_tbl, row.names = FALSE)
  cat("=====================================================================\n")
  cat("\nConfusion matrix, Array vs Sequence:\n"); print(res_array$confusion)
  cat("\nConfusion matrix, Imputed vs Sequence:\n"); print(res_imputed$confusion)
  cat("\nResults written to: ", normalizePath(out_dir), "\n", sep = "")
  
  invisible(list(summary = summary_tbl, array = res_array, imputed = res_imputed))
}

# ---------------------------------------------------------------------------
# Execute when run as a script (Rscript mtdna_concordance.R)
# ---------------------------------------------------------------------------
if (sys.nframe() == 0) {
  run_mtdna_concordance(array_vcf_path, imputed_vcf_path, sequence_vcf_path,
                        mito_chr = mito_chr, out_dir = out_dir)
}

