# =================================================================
# mtDNA HAPLOGROUP concordance: exact-match proportion + Hierarchical
# F-measure (PhyloTree Build 17)
#
# Replaces MCC for haplogroup-level evaluation (kept separately for
# genotype-level concordance, see mtDNA_genotype_concordance_MCC.R).
# Reports two complementary metrics per platform/method:
#   - Proportion of exact matches (with Wilson score CI, and a
#     paired McNemar test + paired bootstrap for Imputed vs Array)
#   - Hierarchical Precision / Recall / F-measure (tree-aware credit
#     for under-resolved / same-branch calls), with paired bootstrap
#
# AUTO-DETECTS cohort structure from the input CSV's column names:
#   - Single-array cohort : Sample_ID, True_Haplo, Array, Imputed
#   - Multi-platform cohort: Sample_ID, True_Haplo,
#         Array_<Platform>, Imputed_<Platform>  (repeated per platform)
#
# Requires 'phylotree17_parents.csv' (PhyloTree Build 17 parent
# table, generated previously) in the working directory.
# =================================================================

## 0. Setup -----------------------------------------------------------
required_pkgs <- c("data.table", "ggplot2", "boot")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(data.table)
library(ggplot2)
library(boot)

## 1. Load data & auto-detect cohort structure ---------------------------
infile <- "1000G_haplogroups.csv"   # <-- set to your actual file path
dt <- fread(infile, na.strings = c("", "NA", "N/A", "."))
stopifnot("True_Haplo" %in% names(dt))

detect_platform_map <- function(dt) {
  cn <- names(dt)
  if (all(c("Array", "Imputed") %in% cn)) {
    return(data.table(Platform = "Overall", Array_col = "Array", Imputed_col = "Imputed"))
  }
  array_cols <- grep("^Array_", cn, value = TRUE)
  imputed_cols <- grep("^Imputed_", cn, value = TRUE)
  platforms_arr <- sub("^Array_", "", array_cols)
  platforms_imp <- sub("^Imputed_", "", imputed_cols)
  common <- intersect(platforms_arr, platforms_imp)
  if (length(common) == 0) {
    stop("Could not detect a valid Array/Imputed column structure. Expected either ",
         "'Array' & 'Imputed' columns, or 'Array_<Platform>' & 'Imputed_<Platform>' pairs.")
  }
  data.table(Platform = common, Array_col = paste0("Array_", common), Imputed_col = paste0("Imputed_", common))
}

platform_map <- detect_platform_map(dt)
is_single_platform <- identical(platform_map$Platform, "Overall")
cat(sprintf("Detected cohort structure: %s (%d platform%s)\n",
            ifelse(is_single_platform, "single-array", "multi-platform"),
            nrow(platform_map), ifelse(nrow(platform_map) == 1, "", "s")))
print(platform_map)

## 2. PhyloTree hierarchical machinery -------------------------------------
tree_tbl <- fread("../../phylotree17_parents.csv", na.strings = "")
parent_map <- setNames(tree_tbl$parent, tree_tbl$haplogroup)
valid_names <- tree_tbl$haplogroup
ROOT <- "mt-MRCA"

match_haplogroup <- function(x, valid = valid_names) {
  if (is.na(x) || x == "") return(NA_character_)
  x <- trimws(x)
  if (x %in% valid) return(x)
  x2 <- sub("[!?].*$", "", x)
  x2 <- sub("\\s*\\(.*\\)\\s*$", "", x2)
  x2 <- trimws(x2)
  if (x2 %in% valid) return(x2)
  cand <- x2
  while (nchar(cand) > 0) {
    if (cand %in% valid) return(cand)
    cand <- substr(cand, 1, nchar(cand) - 1)
  }
  NA_character_
}

ancestor_cache <- new.env()
unmatched_labels <- character(0)

get_ancestors <- function(hg_raw) {
  key <- as.character(hg_raw)
  if (exists(key, envir = ancestor_cache, inherits = FALSE)) {
    return(get(key, envir = ancestor_cache, inherits = FALSE))
  }
  hg <- match_haplogroup(key)
  if (is.na(hg)) {
    unmatched_labels[[length(unmatched_labels) + 1]] <<- key
    assign(key, NA_character_, envir = ancestor_cache)
    return(NA_character_)
  }
  chain <- character(0)
  cur <- hg
  repeat {
    chain <- c(chain, cur)
    nxt <- parent_map[[cur]]
    if (is.null(nxt) || is.na(nxt) || nxt == "") break
    cur <- nxt
  }
  chain <- chain[chain != ROOT]
  assign(key, chain, envir = ancestor_cache)
  chain
}

hier_prf <- function(true_hg, pred_hg) {
  anc_t <- get_ancestors(true_hg)
  anc_p <- get_ancestors(pred_hg)
  if (length(anc_t) == 1 && is.na(anc_t)) return(c(hP = NA_real_, hR = NA_real_, hF = NA_real_))
  if (length(anc_p) == 1 && is.na(anc_p)) return(c(hP = NA_real_, hR = NA_real_, hF = NA_real_))
  if (length(anc_t) == 0 || length(anc_p) == 0) return(c(hP = NA_real_, hR = NA_real_, hF = NA_real_))
  shared <- length(intersect(anc_t, anc_p))
  hP <- shared / length(anc_p)
  hR <- shared / length(anc_t)
  hF <- if (hP + hR == 0) 0 else 2 * hP * hR / (hP + hR)
  c(hP = hP, hR = hR, hF = hF)
}

# Light normalization for EXACT match (whitespace / uncertainty-annotation only -
# no truncation, since that would no longer be an "exact" match)
normalize_for_exact <- function(x) {
  if (is.na(x) || x == "") return(NA_character_)
  x <- trimws(x)
  x <- sub("[!?].*$", "", x)
  x <- sub("\\s*\\(.*\\)\\s*$", "", x)
  trimws(x)
}

## 3. Wilson score CI for a proportion -------------------------------------
wilson_ci <- function(x, n, conf = 0.95) {
  if (n == 0) return(c(low = NA_real_, high = NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  phat <- x / n
  denom <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denom
  half <- (z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2))) / denom
  c(low = max(0, center - half), high = min(1, center + half))
}

## 4. Per-sample scores: exact match + hierarchical P/R/F -------------------
per_sample <- data.table()
for (i in seq_len(nrow(platform_map))) {
  p <- platform_map$Platform[i]
  for (method in c("Array", "Imputed")) {
    col <- if (method == "Array") platform_map$Array_col[i] else platform_map$Imputed_col[i]
    if (!col %in% names(dt)) next
    sub <- dt[!is.na(True_Haplo) & !is.na(get(col)),
              .(Sample_ID, True_Haplo, Call = get(col))]
    if (nrow(sub) == 0) next
    
    sub[, ExactMatch := mapply(function(a, b) {
      na <- normalize_for_exact(a); nb <- normalize_for_exact(b)
      if (is.na(na) || is.na(nb)) NA else (na == nb)
    }, True_Haplo, Call)]
    
    scores <- t(mapply(hier_prf, sub$True_Haplo, sub$Call))
    sub[, `:=`(hP = scores[, "hP"], hR = scores[, "hR"], hF = scores[, "hF"],
               Platform = p, Method = method)]
    per_sample <- rbind(per_sample, sub, fill = TRUE)
  }
}

n_unmatched <- sum(is.na(per_sample$hF))
if (n_unmatched > 0) {
  warning(sprintf(
    "%d sample-comparisons had a haplogroup label not found in PhyloTree Build 17 and were excluded from hF. See 'unmatched_haplogroup_labels.csv'.",
    n_unmatched))
  fwrite(data.table(label = unique(unmatched_labels)), "unmatched_haplogroup_labels.csv")
}
fwrite(per_sample, "haplogroup_per_sample_scores.csv")

## 5. Summary per platform x method -----------------------------------------
set.seed(123)
n_boot <- 1000

boot_mean_ci <- function(x, R = n_boot) {
  b <- boot(data = x, statistic = function(d, i) mean(d[i]), R = R)
  ci <- tryCatch(boot.ci(b, type = "perc")$percent[4:5], error = function(e) c(NA, NA))
  c(CI_low = ci[1], CI_high = ci[2])
}

summary_dt <- per_sample[, {
  em_valid <- ExactMatch[!is.na(ExactMatch)]
  n_em <- length(em_valid)
  x_em <- sum(em_valid)
  wci <- wilson_ci(x_em, n_em)
  
  hf_valid <- hF[!is.na(hF)]
  hci <- if (length(hf_valid) >= 2) boot_mean_ci(hf_valid) else c(CI_low = NA, CI_high = NA)
  
  list(N_EM = n_em, ExactMatch_Prop = x_em / n_em,
       EM_CI_low = wci["low"], EM_CI_high = wci["high"],
       N_hF = length(hf_valid), hP = mean(hP, na.rm = TRUE), hR = mean(hR, na.rm = TRUE),
       hF = mean(hf_valid), hF_CI_low = hci["CI_low"], hF_CI_high = hci["CI_high"])
}, by = .(Platform, Method)]

cat("\n=== Haplogroup concordance summary (exact match + hierarchical F) ===\n")
print(summary_dt)
fwrite(summary_dt, "haplogroup_concordance_summary.csv")

## 6. Paired comparison per platform: Imputed vs Array ------------------------
diff_dt <- data.table()
for (p in platform_map$Platform) {
  sub_p <- per_sample[Platform == p]
  
  ## --- exact match: paired bootstrap + McNemar's test ---
  wide_em <- dcast(sub_p, Sample_ID + True_Haplo ~ Method, value.var = "ExactMatch")
  wide_em <- wide_em[!is.na(Array) & !is.na(Imputed)]
  n_pair <- nrow(wide_em)
  
  em_diff <- NA_real_; em_ci <- c(NA, NA); em_boot_p <- NA_real_; mcnemar_p <- NA_real_
  if (n_pair >= 2) {
    diff_fn <- function(d, i) mean(d$Imputed[i]) - mean(d$Array[i])
    em_diff <- diff_fn(wide_em, seq_len(n_pair))
    b <- boot(data = wide_em, statistic = diff_fn, R = 2000)
    em_ci <- tryCatch(boot.ci(b, type = "perc")$percent[4:5], error = function(e) c(NA, NA))
    em_boot_p <- min(1, 2 * min(mean(b$t <= 0), mean(b$t >= 0)))
    
    tbl <- table(Array = wide_em$Array, Imputed = wide_em$Imputed)
    # ensure 2x2 with both levels present
    tbl <- tbl[c("FALSE", "TRUE")[c("FALSE", "TRUE") %in% rownames(tbl)],
               c("FALSE", "TRUE")[c("FALSE", "TRUE") %in% colnames(tbl)], drop = FALSE]
    if (all(dim(tbl) == c(2, 2))) {
      mcnemar_p <- tryCatch(mcnemar.test(tbl, correct = TRUE)$p.value, error = function(e) NA_real_)
    }
  }
  
  ## --- hierarchical F: paired bootstrap ---
  wide_hf <- dcast(sub_p, Sample_ID + True_Haplo ~ Method, value.var = "hF")
  wide_hf <- wide_hf[!is.na(Array) & !is.na(Imputed)]
  n_hf <- nrow(wide_hf)
  
  hf_diff <- NA_real_; hf_ci <- c(NA, NA); hf_p <- NA_real_
  if (n_hf >= 2) {
    diff_fn2 <- function(d, i) mean(d$Imputed[i]) - mean(d$Array[i])
    hf_diff <- diff_fn2(wide_hf, seq_len(n_hf))
    b2 <- boot(data = wide_hf, statistic = diff_fn2, R = 2000)
    hf_ci <- tryCatch(boot.ci(b2, type = "perc")$percent[4:5], error = function(e) c(NA, NA))
    hf_p <- min(1, 2 * min(mean(b2$t <= 0), mean(b2$t >= 0)))
  }
  
  diff_dt <- rbind(diff_dt, data.table(
    Platform = p,
    N_EM_paired = n_pair, EM_Diff = em_diff, EM_CI_low = em_ci[1], EM_CI_high = em_ci[2],
    EM_boot_p = em_boot_p, EM_McNemar_p = mcnemar_p,
    N_hF_paired = n_hf, hF_Diff = hf_diff, hF_CI_low = hf_ci[1], hF_CI_high = hf_ci[2], hF_p = hf_p
  ))
}

cat("\n=== Imputed - Array: paired differences per platform ===\n")
print(diff_dt)
fwrite(diff_dt, "haplogroup_concordance_diff.csv")

## 7. Combined plot: Exact-match proportion + Hierarchical F ------------------
combined_long <- rbind(
  summary_dt[, .(Platform, Method, N = N_EM, Value = ExactMatch_Prop,
                 CI_low = EM_CI_low, CI_high = EM_CI_high, Metric = "Exact match")],
  summary_dt[, .(Platform, Method, N = N_hF, Value = hF,
                 CI_low = hF_CI_low, CI_high = hF_CI_high, Metric = "Hierarchical F")]
)
combined_long[, Metric := factor(Metric, levels = c("Exact match", "Hierarchical F"))]

p_bars <- ggplot(combined_long, aes(x = Platform, y = Value, fill = Method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position = position_dodge(width = 0.7), width = 0.2) +
  facet_wrap(~Metric, ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(title = "Array vs Imputed mtDNA Haplogroup Concordance",
       subtitle = "Exact-match proportion vs tree-aware Hierarchical F-metric (PhyloTree Build 17)",
       y = "Score", x = if (is_single_platform) NULL else "Genotyping Platform") +
  scale_fill_manual(values = c(Array = "#C2185B", Imputed = "#00796B")) +
  theme(axis.text.x = element_text(angle = if (is_single_platform) 0 else 30, hjust = if (is_single_platform) 0.5 else 1),
        plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"))

plot_width <- if (is_single_platform) 7 else 9
plot_height <- if (is_single_platform) 8 else 9

p_bars

ggsave("haplogroup_concordance_barplot.png", p_bars, width = plot_width, height = plot_height, dpi = 300)

## 8. Markdown report -------------------------------------------------------------
fmt_metric <- function(diff, lo, hi, p) {
  if (is.na(diff)) return("N/A")
  sig_flag <- if (!is.na(p) && p < 0.05) " **(significant)**" else ""
  sprintf("%+.3f [%.3f, %.3f], p=%.3f%s", diff, lo, hi, p, sig_flag)
}

report_lines <- c(
  "# mtDNA Haplogroup Concordance Report",
  "",
  "Array- vs Imputation-derived mtDNA haplogroup calls compared against ground truth, using:",
  "- **Proportion of exact matches** (Wilson score 95% CI; paired McNemar's test and paired",
  "  bootstrap for the Imputed-vs-Array difference)",
  "- **Hierarchical Precision/Recall/F-measure** (Kiritchenko et al. 2006) against the true",
  "  PhyloTree Build 17 topology, crediting under-resolved or same-branch calls proportionally",
  "  to shared ancestry rather than scoring every mismatch as equally wrong",
  ""
)

for (p in platform_map$Platform) {
  s <- summary_dt[Platform == p]
  d <- diff_dt[Platform == p]
  s_arr <- s[Method == "Array"]; s_imp <- s[Method == "Imputed"]
  report_lines <- c(report_lines,
                    sprintf("## %s", p),
                    sprintf("- Exact match: Array = %.3f, Imputed = %.3f, Diff = %s (McNemar p = %s)",
                            s_arr$ExactMatch_Prop, s_imp$ExactMatch_Prop,
                            fmt_metric(d$EM_Diff, d$EM_CI_low, d$EM_CI_high, d$EM_boot_p),
                            ifelse(is.na(d$EM_McNemar_p), "N/A", sprintf("%.3f", d$EM_McNemar_p))),
                    sprintf("- Hierarchical F: Array = %.3f, Imputed = %.3f, Diff = %s",
                            s_arr$hF, s_imp$hF, fmt_metric(d$hF_Diff, d$hF_CI_low, d$hF_CI_high, d$hF_p)),
                    ""
  )
}
writeLines(report_lines, "haplogroup_concordance_report.md")

cat("\nDone. Outputs written:\n",
    " - haplogroup_per_sample_scores.csv\n",
    " - haplogroup_concordance_summary.csv, haplogroup_concordance_diff.csv\n",
    " - haplogroup_concordance_barplot.png, haplogroup_concordance_diff_plot.png\n",
    " - haplogroup_concordance_report.md\n")
