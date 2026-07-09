#!/usr/bin/env Rscript
# ============================================================================
# Venn diagrams of imputed vs. sequenced mtDNA variant overlap, per platform
#
# Expects a CSV with (at minimum) these four pieces of information per row:
#   - platform name            (e.g. "H3Africa")
#   - total imputed variants
#   - total sequenced variants
#   - number of variants common to both
#
# Column names are matched flexibly (case-insensitive, partial match on
# key words) -- see `match_column()` below. Edit `col_*` patterns if your
# headers differ substantially from the defaults.
#
# Requires: VennDiagram, gridExtra  (install.packages(c("VennDiagram","gridExtra")))
# ============================================================================

if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  stop("Package 'VennDiagram' is required. Install with: install.packages('VennDiagram')")
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package 'gridExtra' is required. Install with: install.packages('gridExtra')")
}
suppressMessages(library(VennDiagram))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))

# ---------------------------------------------------------------------------
# 0. USER SETTINGS
# ---------------------------------------------------------------------------
input_csv <- "bcftools_stats.csv"
out_dir   <- "venn_diagrams"

# ---------------------------------------------------------------------------
# 1. Flexible column matching, so small header naming differences don't break
# ---------------------------------------------------------------------------
match_column <- function(df, patterns) {
  hit <- which(sapply(tolower(names(df)), function(nm) {
    any(vapply(patterns, function(p) grepl(p, nm), logical(1)))
  }))
  if (length(hit) == 0) return(NA_character_)
  names(df)[hit[1]]
}

load_counts <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  
  col_platform <- match_column(df, c("platform", "array", "chip", "name"))
  col_imputed  <- match_column(df, c("imput"))
  col_seq      <- match_column(df, c("sequenc", "^seq"))
  col_common   <- match_column(df, c("common", "overlap", "shared", "intersect"))
  
  missing <- c(platform = col_platform, imputed = col_imputed,
               sequenced = col_seq, common = col_common)
  if (any(is.na(missing))) {
    stop("Could not identify these required column(s) in '", path, "': ",
         paste(names(missing)[is.na(missing)], collapse = ", "),
         ". Found columns: ", paste(names(df), collapse = ", "),
         ". Rename your headers or edit match_column() patterns.")
  }
  
  out <- data.frame(
    platform  = as.character(df[[col_platform]]),
    n_imputed = as.numeric(df[[col_imputed]]),
    n_seq     = as.numeric(df[[col_seq]]),
    n_common  = as.numeric(df[[col_common]]),
    stringsAsFactors = FALSE
  )
  
  # sanity checks
  bad <- out$n_common > pmin(out$n_imputed, out$n_seq)
  if (any(bad)) {
    warning("Row(s) where 'common' exceeds imputed and/or sequenced total (impossible for a subset ",
            "count) -- check: ", paste(out$platform[bad], collapse = ", "))
  }
  out
}

# ---------------------------------------------------------------------------
# 2. Draw one pairwise (two-set) Venn diagram for a single platform
# ---------------------------------------------------------------------------
make_venn_grob <- function(platform, n_imputed, n_seq, n_common,
                           fill = c("#00796B", "#283593")) {
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  grob <- VennDiagram::draw.pairwise.venn(
    area1 = n_imputed,
    area2 = n_seq,
    cross.area = n_common,
    category = c("Imputed", "Sequenced"),
    fill = fill,
    alpha = 0.75,
    lty = "blank",
    fontfamily = "sans",
    cex = 1.3,
    cat.cex = 1.1,
    cat.pos = c(-20, 20),
    cat.fontfamily = "sans",
    cat.dist = c(0.055, 0.055),
    ind = FALSE,
    scaled = TRUE
  )
  gTree(children = grob, vp = viewport(width = 0.9, height = 0.85))
}

# ---------------------------------------------------------------------------
# 3. RUN: one PNG per platform + one combined multi-panel PNG
# ---------------------------------------------------------------------------
run_venn_diagrams <- function(input_csv, out_dir = "venn_diagrams") {
  counts <- load_counts(input_csv)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  grobs <- vector("list", nrow(counts))
  names(grobs) <- counts$platform
  failed <- character(0)
  
  for (i in seq_len(nrow(counts))) {
    row <- counts[i, ]
    message("Drawing Venn diagram for: ", row$platform)
    
    g <- tryCatch(
      make_venn_grob(row$platform, row$n_imputed, row$n_seq, row$n_common),
      error = function(e) {
        warning("Skipping '", row$platform, "': ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(g)) { failed <- c(failed, row$platform); next }
    grobs[[i]] <- g
    
    safe_name <- gsub("[^A-Za-z0-9_-]+", "_", row$platform)
    png(file.path(out_dir, paste0("venn_", safe_name, ".png")),
        width = 1400, height = 1400, res = 220)
    grid.newpage()
    grid.draw(g)
    grid.text(row$platform, y = unit(0.97, "npc"), gp = gpar(fontsize = 16, fontface = "bold"))
    dev.off()
  }
  
  grobs <- grobs[!vapply(grobs, is.null, logical(1))]
  counts_ok <- counts[!(counts$platform %in% failed), , drop = FALSE]
  
  # combined multi-panel figure, arranged in a grid with up to 3 columns
  n <- length(grobs)
  if (n == 0) {
    warning("No Venn diagrams could be drawn -- check your input counts.")
    return(invisible(list(counts = counts, combined_path = NA_character_)))
  }
  ncol_panel <- min(3, n)
  nrow_panel <- ceiling(n / ncol_panel)
  
  labelled_grobs <- lapply(seq_len(n), function(i) {
    gridExtra::arrangeGrob(grobs[[i]], top = counts_ok$platform[i])
  })
  
  combined_path <- file.path(out_dir, "venn_all_platforms_combined.png")
  png(combined_path, width = 500 * ncol_panel, height = 550 * nrow_panel, res = 150)
  gridExtra::grid.arrange(grobs = labelled_grobs, ncol = ncol_panel)
  dev.off()
  
  if (length(failed) > 0) {
    message("Skipped ", length(failed), " platform(s) with invalid counts: ",
            paste(failed, collapse = ", "))
  }
  message("\nSaved ", n, " individual Venn diagram(s) and one combined figure to: ",
          normalizePath(out_dir))
  invisible(list(counts = counts, combined_path = combined_path))
}

# ---------------------------------------------------------------------------
# Execute when run as a script (Rscript venn_diagrams.R)
# ---------------------------------------------------------------------------
if (sys.nframe() == 0) {
  run_venn_diagrams(input_csv, out_dir)
}

