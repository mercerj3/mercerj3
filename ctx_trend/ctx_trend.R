#!/usr/bin/env Rscript
#
# C-Telopeptide (CTX) trend by PIN and visit.
#
# Reads a Labcorp study-results export in the "Study Data Entry Template" layout
# and pulls the serum C-Telopeptide column into a tidy long table keyed by
# participant (PIN) and visit, then writes a wide PIN x Visit table and a trend
# plot.
#
# Layout assumed (matches Study_Data_Entry_Template_..._Labcorp_Results.csv):
#     row 1  column headers
#     row 2  units             (C-Telopeptide,Serum -> pg/mL)
#     row 3  reference intervals
#     row 4+ one row per participant / visit / draw
#
# Usage:
#     Rscript ctx_trend.R --csv path/to/Labcorp_Results.csv --outdir out
#     Rscript ctx_trend.R --csv path/to/Labcorp_Results.csv --x days

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

CTX_COL <- "C-Telopeptide,Serum"
UNITS_ROW <- 1L # position within the two non-data header rows
REF_ROW <- 2L

#' Read the template CSV.
#'
#' The two metadata rows under the header are stripped off and returned as
#' named character vectors so the CTX units/reference range can be used in
#' labels.
#'
#' @return list(data = tibble, units = named chr, refs = named chr)
load_labs <- function(csv_path) {
  raw <- read_csv(csv_path, col_types = cols(.default = col_character()),
                  name_repair = "minimal", progress = FALSE)

  units <- unlist(raw[UNITS_ROW, ], use.names = TRUE)
  refs <- unlist(raw[REF_ROW, ], use.names = TRUE)
  data <- raw[-seq_len(REF_ROW), ]

  # Drop the blank filler rows at the bottom of the template.
  data <- data[!is.na(data$PIN) & trimws(data$PIN) != "", ]

  list(data = data, units = units, refs = refs)
}

#' Tidy one-row-per-draw CTX table with baseline-relative columns.
#'
#' Participants 001-006 record the pre/post draw flag in `Pre_Post_Old` and
#' 007+ in `Pre_Post_New`; the two are coalesced into a single `pre_post`.
#' A visit can appear twice for the same PIN (a pre and a post draw), so the
#' natural key is (pin, visit, pre_post).
ctx_by_pin_visit <- function(data) {
  blank_to_na <- function(x) dplyr::na_if(trimws(x), "")

  ctx <- tibble(
    pin = trimws(data$PIN),
    visit = suppressWarnings(as.integer(blank_to_na(data$Visit_Number))),
    pre_post = suppressWarnings(as.integer(coalesce(
      blank_to_na(data$Pre_Post_New),
      blank_to_na(data$Pre_Post_Old)
    ))),
    collection_date = as.Date(blank_to_na(data$Collection_Date)),
    ctx_pg_ml = suppressWarnings(as.numeric(blank_to_na(data[[CTX_COL]])))
  )

  ctx %>%
    # Keep only draws that actually have a CTX result.
    filter(!is.na(ctx_pg_ml)) %>%
    arrange(pin, visit, pre_post) %>%
    mutate(
      draw = case_when(
        pre_post == 0 ~ "pre",
        pre_post == 1 ~ "post",
        TRUE ~ "single"
      ),
      visit_label = ifelse(draw == "single",
                           paste0("V", visit),
                           paste0("V", visit, " (", draw, ")"))
    ) %>%
    # Baseline = each participant's earliest CTX draw.
    group_by(pin) %>%
    mutate(
      baseline_pg_ml = first(ctx_pg_ml),
      change_from_baseline = ctx_pg_ml - baseline_pg_ml,
      pct_change_from_baseline = 100 * change_from_baseline / baseline_pg_ml,
      days_from_baseline = as.integer(collection_date - first(collection_date))
    ) %>%
    ungroup()
}

#' Visit labels ordered by visit number, then pre before post.
visit_label_levels <- function(ctx) {
  ctx %>%
    distinct(visit_label, visit, pre_post) %>%
    arrange(visit, pre_post) %>%
    pull(visit_label)
}

#' PIN x visit matrix of CTX values (rows = PIN, columns = visit label).
ctx_wide <- function(ctx) {
  ctx %>%
    select(pin, visit_label, ctx_pg_ml) %>%
    pivot_wider(names_from = visit_label, values_from = ctx_pg_ml,
                values_fn = first) %>%
    arrange(pin) %>%
    select(pin, any_of(visit_label_levels(ctx)))
}

#' Line plot of CTX over time, one series per PIN.
plot_ctx_trend <- function(ctx, out_path, x = "visit", units = "pg/mL",
                           ref_range = "") {
  x_col <- switch(x, visit = "visit", days = "days_from_baseline")
  x_label <- switch(x, visit = "Visit number", days = "Days from baseline draw")

  p <- ggplot(ctx, aes(x = .data[[x_col]], y = ctx_pg_ml,
                       colour = pin, group = pin))

  # Shade the reference interval when the export supplies one for CTX.
  if (nzchar(ref_range) && grepl("-", ref_range, fixed = TRUE)) {
    bounds <- as.numeric(strsplit(ref_range, "-", fixed = TRUE)[[1]])
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = bounds[1], ymax = bounds[2],
                      fill = "grey85", alpha = 0.6)
  }

  p <- p +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.4) +
    labs(
      title = "C-Telopeptide trend by participant and visit",
      x = x_label,
      y = paste0("C-Telopeptide, serum (", units, ")"),
      colour = "Participant"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  if (x == "visit") {
    p <- p + scale_x_continuous(breaks = sort(unique(ctx$visit)))
  }

  ggsave(out_path, p, width = 9, height = 5.5, dpi = 150)
  out_path
}

#' Minimal --flag value parser so the script needs only the tidyverse.
parse_args <- function(args) {
  opts <- list(csv = NULL, outdir = "out", x = "visit")
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!key %in% names(opts) || i == length(args)) {
      stop("usage: Rscript ctx_trend.R --csv FILE [--outdir DIR] [--x visit|days]",
           call. = FALSE)
    }
    opts[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  if (is.null(opts$csv)) stop("--csv is required", call. = FALSE)
  if (!opts$x %in% c("visit", "days")) stop("--x must be 'visit' or 'days'", call. = FALSE)
  opts
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

  labs <- load_labs(opts$csv)
  ctx <- ctx_by_pin_visit(labs$data)
  wide <- ctx_wide(ctx)

  long_path <- file.path(opts$outdir, "ctx_by_pin_visit.csv")
  wide_path <- file.path(opts$outdir, "ctx_pin_by_visit_wide.csv")
  # Blank rather than "NA" for missing cells, so the CSVs import cleanly elsewhere.
  write_csv(ctx, long_path, na = "")
  write_csv(wide, wide_path, na = "")

  plot_path <- plot_ctx_trend(
    ctx,
    file.path(opts$outdir, "ctx_trend.png"),
    x = opts$x,
    units = if (is.na(labs$units[[CTX_COL]])) "pg/mL" else labs$units[[CTX_COL]],
    ref_range = if (is.na(labs$refs[[CTX_COL]])) "" else labs$refs[[CTX_COL]]
  )

  cat(sprintf("CTX results: %d draws across %d participants\n\n",
              nrow(ctx), n_distinct(ctx$pin)))
  print(as.data.frame(ctx[, c("pin", "visit", "draw", "collection_date",
                              "ctx_pg_ml", "pct_change_from_baseline")]),
        row.names = FALSE)
  cat("\nPIN x visit matrix (pg/mL):\n")
  print(as.data.frame(wide), row.names = FALSE)
  cat(sprintf("\nWrote %s\n      %s\n      %s\n", long_path, wide_path, plot_path))
}

if (sys.nframe() == 0L) {
  main()
}
