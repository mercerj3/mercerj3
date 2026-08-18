#!/usr/bin/env Rscript
#
# C-Telopeptide (CTX) trend by PIN and visit.
#
# Works on a data frame already loaded in your session — named `study` below —
# holding a Labcorp study-results export in the "Study Data Entry Template"
# layout. Pulls the serum C-Telopeptide column into a tidy long table keyed by
# participant (PIN) and visit, then builds a wide PIN x Visit table and a trend
# plot.
#
# Layout assumed (matches Study_Data_Entry_Template_..._Labcorp_Results.csv):
#     row 1  units             (C-Telopeptide,Serum -> pg/mL)
#     row 2  reference intervals
#     row 3+ one row per participant / visit / draw
#
# Those first two rows sit under the header, so they arrive as data rows in
# `study` no matter how it was read; `ctx_by_pin_visit()` drops them, along with
# the blank filler rows at the bottom of the template.
#
# Usage:
#     source("ctx_trend.R")
#     ctx  <- ctx_by_pin_visit(study)
#     wide <- ctx_wide(ctx)
#     plot_ctx_trend(ctx)                 # x = visit number
#     plot_ctx_trend(ctx, x = "days")     # x = days from baseline draw

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

#' Locate a column by name, tolerating however `study` was read in.
#'
#' `read.csv()` turns "C-Telopeptide,Serum" into "C.Telopeptide.Serum" and
#' `janitor::clean_names()` into "c_telopeptide_serum", while `readr::read_csv()`
#' keeps it verbatim. Comparing on alphanumerics only makes all of these match.
find_col <- function(df, name) {
  norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
  hit <- which(norm(names(df)) == norm(name))
  if (length(hit) == 0) {
    stop(sprintf("column '%s' not found in the data frame", name), call. = FALSE)
  }
  names(df)[hit[1]]
}

#' Pull a column as trimmed character, blanks as NA.
#'
#' Columns may already be typed (numeric/Date) if the reader guessed, so
#' everything is routed through character first for one coercion path.
col_chr <- function(df, name) {
  x <- trimws(as.character(df[[find_col(df, name)]]))
  dplyr::na_if(x, "")
}

#' Tidy one-row-per-draw CTX table with baseline-relative columns.
#'
#' Participants 001-006 record the pre/post draw flag in `Pre_Post_Old` and
#' 007+ in `Pre_Post_New`; the two are coalesced into a single `pre_post`.
#' A visit can appear twice for the same PIN (a pre and a post draw), so the
#' natural key is (pin, visit, pre_post).
#'
#' @param study data frame of the export, header row already used as names
ctx_by_pin_visit <- function(study) {
  # Drop the units / reference-interval rows and the blank filler rows first, so
  # the header text ("yyyy-mm-dd", "pg/mL") never reaches the type coercions.
  pin <- col_chr(study, "PIN")
  study <- study[!is.na(pin) & !grepl("^(ParticipantID|ReferenceIntervals)", pin), ]

  ctx <- tibble(
    pin = col_chr(study, "PIN"),
    visit = suppressWarnings(as.integer(col_chr(study, "Visit_Number"))),
    pre_post = suppressWarnings(as.integer(coalesce(
      col_chr(study, "Pre_Post_New"),
      col_chr(study, "Pre_Post_Old")
    ))),
    collection_date = as.Date(col_chr(study, "Collection_Date"), format = "%Y-%m-%d"),
    ctx_pg_ml = suppressWarnings(as.numeric(col_chr(study, "C-Telopeptide,Serum")))
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

#' Units / reference interval for a lab column, read off the template's row 1-2.
#'
#' Returns "" when the export leaves the cell blank, as it does for CTX.
lab_metadata <- function(study, name = "C-Telopeptide,Serum") {
  pin <- col_chr(study, "PIN")
  value_at <- function(prefix) {
    row <- which(grepl(prefix, pin))[1]
    if (is.na(row)) return("")
    val <- as.character(study[[find_col(study, name)]][row])
    if (is.na(val)) "" else trimws(val)
  }
  list(units = value_at("^ParticipantID"), ref_range = value_at("^ReferenceIntervals"))
}

#' Line plot of CTX over time, one series per PIN.
plot_ctx_trend <- function(ctx, x = c("visit", "days"), units = "pg/mL",
                           ref_range = "") {
  x <- match.arg(x)
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

  p
}

# ---------------------------------------------------------------------------
# Run on `study`
# ---------------------------------------------------------------------------
# Guarded so `source("ctx_trend.R")` just defines the functions when `study`
# is not (yet) in the session.

if (exists("study")) {
  ctx <- ctx_by_pin_visit(study)
  wide <- ctx_wide(ctx)
  meta <- lab_metadata(study)

  ctx_plot <- plot_ctx_trend(
    ctx,
    x = "visit",
    units = if (nzchar(meta$units)) meta$units else "pg/mL",
    ref_range = meta$ref_range
  )

  cat(sprintf("CTX results: %d draws across %d participants\n\n",
              nrow(ctx), n_distinct(ctx$pin)))
  print(as.data.frame(ctx[, c("pin", "visit", "draw", "collection_date",
                              "ctx_pg_ml", "pct_change_from_baseline")]),
        row.names = FALSE)
  cat("\nPIN x visit matrix (pg/mL):\n")
  print(as.data.frame(wide), row.names = FALSE)

  print(ctx_plot)

  # To save:
  # write.csv(ctx, "ctx_by_pin_visit.csv", row.names = FALSE, na = "")
  # write.csv(wide, "ctx_pin_by_visit_wide.csv", row.names = FALSE, na = "")
  # ggsave("ctx_trend.png", ctx_plot, width = 9, height = 5.5, dpi = 150)
}
