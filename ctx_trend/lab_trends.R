#!/usr/bin/env Rscript
#
# Per-lab trend graphs, one panel per participant, on shared axes.
#
# Works on a data frame already loaded in your session — named `labcorp` below —
# holding the Labcorp study results (one row per participant / visit / draw,
# one column per lab). For every lab it builds a single figure faceted by
# participant, so a lab can be read across people at a glance.
#
# Axes are held identical:
#   * y  — one scale per lab, shared by every participant panel, spanning that
#          lab's full observed range, so panels are directly comparable.
#   * x  — the same visit range on every panel of every lab, taken from all
#          visits observed anywhere in the data.
#
# Pre/post draws are recorded two different ways, and both are handled:
#   * PINs 001-006 ("old" scheme, Pre_Post_Old) — one draw per visit, and the
#     flag labels the whole visit: visits 4/9/14/19 are pre, 8/13/18/23 post.
#   * PINs 007+ ("new" scheme, Pre_Post_New) — visits 8/13/18/23 are drawn
#     twice, a pre and a post within the same visit.
#   Under the new scheme the two draws share an x position, so the post point
#   is nudged right by `nudge_post` to keep both visible; the segment joining
#   them is the within-visit change.
#
# Usage:
#     source("lab_trends.R")
#     long <- labs_long(labcorp)
#
#     plot_lab_trend(long, "Hemoglobin")     # one lab
#     plots <- plot_all_labs(long)           # named list, one ggplot per lab
#     save_lab_trends(long, "lab_trends")    # multi-page PDF + one PNG per lab

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Columns that identify the draw rather than report a result.
ID_COLS <- c("PIN", "Visit_Number", "Pre_Post_Old", "Pre_Post_New",
             "Collection_Date", "Ahuja_ID", "Specimen_ID", "ControlID", "Notes")

# Units and male reference intervals, taken from the two header rows of the
# CSV version of this template (the xlsx export drops them). Names are matched
# loosely, so they survive read.csv()/janitor renaming.
LAB_UNITS <- c(
  "AbsoluteCD4Helper" = "/uL", "%CD4Pos.Lymph." = "%", "Abs.CD8Suppressor" = "/uL",
  "%CD8Pos.LympH" = "%", "WBC" = "x10E3/uL", "RBC" = "x10E6/uL",
  "Hemoglobin" = "g/dL", "Hematocrit" = "%", "MCV" = "fL", "MCH" = "pg",
  "MCHC" = "g/dL", "RDW" = "%", "Platelets" = "x10E3/uL", "Neutrophils" = "%",
  "Lymphs" = "%", "Monocytes" = "%", "Eos" = "%", "Basos" = "%",
  "Neutrophils(Absolute)" = "x10E3/uL", "Lymphs(Absolute)" = "x10E3/uL",
  "Monocytes(Absolute)" = "x10E3/uL", "Eos(Absolute)" = "x10E3/uL",
  "Baso(Absolute)" = "x10E3/uL", "ImmatureGranulocytes" = "%",
  "ImmatureGrans(Abs)" = "x10E3/uL", "Glucose" = "mg/dL", "BUN" = "mg/dL",
  "Creatinine" = "mg/dL", "eGFR" = "mL/min/1.73", "Sodium" = "mmol/L",
  "Potassium" = "mmol/L", "Chloride" = "mmol/L", "CarbonDioxide,Total" = "mmol/L",
  "Calcium" = "mg/dL", "Protein,Total" = "g/dL", "Albumin" = "g/dL",
  "Globulin,Total" = "g/dL", "Bilirubin,Total" = "mg/dL",
  "AlkalinePhosphatase" = "IU/L", "AST(SGOT)" = "IU/L", "ALT(SGPT)" = "IU/L",
  "Cholesterol,Total" = "mg/dL", "Triglycerides" = "mg/dL",
  "HDLCholesterol" = "mg/dL", "VLDLCholesterolCal" = "mg/dL",
  "LDLCholCalc(NIH)" = "mg/dL", "Testosterone" = "ng/dL",
  "VitaminD,25-Hydroxy" = "ng/mL", "C-ReactiveProtein,Cardiac" = "mg/L",
  "Iron" = "ug/dL", "CreatineKinase,Total" = "U/L", "Ferritin" = "ng/mL",
  "Myoglobin,Serum" = "ng/mL", "Prothrombin_Time" = "sec", "aPTT" = "sec",
  "HemoglobinA1c" = "%", "C-Telopeptide,Serum" = "pg/mL", "LDH" = "IU/L",
  "(LD)Fraction1" = "%", "(LD)Fraction2" = "%", "(LD)Fraction3" = "%",
  "(LD)Fraction4" = "%", "(LD)Fraction5" = "%",
  "ImmunoglobulinG,Qn,Serum" = "mg/dL", "ImmunoglobulinA,Qn,Serum" = "mg/dL",
  "ImmunoglobulinM,Qn,Serum" = "mg/dL", "ImmunoglobulinE,Total" = "IU/mL"
)

LAB_REF_RANGES <- c(
  "AbsoluteCD4Helper" = "359-1519", "%CD4Pos.Lymph." = "30.8-58.5",
  "Abs.CD8Suppressor" = "109-897", "%CD8Pos.LympH" = "12.0-35.5",
  "CD4/CD8Ratio" = "0.92-3.72", "WBC" = "3.4-10.8", "RBC" = "4.14-5.80",
  "Hemoglobin" = "13.0-17.7", "Hematocrit" = "37.5-51.0", "MCV" = "79-97",
  "MCH" = "26.6-33.0", "MCHC" = "31.5-35.7", "RDW" = "11.6-15.4",
  "Platelets" = "150-450", "Neutrophils(Absolute)" = "1.4-7.0",
  "Lymphs(Absolute)" = "0.7-3.1", "Monocytes(Absolute)" = "0.1-0.9",
  "Eos(Absolute)" = "0.0-0.4", "Baso(Absolute)" = "0.0-0.2",
  "ImmatureGrans(Abs)" = "0.0-0.1", "Glucose" = "70-99", "BUN" = "6-20",
  "Creatinine" = "0.76-1.27", "eGFR" = ">59", "BUN/CreatinineRatio" = "9-20",
  "Sodium" = "134-144", "Potassium" = "3.5-5.2", "Chloride" = "96-106",
  "CarbonDioxide,Total" = "20-29", "Calcium" = "8.7-10.2",
  "Protein,Total" = "6.0-8.5", "Albumin" = "4.1-5.1", "Globulin,Total" = "1.5-4.5",
  "Bilirubin,Total" = "0.0-1.2", "AlkalinePhosphatase" = "44-121",
  "AST(SGOT)" = "0-40", "ALT(SGPT)" = "0-44", "Cholesterol,Total" = "100-199",
  "Triglycerides" = "0-149", "HDLCholesterol" = ">39",
  "VLDLCholesterolCal" = "5-40", "LDLCholCalc(NIH)" = "0-99",
  "Testosterone" = "264-916", "VitaminD,25-Hydroxy" = "30.0-100.0",
  "C-ReactiveProtein,Cardiac" = "0.00-3.00", "Iron" = "38-169",
  "CreatineKinase,Total" = "49-439", "Ferritin" = "30-400",
  "Myoglobin,Serum" = "28-72", "INR" = "0.9-1.2", "Prothrombin_Time" = "9.1-12.0",
  "aPTT" = "24-33", "HemoglobinA1c" = "4.8-5.6", "LDH" = "121-224",
  "(LD)Fraction1" = "17-32", "(LD)Fraction2" = "25-40", "(LD)Fraction3" = "17-27",
  "(LD)Fraction4" = "5-13", "(LD)Fraction5" = "4-20",
  "ImmunoglobulinG,Qn,Serum" = "603-1613", "ImmunoglobulinA,Qn,Serum" = "90-386",
  "ImmunoglobulinM,Qn,Serum" = "20-172", "ImmunoglobulinE,Total" = "6-495"
)

# ---------------------------------------------------------------------------
# Reading whatever `labcorp` looks like
# ---------------------------------------------------------------------------

#' Compare names on alphanumerics only.
#'
#' `read.csv()` turns "Bilirubin,Total" into "Bilirubin.Total" and
#' `janitor::clean_names()` into "bilirubin_total", while readxl/readr keep it
#' verbatim. Normalising makes all of them match.
norm_name <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

#' Locate a column by name, tolerating however `labcorp` was read in.
find_col <- function(df, name, required = TRUE) {
  hit <- which(norm_name(names(df)) == norm_name(name))
  if (length(hit) == 0) {
    if (required) stop(sprintf("column '%s' not found in the data frame", name), call. = FALSE)
    return(NA_character_)
  }
  names(df)[hit[1]]
}

#' Pull a column as trimmed character, blanks as NA.
#'
#' Everything is routed through character so one coercion path covers every
#' reader — Visit_Number in particular arrives mixed numeric/text for PIN 009.
col_chr <- function(df, name, required = TRUE) {
  col <- find_col(df, name, required = required)
  if (is.na(col)) return(rep(NA_character_, nrow(df)))
  dplyr::na_if(trimws(as.character(df[[col]])), "")
}

#' Pull a column as Date, whether it arrived as Date, POSIXct, or text.
col_date <- function(df, name) {
  x <- df[[find_col(df, name)]]
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  as.Date(trimws(as.character(x)), format = "%Y-%m-%d")
}

#' Look a lab up in a reference table by normalised name.
lookup <- function(table, lab) {
  hit <- which(norm_name(names(table)) == norm_name(lab))
  if (length(hit) == 0) "" else unname(table[hit[1]])
}

#' Parse a lab result, keeping censored values.
#'
#' Results come through mostly numeric, but the export carries a few
#' below-detection values ("<0.2" in Bilirubin) and lab comment codes ("WF" in
#' the LD fractions). Censored values are kept at their detection limit and
#' flagged; comment codes become NA and drop out.
parse_result <- function(x) {
  chr <- trimws(as.character(x))
  chr[chr == ""] <- NA_character_
  censored <- grepl("^[<>]", chr) & !is.na(chr)
  value <- suppressWarnings(as.numeric(sub("^[<>=]+", "", chr)))
  list(value = value, censored = censored & !is.na(value))
}

# ---------------------------------------------------------------------------
# Reshaping
# ---------------------------------------------------------------------------

#' Which columns hold lab results?
#'
#' Everything that is not an identifier and parses to at least one number.
lab_columns <- function(labcorp) {
  candidates <- setdiff(names(labcorp), sapply(ID_COLS, find_col, df = labcorp, required = FALSE))
  keep <- vapply(candidates, function(nm) any(!is.na(parse_result(labcorp[[nm]])$value)), logical(1))
  candidates[keep]
}

#' Tidy one row per participant / visit / draw / lab.
#'
#' @param labcorp data frame of the export, header row already used as names
#' @param labs    optional character vector of labs to keep (default: all)
#' @return tibble with pin, visit, pre_post, draw, scheme, paired_visit,
#'   collection_date, lab, value, censored
labs_long <- function(labcorp, labs = NULL) {
  pin <- col_chr(labcorp, "PIN")

  # Drop the units / reference-interval rows if this copy still carries them,
  # plus any blank filler rows at the bottom of the template.
  keep <- !is.na(pin) & !grepl("^(ParticipantID|ReferenceIntervals)", pin)
  labcorp <- labcorp[keep, , drop = FALSE]

  pre_post_old <- suppressWarnings(as.integer(col_chr(labcorp, "Pre_Post_Old", required = FALSE)))
  pre_post_new <- suppressWarnings(as.integer(col_chr(labcorp, "Pre_Post_New", required = FALSE)))

  draws <- tibble(
    row_id = seq_len(nrow(labcorp)),
    pin = col_chr(labcorp, "PIN"),
    visit = suppressWarnings(as.numeric(col_chr(labcorp, "Visit_Number"))),
    pre_post = coalesce(pre_post_new, pre_post_old),
    # Which column the flag came from — the schemes mean different things.
    scheme = case_when(!is.na(pre_post_new) ~ "new", !is.na(pre_post_old) ~ "old",
                       TRUE ~ NA_character_),
    collection_date = col_date(labcorp, "Collection_Date")
  ) %>%
    mutate(
      draw = case_when(pre_post == 0 ~ "pre", pre_post == 1 ~ "post", TRUE ~ "single"),
      draw = factor(draw, levels = c("pre", "post", "single"))
    ) %>%
    # A participant's scheme is a property of the participant, not of one row:
    # carry it across the visits where no flag was recorded.
    group_by(pin) %>%
    mutate(scheme = if (all(is.na(scheme))) NA_character_ else
      names(sort(table(scheme), decreasing = TRUE))[1]) %>%
    # Under the new scheme a visit holds both a pre and a post draw; under the
    # old one the flag labels the whole visit, one draw each.
    group_by(pin, visit) %>%
    mutate(paired_visit = n() > 1) %>%
    ungroup()

  labs <- labs %||% lab_columns(labcorp)
  missing <- setdiff(labs, names(labcorp))
  if (length(missing)) {
    stop("lab column(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  values <- lapply(labs, function(nm) {
    parsed <- parse_result(labcorp[[nm]])
    tibble(row_id = seq_len(nrow(labcorp)),
           # One header ("ImmunoglobulinG,Qn,\nSerum") wraps across two lines;
           # the stray newline breaks text rendering in the PDF device.
           lab = gsub("\\s+", " ", trimws(nm)),
           value = parsed$value, censored = parsed$censored)
  })

  draws %>%
    inner_join(bind_rows(values), by = "row_id", relationship = "one-to-many") %>%
    filter(!is.na(value), !is.na(pin)) %>%
    # Protocol order, not date order: one post draw (PIN 007 visit 23) is dated
    # two days before its own pre draw, so dates would reverse that pair.
    arrange(lab, pin, visit, pre_post, collection_date) %>%
    select(-row_id)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

#' Parse "70-99" or ">59" into numeric bounds; NULL when there is no range.
parse_ref_range <- function(txt) {
  if (!nzchar(txt)) return(NULL)
  if (grepl("^>", txt)) return(c(as.numeric(sub("^>=?", "", txt)), Inf))
  if (grepl("^<", txt)) return(c(-Inf, as.numeric(sub("^<=?", "", txt))))
  bounds <- suppressWarnings(as.numeric(strsplit(txt, "-", fixed = TRUE)[[1]]))
  if (length(bounds) != 2 || anyNA(bounds)) NULL else bounds
}

#' Trend graph for one lab: a panel per participant, all on the same axes.
#'
#' @param long       output of `labs_long()`
#' @param lab        lab name (matched loosely against the column name)
#' @param x_limits   shared visit range; defaults to every visit in `long`, so
#'   the x axis is identical across labs as well as across participants
#' @param pins       participants to panel; defaults to every PIN in `long`, so
#'   the panel layout is identical from lab to lab
#' @param ref_band   shade the male reference interval where one is known
#' @param nudge_post offset for the post draw of a paired visit, so the two
#'   draws at one visit stay distinguishable (0 to stack them exactly)
#' @param ncol       panel columns
plot_lab_trend <- function(long, lab, x_limits = NULL, pins = NULL,
                           ref_band = TRUE, nudge_post = 0.35, ncol = 3) {
  lab_name <- unique(long$lab[norm_name(long$lab) == norm_name(lab)])
  if (length(lab_name) == 0) stop("no data for lab '", lab, "'", call. = FALSE)

  pins <- pins %||% sort(unique(long$pin))
  x_limits <- x_limits %||% range(long$visit, na.rm = TRUE)

  dat <- long %>%
    filter(lab == lab_name, pin %in% pins) %>%
    mutate(
      pin = factor(pin, levels = pins),
      # Separate the two draws of a paired visit; single draws stay on the tick.
      x = visit + ifelse(paired_visit & draw == "post", nudge_post, 0)
    )
  if (nrow(dat) == 0) stop("no data for lab '", lab, "'", call. = FALSE)

  units <- lookup(LAB_UNITS, lab_name)
  y_label <- if (nzchar(units)) paste0(lab_name, " (", units, ")") else lab_name
  # Pad the shared range so a participant sitting at the lab's min or max is not
  # clipped in half by coord_cartesian().
  y_limits <- range(dat$value, na.rm = TRUE)
  pad <- if (diff(y_limits) > 0) 0.04 * diff(y_limits) else max(abs(y_limits[1]) * 0.05, 0.5)
  y_limits <- y_limits + c(-pad, pad)

  p <- ggplot(dat, aes(x = x, y = value))

  # Reference interval, clipped by coord_cartesian below so it never widens the
  # shared y axis.
  ref <- if (ref_band) parse_ref_range(lookup(LAB_REF_RANGES, lab_name)) else NULL
  if (!is.null(ref)) {
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = ref[1], ymax = ref[2], fill = "grey88", alpha = 0.55)
  }

  # Participants 010-014 have a single draw; a one-point line warns and draws
  # nothing, so the line layer only sees panels that can actually show a trend.
  joined <- dat %>% group_by(pin) %>% filter(n() > 1) %>% ungroup()

  p +
    geom_line(data = joined, aes(group = pin), colour = "grey45", linewidth = 0.5) +
    geom_point(aes(colour = draw, shape = censored), size = 2.1) +
    facet_wrap(~pin, ncol = ncol, drop = FALSE,
               labeller = labeller(pin = function(x) paste("PIN", x))) +
    # A tick at every visit, but adjacent visits (3/4, 8/9, ...) collide when
    # labelled, so let the guide drop labels it cannot fit.
    scale_x_continuous(breaks = sort(unique(long$visit)),
                       guide = guide_axis(check.overlap = TRUE)) +
    scale_colour_manual(values = c(pre = "#2c7fb8", post = "#d95f02", single = "grey30"),
                        drop = FALSE, name = "Draw") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
    # Fixed on both axes: every panel of every lab shares the same window.
    coord_cartesian(xlim = x_limits + c(0, nudge_post), ylim = y_limits) +
    labs(
      title = paste0(lab_name, " by visit"),
      subtitle = "One panel per participant; axes shared across panels",
      x = "Visit number", y = y_label
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          axis.text.x = element_text(size = 7),
          strip.text = element_text(face = "bold"))
}

#' Every lab in `long`, as a named list of ggplots on a common x axis.
plot_all_labs <- function(long, ...) {
  labs <- unique(long$lab)
  x_limits <- range(long$visit, na.rm = TRUE)
  pins <- sort(unique(long$pin))
  setNames(
    lapply(labs, function(nm) plot_lab_trend(long, nm, x_limits = x_limits, pins = pins, ...)),
    labs
  )
}

#' Write the trend graphs to disk: one multi-page PDF, plus a PNG per lab.
#'
#' The PDF is usually the useful artefact — one lab per page, flip through them.
save_lab_trends <- function(long, outdir = "lab_trends", png = TRUE,
                            width = 12, height = 7.5, dpi = 150, ...) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  plots <- plot_all_labs(long, ...)

  pdf_path <- file.path(outdir, "lab_trends.pdf")
  grDevices::pdf(pdf_path, width = width, height = height)
  for (p in plots) print(p)
  grDevices::dev.off()

  if (png) {
    for (nm in names(plots)) {
      safe <- gsub("[^A-Za-z0-9]+", "_", nm)
      ggsave(file.path(outdir, paste0(safe, ".png")), plots[[nm]],
             width = width, height = height, dpi = dpi)
    }
  }
  cat(sprintf("Wrote %d lab graphs to %s/\n", length(plots), outdir))
  invisible(plots)
}

# ---------------------------------------------------------------------------
# Run on `labcorp`
# ---------------------------------------------------------------------------
# Guarded so `source("lab_trends.R")` just defines the functions when `labcorp`
# is not (yet) in the session.

if (exists("labcorp")) {
  long <- labs_long(labcorp)

  cat(sprintf("%d results across %d labs, %d participants, visits %s-%s\n",
              nrow(long), dplyr::n_distinct(long$lab), dplyr::n_distinct(long$pin),
              min(long$visit, na.rm = TRUE), max(long$visit, na.rm = TRUE)))

  # One figure per lab, panelled by participant:
  #   plot_lab_trend(long, "Hemoglobin")
  #   plots <- plot_all_labs(long)
  #   save_lab_trends(long, "lab_trends")
}
