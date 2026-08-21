# =============================================================================
# Ionized calcium report
# -----------------------------------------------------------------------------
# Builds a multi-page landscape PDF from a blood-gas export:
#
#   1  Scope and methods
#   2  Data and measurement checks
#   3  Pre/post pair table
#   4  Raw ionized calcium, per participant
#   5  Raw ionized calcium, all participants
#   6  Estimated plasma-volume change (two estimators)
#   7  Post-exercise iCa: raw and plasma-volume adjusted
#   8  Post minus pre iCa across the three definitions
#   9  Haemoconcentration markers (ctHb, Hct)
#  10  Sodium (cNa+)
#  11  Chloride (cCl-)
#  12  Dyshaemoglobin fractions (FCOHb, FMetHb)
#      (+ optional pH pages when PH_MODE = "sensitivity")
#
# Method:
#   %dPV (Dill-Costill) = 100*[(Hb_pre/Hb_post)*((1-Hct_post)/(1-Hct_pre)) - 1]
#   %dPV (Hb only)      = 100*[(Hb_pre/Hb_post) - 1]
#   iCa_PVadj,post      = iCa_post * (1 + %dPV/100)
#
# Raw iCa is the concentration relevant to calcium-sensing receptor signalling.
# PV-adjusted iCa estimates whether vascular calcium CONTENT changed
# independently of haemoconcentration. Both are reported; they answer
# different questions.
# =============================================================================


# =============================================================================
# 1. CONFIGURATION  -- edit this block only
# =============================================================================

DATA_FILE  <- "icalpin14.xlsx"                 # path to the workbook
SHEET      <- "Blood_Gas"                      # sheet name or index
OUT_PDF    <- "ionized_calcium_report.pdf"     # output file

# Study structure: which visit belongs to which week block.
# Each block is drawn as a shaded band behind every panel.
WEEK_MAP <- c(
  "4"  = "Control", "8"  = "Control",
  "9"  = "Heat 1",  "13" = "Heat 1",
  "14" = "Heat 2",  "18" = "Heat 2",
  "19" = "Heat 3",  "23" = "Heat 3"
)
WEEK_LEVELS <- c("Control", "Heat 1", "Heat 2", "Heat 3")

# pH handling.
#   "exclude"     - pH is treated as unusable: no pH-normalised iCa, no pH
#                   covariate, no acid-base commentary. The scope page says so
#                   explicitly so the gap is on the record.
#   "sensitivity" - adds a pH page and a pH-normalised iCa page, presented as a
#                   sensitivity analysis, NOT as the primary readout.
PH_MODE     <- "exclude"
PH_REF      <- 7.40    # reference pH for normalisation
PH_SLOPE    <- 0.05    # fractional change in iCa per 0.1 pH unit (~5%)

# Repeated-measures model. Left off by default: with a handful of pairs across
# a couple of participants a random-intercept model is not estimable, and a
# printed table would imply more than the data support.
FIT_MODEL   <- FALSE

# Analytes carried through the report. Everything else on the blood-gas panel
# is dropped rather than reported without comment.
MEASURES_USED <- c("Ionized_Calcium", "cNa+", "cCl-", "ctHb", "Hct",
                   "FCOHb", "FMetHb")

FCOHB_FLAG  <- 3.0     # flag carboxyhaemoglobin at or above this % of ctHb
DPV_EXTREME <- 15      # flag |%dPV| at or above this in the checks page

PAGE_W <- 11; PAGE_H <- 8.5    # US letter, landscape


# =============================================================================
# 2. PACKAGES
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(grid)
})


# =============================================================================
# 3. LOAD AND CLEAN
# =============================================================================

# Read everything as text. The export mixes numeric and character storage in
# the same column (1.09 stored as a number, "1.14" as text, ".5" as text), so
# letting readxl guess produces a column of NAs and list-columns.
raw <- readxl::read_excel(DATA_FILE, sheet = SHEET, col_types = "text")

# Column lookup that tolerates the small naming variations between exports
# (cCl- vs cCl−, "Date (yyyy-mm-dd)" vs "Date").
pick <- function(df, ..., required = TRUE) {
  wanted <- c(...)
  norm <- function(x) tolower(gsub("[^a-z0-9]", "", tolower(x)))
  hit <- match(norm(wanted), norm(names(df)))
  hit <- hit[!is.na(hit)]
  if (!length(hit)) {
    if (required) stop("Column not found: ", paste(wanted, collapse = " / "))
    return(NULL)
  }
  names(df)[hit[1]]
}

num <- function(x) suppressWarnings(as.numeric(trimws(x)))

parse_excel_date <- function(x) {
  x <- trimws(x)
  serial <- suppressWarnings(as.numeric(x))
  out <- as.Date(rep(NA_real_, length(x)), origin = "1970-01-01")
  ok <- !is.na(serial)
  out[ok] <- as.Date(serial[ok], origin = "1899-12-30")
  if (any(!ok & !is.na(x))) {
    txt <- suppressWarnings(as.Date(substr(x[!ok & !is.na(x)], 1, 10)))
    out[!ok & !is.na(x)] <- txt
  }
  out
}

dat <- tibble::tibble(
  pin            = trimws(raw[[pick(raw, "PIN")]]),
  visit          = num(raw[[pick(raw, "Visit_Number", "Visit")]]),
  phase_raw      = tolower(trimws(raw[[pick(raw, "Pre_Post", "Phase")]])),
  collection_date= parse_excel_date(raw[[pick(raw, "Date (yyyy-mm-dd)", "Date")]]),
  ica            = num(raw[[pick(raw, "Ionized_Calcium")]]),
  ph             = num(raw[[pick(raw, "pH")]]),
  hb             = num(raw[[pick(raw, "ctHb")]]),
  hct            = num(raw[[pick(raw, "Hct")]]),
  na             = num(raw[[pick(raw, "cNa+")]]),
  cl             = num(raw[[pick(raw, "cCl-", "cCl−")]]),
  fcohb          = num(raw[[pick(raw, "FCOHb")]]),
  fmethb         = num(raw[[pick(raw, "FMetHb")]]),
  notes          = if (!is.null(pick(raw, "Notes", required = FALSE)))
                     trimws(raw[[pick(raw, "Notes")]]) else NA_character_
) %>%
  filter(!is.na(pin), !is.na(visit)) %>%
  mutate(
    phase = factor(ifelse(phase_raw %in% c("pre", "pre-exercise", "preexercise"),
                          "Pre-exercise", "Post-exercise"),
                   levels = c("Pre-exercise", "Post-exercise")),
    # PIN 09 / PIN 12 / PIN 14 - zero-padded so panels sort naturally
    pin_label = factor(sprintf("PIN %02d", num(pin)),
                       levels = sprintf("PIN %02d", sort(unique(num(pin))))),
    week = factor(unname(WEEK_MAP[as.character(visit)]), levels = WEEK_LEVELS)
  ) %>%
  select(-phase_raw) %>%
  arrange(pin_label, visit, phase)

# Hct enters Dill-Costill as a FRACTION. The export reports percent.
HCT_IS_PERCENT <- isTRUE(median(dat$hct, na.rm = TRUE) > 1)
dat$hct_frac <- if (HCT_IS_PERCENT) dat$hct / 100 else dat$hct

# Visits actually present, used as the shared x breaks on every panel.
VISITS <- sort(unique(dat$visit))

if (any(is.na(dat$week))) {
  warning("Visits with no week assignment in WEEK_MAP: ",
          paste(sort(unique(dat$visit[is.na(dat$week)])), collapse = ", "))
}


# =============================================================================
# 4. DERIVED QUANTITIES
# =============================================================================

dill_costill <- function(hb_pre, hb_post, hct_pre, hct_post) {
  100 * ((hb_pre / hb_post) * ((1 - hct_post) / (1 - hct_pre)) - 1)
}
hb_only <- function(hb_pre, hb_post) {
  100 * (hb_pre / hb_post - 1)
}
pv_adjust <- function(ica_post, dpv) {
  ica_post * (1 + dpv / 100)
}
# iCa normalised to PH_REF: measured iCa rises ~5% per 0.1 unit fall in pH.
ph_normalise <- function(ica, ph) {
  ica * (1 + (PH_SLOPE / 0.1) * (ph - PH_REF))
}

wide <- dat %>%
  select(pin, pin_label, visit, week, phase, collection_date,
         ica, ph, hb, hct, hct_frac, na, cl) %>%
  pivot_wider(
    id_cols     = c(pin, pin_label, visit, week),
    names_from  = phase,
    values_from = c(collection_date, ica, ph, hb, hct, hct_frac, na, cl),
    names_sep   = "."
  )

nm <- function(x, phase) {
  col <- paste0(x, ".", phase)
  if (col %in% names(wide)) wide[[col]] else rep(NA_real_, nrow(wide))
}

pairs_tbl <- wide %>%
  mutate(
    ica_pre   = nm("ica",  "Pre-exercise"),  ica_post  = nm("ica",  "Post-exercise"),
    ph_pre    = nm("ph",   "Pre-exercise"),  ph_post   = nm("ph",   "Post-exercise"),
    hb_pre    = nm("hb",   "Pre-exercise"),  hb_post   = nm("hb",   "Post-exercise"),
    hct_pre   = nm("hct",  "Pre-exercise"),  hct_post  = nm("hct",  "Post-exercise"),
    hctf_pre  = nm("hct_frac", "Pre-exercise"),
    hctf_post = nm("hct_frac", "Post-exercise"),
    na_pre    = nm("na",   "Pre-exercise"),  na_post   = nm("na",   "Post-exercise"),
    cl_pre    = nm("cl",   "Pre-exercise"),  cl_post   = nm("cl",   "Post-exercise")
  ) %>%
  filter(!is.na(ica_pre) & !is.na(ica_post)) %>%
  mutate(
    dpv_dc      = dill_costill(hb_pre, hb_post, hctf_pre, hctf_post),
    dpv_hb      = hb_only(hb_pre, hb_post),
    pvadj_dc    = pv_adjust(ica_post, dpv_dc),
    pvadj_hb    = pv_adjust(ica_post, dpv_hb),
    d_raw       = ica_post - ica_pre,
    d_pvadj_dc  = pvadj_dc - ica_pre,
    d_pvadj_hb  = pvadj_hb - ica_pre,
    d_na        = na_post - na_pre,
    d_cl        = cl_post - cl_pre,
    d_ph        = ph_post - ph_pre
  ) %>%
  arrange(pin_label, visit)

# Shared y-limits for every iCa panel: raw and both adjusted series on one
# scale, padded by 6% of the range, so the pages are directly comparable.
ica_all <- c(dat$ica, pairs_tbl$pvadj_dc, pairs_tbl$pvadj_hb)
if (identical(PH_MODE, "sensitivity")) {
  ica_all <- c(ica_all, ph_normalise(dat$ica, dat$ph))
}
ica_all  <- ica_all[is.finite(ica_all)]
ICA_YLIM <- range(ica_all) + c(-1, 1) * 0.06 * diff(range(ica_all))

N_ICA        <- sum(!is.na(dat$ica))
N_PAIRS      <- nrow(pairs_tbl)
N_PAIRS_DPV  <- sum(!is.na(pairs_tbl$dpv_dc))
POST_VISITS  <- sort(unique(dat$visit[dat$phase == "Post-exercise"]))


# =============================================================================
# 5. LOOK AND FEEL
# =============================================================================

WEEK_COLS <- c("Control" = "#D9D9D9", "Heat 1" = "#FDD49E",
               "Heat 2"  = "#FC8D59", "Heat 3" = "#D7301F")
DRAW_COLS <- c("Pre-exercise" = "#2C7FB8", "Post-exercise" = "#D95F02")
EST_COLS  <- c("Dill-Costill (Hb + Hct)" = "#762A83",
               "Haemoglobin only"        = "#1B7837")
POST_COLS <- c("Post: raw"                          = "#D95F02",
               "Post: PV-adjusted (Dill-Costill)"   = "#762A83",
               "Post: PV-adjusted (Hb only)"        = "#1B7837")
POST_SHAPES <- c("Post: raw" = 16,
                 "Post: PV-adjusted (Dill-Costill)" = 15,
                 "Post: PV-adjusted (Hb only)" = 17)

dark2 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
           "#66A61E", "#E6AB02", "#A6761D", "#666666")
PIN_COLS <- setNames(rep_len(dark2, nlevels(dat$pin_label)),
                     levels(dat$pin_label))

theme_report <- function() {
  theme_bw(base_size = 11, base_family = "sans") +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey92", linewidth = 0.25),
      panel.border      = element_rect(colour = "grey20", fill = NA, linewidth = 0.5),
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.background   = element_rect(fill = "white", colour = NA),
      strip.background  = element_rect(fill = "grey92", colour = "grey20", linewidth = 0.5),
      strip.text        = element_text(size = 9, colour = "black",
                                       margin = margin(4, 4, 4, 4)),
      axis.text         = element_text(size = 9, colour = "grey20"),
      axis.title        = element_text(size = 11),
      axis.ticks        = element_line(colour = "grey20", linewidth = 0.4),
      plot.title        = element_text(size = 13, hjust = 0, margin = margin(b = 3)),
      plot.subtitle     = element_text(size = 11, hjust = 0, colour = "grey25",
                                       margin = margin(b = 8)),
      plot.caption      = element_text(size = 8, hjust = 0, colour = "grey30",
                                       margin = margin(t = 8)),
      legend.position   = "bottom",
      legend.title      = element_text(size = 11),
      legend.text       = element_text(size = 9),
      legend.key        = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.margin     = margin(2, 2, 2, 2)
    )
}

# Week bands sit behind everything, spanning the visits assigned to each week.
week_bands <- dat %>%
  filter(!is.na(week)) %>%
  group_by(week) %>%
  summarise(xmin = min(visit) - 0.5, xmax = max(visit) + 0.5, .groups = "drop")

layer_week_bands <- function(alpha = 0.35) {
  list(
    geom_rect(data = week_bands, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = week),
              alpha = alpha),
    scale_fill_manual(name = "Week", values = WEEK_COLS, drop = FALSE)
  )
}

scale_visit <- function() {
  # limits pin the full visit set on every page, including analytes that are
  # missing at the first or last visit; oob_keep lets the week bands run past
  # the outermost visit to the panel edge instead of being dropped.
  scale_x_continuous(name = "Visit number", breaks = VISITS,
                     limits = range(VISITS), oob = scales::oob_keep,
                     expand = expansion(mult = 0.05))
}

wrap_txt <- function(..., width = 150) {
  paste(strwrap(paste0(...), width = width), collapse = "\n")
}


# =============================================================================
# 6. TEXT-PAGE RENDERER
#    Headings in bold sans; body in monospace so that printed tables and the
#    formulae line up. Positions are in points from the top-left of the page.
# =============================================================================

pt_u  <- function(x) unit(x / 72, "in")
top_u <- function(x) unit(1, "npc") - pt_u(x)

text_page <- function(heading, lines,
                      size = 8, first_top = 78, lead = 15.15,
                      x_left = 54.1, indent = 9.6, fit = FALSE) {
  grid.newpage()
  grid.text(heading, x = pt_u(19.9), y = top_u(19), just = c("left", "top"),
            gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = 13))

  n <- length(lines)
  if (fit && n > 1) {
    avail <- (PAGE_H * 72) - 36 - first_top
    lead  <- min(lead, avail / (n - 1))
  }
  for (i in seq_along(lines)) {
    ln <- lines[i]
    if (!nzchar(sub("\\s+$", "", ln))) next
    ind <- grepl("^  ", ln)
    grid.text(sub("^  ", "", ln),
              x = pt_u(if (ind) x_left + indent else x_left),
              y = top_u(first_top + (i - 1) * lead),
              just = c("left", "top"),
              gp = gpar(fontfamily = "mono", fontsize = size))
  }
  invisible(NULL)
}

# Body text for the scope page is authored at a fixed column count so it reads
# the same regardless of the data.
body <- function(..., width = 74) paste0("  ", strwrap(paste0(...), width = width))

plural <- function(n, one, many) if (n == 1) one else many

# Render a data frame as monospace lines, indented one level.
mono_table <- function(df, digits = NULL) {
  if (!nrow(df)) return("  (none)")
  if (!is.null(digits)) df <- as.data.frame(df)
  paste0("  ", capture.output(print(as.data.frame(df), row.names = FALSE)))
}


# =============================================================================
# 7. PAGES 1-3: SCOPE, CHECKS, PAIR TABLE
# =============================================================================

page_scope <- function() {
  ph_block <- if (identical(PH_MODE, "exclude")) {
    c("pH excluded",
      body("pH was unreliable in this run, so no pH-normalised iCa, no pH ",
           "covariate, and no acid-base commentary appear here. Ionized ",
           "calcium is sensitive to pH through albumin binding, so the ",
           "acid-base component of the pre-to-post change is UNADDRESSED ",
           "rather than controlled. Differences reported below may contain ",
           "an acid-base contribution that cannot be separated with the ",
           "available data."))
  } else {
    c("pH: sensitivity analysis only",
      body("Measured iCa is reported as the primary readout. iCa normalised ",
           "to pH ", sprintf("%.2f", PH_REF), " (", format(PH_SLOPE * 100),
           "% per 0.1 unit) is shown separately as a sensitivity analysis. ",
           "The in-vivo pH shift is part of the exercise response, so the ",
           "normalised series is not treated as the physiological value."))
  }

  model_block <- if (FIT_MODEL) {
    c("Model", body("Repeated-measures model fitted; see the model page."))
  } else {
    c("Not run",
      body("Repeated-measures model (FIT_MODEL = FALSE). With ", N_PAIRS,
           plural(N_PAIRS, " pair", " pairs"), " across ",
           nlevels(dat$pin_label),
           plural(nlevels(dat$pin_label), " participant", " participants"),
           " it is not estimable."))
  }

  lines <- c(
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    paste0("Participants: ",
           paste(as.numeric(sub("PIN ", "", levels(dat$pin_label))), collapse = ", "),
           " | iCa values: ", N_ICA,
           " | pre/post pairs: ", N_PAIRS,
           " | pairs with %dPV: ", N_PAIRS_DPV),
    "",
    "Measures used",
    body(paste(MEASURES_USED, collapse = ", "), "."),
    body("All other channels on the blood-gas panel are excluded from this report."),
    "",
    ph_block,
    "",
    "Quantities computed",
    "  %dPV (Dill-Costill) = 100*[(Hb_pre/Hb_post)*((1-Hct_post)/(1-Hct_pre)) - 1]",
    "  %dPV (Hb only)      = 100*[(Hb_pre/Hb_post) - 1]",
    "  iCa_PVadj,post      = iCa_post * (1 + %dPV/100)",
    "",
    "Interpretation",
    body("Raw iCa is the concentration relevant to calcium-sensing receptor ",
         "signalling. PV-adjusted iCa estimates whether vascular calcium ",
         "content changed independently of haemoconcentration. They answer ",
         "different questions and both are reported."),
    "",
    "Axis convention",
    sprintf("  All iCa panels share y = %.3f to %.3f mmol/L; x breaks are the full visit set.",
            ICA_YLIM[1], ICA_YLIM[2]),
    body("Every other analyte uses one fixed y-scale across all panels on its page."),
    "",
    model_block
  )
  text_page("Ionized calcium report - scope and methods", lines)
}

page_checks <- function() {
  # -- 1. is Hct independent of ctHb, as Dill-Costill assumes?
  ratio <- with(dat[!is.na(dat$hct) & !is.na(dat$hb) & dat$hb > 0, ], hct / hb)
  ratio_mean <- mean(ratio); ratio_cv <- 100 * sd(ratio) / ratio_mean
  infl <- with(pairs_tbl[is.finite(pairs_tbl$dpv_dc) & is.finite(pairs_tbl$dpv_hb) &
                           pairs_tbl$dpv_hb != 0, ],
               mean(dpv_dc / dpv_hb))

  # -- 2. reported Hb resolution against the effect being estimated
  hb_mean <- mean(dat$hb, na.rm = TRUE)
  hb_step <- 100 * 0.1 / hb_mean

  # -- 4. dyshaemoglobin fractions, components of the ctHb the adjustment uses
  flagged <- dat %>%
    filter(!is.na(fcohb), fcohb >= FCOHB_FLAG) %>%
    transmute(PIN = as.numeric(pin), Visit_Number = visit,
              phase = as.character(phase), FCOHb = fcohb, FMetHb = fmethb)

  # -- 5. pairs with iCa but no Hb/Hct, so no plasma-volume estimate
  no_dpv <- pairs_tbl %>%
    filter(is.na(dpv_dc)) %>%
    transmute(PIN = as.numeric(pin), Visit_Number = visit,
              iCa_pre = ica_pre, iCa_post = ica_post)

  # -- 6. rows carrying no iCa result at all
  no_ica <- dat %>%
    filter(is.na(ica)) %>%
    transmute(PIN = as.numeric(pin), Visit_Number = visit,
              phase = as.character(phase),
              collection_date = as.character(collection_date))

  # -- 7. pre/post structure
  post_ok <- dat %>% filter(phase == "Post-exercise") %>%
    left_join(dat %>% filter(phase == "Pre-exercise") %>%
                transmute(pin, visit, has_pre = TRUE),
              by = c("pin", "visit"))
  complete_txt <- if (all(!is.na(post_ok$has_pre)))
    "All designated visits have a complete pair." else
    "Some post draws have no matching pre draw."

  extreme <- pairs_tbl %>% filter(abs(dpv_dc) >= DPV_EXTREME)

  lines <- c(
    "1. Hct is not independent of ctHb in this file.",
    body(sprintf("Hct / ctHb = %.3f, CV %.2f%% across %d samples.",
                 ratio_mean, ratio_cv, length(ratio))),
    body("A separately measured haematocrit would vary with MCHC between ",
         "people and states; a near-constant ratio indicates the analyzer ",
         "derives Hct from ctHb. Dill-Costill assumes two independent ",
         "measurements. Here the Hct term adds no information and amplifies ",
         "the Hb signal: %dPV(DC) averages ", sprintf("%.2fx", infl),
         " the Hb-only estimate. Both are plotted; neither is a clean ",
         "paired-marker estimate. The ratio is measured; the derivation ",
         "mechanism is inference."),
    "",
    "2. Reported Hb resolution vs the effect size.",
    body(sprintf("ctHb is reported to 0.1 g/dL; at the mean of %.1f g/dL one digit is %.1f%% of the value.",
                 hb_mean, hb_step)),
    body("One rounding step therefore moves %dPV by a few percent, against ",
         "an expected plasma-volume effect on the order of 10-14%."),
    "",
    "3. Observed %dPV range.",
    body(sprintf("Dill-Costill: %.1f to %.1f%%",
                 min(pairs_tbl$dpv_dc, na.rm = TRUE), max(pairs_tbl$dpv_dc, na.rm = TRUE))),
    body(sprintf("Hb-only     : %.1f to %.1f%%",
                 min(pairs_tbl$dpv_hb, na.rm = TRUE), max(pairs_tbl$dpv_hb, na.rm = TRUE))),
    body(if (nrow(extreme))
           paste0(nrow(extreme), " of ", N_PAIRS_DPV,
                  plural(nrow(extreme), " pair exceeds", " pairs exceed"),
                  " the 10-14% range typically reported after ",
                  "moderate-to-high-intensity exercise.")
         else "All pairs fall within the 10-14% range typically reported after exercise."),
    "",
    "4. Dyshaemoglobin fractions (components of ctHb).",
    body(sprintf("FCOHb %.1f to %.1f%% | FMetHb %.1f to %.1f%%",
                 min(dat$fcohb, na.rm = TRUE), max(dat$fcohb, na.rm = TRUE),
                 min(dat$fmethb, na.rm = TRUE), max(dat$fmethb, na.rm = TRUE))),
    if (nrow(flagged)) mono_table(flagged) else
      body(sprintf("No sample reaches FCOHb %.1f%%.", FCOHB_FLAG)),
    "",
    "5. Pairs without Hb/Hct, so no %dPV:",
    mono_table(no_dpv),
    "",
    "6. Rows with no iCa result:",
    mono_table(no_ica),
    "",
    "7. Pre/post structure:",
    body(paste0("All post draws fall on visits ",
                paste(POST_VISITS, collapse = "/"), ".")),
    body(complete_txt)
  )
  text_page("Data and measurement checks", lines, fit = TRUE)
}

page_pair_table <- function() {
  tab <- pairs_tbl %>%
    transmute(
      PIN      = as.numeric(pin),
      Visit    = visit,
      iCa_pre  = sprintf("%.2f", ica_pre),
      iCa_post = sprintf("%.2f", ica_post),
      Hb_pre   = fnum(hb_pre, 1),
      Hb_post  = fnum(hb_post, 1),
      Hct_pre  = fnum(hct_pre, 1),
      Hct_post = fnum(hct_post, 1),
      dPV_DC   = fnum(dpv_dc, 1),
      dPV_Hb   = fnum(dpv_hb, 1),
      PVadj    = fnum(pvadj_dc, 3),
      d_raw    = sprintf("%.2f", d_raw),
      d_PVadj  = fnum(d_pvadj_dc, 3),
      d_Na     = fnum(d_na, 0),
      d_Cl     = fnum(d_cl, 0)
    )
  if (identical(PH_MODE, "sensitivity")) {
    tab$pH_pre  <- fnum(pairs_tbl$ph_pre, 3)
    tab$pH_post <- fnum(pairs_tbl$ph_post, 3)
    tab$d_pH    <- fnum(pairs_tbl$d_ph, 3)
  }

  old <- options(width = 96); on.exit(options(old), add = TRUE)
  lines <- capture.output(print(as.data.frame(tab), row.names = FALSE))
  text_page("Pre/post pair table", lines,
            size = 7, first_top = 89.7, lead = 26.3,
            x_left = 58.3, indent = 0, fit = TRUE)
}

# NA-safe fixed-decimal formatting: keeps "NA" in the table rather than "  NA"
# turning into a numeric coercion warning.
fnum <- function(x, d) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))


# =============================================================================
# 8. PLOT PAGES
# =============================================================================

# ---- page 4: raw iCa, per participant ---------------------------------------
plot_ica_facet <- function() {
  pre_series <- dat %>% filter(phase == "Pre-exercise", !is.na(ica))
  links <- dat %>%
    semi_join(pairs_tbl %>% select(pin, visit), by = c("pin", "visit")) %>%
    filter(!is.na(ica))

  ggplot(dat %>% filter(!is.na(ica)), aes(visit, ica)) +
    layer_week_bands() +
    geom_line(data = pre_series, aes(group = pin_label),
              colour = "grey35", linewidth = 0.4) +
    geom_line(data = links, aes(group = interaction(pin_label, visit)),
              colour = "grey35", linewidth = 0.4, linetype = "22") +
    geom_point(colour = "white", size = 3.2) +
    geom_point(aes(colour = phase), size = 2.2) +
    facet_wrap(~pin_label, nrow = 1) +
    scale_visit() +
    scale_y_continuous(name = "Ionized calcium (mmol/L)",
                       limits = ICA_YLIM, expand = c(0, 0)) +
    scale_colour_manual(name = "Draw", values = DRAW_COLS, drop = FALSE) +
    labs(title    = "Raw ionized calcium (as measured)",
         subtitle = "Primary readout: the concentration seen by the calcium-sensing receptor",
         caption  = wrap_txt("Solid line: pre-exercise series. Dashed: same-day ",
                             "pre-to-post pair. Y-axis is shared with the PV-adjusted page.")) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
    theme_report()
}

# ---- page 5: raw iCa, all participants on one panel -------------------------
plot_ica_pooled <- function() {
  pre_series <- dat %>% filter(phase == "Pre-exercise", !is.na(ica))
  links <- dat %>%
    semi_join(pairs_tbl %>% select(pin, visit), by = c("pin", "visit")) %>%
    filter(!is.na(ica))

  ggplot(dat %>% filter(!is.na(ica)), aes(visit, ica, colour = pin_label)) +
    layer_week_bands() +
    geom_line(data = pre_series, aes(group = pin_label), linewidth = 0.4) +
    geom_line(data = links, aes(group = interaction(pin_label, visit)),
              linewidth = 0.4, linetype = "22") +
    geom_point(aes(shape = phase), size = 2.2) +
    scale_visit() +
    scale_y_continuous(name = "Ionized calcium (mmol/L)",
                       limits = ICA_YLIM, expand = c(0, 0)) +
    scale_colour_manual(name = "Participant", values = PIN_COLS) +
    scale_shape_manual(name = "Draw", values = c("Pre-exercise" = 16,
                                                 "Post-exercise" = 17)) +
    labs(title   = "Raw ionized calcium, all participants",
         caption = "Same y-limits as the per-participant page.") +
    guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2),
           colour = guide_legend(order = 3)) +
    theme_report()
}

# ---- page 6: plasma-volume change, two estimators ---------------------------
plot_dpv <- function() {
  dpv_long <- pairs_tbl %>%
    select(pin_label, visit, week,
           `Dill-Costill (Hb + Hct)` = dpv_dc,
           `Haemoglobin only`        = dpv_hb) %>%
    pivot_longer(-c(pin_label, visit, week),
                 names_to = "estimator", values_to = "dpv") %>%
    filter(!is.na(dpv)) %>%
    mutate(estimator = factor(estimator, levels = names(EST_COLS)))

  ggplot(dpv_long, aes(visit, dpv, colour = estimator)) +
    layer_week_bands() +
    geom_hline(yintercept = 0, colour = "grey40", linetype = "22", linewidth = 0.3) +
    geom_line(aes(group = estimator), linewidth = 0.4) +
    geom_point(size = 2.2) +
    facet_wrap(~pin_label, nrow = 1) +
    scale_visit() +
    scale_y_continuous(name = "%dPV (negative = contraction)") +
    scale_colour_manual(name = NULL, values = EST_COLS) +
    labs(title   = "Estimated plasma-volume change, pre to post (%)",
         caption = wrap_txt("Two estimators shown because Hct in this file is a ",
                            "near-constant multiple of ctHb (see checks page): the ",
                            "Dill-Costill Hct term carries no independent information ",
                            "here and inflates the magnitude relative to the Hb-only estimate.")) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
    theme_report()
}

# ---- page 7: post-exercise iCa, raw vs PV-adjusted --------------------------
plot_post_adjusted <- function() {
  post_long <- pairs_tbl %>%
    select(pin_label, visit, week,
           `Post: raw`                        = ica_post,
           `Post: PV-adjusted (Dill-Costill)` = pvadj_dc,
           `Post: PV-adjusted (Hb only)`      = pvadj_hb) %>%
    pivot_longer(-c(pin_label, visit, week),
                 names_to = "series", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(series = factor(series, levels = names(POST_COLS)))

  ggplot(post_long, aes(visit, value, colour = series, shape = series)) +
    layer_week_bands() +
    geom_point(data = pairs_tbl, inherit.aes = FALSE,
               aes(visit, ica_pre), shape = 4, size = 2.4,
               colour = "grey25", stroke = 0.6) +
    geom_point(size = 2.2) +
    facet_wrap(~pin_label, nrow = 1) +
    scale_visit() +
    scale_y_continuous(name = "Ionized calcium (mmol/L)",
                       limits = ICA_YLIM, expand = c(0, 0)) +
    scale_colour_manual(name = NULL, values = POST_COLS) +
    scale_shape_manual(name = NULL, values = POST_SHAPES) +
    labs(title    = "Post-exercise ionized calcium: raw and plasma-volume adjusted",
         subtitle = "Grey crosses are the same-visit pre-exercise value",
         caption  = wrap_txt("PV-adjusted values estimate vascular calcium CONTENT, ",
                             "not the concentration sensed in vivo. Neither version is ",
                             "corrected for acid-base status, which is ",
                             if (identical(PH_MODE, "exclude"))
                               "not addressed in this report."
                             else "reported separately as a sensitivity analysis.")) +
    guides(fill = guide_legend(order = 1),
           colour = guide_legend(order = 2), shape = guide_legend(order = 2)) +
    theme_report()
}

# ---- page 8: post minus pre across the three definitions --------------------
plot_deltas <- function() {
  d_long <- pairs_tbl %>%
    select(pin_label, visit,
           `Raw`                             = d_raw,
           `PV-adjusted\n(Dill-Costill)`     = d_pvadj_dc,
           `PV-adjusted\n(Hb only)`          = d_pvadj_hb) %>%
    pivot_longer(-c(pin_label, visit),
                 names_to = "definition", values_to = "delta") %>%
    filter(!is.na(delta)) %>%
    mutate(definition = factor(definition,
                               levels = c("Raw", "PV-adjusted\n(Dill-Costill)",
                                          "PV-adjusted\n(Hb only)")),
           visit_f = factor(visit, levels = POST_VISITS))

  ggplot(d_long, aes(definition, delta,
                     group = interaction(pin_label, visit),
                     colour = pin_label)) +
    geom_hline(yintercept = 0, colour = "grey40", linetype = "22", linewidth = 0.3) +
    geom_line(linewidth = 0.4) +
    geom_point(aes(shape = visit_f), size = 2.4) +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "Post - pre iCa (mmol/L)") +
    scale_colour_manual(name = "Participant", values = PIN_COLS) +
    scale_shape_manual(name = "Visit",
                       values = setNames(rep_len(c(16, 17, 15, 18, 8, 7, 3, 4),
                                                 length(POST_VISITS)),
                                         as.character(POST_VISITS)),
                       drop = TRUE) +
    labs(title   = "Post minus pre ionized calcium, raw and PV-adjusted",
         caption = wrap_txt("Each line is one pre/post pair carried across the ",
                            "definitions. Sign changes indicate the direction of the ",
                            "effect is adjustment-dependent.")) +
    theme_report()
}

# ---- page 9: the inputs to the adjustment -----------------------------------
plot_hb_hct <- function() {
  d <- dat %>%
    select(pin_label, visit, phase, `ctHb (g/dL)` = hb, `Hct (%)` = hct) %>%
    pivot_longer(-c(pin_label, visit, phase),
                 names_to = "measure", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(measure = factor(measure, levels = c("ctHb (g/dL)", "Hct (%)")))

  ggplot(d, aes(visit, value, colour = phase)) +
    geom_line(data = d %>% filter(phase == "Pre-exercise"),
              aes(group = interaction(pin_label, measure)),
              colour = "grey35", linewidth = 0.4) +
    geom_point(size = 2.2) +
    facet_grid(measure ~ pin_label, scales = "free_y", switch = NULL) +
    scale_visit() +
    scale_y_continuous(name = NULL) +
    scale_colour_manual(name = "Draw", values = DRAW_COLS, drop = FALSE) +
    labs(title   = "Haemoconcentration markers: the inputs to the PV adjustment",
         caption = wrap_txt("Rows use separate scales because the units differ; ",
                            "within a row every participant panel shares one scale. ",
                            "The two rows are near-identical in shape because Hct ",
                            "tracks ctHb almost exactly in this file.")) +
    theme_report()
}

# ---- pages 10-11: sodium and chloride ---------------------------------------
plot_analyte <- function(var, title, subtitle = NULL, ylab) {
  d <- dat %>%
    mutate(value = .data[[var]]) %>%
    filter(!is.na(value))
  links <- d %>% semi_join(pairs_tbl %>% select(pin, visit), by = c("pin", "visit"))

  ggplot(d, aes(visit, value)) +
    layer_week_bands() +
    geom_line(data = d %>% filter(phase == "Pre-exercise"),
              aes(group = pin_label), colour = "grey35", linewidth = 0.4) +
    geom_line(data = links, aes(group = interaction(pin_label, visit)),
              colour = "grey35", linewidth = 0.4, linetype = "22") +
    geom_point(colour = "white", size = 3.2) +
    geom_point(aes(colour = phase), size = 2.2) +
    facet_wrap(~pin_label, nrow = 1) +
    scale_visit() +
    scale_y_continuous(name = ylab) +
    scale_colour_manual(name = "Draw", values = DRAW_COLS, drop = FALSE) +
    labs(title = title, subtitle = subtitle,
         caption = wrap_txt("Solid line: pre-exercise series. Dashed: same-day ",
                            "pre-to-post pair. All panels on this page share ",
                            "identical x and y limits.")) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
    theme_report()
}

# ---- page 12: dyshaemoglobin fractions --------------------------------------
plot_dyshb <- function() {
  d <- dat %>%
    select(pin_label, visit, phase, FCOHb = fcohb, FMetHb = fmethb) %>%
    pivot_longer(-c(pin_label, visit, phase),
                 names_to = "measure", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(measure = factor(measure, levels = c("FCOHb", "FMetHb")))

  ggplot(d, aes(visit, value, colour = phase)) +
    geom_line(data = d %>% filter(phase == "Pre-exercise"),
              aes(group = interaction(pin_label, measure)),
              colour = "grey35", linewidth = 0.4) +
    geom_point(size = 2.2) +
    facet_grid(measure ~ pin_label, scales = "free_y") +
    scale_visit() +
    scale_y_continuous(name = "% of total haemoglobin") +
    scale_colour_manual(name = "Draw", values = DRAW_COLS, drop = FALSE) +
    labs(title   = "Dyshaemoglobin fractions (FCOHb, FMetHb)",
         caption = wrap_txt("Shown because both are components of the ctHb on which ",
                            "the plasma-volume adjustment is built. Within a row, all ",
                            "participant panels share one scale.")) +
    theme_report()
}

# ---- optional pH pages ------------------------------------------------------
plot_ph <- function() {
  plot_analyte("ph", "Whole-blood pH by visit",
               "Reported for context; the acid-base response is part of the exercise effect",
               "pH")
}

plot_ica_ph <- function() {
  d <- dat %>%
    filter(!is.na(ica), !is.na(ph)) %>%
    mutate(`Measured`                        = ica,
           `Normalised to pH 7.40`           = ph_normalise(ica, ph)) %>%
    select(pin_label, visit, phase, Measured, `Normalised to pH 7.40`) %>%
    pivot_longer(-c(pin_label, visit, phase),
                 names_to = "series", values_to = "value") %>%
    mutate(series = factor(series, levels = c("Measured", "Normalised to pH 7.40")))

  ggplot(d, aes(visit, value, colour = series, shape = series)) +
    layer_week_bands() +
    geom_point(size = 2.2) +
    facet_wrap(~pin_label, nrow = 1) +
    scale_visit() +
    scale_y_continuous(name = "Ionized calcium (mmol/L)",
                       limits = ICA_YLIM, expand = c(0, 0)) +
    scale_colour_manual(name = NULL, values = c("Measured" = "#D95F02",
                                                "Normalised to pH 7.40" = "#2C7FB8")) +
    scale_shape_manual(name = NULL, values = c("Measured" = 16,
                                               "Normalised to pH 7.40" = 17)) +
    labs(title    = "Sensitivity analysis: iCa normalised to pH 7.40",
         subtitle = sprintf("%.0f%% change in iCa per 0.1 pH unit", PH_SLOPE * 100),
         caption  = wrap_txt("Not the primary readout. The in-vivo pH shift is part of ",
                             "the exercise response, so the measured value is what the ",
                             "calcium-sensing receptor experienced.")) +
    guides(fill = guide_legend(order = 1),
           colour = guide_legend(order = 2), shape = guide_legend(order = 2)) +
    theme_report()
}

# ---- optional model page ----------------------------------------------------
page_model <- function() {
  ok <- requireNamespace("lme4", quietly = TRUE)
  if (!ok) {
    text_page("Repeated-measures model",
              c("Not run", body("Package lme4 is not installed.")))
    return(invisible(NULL))
  }
  d <- dat %>% filter(!is.na(ica))
  fit <- try(lme4::lmer(ica ~ phase + (1 | pin), data = d), silent = TRUE)
  lines <- if (inherits(fit, "try-error"))
    c("Not estimable", body(as.character(fit)))
  else
    paste0("  ", capture.output(summary(fit)))
  text_page("Repeated-measures model: iCa ~ phase + (1 | PIN)", lines,
            size = 7, first_top = 78, lead = 12, fit = TRUE)
}


# =============================================================================
# 9. RENDER
# =============================================================================

pdf(OUT_PDF, width = PAGE_W, height = PAGE_H, onefile = TRUE, paper = "special")

page_scope()
page_checks()
page_pair_table()

print(plot_ica_facet())
print(plot_ica_pooled())
print(plot_dpv())
print(plot_post_adjusted())
print(plot_deltas())
print(plot_hb_hct())
print(plot_analyte("na", "Sodium (cNa+) by visit",
                   "Rises with plasma-volume contraction, but is also actively regulated -- not a clean volume marker",
                   "cNa+ (mmol/L)"))
print(plot_analyte("cl", "Chloride (cCl-) by visit", NULL, "cCl- (mmol/L)"))
print(plot_dyshb())

if (identical(PH_MODE, "sensitivity")) {
  print(plot_ph())
  print(plot_ica_ph())
}
if (FIT_MODEL) page_model()

invisible(dev.off())

message("Wrote ", normalizePath(OUT_PDF))
