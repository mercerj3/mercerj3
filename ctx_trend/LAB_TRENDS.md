# Per-lab trend graphs

`lab_trends.R` builds one trend figure per lab from a `labcorp` data frame already
loaded in your session — a panel per participant, all panels on the same axes.

## Usage

```r
install.packages(c("dplyr", "tidyr", "ggplot2"))

labcorp <- readxl::read_excel("labcorp.xlsx")   # or read.csv(), read_csv(), ...
source("lab_trends.R")

long <- labs_long(labcorp)                      # tidy: one row per draw per lab

plot_lab_trend(long, "Hemoglobin")              # a single lab
plots <- plot_all_labs(long)                    # named list, one ggplot per lab
save_lab_trends(long, "lab_trends")             # multi-page PDF + a PNG per lab
```

On the current export this yields **7,234 results across 74 labs, 14 participants,
visits 1–25**. `save_lab_trends()` writes `lab_trends.pdf` with one lab per page —
usually the artefact you want, since 74 separate PNGs are awkward to flip through.

## Axes are held identical

- **y** — one scale per lab, shared by every participant panel, spanning that lab's
  full observed range across all participants. A panel is therefore directly
  comparable to the one beside it.
- **x** — the same visit range on every panel of every lab, taken from all visits
  observed anywhere in the data, so figures line up lab-to-lab as well.
- Limits are applied with `coord_cartesian()`, which clips rather than drops, and the
  y range is padded 4% so a participant sitting at the lab's min or max is not cut in
  half.

## Pre/post draws — two schemes

The protocol changed partway through, and the two schemes mean different things:

| | PINs 001–006 (`Pre_Post_Old`) | PINs 007+ (`Pre_Post_New`) |
| --- | --- | --- |
| draws per visit | one | two on visits 8, 13, 18, 23 |
| what the flag marks | the **whole visit** — 4/9/14/19 pre, 8/13/18/23 post | the **draw within** a visit |

Both are coalesced into one `pre_post` column, with `scheme` recording which column
the participant's flag came from and `paired_visit` marking the visits that hold two
draws. Because a paired visit puts both draws at the same visit number, the post point
is nudged right by `nudge_post` (default `0.35`) so the pair stays legible — the short
segment joining them is the within-visit change. Set `nudge_post = 0` to stack them.

## Data handling

- **Column names** are matched on alphanumerics only, so `Bilirubin,Total`,
  `Bilirubin.Total` (`read.csv`), and `bilirubin_total` (`janitor`) all resolve.
- **Lab detection** is automatic: every non-identifier column that parses to at least
  one number. `LDIsoenzymes:` is empty and drops out; `Notes` is excluded by name.
- **`Visit_Number` arrives mixed** numeric/text (PIN 009's later visits are stored as
  text), so everything is coerced through character.
- **Censored results** (`<0.2`, ten of them in Bilirubin) are kept at the detection
  limit and flagged in `censored`, drawn as hollow points.
- **Lab comment codes** (`WF` in the LD fractions) parse to NA and drop out.
- **Metadata rows** are dropped if present. The xlsx has none; the CSV version of this
  template carries units and reference intervals in its first two rows.
- **Row order is protocol order** (visit, then pre/post), not collection date — see the
  data note below.
- **Units and male reference intervals** are held in `LAB_UNITS` and `LAB_REF_RANGES`,
  lifted from the CSV template's header rows since the xlsx drops them. Units label the
  y axis; the reference interval is shaded, clipped so it never widens the shared axis.
  Turn the band off with `ref_band = FALSE`.

## Data note

PIN 007, visit 23: the **post** draw is dated 2026-04-22, two days *before* its own pre
draw (2026-04-24). Every other paired visit has both draws on one date. Ordering by
collection date would reverse that pre/post pair, so the tidy table and the plots order
by visit and pre/post instead. Worth checking against source records.

## Functions

- `labs_long(labcorp, labs = NULL)` → tidy frame: `pin`, `visit`, `pre_post`, `draw`,
  `scheme`, `paired_visit`, `collection_date`, `lab`, `value`, `censored`
- `lab_columns(labcorp)` → the lab columns detected
- `plot_lab_trend(long, lab, ...)` → one ggplot
- `plot_all_labs(long, ...)` → named list of ggplots on a common x axis
- `save_lab_trends(long, outdir, png = TRUE, ...)` → writes the PDF and PNGs

Sourcing `lab_trends.R` without a `labcorp` object just defines the functions.
