# C-Telopeptide (CTX) trend by PIN and visit

`ctx_trend.py` and `ctx_trend.R` pull the serum C-Telopeptide column out of a Labcorp
study-results export in the "Study Data Entry Template" layout and produce the tables
and plot needed for a CTX trend.

The R version works on a data frame already loaded in your session, named `study`;
the Python version is a command-line script that reads the CSV itself.

## Input layout

Both scripts expect the template as exported:

| row | contents |
| --- | --- |
| 1 | column headers (`PIN`, `Visit_Number`, `Pre_Post_Old`, `Pre_Post_New`, `Collection_Date`, … `C-Telopeptide,Serum`, …) |
| 2 | units — CTX is `pg/mL` |
| 3 | reference intervals (blank for CTX in the current export) |
| 4+ | one row per participant / visit / draw; trailing blank rows are dropped |

Two details the parser handles explicitly:

- **Two pre/post columns.** PINs 001–006 use `Pre_Post_Old`, 007+ use `Pre_Post_New`
  (both `0=pre`, `1=post`). They are coalesced into a single `pre_post` field, so the
  natural key is `(pin, visit, pre_post)` — a visit can appear twice for one PIN.
- **Sparse CTX.** Most rows have no CTX result. Only rows with a value are kept, so the
  trend is built from actual draws rather than empty cells.

## Usage

### R — operates on `study`

Load the export however you normally would, then source the script:

```r
install.packages(c("dplyr", "tidyr", "ggplot2"))

study <- read.csv("path/to/Labcorp_Results.csv")   # or readr::read_csv(), etc.
source("ctx_trend.R")
```

Sourcing with `study` in the session prints the tidy table and the PIN x visit matrix,
and draws the plot. The objects it leaves behind are `ctx` (long), `wide`, and
`ctx_plot`. To work with them directly instead:

```r
ctx  <- ctx_by_pin_visit(study)
wide <- ctx_wide(ctx)
plot_ctx_trend(ctx)                 # x = visit number
plot_ctx_trend(ctx, x = "days")     # x = days from baseline draw
```

The script does not care how `study` was read in:

- **Column names** are matched on alphanumerics only, so `C-Telopeptide,Serum`
  (`readr::read_csv`), `C.Telopeptide.Serum` (`read.csv`), and `c_telopeptide_serum`
  (`janitor::clean_names`) all resolve.
- **The two metadata rows** are dropped if present — they arrive as data rows in `study`
  under any reader — and it is equally fine if you already removed them.
- **Column types** are routed through character first, so a reader that guessed types
  changes nothing.

### Python — command line

```bash
pip install pandas matplotlib
python ctx_trend.py --csv path/to/Labcorp_Results.csv --outdir out
python ctx_trend.py --csv path/to/Labcorp_Results.csv --outdir out --x days
```

`visit` (default) plots against visit number; `days` plots against days from each
participant's baseline draw, which spaces the points by real elapsed time.

## Outputs

- **tidy, one row per draw** — `pin`, `visit`, `pre_post`, `draw`, `collection_date`,
  `ctx_pg_ml`, `baseline_pg_ml`, `change_from_baseline`, `pct_change_from_baseline`,
  `days_from_baseline`
- **PIN × visit matrix** of CTX values
- **trend plot** — one line per participant

Python writes these to `--outdir`; R returns them as `ctx`, `wide`, and `ctx_plot`, with
commented `write.csv()` / `ggsave()` lines at the bottom of the script if you want files.

Baseline is each participant's earliest CTX draw, not a fixed visit number, since
participants entered CTX collection at different visits.

## Functions

- `ctx_by_pin_visit(study)` → tidy long CTX frame (R) / `ctx_by_pin_visit(data)` (Python)
- `ctx_wide(ctx)` → PIN × visit matrix, columns ordered by visit then pre/post
- `plot_ctx_trend(ctx, ...)` → a ggplot object (R) / saved figure path (Python)
- `lab_metadata(study, name)` → units and reference interval read off the template's
  metadata rows (R); `load_labs(csv_path)` returns the same, plus the data (Python)

Sourcing `ctx_trend.R` without a `study` object in the session just defines the
functions, so it is safe to source from an Rmd or another script.

## Note on data

No participant data is committed to this repository — load it locally into `study` (R) or
pass the export path with `--csv` (Python).
