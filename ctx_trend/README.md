# C-Telopeptide (CTX) trend by PIN and visit

`ctx_trend.py` pulls the serum C-Telopeptide column out of a Labcorp study-results
export in the "Study Data Entry Template" layout and produces the tables and plot
needed for a CTX trend.

## Input layout

The script expects the template as exported:

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

```bash
pip install pandas matplotlib
python ctx_trend.py --csv path/to/Labcorp_Results.csv --outdir out
python ctx_trend.py --csv path/to/Labcorp_Results.csv --outdir out --x days
```

`--x visit` (default) plots against visit number; `--x days` plots against days from
each participant's baseline draw, which spaces the points by real elapsed time.

## Outputs

- `out/ctx_by_pin_visit.csv` — tidy, one row per draw: `pin`, `visit`, `pre_post`,
  `draw`, `collection_date`, `ctx_pg_ml`, `baseline_pg_ml`, `change_from_baseline`,
  `pct_change_from_baseline`, `days_from_baseline`
- `out/ctx_pin_by_visit_wide.csv` — PIN × visit matrix of CTX values
- `out/ctx_trend.png` — one line per participant

Baseline is each participant's earliest CTX draw, not a fixed visit number, since
participants entered CTX collection at different visits.

## Functions

- `load_labs(csv_path)` → `(data, units, refs)` — strips the two metadata rows and returns
  them as `{column: value}` dicts
- `ctx_by_pin_visit(data)` → tidy long CTX frame
- `ctx_wide(ctx)` → PIN × visit matrix, columns ordered by visit then pre/post
- `plot_ctx_trend(ctx, out_path, x=..., units=..., ref_range=...)` → saved figure path

## Note on data

No participant data is committed to this repository — pass the export path with `--csv`.
