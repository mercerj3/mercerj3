"""C-Telopeptide (CTX) trend by PIN and visit.

Reads a Labcorp study-results export in the "Study Data Entry Template" layout
and pulls the serum C-Telopeptide column into a tidy long table keyed by
participant (PIN) and visit, then writes a wide PIN x Visit table and a trend
plot.

Layout assumed (matches Study_Data_Entry_Template_..._Labcorp_Results.csv):
    row 1  column headers
    row 2  units             (C-Telopeptide,Serum -> pg/mL)
    row 3  reference intervals
    row 4+ one row per participant / visit / draw

Usage:
    python ctx_trend.py --csv path/to/Labcorp_Results.csv --outdir out
    python ctx_trend.py --csv path/to/Labcorp_Results.csv --x days
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

CTX_COL = "C-Telopeptide,Serum"
UNITS_ROW = 0        # 0-based position within the two non-data header rows
REF_ROW = 1
ID_COLS = ["PIN", "Visit_Number", "Pre_Post_Old", "Pre_Post_New", "Collection_Date"]


def load_labs(csv_path: str | Path) -> tuple[pd.DataFrame, dict[str, str], dict[str, str]]:
    """Return (data, units, reference_intervals) from the template CSV.

    The two metadata rows under the header are stripped off and returned as
    {column -> value} dicts so the CTX units/ref range can be used in labels.
    """
    raw = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    units = raw.iloc[UNITS_ROW].to_dict()
    refs = raw.iloc[REF_ROW].to_dict()
    data = raw.iloc[REF_ROW + 1 :].copy()

    # Drop the blank filler rows at the bottom of the template.
    data = data[data["PIN"].str.strip() != ""]
    return data.reset_index(drop=True), units, refs


def ctx_by_pin_visit(data: pd.DataFrame) -> pd.DataFrame:
    """Tidy one-row-per-draw CTX table with baseline-relative columns.

    Participants 001-006 record the pre/post draw flag in ``Pre_Post_Old`` and
    007+ in ``Pre_Post_New``; the two are coalesced into a single ``pre_post``.
    A visit can appear twice for the same PIN (a pre and a post draw), so the
    natural key is (pin, visit, pre_post).
    """
    ctx = pd.DataFrame(
        {
            "pin": data["PIN"].str.strip(),
            "visit": pd.to_numeric(data["Visit_Number"], errors="coerce").astype("Int64"),
            "pre_post": pd.to_numeric(
                data["Pre_Post_New"].where(data["Pre_Post_New"].str.strip() != "", data["Pre_Post_Old"]),
                errors="coerce",
            ).astype("Int64"),
            "collection_date": pd.to_datetime(data["Collection_Date"], errors="coerce"),
            "ctx_pg_ml": pd.to_numeric(data[CTX_COL], errors="coerce"),
        }
    )

    # Keep only draws that actually have a CTX result.
    ctx = ctx.dropna(subset=["ctx_pg_ml"]).sort_values(["pin", "visit", "pre_post"])

    ctx["draw"] = ctx["pre_post"].map({0: "pre", 1: "post"}).fillna("single")
    ctx["visit_label"] = [
        f"V{v}" if d == "single" else f"V{v} ({d})" for v, d in zip(ctx["visit"], ctx["draw"])
    ]

    # Baseline = each participant's earliest CTX draw.
    baseline_value = ctx.groupby("pin")["ctx_pg_ml"].transform("first")
    baseline_date = ctx.groupby("pin")["collection_date"].transform("first")
    ctx["baseline_pg_ml"] = baseline_value
    ctx["change_from_baseline"] = ctx["ctx_pg_ml"] - baseline_value
    ctx["pct_change_from_baseline"] = 100 * ctx["change_from_baseline"] / baseline_value
    ctx["days_from_baseline"] = (ctx["collection_date"] - baseline_date).dt.days

    return ctx.reset_index(drop=True)


def ctx_wide(ctx: pd.DataFrame) -> pd.DataFrame:
    """PIN x visit matrix of CTX values (rows = PIN, columns = visit label)."""
    wide = ctx.pivot_table(
        index="pin", columns="visit_label", values="ctx_pg_ml", aggfunc="first"
    )
    # Order columns by visit number, then pre before post.
    order = (
        ctx[["visit_label", "visit", "pre_post"]]
        .drop_duplicates()
        .sort_values(["visit", "pre_post"])["visit_label"]
    )
    return wide.reindex(columns=[c for c in order if c in wide.columns])


def plot_ctx_trend(
    ctx: pd.DataFrame,
    out_path: str | Path,
    x: str = "visit",
    units: str = "pg/mL",
    ref_range: str = "",
) -> Path:
    """Line plot of CTX over time, one series per PIN."""
    x_col = {"visit": "visit", "days": "days_from_baseline"}[x]
    x_label = {"visit": "Visit number", "days": "Days from baseline draw"}[x]

    fig, ax = plt.subplots(figsize=(9, 5.5))
    for pin, grp in ctx.groupby("pin"):
        grp = grp.sort_values(x_col)
        ax.plot(grp[x_col], grp["ctx_pg_ml"], marker="o", linewidth=1.8, label=f"PIN {pin}")

    if ref_range and "-" in ref_range:
        low, high = (float(v) for v in ref_range.split("-", 1))
        ax.axhspan(low, high, color="0.85", zorder=0, label=f"Reference {ref_range} {units}")

    ax.set_xlabel(x_label)
    ax.set_ylabel(f"C-Telopeptide, serum ({units})")
    ax.set_title("C-Telopeptide trend by participant and visit")
    if x == "visit":
        ax.set_xticks(sorted(ctx["visit"].dropna().unique()))
    ax.grid(axis="y", alpha=0.3)
    ax.legend(title="Participant", frameon=False)
    fig.tight_layout()

    out_path = Path(out_path)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", required=True, help="Labcorp results CSV in the study template layout")
    parser.add_argument("--outdir", default="out", help="directory for the CSV/PNG outputs")
    parser.add_argument("--x", choices=["visit", "days"], default="visit", help="x-axis for the trend plot")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    data, units, refs = load_labs(args.csv)
    ctx = ctx_by_pin_visit(data)
    wide = ctx_wide(ctx)

    long_path = outdir / "ctx_by_pin_visit.csv"
    wide_path = outdir / "ctx_pin_by_visit_wide.csv"
    ctx.to_csv(long_path, index=False)
    wide.to_csv(wide_path)

    plot_path = plot_ctx_trend(
        ctx,
        outdir / "ctx_trend.png",
        x=args.x,
        units=units.get(CTX_COL, "pg/mL"),
        ref_range=refs.get(CTX_COL, ""),
    )

    print(f"CTX results: {len(ctx)} draws across {ctx['pin'].nunique()} participants\n")
    print(
        ctx[["pin", "visit", "draw", "collection_date", "ctx_pg_ml", "pct_change_from_baseline"]].to_string(
            index=False
        )
    )
    print("\nPIN x visit matrix (pg/mL):")
    print(wide.to_string())
    print(f"\nWrote {long_path}\n      {wide_path}\n      {plot_path}")


if __name__ == "__main__":
    main()
