"""Helper functions for participant-level spectral trend analysis.

These functions do NOT touch the EEG signal processing. All of the filtering,
epoching, artifact rejection and multitaper spectral estimation happens in
``eeg_report.py`` (the pipeline written by the lab). This module only reads the
two CSV files that ``eeg_report.generate_reports(kind="export")`` writes out,
and turns them into per-participant trends, figures and statistics.

Read top to bottom - each function is short and does one thing.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats

# ---------------------------------------------------------------- constants
# Band definitions and colours are copied from eeg_report.py so that every
# figure in this analysis matches the lab's existing reports.
BANDS = ["delta", "theta", "alpha", "beta"]
BAND_HZ = {"delta": (0.5, 4), "theta": (4, 8), "alpha": (8, 13), "beta": (13, 30)}
BCOL = {"delta": "#378ADD", "theta": "#1D9E75", "alpha": "#D85A30", "beta": "#7F77DD"}

INK, MUTED, GRIDC = "#22262b", "#6b7280", "#d8dbe0"
MIN_EPOCHS = 20          # same quality threshold the lab's reports use (THR)

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 160, "savefig.bbox": "tight",
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.titlesize": 10, "axes.titleweight": "bold", "axes.labelsize": 9,
    "axes.edgecolor": GRIDC, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.7, "grid.alpha": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.fontsize": 8,
})


# ---------------------------------------------------------------- 1. loading
def load(folder="eeg_reports_output"):
    """Read the two CSVs written by the export step.

    Returns (sessions, epochs):
      sessions - one row per participant x visit x phase x condition x region x band
      epochs   - one row per single clean epoch (frontal No-Go only)
    """
    s = pd.read_csv(os.path.join(folder, "band_power_session.csv"))
    e_path = os.path.join(folder, "band_power_epochs.csv")
    e = pd.read_csv(e_path) if os.path.exists(e_path) else None
    for df in (s, e):
        if df is not None:
            df["participant"] = df["participant"].astype(int)
            df["visit"] = df["visit"].astype("Int64")
    return s, e


def inventory(s):
    """One row per participant: what was actually recorded and how clean it was.

    Run this FIRST. It tells you which participants and visits you really have
    before you interpret any trend.
    """
    f = s[(s.region == "frontal") & (s.condition == "NoGo") & (s.band == "delta")]
    out = []
    for pid, g in f.groupby("participant"):
        visits = sorted(g.visit.dropna().unique())
        have = g.pivot_table(index="visit", columns="phase", values="n_epochs", aggfunc="first")
        have = have.reindex(columns=["pre", "post"])
        paired = have.dropna()
        usable = paired[(paired.pre >= MIN_EPOCHS) & (paired.post >= MIN_EPOCHS)]
        out.append(dict(
            participant=pid,
            visits_recorded=len(visits),
            visit_range=f"V{int(min(visits))}-V{int(max(visits))}" if visits else "-",
            visits_with_pre_and_post=len(paired),
            visits_usable=len(usable),
            median_clean_epochs=float(np.nanmedian(g.n_epochs)) if len(g) else np.nan,
            median_pct_rejected=round(float(np.nanmedian(g.pct_rejected)), 1) if len(g) else np.nan,
        ))
    return pd.DataFrame(out).sort_values("participant").reset_index(drop=True)


# ---------------------------------------------------------------- 2. shaping
def prepare(s, region="frontal", condition="NoGo", min_epochs=MIN_EPOCHS):
    """Filter to one region/condition, drop low-quality sessions, add relative power.

    relative power = this band's share of the total 0.5-30 Hz power, in percent.
    It is computed only for sessions where all four bands survived the quality
    filter, so the four shares always add up to 100.
    """
    t = s[(s.region == region) & (s.condition == condition)].copy()
    t = t[t.n_epochs >= min_epochs]
    t = t[t.band.isin(BANDS)]

    key = ["participant", "visit", "phase"]
    complete = t.groupby(key).band.transform("nunique") == len(BANDS)
    t = t[complete].copy()

    tot = t.groupby(key).abs_power_uv2.transform("sum")
    t["rel_power_pct"] = np.where(tot > 0, t.abs_power_uv2 / tot * 100, np.nan)
    t["log2_abs"] = np.log2(t.abs_power_uv2.where(t.abs_power_uv2 > 0))
    t["band"] = pd.Categorical(t.band, categories=BANDS, ordered=True)
    return t.sort_values(["participant", "visit", "phase", "band"]).reset_index(drop=True)


def acute(t, value="abs_power_uv2"):
    """POST-vs-PRE change at each visit, per participant and band.

    log2_ratio = log2(POST / PRE).  0 means no change, +1 means POST is double
    PRE, -1 means POST is half PRE. Using log2 makes an increase and the matching
    decrease the same size, which raw differences do not.
    """
    w = t.pivot_table(index=["participant", "visit", "band"], columns="phase",
                      values=value, observed=True)
    w = w.reindex(columns=["pre", "post"]).dropna()
    w = w[(w.pre > 0) & (w.post > 0)]
    w = w.assign(log2_ratio=np.log2(w.post) - np.log2(w.pre),
                 percent_change=(w.post / w.pre - 1) * 100)
    return w.reset_index()


# ---------------------------------------------------------------- 3. figures
def _style_visit_axis(ax, visits):
    if len(visits):
        ax.set_xticks(list(visits))
    ax.tick_params(length=0)


def _log_y(ax):
    """Log y-axis with ordinary numbers on the ticks.

    Band power spans orders of magnitude and changes multiplicatively, so a log
    axis is the honest one. Matplotlib's default log labels are scientific
    notation ("2.5 x 10^0"), which is unreadable at these ranges - so the tick
    labels are forced back to plain numbers.
    """
    ax.set_yscale("log")
    for axis in (ax.yaxis,):
        axis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        axis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)


def _phase_legend(fig, color="#6b7280", ncol=2, y=0.0):
    """One legend for the whole figure instead of the same key in every panel."""
    handles = [
        matplotlib.lines.Line2D([], [], color=color, ls="--", marker="o", ms=6, lw=1.8,
                                markerfacecolor="white", markeredgecolor=color,
                                markeredgewidth=1.6, label="PRE (before exercise)"),
        matplotlib.lines.Line2D([], [], color=color, ls="-", marker="s", ms=6, lw=1.8,
                                markerfacecolor=color, markeredgecolor=color,
                                markeredgewidth=1.6, label="POST (after exercise)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=ncol, bbox_to_anchor=(0.5, y),
               frameon=False, fontsize=9)


def plot_participant_trends(t, pid, value="abs_power_uv2", save_dir=None):
    """One participant: each band's power at every visit, PRE vs POST.

    Four panels (delta / theta / alpha / beta). Dashed line = PRE, solid = POST,
    matching the lab's existing reports. Identity is carried by line style as
    well as colour, so the figure still reads in greyscale.
    """
    g = t[t.participant == pid]
    if not len(g):
        print(f"P{pid}: no usable sessions after quality filtering"); return None
    rel = value == "rel_power_pct"
    ylab = "share of total power (%)" if rel else "power (uV^2)"
    visits = sorted(g.visit.dropna().unique())

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 6.6), sharex=True)
    for ax, b in zip(axes.ravel(), BANDS):
        sb = g[g.band == b]
        for phase, ls, marker in [("pre", "--", "o"), ("post", "-", "s")]:
            p = sb[sb.phase == phase].sort_values("visit")
            if not len(p):
                continue
            ax.plot(p.visit, p[value], ls, marker=marker, ms=6, lw=1.8,
                    color=BCOL[b], markerfacecolor="white" if phase == "pre" else BCOL[b],
                    markeredgecolor=BCOL[b], markeredgewidth=1.6, label=phase.upper())
        lo, hi = BAND_HZ[b]
        ax.set_title(f"{b}  ({lo}-{hi} Hz)", color=BCOL[b])
        ax.set_ylabel(ylab)
        if not rel:
            _log_y(ax)                      # band power is strongly right-skewed
        _style_visit_axis(ax, visits)
    for ax in axes[-1]:
        ax.set_xlabel("visit")
    fig.suptitle(f"Participant {pid} - frontal No-Go band power across visits"
                 f"{' (relative)' if rel else ' (absolute, log scale)'}",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    _phase_legend(fig, y=0.005)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        fig.savefig(os.path.join(save_dir, f"P{pid}_trends_{'rel' if rel else 'abs'}.png"))
    return fig


def plot_band_across_participants(t, band, value="abs_power_uv2", save_dir=None):
    """One band, every participant - a small multiple per participant.

    Small multiples rather than nine overlaid lines: nine cycled colours are not
    distinguishable, and each participant's power is on a different scale anyway.
    """
    pids = sorted(t.participant.unique())
    if not pids:
        print("nothing to plot"); return None
    rel = value == "rel_power_pct"
    ncol = min(3, len(pids)); nrow = int(np.ceil(len(pids) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 2.8 * nrow), squeeze=False)
    for ax, pid in zip(axes.ravel(), pids):
        g = t[(t.participant == pid) & (t.band == band)]
        for phase, ls, marker in [("pre", "--", "o"), ("post", "-", "s")]:
            p = g[g.phase == phase].sort_values("visit")
            if not len(p):
                continue
            ax.plot(p.visit, p[value], ls, marker=marker, ms=5, lw=1.6,
                    color=BCOL[band], markerfacecolor="white" if phase == "pre" else BCOL[band],
                    markeredgecolor=BCOL[band], markeredgewidth=1.4)
        ax.set_title(f"P{pid}", fontsize=9)
        ax.tick_params(labelsize=7, length=0)
        if not rel:
            _log_y(ax)
    for ax in axes.ravel()[len(pids):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("visit", fontsize=8)
    for r in range(nrow):
        axes[r, 0].set_ylabel("share (%)" if rel else "power (uV^2)", fontsize=8)
    lo, hi = BAND_HZ[band]
    fig.suptitle(f"{band} ({lo}-{hi} Hz) - frontal No-Go band power across visits"
                 f"{' (relative)' if rel else ' (absolute, log scale)'}\n"
                 f"each panel has its own y-scale - compare the shape of a trend, not its height",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    _phase_legend(fig, color=BCOL[band], y=0.005)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        fig.savefig(os.path.join(save_dir, f"band_{band}_{'rel' if rel else 'abs'}_all.png"))
    return fig


def plot_acute_trends(a, pid, save_dir=None):
    """One participant: does the POST-minus-PRE effect itself drift across visits?

    Each point is one visit. The line is a least-squares fit, shown only as a
    visual guide - the reported statistic is Spearman's rho in trend_table().
    """
    g = a[a.participant == pid]
    if not len(g):
        print(f"P{pid}: no paired PRE+POST visits"); return None
    visits = sorted(g.visit.dropna().unique())
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 6.0), sharex=True, sharey=True)
    for ax, b in zip(axes.ravel(), BANDS):
        sb = g[g.band == b].sort_values("visit")
        ax.axhline(0, color=MUTED, lw=1, ls=":")
        ax.plot(sb.visit, sb.log2_ratio, "o", ms=7, color=BCOL[b],
                markeredgecolor="white", markeredgewidth=1.4, label="visit")
        if len(sb) >= 2 and sb.visit.nunique() > 1:
            x = sb.visit.astype(float).to_numpy(); y = sb.log2_ratio.to_numpy()
            m, c = np.polyfit(x, y, 1)
            xs = np.linspace(x.min(), x.max(), 20)
            ax.plot(xs, m * xs + c, "-", lw=1.8, color=BCOL[b], alpha=0.55)
        lo, hi = BAND_HZ[b]
        ax.set_title(f"{b}  ({lo}-{hi} Hz)", color=BCOL[b])
        _style_visit_axis(ax, visits)
    for ax in axes[:, 0]:
        ax.set_ylabel("log2(POST / PRE)")
    for ax in axes[-1]:
        ax.set_xlabel("visit")
    fig.suptitle(f"Participant {pid} - acute exercise effect across visits (frontal No-Go)\n"
                 f"each dot is one visit; line is a least-squares guide.  "
                 f"above 0 = power rose after exercise, below 0 = power fell",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        fig.savefig(os.path.join(save_dir, f"P{pid}_acute.png"))
    return fig


# ---------------------------------------------------------------- 4. statistics
def _holm(pvals):
    """Holm-Bonferroni step-down correction for multiple testing."""
    p = np.asarray(pvals, float)
    ok = np.isfinite(p)
    adj = np.full(p.shape, np.nan)
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    m = len(order)
    running = 0.0
    for rank, i in enumerate(order):
        val = (m - rank) * p[i]
        running = max(running, val)
        adj[i] = min(1.0, running)
    return adj


def trend_table(a, min_visits=5):
    """Per participant x band, the two questions you actually want answered.

    1. Is there an acute effect at all?  Wilcoxon signed-rank on log2(POST/PRE)
       against zero, across visits. Non-parametric, so it does not assume the
       effect is normally distributed.
    2. Does that effect trend across the study?  Spearman's rho between visit
       number and log2(POST/PRE). Non-parametric, so a single odd visit cannot
       drag the trend around the way a Pearson correlation would.

    p-values are Holm-corrected across all participant x band tests in the table.
    """
    rows = []
    for (pid, b), g in a.groupby(["participant", "band"], observed=True):
        g = g.sort_values("visit")
        y = g.log2_ratio.to_numpy(float)
        x = g.visit.astype(float).to_numpy()
        n = len(y)
        rec = dict(participant=pid, band=str(b), n_visits=n,
                   median_log2_ratio=float(np.median(y)) if n else np.nan,
                   median_percent_change=float(np.median(g.percent_change)) if n else np.nan,
                   wilcoxon_p=np.nan, spearman_rho=np.nan, spearman_p=np.nan)
        if n >= min_visits and np.ptp(y) > 0:
            rec["wilcoxon_p"] = float(stats.wilcoxon(y).pvalue)
        if n >= min_visits and np.ptp(x) > 0 and np.ptp(y) > 0:
            r = stats.spearmanr(x, y)
            rec["spearman_rho"] = float(r.statistic); rec["spearman_p"] = float(r.pvalue)
        rows.append(rec)
    d = pd.DataFrame(rows)
    if len(d):
        d["wilcoxon_p_holm"] = _holm(d.wilcoxon_p)
        d["spearman_p_holm"] = _holm(d.spearman_p)
        d["band"] = pd.Categorical(d.band, categories=BANDS, ordered=True)
        d = d.sort_values(["participant", "band"]).reset_index(drop=True)
    return d


def friedman_across_visits(t, band, phase="pre", value="abs_power_uv2", min_participants=3):
    """Group-level Friedman test: does band power differ across visits?

    Blocks are participants, treatments are visits. Friedman needs a COMPLETE
    block design, so this keeps only the visits every remaining participant has,
    and reports exactly what was kept. Read `n_participants` and `visits` before
    reading the p-value - with dropout the test often runs on very few blocks.
    """
    g = t[(t.band == band) & (t.phase == phase)]
    m = g.pivot_table(index="participant", columns="visit", values=value, observed=True)
    best = None
    for min_v in range(m.shape[1], 2, -1):
        for cols in [m.columns[i:i + min_v] for i in range(m.shape[1] - min_v + 1)]:
            block = m[list(cols)].dropna()
            if len(block) >= min_participants:
                best = block; break
        if best is not None:
            break
    if best is None or best.shape[1] < 3:
        return dict(band=band, phase=phase, n_participants=0, n_visits=0, visits=[],
                    statistic=np.nan, p=np.nan,
                    note="no complete block of >=3 visits shared by >=3 participants")
    stat, p = stats.friedmanchisquare(*[best[c].to_numpy() for c in best.columns])
    return dict(band=band, phase=phase, n_participants=int(best.shape[0]),
                n_visits=int(best.shape[1]), visits=[int(c) for c in best.columns],
                statistic=float(stat), p=float(p),
                note="blocks = participants, treatments = visits (complete cases only)")


def wilcoxon_pre_post_group(t, band, value="abs_power_uv2"):
    """Group-level PRE vs POST: one value per participant (median across visits).

    With only two related conditions a Friedman test reduces to a sign test, so
    the Wilcoxon signed-rank test is the correct choice here.
    """
    g = t[t.band == band]
    m = g.pivot_table(index="participant", columns="phase", values=value, aggfunc="median", observed=True)
    m = m.reindex(columns=["pre", "post"]).dropna()
    if len(m) < 5:
        return dict(band=band, n_participants=len(m), statistic=np.nan, p=np.nan,
                    note="fewer than 5 participants with both PRE and POST - not tested")
    r = stats.wilcoxon(m.post.to_numpy(), m.pre.to_numpy())
    return dict(band=band, n_participants=int(len(m)),
                median_pre=float(m.pre.median()), median_post=float(m.post.median()),
                statistic=float(r.statistic), p=float(r.pvalue),
                note="paired across participants, median over each participant's visits")
