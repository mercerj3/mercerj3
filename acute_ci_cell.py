# ---- Acute effect WITH uncertainty, one PNG per band ------------------------
# Needs `trends`, `epochs` and `FIG_DIR`, which you already have.
#
# WHERE THE UNCERTAINTY COMES FROM
#   Each session's band power is an average over clean epochs. Those epochs vary,
#   so the average has a margin of error. This resamples the epochs (a bootstrap:
#   draw n epochs at random WITH replacement, recompute the average, repeat 2000
#   times) and reports the middle 95% of the resulting log2(POST/PRE) values.
#   No bell-curve assumption is made, which matters because band power is skewed.
#
# WHAT IT DOES *NOT* COVER  <-- read this before quoting any interval
#   This is WITHIN-session precision only: how well that day's power was measured
#   given epoch-to-epoch noise. It does NOT include day-to-day variation from
#   electrode placement, gel, hair or hydration, which is usually much larger.
#   So these bars are optimistic. Treat a bar excluding 0 as "measured precisely
#   that day", NOT as "a real effect" - for that, use the visit-to-visit tests in
#   Step 9, which do capture between-session variability.
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np, pandas as pd, os

BANDS = ["delta", "theta", "alpha", "beta"]
BCOL = {"delta": "#378ADD", "theta": "#1D9E75", "alpha": "#D85A30", "beta": "#7F77DD"}
BAND_HZ = {"delta": "0.5-4 Hz", "theta": "4-8 Hz", "alpha": "8-13 Hz", "beta": "13-30 Hz"}
MUTED = "#8a8f98"


def acute_with_ci(trends, epochs, n_boot=2000, seed=0):
    """log2(POST/PRE) per participant x visit x band, with a bootstrap 95% CI."""
    rng = np.random.default_rng(seed)
    keep = set(map(tuple, trends[["participant", "visit", "phase"]].drop_duplicates().values))
    e = epochs[[tuple(r) in keep for r in
                epochs[["participant", "visit", "phase"]].values]]
    out = []
    for (pid, visit, band), g in e.groupby(["participant", "visit", "band"]):
        pre = g.loc[g.phase == "pre", "abs_power_uv2"].to_numpy(float)
        post = g.loc[g.phase == "post", "abs_power_uv2"].to_numpy(float)
        pre = pre[np.isfinite(pre) & (pre > 0)]
        post = post[np.isfinite(post) & (post > 0)]
        if len(pre) < 2 or len(post) < 2:
            continue
        point = np.log2(post.mean()) - np.log2(pre.mean())
        bp = pre[rng.integers(0, len(pre), size=(n_boot, len(pre)))].mean(1)
        bq = post[rng.integers(0, len(post), size=(n_boot, len(post)))].mean(1)
        draws = np.log2(bq) - np.log2(bp)
        lo, hi = np.percentile(draws, [2.5, 97.5])
        out.append(dict(participant=int(pid), visit=int(visit), band=str(band),
                        log2_ratio=float(point), ci_lo=float(lo), ci_hi=float(hi),
                        ci_width=float(hi - lo), n_pre=len(pre), n_post=len(post)))
    d = pd.DataFrame(out)
    if len(d):
        d["band"] = pd.Categorical(d.band, categories=BANDS, ordered=True)
        d = d.sort_values(["participant", "band", "visit"]).reset_index(drop=True)
    return d


def plot_acute_band(a_ci, band, save_dir=None, max_ticks=8, show_ci=True, clip=None):
    """One band, every participant, one PNG. Error bars = bootstrap 95% CI."""
    pids = sorted(a_ci.participant.unique())
    if not pids:
        print("nothing to plot"); return None
    col = BCOL[band]
    ncol = min(3, len(pids)); nrow = int(np.ceil(len(pids) / ncol))
    vis = sorted(int(v) for v in a_ci.visit.dropna().unique())
    step = max(1, int(np.ceil(len(vis) / max_ticks)))
    ticks = vis[::step]

    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.0 * nrow),
                             squeeze=False, sharex=True, sharey=True)
    for ax, pid in zip(axes.ravel(), pids):
        d = a_ci[(a_ci.participant == pid) & (a_ci.band == band)].sort_values("visit")
        ax.axhline(0, color=MUTED, lw=1, ls=":")
        if len(d):
            if show_ci:
                ax.errorbar(d.visit, d.log2_ratio,
                            yerr=[d.log2_ratio - d.ci_lo, d.ci_hi - d.log2_ratio],
                            fmt="none", ecolor=col, elinewidth=1.2, capsize=2.5, alpha=.75)
            ax.plot(d.visit, d.log2_ratio, "o", ms=5, color=col,
                    markeredgecolor="white", markeredgewidth=1.0)
            if len(d) >= 2 and d.visit.nunique() > 1:
                x = d.visit.astype(float).to_numpy(); y = d.log2_ratio.to_numpy()
                m, b = np.polyfit(x, y, 1)
                xs = np.linspace(x.min(), x.max(), 20)
                ax.plot(xs, m * xs + b, "-", lw=1.8, color=col, alpha=.55)
        else:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center",
                    va="center", fontsize=8, color=MUTED)
        ax.set_title(f"P{pid}", fontsize=10, fontweight="bold")
        ax.set_xticks(ticks)
        ax.set_xlim(min(vis) - 0.7, max(vis) + 0.7)
        if clip:
            ax.set_ylim(-clip, clip)
        ax.tick_params(labelsize=8, length=0)
        ax.grid(alpha=.4)
    for ax in axes.ravel()[len(pids):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("visit", fontsize=9)
    for r in range(nrow):
        axes[r, 0].set_ylabel("log2(POST / PRE)", fontsize=9)

    fig.suptitle(f"{band} ({BAND_HZ[band]}) - acute exercise effect across visits, frontal No-Go\n"
                 f"above 0 = power rose after exercise, below 0 = it fell"
                 f"{'  |  bars = bootstrap 95% CI across epochs (within-session only)' if show_ci else ''}",
                 fontsize=11.5, fontweight="bold")
    fig.legend(handles=[Line2D([], [], ls=":", color=MUTED, label="no change (0)")],
               loc="lower center", frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.004))
    fig.tight_layout(rect=(0, 0.04, 1, 1))

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        p = os.path.join(save_dir, f"acute_{band}_all.png")
        fig.savefig(p, dpi=160, facecolor="white")
        print("saved", p)
    return fig


# ---- run it -----------------------------------------------------------------
acute_ci = acute_with_ci(trends, epochs)
acute_ci.to_csv(os.path.join(OUT_FOLDER, "acute_with_ci.csv"), index=False)
print(f"{len(acute_ci)} participant x visit x band estimates, each with a 95% CI")
print(f"median CI width: {acute_ci.ci_width.median():.3f} log2 units")

for band in BANDS:
    plot_acute_band(acute_ci, band, save_dir=FIG_DIR)
