# ---- One band, all participants - every panel on the SAME x-axis ------------
# Paste-and-run. Uses only `trends` and `FIG_DIR`, which already exist.
# Nothing to install, no file to replace, no kernel restart.
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np, os

BCOL = {"delta": "#378ADD", "theta": "#1D9E75", "alpha": "#D85A30", "beta": "#7F77DD"}
BAND_HZ = {"delta": "0.5-4 Hz", "theta": "4-8 Hz", "alpha": "8-13 Hz", "beta": "13-30 Hz"}


def plot_band_grid(trends, band, value="rel_power_pct", save_dir=None, max_ticks=8,
                   share_y=False):
    """All participants for one band. Every panel gets an identical x-axis."""
    d_all = trends[trends.phase.isin(["pre", "post"])]
    pids = sorted(d_all.participant.unique())
    if not pids:
        print("nothing to plot"); return None
    rel = value == "rel_power_pct"
    col = BCOL[band]
    ncol = min(3, len(pids)); nrow = int(np.ceil(len(pids) / ncol))

    # ONE visit range and ONE set of whole-number ticks, shared by every panel
    vis = sorted(int(v) for v in d_all.visit.dropna().unique())
    step = max(1, int(np.ceil(len(vis) / max_ticks)))
    ticks = vis[::step]

    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.0 * nrow),
                             squeeze=False, sharex=True, sharey=share_y)
    for ax, pid in zip(axes.ravel(), pids):
        d = d_all[(d_all.participant == pid) & (d_all.band == band)]
        for phase, ls, mk, mfc in [("pre", "--", "o", "white"), ("post", "-", "s", col)]:
            g = d[d.phase == phase].sort_values("visit")
            if len(g):
                ax.plot(g.visit, g[value], ls=ls, marker=mk, ms=5, lw=1.6, color=col,
                        markerfacecolor=mfc, markeredgecolor=col, markeredgewidth=1.4)
        if not len(d):
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center",
                    va="center", fontsize=8, color="#8a8f98")
        ax.set_title(f"P{pid}", fontsize=10, fontweight="bold")
        ax.set_xticks(ticks)                      # whole visit numbers, never 7.5
        ax.set_xlim(min(vis) - 0.7, max(vis) + 0.7)
        if not rel:
            ax.set_yscale("log")
        ax.tick_params(labelsize=8, length=0)
        ax.grid(alpha=.45)

    for ax in axes.ravel()[len(pids):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("visit", fontsize=9)
    for r in range(nrow):
        axes[r, 0].set_ylabel("share of total (%)" if rel else "power (uV^2)", fontsize=9)

    ynote = "all panels share one y-scale" if share_y else "each panel has its own y-scale"
    fig.suptitle(f"{band} ({BAND_HZ[band]}) - frontal No-Go, "
                 f"{'relative' if rel else 'absolute'} band power across visits\n"
                 f"same visit axis in every panel; {ynote}",
                 fontsize=12, fontweight="bold")
    handles = [Line2D([], [], color=col, ls="--", marker="o", ms=6, lw=1.6,
                      markerfacecolor="white", markeredgecolor=col, label="PRE (before exercise)"),
               Line2D([], [], color=col, ls="-", marker="s", ms=6, lw=1.6,
                      markerfacecolor=col, markeredgecolor=col, label="POST (after exercise)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, 0.005))
    fig.tight_layout(rect=(0, 0.05, 1, 1))

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        tag = "rel" if rel else "abs"
        p = os.path.join(save_dir, f"band_{band}_{tag}_all.png")
        fig.savefig(p, dpi=160, facecolor="white")
        print("saved", p)
    return fig


for band in ["delta", "theta", "alpha", "beta"]:
    plot_band_grid(trends, band, value="rel_power_pct", save_dir=FIG_DIR)
    plot_band_grid(trends, band, value="abs_power_uv2", save_dir=FIG_DIR)
