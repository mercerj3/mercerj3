# ---- Acute exercise effect: ALL participants and ALL bands in ONE png -------
# Replaces the per-participant loop in Step 7. Paste-and-run.
#   rows    = band (delta / theta / alpha / beta)
#   columns = participant
#   each dot = one visit;  y = log2(POST / PRE)
#   above 0 = power rose after exercise, below 0 = it fell
# Every panel shares the same x-axis, and panels in the same ROW share a y-axis
# so participants are directly comparable within a band.
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np, os

BANDS = ["delta", "theta", "alpha", "beta"]
BCOL = {"delta": "#378ADD", "theta": "#1D9E75", "alpha": "#D85A30", "beta": "#7F77DD"}
BAND_HZ = {"delta": "0.5-4 Hz", "theta": "4-8 Hz", "alpha": "8-13 Hz", "beta": "13-30 Hz"}
MUTED = "#8a8f98"


def plot_acute_grid(acute, save_dir=None, max_ticks=6, trend_line=True, clip=None):
    """All participants x all bands on one figure.

    clip: optional y-limit, e.g. clip=3 shows only -3..+3 (a 8x change either
          way) so one extreme visit cannot flatten every other panel.
    """
    pids = sorted(acute.participant.unique())
    if not pids:
        print("nothing to plot"); return None
    n = len(pids)

    # one visit range + whole-number ticks for every panel
    vis = sorted(int(v) for v in acute.visit.dropna().unique())
    step = max(1, int(np.ceil(len(vis) / max_ticks)))
    ticks = vis[::step]

    fig, axes = plt.subplots(len(BANDS), n, figsize=(2.6 * n + 1.0, 2.5 * len(BANDS) + 1.0),
                             squeeze=False, sharex=True, sharey="row")
    for r, band in enumerate(BANDS):
        col = BCOL[band]
        for c, pid in enumerate(pids):
            ax = axes[r][c]
            d = acute[(acute.participant == pid) & (acute.band == band)].sort_values("visit")
            ax.axhline(0, color=MUTED, lw=1, ls=":")
            if len(d):
                ax.plot(d.visit, d.log2_ratio, "o", ms=5, color=col,
                        markeredgecolor="white", markeredgewidth=1.0)
                if trend_line and len(d) >= 2 and d.visit.nunique() > 1:
                    x = d.visit.astype(float).to_numpy(); y = d.log2_ratio.to_numpy()
                    m, b = np.polyfit(x, y, 1)
                    xs = np.linspace(x.min(), x.max(), 20)
                    ax.plot(xs, m * xs + b, "-", lw=1.8, color=col, alpha=0.55)
            else:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center",
                        va="center", fontsize=7, color=MUTED)
            ax.set_xticks(ticks)
            ax.set_xlim(min(vis) - 0.7, max(vis) + 0.7)
            if clip:
                ax.set_ylim(-clip, clip)
            ax.tick_params(labelsize=7.5, length=0)
            ax.grid(alpha=.4)
            if r == 0:
                ax.set_title(f"P{pid}", fontsize=10, fontweight="bold")
            if c == 0:
                ax.set_ylabel(f"{band}\n{BAND_HZ[band]}\nlog2(POST/PRE)",
                              fontsize=8.5, color=col, fontweight="bold")
            if r == len(BANDS) - 1:
                ax.set_xlabel("visit", fontsize=8.5)

    fig.suptitle("Acute exercise effect across visits - frontal No-Go, all participants\n"
                 "each dot is one visit; line is a least-squares guide.  "
                 "above 0 = power rose after exercise, below 0 = it fell",
                 fontsize=12, fontweight="bold")
    fig.legend(handles=[Line2D([], [], ls=":", color=MUTED, label="no change (0)")],
               loc="lower center", ncol=1, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, 0.004))
    fig.tight_layout(rect=(0, 0.035, 1, 1))

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        p = os.path.join(save_dir, "acute_all_participants.png")
        fig.savefig(p, dpi=160, facecolor="white")
        print("saved", p)
    return fig


plot_acute_grid(acute, save_dir=FIG_DIR)
