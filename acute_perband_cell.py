# ---- Acute effect: one PNG per band, all participants -----------------------
# Needs `acute` and `FIG_DIR`, which you already have. No CIs, no new data.
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np, os

BANDS = ["delta", "theta", "alpha", "beta"]
BCOL = {"delta": "#378ADD", "theta": "#1D9E75", "alpha": "#D85A30", "beta": "#7F77DD"}
BAND_HZ = {"delta": "0.5-4 Hz", "theta": "4-8 Hz", "alpha": "8-13 Hz", "beta": "13-30 Hz"}
MUTED = "#8a8f98"


def plot_acute_band(acute, band, save_dir=None, max_ticks=8, clip=None):
    """One band, every participant, one PNG. Same x-axis and y-axis on every panel."""
    pids = sorted(acute.participant.unique())
    if not pids:
        print("nothing to plot"); return None
    col = BCOL[band]
    ncol = min(3, len(pids)); nrow = int(np.ceil(len(pids) / ncol))
    vis = sorted(int(v) for v in acute.visit.dropna().unique())
    step = max(1, int(np.ceil(len(vis) / max_ticks)))
    ticks = vis[::step]

    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.0 * nrow),
                             squeeze=False, sharex=True, sharey=True)
    for ax, pid in zip(axes.ravel(), pids):
        d = acute[(acute.participant == pid) & (acute.band == band)].sort_values("visit")
        ax.axhline(0, color=MUTED, lw=1, ls=":")
        if len(d):
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
                 f"each dot is one visit; line is a least-squares guide.  "
                 f"above 0 = power rose after exercise, below 0 = it fell",
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


for band in BANDS:
    plot_acute_band(acute, band, save_dir=FIG_DIR)
