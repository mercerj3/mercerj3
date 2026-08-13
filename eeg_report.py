"""EEG Go/No-Go report generator (self-contained, EEG only).

Point it at a folder of raw recordings (any layout - it searches recursively for files
named  PIN <id> V<k> <PRE|POST> ... Go No-Go ... .csv / .csv.zip) and it computes
everything from raw. Three report types:

    import eeg_report
    eeg_report.generate_reports()                 # asks: daily / weekly / complete / all
    eeg_report.generate_reports(kind="daily",  participant=9, visit=13)
    eeg_report.generate_reports(kind="daily",  participant="all", visit="all")  # every participant + visit
    eeg_report.generate_reports(kind="weekly", participant=9, week=2)  # week 1-4: Control/Heat1-3
    eeg_report.generate_reports(kind="weekly", participant="all", week="all")
    eeg_report.generate_reports(kind="complete")                    # group (all participants), asks scope
    eeg_report.generate_reports(kind="complete", participant="all") # group, all participants
    eeg_report.generate_reports(kind="complete", participant=1)     # individual complete report
    eeg_report.generate_reports(kind="all")       # daily + weekly + complete for everything

Outputs land in eeg_reports_output/reports/ with clickable download links in a notebook.
"""
import os, re, glob, zipfile, warnings, sys, subprocess, pickle, hashlib, textwrap


def _ensure_packages():
    for pip_name, module in [("numpy", "numpy"), ("pandas", "pandas"), ("scipy", "scipy"), ("matplotlib", "matplotlib")]:
        try:
            __import__(module)
        except ImportError:
            try:
                print(f"Installing '{pip_name}' (one-time)...")
                subprocess.check_call([sys.executable, "-m", "pip", "install", "--quiet", pip_name])
            except Exception as e:
                print(f"Could not auto-install '{pip_name}': {e}")


_ensure_packages()
import numpy as np, pandas as pd  # noqa: E402
from scipy import signal, stats  # noqa: E402
from scipy.signal.windows import dpss  # noqa: E402
import matplotlib; matplotlib.use("Agg")  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
warnings.filterwarnings("ignore")

# ---- plain-language text shown in the reports ----
ERP_TEXT = ("What is an ERP?  An ERP (event-related potential) is the brain's average electrical "
            "response in the roughly 1 second right after a stimulus appears. Because a single trial is buried "
            "in noise, we line up many trials to the moment the stimulus appears and average them; the "
            "shared, time-locked response survives and the random noise cancels out. The wiggles (called "
            "N2, P3, etc.) reflect how the brain detected and acted on the stimulus.")
NOGO_TEXT = ("Why No-Go and not Go?  In this task, 'Go' trials just ask for a quick button press, while "
             "'No-Go' trials ask the person to STOP themselves from pressing - which requires active "
             "inhibitory/attentional control. That effortful control is the cognitive process this study "
             "cares about, so the main analyses use the No-Go response (Go is shown alongside for comparison).")
DESC = {
    "signal": "Raw EEG straight off the sensor (top), the full cleaned/filtered recording (middle; still contains movement bursts), and a 10-second zoom on a clean stretch (bottom) - a quality check.",
    "spectrogram": "How EEG power at each frequency (vertical) changes across the whole recording (horizontal); warmer = more power.",
    "psd": "Average signal strength at each frequency for this recording.",
    "imu": "How much the head moved during each visit (higher = more movement) - a data-quality context measure.",
    "qc": "How much data survived cleaning: share of 1-second snippets removed as noise, and how many clean snippets remained.",
    "erp": ("The brain's average response in roughly the 1 second after each stimulus, averaged over many trials so "
            "random noise cancels out. No-Go trials (holding back a button press) engage inhibitory control - the process "
            "of interest here; Go trials are shown alongside for comparison."),
    "tfr": "How brain-wave power rises/falls moment-to-moment after the stimulus, by frequency band.",
    "bandpower_session": "How each brain-wave band's power changes across the recording.",
    "bandbars": "Absolute power in each brain-wave band before vs after exercise (error bars = 95% confidence interval across clean epochs).",
    "bandbars_rel": "Each band's SHARE of total power (%) before vs after exercise - shows how the balance between bands shifts, independent of overall signal size (error bars = 95% CI across clean epochs).",
    "group_bandabs": "Absolute power in each band, averaged across participants, before vs after exercise (error bars = 95% CI across participants).",
    "group_bandrel": "Each band's SHARE of total power (%), averaged across participants, before vs after exercise (error bars = 95% CI across participants).",
    "eeg_levels": "Overall (all bands combined) and each band on its own graph - EEG power before (dashed) and after (solid) exercise at each visit.",
    "eeg_levels_rel": "Each band's SHARE of total power (%) at each visit, before (dashed) vs after (solid) exercise - shows how the balance between bands shifts across the week, independent of overall signal size.",
    "acute": "How the before-vs-after-exercise change in each band shifts across visits (the 'adaptiveness' trend). No-Go is the response of interest; Go is shown alongside for comparison.",
    "rejection_time": "Share of EEG removed as noise at each visit - a data-quality check across the week.",
    "adaptiveness_bar": "Each participant's overall adaptiveness: did their acute response trend up or down over the study.",
    "group_prepost": "Average before-vs-after-exercise response across all participants, with confidence intervals.",
    "classify": ("Each participant's acute effect (POST minus PRE) at every visit. Responder type is from a Hotelling T-squared "
                 "test across the four bands (frontal No-Go): rapid responder = power rose significantly (p<.05, net positive); "
                 "decreaser = fell significantly; nonresponder = no significant joint change. Faded lines are individuals; thick "
                 "line is the group mean with 95% CI."),
    "erp_byclass": "The average PRE and POST brain response (frontal No-Go), shown separately for each responder type, with 95% CI across participants.",
    "over_visits": "Average frontal EEG power before (PRE) and after (POST) exercise at each visit - overall and split by responder type.",
}
WEEKS = [("Control", 4, 8), ("Heat week 1", 9, 13), ("Heat week 2", 14, 18), ("Heat week 3", 19, 23)]


def generate_reports(kind=None, participant=None, visit=None, week=None, base_dir=None, out_dir="eeg_reports_output",
                     show_links=True, zip_name="eeg_reports.zip", _max_files=None, ask=True):
    FS, THR = 250, 20
    NPRE, NPOST = 50, 200; EP = NPRE + NPOST; ART = 100.0; NFFT = 512; NW, KT = 2.0, 3
    CH = ["F3", "F4", "C3", "Cz", "C4", "P3", "P4"]
    RIDX = {"frontal": [0, 1], "central": [2, 3, 4], "posterior": [5, 6]}
    ACC = ["AccelX", "AccelY", "AccelZ"]; CONDS = {"Go": 5, "NoGo": 6}
    BANDS = ["delta", "theta", "alpha", "beta"]; BAND_HZ = {"delta": (0.5, 4), "theta": (4, 8), "alpha": (8, 13), "beta": (13, 30)}
    BCOL = {"delta": "#378ADD", "theta": "#1D9E75", "alpha": "#D85A30", "beta": "#7F77DD"}
    KORD = ["rapid responder", "nonresponder", "decreaser"]
    KCOL = {"rapid responder": "#1a7f37", "nonresponder": "#8a8f98", "decreaser": "#c0392b"}
    NAVY, GREY, RULE, TOP = "#22406b", "#8a8f98", "#cfd3d8", 0.88
    b_bp, a_bp = signal.butter(4, [0.1 / (FS / 2), 30 / (FS / 2)], btype="band")
    b_n, a_n = signal.iirnotch(60.0, 30.0, FS)
    TAP = dpss(EP, NW, KT, sym=False); TAP = TAP / np.sqrt((TAP ** 2).sum(1, keepdims=True))
    FREQS = np.fft.rfftfreq(NFFT, 1 / FS); DF = FS / NFFT
    BM = {b: ((FREQS >= lo) & (FREQS <= hi)) for b, (lo, hi) in BAND_HZ.items()}
    tms = (np.arange(EP) - NPRE) / FS * 1000.0
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.titlesize": 10, "axes.titleweight": "bold",
                         "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8})

    # ---------- locate recordings ----------
    def _scan(folder):
        fs = []
        for pat in ("**/*Go No-Go*.csv", "**/*Go No-Go*.csv.zip", "*Go No-Go*.csv", "*Go No-Go*.csv.zip"):
            fs += glob.glob(os.path.join(folder, pat), recursive=True)
        return sorted(set(f for f in fs if re.search(r"PIN\s*\d+", os.path.basename(f), re.I)))

    if ask and base_dir is None:
        for _ in range(3):
            try:
                resp = input("Folder where you downloaded the recordings [Enter = current folder]: ").strip().strip('"').strip("'")
            except Exception:
                resp = ""
            cand = os.path.abspath(os.path.expanduser(resp)) if resp else os.getcwd()
            if _scan(cand) or not resp:
                base_dir = cand if resp else None; break
            print(f"  No 'Go No-Go' files under: {cand} - try again.")
    base = os.path.abspath(base_dir) if base_dir else os.getcwd()
    print(f"Scanning for recordings under: {base}", flush=True)
    ALL = _scan(base)
    if not ALL:
        print(f"No Go/No-Go recordings found under {base!r}."); return []
    REP = os.path.join(os.path.abspath(out_dir), "reports"); os.makedirs(REP, exist_ok=True)
    print(f"Found {len(ALL)} recording file(s). Reports will be saved in: {REP}", flush=True)

    def meta(fn):
        up = fn.upper(); p = re.search(r"PIN\s*(\d+)", up); v = re.search(r"\bV(\d+)", up)
        ph = "post" if "POST" in up else ("pre" if "PRE" in up else "na")
        return (int(p.group(1)) if p else None, int(v.group(1)) if v else None, ph)

    def eeg_files(pid):
        fs = [f for f in ALL if meta(os.path.basename(f))[0] == pid]
        seen = {}
        for f in sorted(fs, key=lambda x: x.lower().endswith(".zip")):
            key = meta(os.path.basename(f))[1:]
            if key not in seen: seen[key] = f
        out = sorted(seen.values())
        return out[:_max_files] if _max_files else out

    all_pids = sorted({meta(os.path.basename(f))[0] for f in ALL if meta(os.path.basename(f))[0] is not None})

    def read_cols(path, want):
        try:
            if path.endswith(".zip"):
                with zipfile.ZipFile(path) as z:
                    df = pd.read_csv(z.open(z.namelist()[0]), usecols=lambda c: c.strip().lower().startswith("time") or c in want)
            else:
                df = pd.read_csv(path, usecols=lambda c: c.strip().lower().startswith("time") or c in want)
        except Exception:
            return None
        for c in df.columns: df[c] = pd.to_numeric(df[c], errors="coerce")
        return df

    def mt(ep):
        acc = None
        for k in range(KT):
            Xk = np.fft.rfft(TAP[k][:, None] * ep, n=NFFT, axis=0); acc = np.abs(Xk) ** 2 if acc is None else acc + np.abs(Xk) ** 2
        return acc / KT

    def process(path, want_cont=False):
        # SAME methodology as the graph pipeline. For EACH region independently (a missing/NaN channel in
        # one region -- e.g. parietal -- never drops another region):
        #   SPECTRAL band power: per-channel multitaper power, then AVERAGE the region's channels (pool
        #     AFTER the FFT); reject an epoch if any region channel ptp>100 uV; MEAN across clean epochs.
        #     (Pooling channels in the time domain BEFORE the FFT keeps only phase-locked/evoked power and
        #      is wrong for band power, so channel pooling for power is done after the transform.)
        #   ERP (N200/P300): pool the region's channels in the TIME domain FIRST, pre-stim baseline, reject
        #     on the pooled waveform ptp>100 uV, MEAN across clean epochs.
        df = read_cols(path, ["Marker"] + CH + (ACC if want_cont else []))
        if df is None or "Marker" not in df or not any(c in df for c in CH): return None
        tcol = [c for c in df.columns if c.lower().startswith("time")][0]
        pin, visit, phase = meta(os.path.basename(path))
        mk_all = df["Marker"].fillna(0).to_numpy(); t_all = df[tcol].to_numpy(float)
        rows = []; erp_reg = {r: {c: [] for c in CONDS} for r in RIDX}; band_epochs = {}
        cont = {"raw_front": None, "filt_front": None, "acc": {c: None for c in ACC}, "dur_min": np.nan}
        for region, idx in RIDX.items():
            chans = [CH[i] for i in idx if CH[i] in df]
            if not chans: continue
            good = df[chans].notna().all(axis=1).to_numpy() & np.isfinite(t_all)
            if int(good.sum()) < EP * 2: continue
            ts = t_all[good]; X = df[chans].to_numpy(float)[good]; mk = mk_all[good]
            o = np.argsort(ts, kind="stable"); ts, X, mk = ts[o], X[o], mk[o]
            filt = np.empty_like(X)
            for k in range(X.shape[1]): filt[:, k] = signal.filtfilt(b_n, a_n, signal.filtfilt(b_bp, a_bp, X[:, k]))
            pooled = filt.mean(1)                                       # pool region channels in TIME domain first
            if want_cont and region == "frontal":
                cont["raw_front"] = X.mean(1); cont["filt_front"] = pooled; cont["dur_min"] = len(pooled) / FS / 60.0
                for c in ACC:
                    if c in df: cont["acc"][c] = df[c].to_numpy(float)[good][o]
            erp_by = {c: [] for c in CONDS}; spec_by = {c: [] for c in CONDS}; n_total = {c: 0 for c in CONDS}
            for code, t in zip(mk, ts):
                if code not in CONDS.values(): continue
                cond = "Go" if code == 5 else "NoGo"; s = int(np.searchsorted(ts, t)); a, bb = s - NPRE, s + NPOST
                if a < 0 or bb > len(filt): continue
                n_total[cond] += 1
                segc = filt[a:bb, :] - filt[a:bb, :][:NPRE].mean(0, keepdims=True)   # per-channel, pre-stim baseline
                if np.all(np.ptp(segc, 0) <= ART): spec_by[cond].append(mt(segc))    # PER-CHANNEL power (pool AFTER)
                segp = pooled[a:bb] - pooled[a:bb][:NPRE].mean()                      # pooled region waveform for ERP
                if np.ptp(segp) <= ART: erp_by[cond].append(segp)
            for cond in CONDS:
                erp_reg[region][cond] = erp_by[cond]                                 # ERP (pool FIRST), kept per region
                ne = len(spec_by[cond]); nt = n_total[cond]
                if ne == 0: continue
                stack = np.stack(spec_by[cond])                                       # ne x freq x nch
                Sf = stack.mean(0).mean(1); psd = Sf / FS; psd = psd.copy(); psd[1:-1] *= 2   # MEAN epochs, AVG channels
                for b, mask in BM.items():
                    rows.append(dict(visit=visit, phase=phase, condition=cond, region=region, band=b,
                                     abs_power_uv2=float(np.sum(psd[mask]) * DF), n_epochs=ne,
                                     pct_rejected=100 * (nt - ne) / nt if nt else np.nan))
                if region == "frontal" and cond == "NoGo":                            # per-epoch band power -> daily error bars
                    pe = stack.mean(2) / FS; pe = pe.copy(); pe[:, 1:-1] *= 2
                    band_epochs = {b: (pe[:, mask].sum(1) * DF) for b, mask in BM.items()}
        if not rows: return None
        out = dict(visit=visit, phase=phase, rows=rows, band_epochs=band_epochs,
                   nogo=erp_reg["frontal"]["NoGo"], go=erp_reg["frontal"]["Go"],
                   nogo_par=erp_reg["posterior"]["NoGo"], go_par=erp_reg["posterior"]["Go"])
        if want_cont:
            out.update(raw_front=cont["raw_front"], filt_front=cont["filt_front"], acc=cont["acc"], dur_min=cont["dur_min"])
        return out

    # ---------- shared plotting ----------
    def styled(section, sup=None):
        fig = plt.figure(figsize=(8.5, 11))
        if sup: fig.text(0.07, 0.955, sup, fontsize=10, color=NAVY, fontweight="bold")
        fig.text(0.07, 0.928, section, fontsize=14, color=NAVY, fontweight="bold")
        fig.add_artist(Line2D([0.07, 0.93], [0.915, 0.915], color=RULE, lw=1.1, transform=fig.transFigure)); return fig

    def caption(fig, key):
        txt = "\n".join(textwrap.wrap(DESC.get(key, ""), 118))  # hard-wrap so it never runs off the page
        fig.text(0.07, 0.052, txt, fontsize=8, color=GREY, va="top")

    def gsp(fig, nr, nc, **kw):
        dd = dict(top=TOP, bottom=0.11, left=0.11, right=0.93); dd.update(kw); return fig.add_gridspec(nr, nc, **dd)

    def dec(y, target=3000): return max(1, len(np.asarray(y)) // target)

    def fitln(ax, x, y, c):
        x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(x) & np.isfinite(y)
        if m.sum() >= 2 and np.ptp(x[m]) > 0:
            b, a = np.polyfit(x[m], y[m], 1); xs = np.linspace(x[m].min(), x[m].max(), 30); ax.plot(xs, a + b * xs, "--", color=c, lw=1.4, alpha=.85)

    def erp_mean(eps):
        return np.mean(eps, 0) if eps else np.full(EP, np.nan)

    def tfr_map(eps):
        if not eps: return None, None, None
        sp = [signal.spectrogram(e, fs=FS, nperseg=48, noverlap=44, nfft=128, mode="psd")[2] for e in eps]
        fq, t, _ = signal.spectrogram(eps[0], fs=FS, nperseg=48, noverlap=44, nfft=128, mode="psd")
        db = 10 * np.log10(np.mean(sp, 0) + 1e-12); base = db[:, t < NPRE / FS].mean(1, keepdims=True)
        return fq, (t - NPRE / FS) * 1000, db - base

    def page_intro(title_line, sub):
        fig = plt.figure(figsize=(8.5, 11))
        fig.add_artist(plt.Circle((0.085, 0.80), 0.015, color=NAVY, transform=fig.transFigure))
        fig.text(0.115, 0.79, title_line, fontsize=21, fontweight="bold", color="#2b2b2b")
        fig.add_artist(Line2D([0.07, 0.93], [0.758, 0.758], color=RULE, lw=1.2, transform=fig.transFigure))
        fig.text(0.09, 0.72, sub, fontsize=11, color=GREY, wrap=True)
        return fig

    def draw_signal(fig, rawf, filtf, ses):  # raw (entire) + cleaned (entire) + 10-second clean zoom
        tmin = np.arange(len(filtf)) / FS / 60.0; s = dec(filtf)
        g = fig.add_gridspec(3, 1, top=0.86, bottom=0.10, left=0.12, right=0.92, hspace=.55)
        a0 = fig.add_subplot(g[0]); a0.plot(tmin[::s], rawf[::s], color="k", lw=.4, rasterized=True)
        a0.set_ylabel("raw [uV]"); a0.set_title(f"Raw EEG - entire recording ({ses})"); a0.grid(alpha=.25)
        a1 = fig.add_subplot(g[1], sharex=a0); a1.plot(tmin[::s], filtf[::s], color="#185FA5", lw=.4, rasterized=True)
        a1.set_ylim(-120, 120); a1.set_ylabel("clean [uV]"); a1.set_xlabel("time (min)")
        a1.set_title("Cleaned EEG - entire recording (continuous; still contains movement bursts)"); a1.grid(alpha=.25)
        w = int(10 * FS); best, bp = None, None
        for s0 in range(0, max(1, len(filtf) - w), 5 * FS):
            seg = filtf[s0:s0 + w]; sd = seg.std(); p2p = np.ptp(seg)
            if sd < 3 or sd > 60 or p2p > 300: continue
            if bp is None or p2p < bp: best, bp = s0, p2p
        if best is None: best = max(0, len(filtf) // 2 - w // 2)
        seg = filtf[best:best + w]; tsec = np.arange(len(seg)) / FS; m = max(20, np.percentile(np.abs(seg), 99) * 1.4)
        a2 = fig.add_subplot(g[2]); a2.plot(tsec, seg, color="#1a7f37", lw=.8); a2.set_ylim(-m, m)
        a2.set_ylabel("clean [uV]"); a2.set_xlabel("time (s)"); a2.set_title(f"10-second clean zoom (@ {best/FS/60:.1f} min) - what good EEG looks like"); a2.grid(alpha=.25)

    def spectrogram_ax(ax, filtf):
        f, t, Sxx = signal.spectrogram(filtf, fs=FS, nperseg=512, noverlap=384); db = 10 * np.log10(Sxx + 1e-12); sel = f <= 35
        im = ax.pcolormesh(t / 60.0, f[sel], db[sel], cmap="turbo", shading="auto", rasterized=True,
                           vmin=np.percentile(db[sel], 5), vmax=np.percentile(db[sel], 99))
        ax.set_ylabel("frequency [Hz]"); return im

    def erp_ax(ax, nogo, go, title):
        ax.axhline(0, color="grey", lw=.6); ax.axvline(0, color="grey", ls=":", lw=.8)
        for eps, col, lab in [(go, "#888", "Go"), (nogo, "#c0392b", "No-Go")]:
            y = erp_mean(eps)
            if np.isfinite(y).any(): ax.plot(tms, y, color=col, lw=1.6, label=f"{lab} (n={len(eps)})")
        ax.set_title(title); ax.set_ylabel("frontal amp (uV)"); ax.legend(fontsize=8); ax.grid(alpha=.3)

    def tfr_ax(ax, nogo):
        f, t, M = tfr_map(nogo)
        if M is None or not M.size: ax.text(.5, .5, "no data", ha="center", transform=ax.transAxes); return None
        sel = f <= 30; vmax = np.nanpercentile(np.abs(M[sel]), 98) or 3
        im = ax.pcolormesh(t, f[sel], M[sel], cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto", rasterized=True)
        ax.axvline(0, color="k", ls=":", lw=.8); ax.set_ylabel("freq (Hz)"); return im

    def bandpower_over(ax, filtf, color_by_band=True):
        win, step = int(2 * FS), int(FS / 2); centers = []; bp = {b: [] for b in BANDS}
        for st in range(0, max(1, len(filtf) - win), step):
            seg = filtf[st:st + win] * np.hanning(win); fr = np.fft.rfftfreq(win, 1 / FS); P = np.abs(np.fft.rfft(seg)) ** 2
            centers.append((st + win / 2) / FS / 60.0)
            for b, (lo, hi) in BAND_HZ.items(): bp[b].append(P[(fr >= lo) & (fr <= hi)].sum())
        for b in BANDS: ax.plot(centers, 10 * np.log10(np.array(bp[b]) + 1e-12), color=BCOL[b], lw=1.1, label=b)
        ax.set_ylabel("band power (dB)"); ax.legend(fontsize=7, ncol=4); ax.grid(alpha=.3)

    def bandpower_frontal_nogo(rows, phase):
        out = {}
        for r in rows:
            if r["phase"] == phase and r["condition"] == "NoGo" and r["region"] == "frontal":
                out[r["band"]] = r["abs_power_uv2"]
        return out

    def qc_frontal_nogo(rows, phase):
        for r in rows:
            if r["phase"] == phase and r["condition"] == "NoGo" and r["region"] == "frontal":
                return r["n_epochs"], r["pct_rejected"]
        return 0, 100.0

    # ==================================================================== DAILY
    def _nodata(ax, msg):
        ax.text(0.5, 0.5, msg, ha="center", va="center", color=GREY, fontsize=9, transform=ax.transAxes)

    def build_daily(pid, vis):
        print(f"Building daily report for P{pid} V{vis}...", flush=True)
        pre = next((f for f in eeg_files(pid) if meta(os.path.basename(f))[1] == vis and meta(os.path.basename(f))[2] == "pre"), None)
        post = next((f for f in eeg_files(pid) if meta(os.path.basename(f))[1] == vis and meta(os.path.basename(f))[2] == "post"), None)
        P = process(pre, want_cont=True) if pre else None
        Q = process(post, want_cont=True) if post else None
        if P is None and Q is None: print(f"P{pid} V{vis}: no usable recording"); return None

        def status(S, has_file):  # present / unclean / missing
            if not has_file: return "missing"
            if S is None or len(S["nogo"]) == 0: return "unclean"
            return "present"
        st = {"PRE": status(P, pre is not None), "POST": status(Q, post is not None)}
        rows = (P["rows"] if P else []) + (Q["rows"] if Q else [])
        path = os.path.join(REP, f"P{pid}_V{vis}_daily.pdf")
        with PdfPages(path) as pdf:
            wk = next((n for n, lo, hi in WEEKS if lo <= vis <= hi), "-")
            note = "  |  ".join(f"{k}: {st[k]}" for k in ("PRE", "POST"))
            pdf.savefig(page_intro(f"Participant {pid} - Visit {vis} (daily)", f"Study week: {wk}     ({note})")); plt.close()
            # signal (uses whichever phase has the most clean epochs, else whichever recorded)
            best = max([S for S in (P, Q) if S], key=lambda S: len(S["nogo"]), default=None)
            fig = styled("EEG signal", "Signal Overview")
            if best is not None: draw_signal(fig, best["raw_front"], best["filt_front"], f"V{vis} {best['phase']}")
            else: _nodata(fig.add_axes([0.1, 0.4, 0.8, 0.3]), "No continuous EEG for this visit")
            caption(fig, "signal"); pdf.savefig(fig); plt.close()
            # spectrogram PRE & POST
            fig = styled("EEG spectrogram - PRE vs POST", "Signal Overview"); g = gsp(fig, 2, 1, hspace=.3, right=0.88); im = None
            for i, (lab, S) in enumerate([("PRE", P), ("POST", Q)]):
                ax = fig.add_subplot(g[i]); ax.set_title(lab)
                if S: im = spectrogram_ax(ax, S["filt_front"])
                else: _nodata(ax, f"{lab}: {st[lab]}")
            fig.axes[-1].set_xlabel("time (min)")
            if im is not None: fig.colorbar(im, cax=fig.add_axes([0.90, 0.15, 0.02, 0.7]), label="dB")
            caption(fig, "spectrogram"); pdf.savefig(fig); plt.close()
            # QC bars PRE & POST
            fig = styled("Data quality (QC) - PRE vs POST", "Electrode Performance"); g = gsp(fig, 1, 2, wspace=.35, bottom=0.16)
            ax0 = fig.add_subplot(g[0]); ax1 = fig.add_subplot(g[1])
            for i, (lab, S) in enumerate([("PRE", P), ("POST", Q)]):
                ne, pct = qc_frontal_nogo(rows, lab.lower())
                if S: ax0.bar(i, pct, color="#c0392b"); ax1.bar(i, ne, color=("#1a7f37" if ne >= THR else "#c0392b")); ax1.text(i, ne, str(ne), ha="center", va="bottom", fontsize=8)
                else: ax1.text(i, 1, st[lab], ha="center", va="bottom", color=GREY, fontsize=8)
            ax0.set_xticks([0, 1]); ax0.set_xticklabels(["PRE", "POST"]); ax0.set_ylabel("% epochs rejected"); ax0.set_title("Artifact rejection"); ax0.grid(alpha=.3, axis="y")
            ax1.set_xticks([0, 1]); ax1.set_xticklabels(["PRE", "POST"]); ax1.set_ylabel("clean epochs"); ax1.axhline(THR, color="k", ls="--", lw=1); ax1.set_title("Epoch count (No-Go, frontal)"); ax1.grid(alpha=.3, axis="y")
            caption(fig, "qc"); pdf.savefig(fig); plt.close()
            # ERP PRE & POST
            fig = styled("ERP (stimulus-locked) - PRE vs POST", "Task Analysis"); g = gsp(fig, 2, 1, hspace=.3); a0 = None
            for i, (lab, S) in enumerate([("PRE", P), ("POST", Q)]):
                ax = fig.add_subplot(g[i], sharex=a0, sharey=a0) if a0 else fig.add_subplot(g[i]); a0 = a0 or ax
                if S and len(S["nogo"]): erp_ax(ax, S["nogo"], S["go"], lab)
                else: ax.set_title(lab); _nodata(ax, f"{lab}: {st[lab]}")
            a0.set_xlabel("time from stimulus (ms)"); caption(fig, "erp"); pdf.savefig(fig); plt.close()
            # TFR PRE & POST
            fig = styled("Time-frequency power (No-Go) - PRE vs POST", "Task Analysis"); g = gsp(fig, 2, 1, hspace=.3, right=0.88); im = None
            for i, (lab, S) in enumerate([("PRE", P), ("POST", Q)]):
                ax = fig.add_subplot(g[i]); ax.set_title(lab)
                if S and len(S["nogo"]): im = tfr_ax(ax, S["nogo"]) or im
                else: _nodata(ax, f"{lab}: {st[lab]}")
            fig.axes[-1].set_xlabel("time (ms)")
            if im is not None: fig.colorbar(im, cax=fig.add_axes([0.90, 0.15, 0.02, 0.7]), label="power vs baseline (dB)")
            caption(fig, "tfr"); pdf.savefig(fig); plt.close()
            # band power over session PRE & POST
            fig = styled("Band power over the session - PRE vs POST", "Task Analysis"); g = gsp(fig, 2, 1, hspace=.3); a0 = None
            for i, (lab, S) in enumerate([("PRE", P), ("POST", Q)]):
                ax = fig.add_subplot(g[i], sharex=a0) if a0 else fig.add_subplot(g[i]); a0 = a0 or ax; ax.set_title(lab)
                if S: bandpower_over(ax, S["filt_front"])
                else: _nodata(ax, f"{lab}: {st[lab]}")
            a0.set_xlabel("time (min)"); caption(fig, "bandpower_session"); pdf.savefig(fig); plt.close()
            # 4 band bars PRE vs POST
            band_bars_pages(pdf, P["band_epochs"] if P else {}, Q["band_epochs"] if Q else {}, st)
        print(f"wrote {os.path.basename(path)}  (PRE {st['PRE']} | POST {st['POST']})")
        return path

    # ==================================================================== WEEKLY (one week)
    def build_weekly(pid, widx=0):
        wname, lo, hi = WEEKS[widx]
        print(f"Building weekly report for P{pid} - {wname}...", flush=True)
        sess = {}
        for f in eeg_files(pid):
            _, v, ph = meta(os.path.basename(f))
            if v is None or not (lo <= v <= hi): continue
            r = process(f)
            if r: sess[(r["visit"], r["phase"])] = r
        if not sess: print(f"P{pid}: no usable EEG in {wname} (V{lo}-V{hi})"); return None
        rows = [x for r in sess.values() for x in r["rows"]]
        d = pd.DataFrame(rows); vs = sorted({v for (v, ph) in sess})
        path = os.path.join(REP, f"P{pid}_{wname.replace(' ', '')}_weekly.pdf")
        with PdfPages(path) as pdf:
            pdf.savefig(page_intro(f"Participant {pid} - {wname} (weekly)", f"Visits V{lo}-V{hi}    ({len(vs)} visit(s) recorded)")); plt.close()
            # 1) artifact rejection over the week's visits
            fig = styled(f"Artifact rejection over visits - {wname}", "Data quality"); ax = fig.add_axes([0.12, 0.16, 0.8, 0.68])
            g = d[d.condition == "NoGo"]
            for r, c in [("frontal", "#D85A30"), ("central", "#1D9E75"), ("posterior", "#378ADD")]:
                s = g[g.region == r].groupby("visit").pct_rejected.mean()
                if len(s): ax.plot(s.index, s.values, "-o", color=c, ms=5, label=r)
            ax.set_xlabel("visit"); ax.set_ylabel("% epochs rejected"); ax.set_ylim(-2, 102)
            ax.set_xticks(vs); ax.set_xlim(lo - 0.4, hi + 0.4); ax.legend(fontsize=8, ncol=3); ax.grid(alpha=.3)
            caption(fig, "rejection_time"); pdf.savefig(fig); plt.close()
            # 2) EEG band power over the week - ABSOLUTE page + RELATIVE page
            sub = d[(d.condition == "NoGo") & (d.region == "frontal")]
            pv = sub.groupby(["visit", "phase", "band"]).abs_power_uv2.mean().reset_index()
            tot = pv.groupby(["visit", "phase"]).abs_power_uv2.sum().rename("tot").reset_index()
            pv = pv.merge(tot, on=["visit", "phase"]); pv["rel"] = np.where(pv.tot > 0, pv.abs_power_uv2 / pv.tot * 100, np.nan)
            # 2a) absolute
            fig = styled(f"EEG power levels over {wname} (absolute)", "Study Trajectories"); g2 = gsp(fig, 3, 2, hspace=.6, wspace=.3, bottom=0.09)
            axA = fig.add_subplot(g2[0, :])
            for ph, col, ls in [("pre", "#7f8c8d", "--"), ("post", "#22406b", "-")]:
                s = sub[sub.phase == ph].groupby("visit").abs_power_uv2.sum()
                if len(s): axA.plot(s.index, s.values, ls, marker="o", ms=5, color=col, label=ph.upper())
            axA.set_title("Overall power (all bands combined, frontal No-Go)"); axA.set_ylabel("power (uV^2)"); axA.set_xlabel("visit")
            axA.set_xticks(vs); axA.set_xlim(lo - 0.4, hi + 0.4); axA.legend(fontsize=8); axA.grid(alpha=.3)
            for i, b in enumerate(BANDS):
                ax = fig.add_subplot(g2[1 + i // 2, i % 2])
                for ph, ls, a in [("post", "-", 1.0), ("pre", "--", 0.6)]:
                    s = pv[(pv.phase == ph) & (pv.band == b)].sort_values("visit")
                    if len(s): ax.plot(s.visit, s.abs_power_uv2, ls, marker="o", ms=3, lw=1.3, color=BCOL[b], alpha=a, label=ph.upper())
                ax.set_title(f"{b}  (solid = POST, dashed = PRE)", fontsize=9); ax.set_xlabel("visit"); ax.set_ylabel("power (uV^2)"); ax.set_xticks(vs); ax.grid(alpha=.3)
                if i == 0: ax.legend(fontsize=7)
            caption(fig, "eeg_levels"); pdf.savefig(fig); plt.close()
            # 2b) relative (% of total)
            fig = styled(f"EEG band share over {wname} (relative, % of total)", "Study Trajectories"); g2 = gsp(fig, 2, 2, hspace=.45, wspace=.3, bottom=0.12)
            for i, b in enumerate(BANDS):
                ax = fig.add_subplot(g2[i // 2, i % 2])
                for ph, ls, a in [("post", "-", 1.0), ("pre", "--", 0.6)]:
                    s = pv[(pv.phase == ph) & (pv.band == b)].sort_values("visit")
                    if len(s): ax.plot(s.visit, s.rel, ls, marker="o", ms=3, lw=1.3, color=BCOL[b], alpha=a, label=ph.upper())
                ax.set_title(f"{b}  (solid = POST, dashed = PRE)", fontsize=9); ax.set_xlabel("visit"); ax.set_ylabel("% of total power"); ax.set_xticks(vs); ax.grid(alpha=.3)
                if i == 0: ax.legend(fontsize=7)
            caption(fig, "eeg_levels_rel"); pdf.savefig(fig); plt.close()
            # 3) acute-effect trajectory - No-Go AND Go
            fig = styled(f"Acute-effect trajectory - {wname}", "Study Trajectories"); g2 = gsp(fig, 2, 3, hspace=.5, wspace=.36, bottom=0.13)
            base = d[(d.region == "frontal") & (d.n_epochs >= THR)]

            def acute(cond):
                s = base[base.condition == cond]
                w = s.pivot_table(index=["visit", "band"], columns="phase", values="abs_power_uv2").reindex(columns=["pre", "post"]).dropna(subset=["pre", "post"]) if len(s) else pd.DataFrame()
                return w.assign(l2=np.log2(w.post) - np.log2(w.pre)).reset_index() if len(w) else pd.DataFrame()
            wN, wG = acute("NoGo"), acute("Go")
            if len(wN) or len(wG):
                for i, b in enumerate(BANDS):
                    ax = fig.add_subplot(g2[i // 3, i % 3]); ax.axhline(0, color="grey", ls=":", lw=.7)
                    for wdat, col, mk, lab in [(wN, BCOL[b], "o", "No-Go"), (wG, "#9aa0a6", "s", "Go")]:
                        sb = wdat[wdat.band == b].sort_values("visit") if len(wdat) else pd.DataFrame()
                        if len(sb): ax.scatter(sb.visit, sb.l2, color=col, s=16, marker=mk, zorder=3, label=lab); fitln(ax, sb.visit.values, sb.l2.values, col)
                    ax.set_title(b); ax.set_xlabel("visit"); ax.set_ylabel("log2(POST/PRE)"); ax.set_xticks(vs); ax.grid(alpha=.3)
                    if i == 0: ax.legend(fontsize=7)
                axm = fig.add_subplot(g2[1, 1]); axm.axhline(0, color="grey", ls=":", lw=.7)
                for wdat, col, mk, lab in [(wN, "#333", "o", "No-Go"), (wG, "#9aa0a6", "s", "Go")]:
                    if len(wdat):
                        gm = wdat.groupby("visit", as_index=False).l2.mean(); axm.scatter(gm.visit, gm.l2, color=col, s=20, marker=mk, zorder=3, label=lab); fitln(axm, gm.visit.values, gm.l2.values, col)
                axm.set_title("4-band mean"); axm.set_xlabel("visit"); axm.set_ylabel("log2(POST/PRE)"); axm.set_xticks(vs); axm.grid(alpha=.3); axm.legend(fontsize=7)
            else:
                fig.text(0.5, 0.5, "Not enough paired PRE+POST sessions (>= 20 clean epochs each) to plot adaptiveness.", ha="center", color=GREY, fontsize=11)
            caption(fig, "acute"); pdf.savefig(fig); plt.close()
            # 4) per-week ERPs: each session PRE + POST, then aggregate
            fig = styled(f"ERPs - {wname} (each session + aggregate)", "Task Analysis")
            nrow = len(vs) + 1; g = gsp(fig, nrow, 2, hspace=.55, wspace=.25, bottom=0.09); a0 = None
            aggP, aggQ = [], []
            for ri, v in enumerate(vs):
                for ci, ph in enumerate(["pre", "post"]):
                    ax = fig.add_subplot(g[ri, ci], sharex=a0, sharey=a0) if a0 else fig.add_subplot(g[ri, ci]); a0 = a0 or ax
                    S = sess.get((v, ph)); ax.axhline(0, color="grey", lw=.5)
                    y = erp_mean(S["nogo"]) if S else np.full(EP, np.nan)
                    if S: (aggP if ph == "pre" else aggQ).extend(S["nogo"])
                    if np.isfinite(y).any(): ax.plot(tms, y, color="#c0392b", lw=1.2)
                    else: ax.text(0.5, 0.5, "no clean\nNo-Go epochs", ha="center", va="center", transform=ax.transAxes, color=GREY, fontsize=6)
                    ax.set_title(f"V{v} {ph.upper()}", fontsize=8); ax.tick_params(labelsize=6)
            for ci, agg in enumerate([aggP, aggQ]):
                ax = fig.add_subplot(g[nrow - 1, ci], sharex=a0, sharey=a0); ax.axhline(0, color="grey", lw=.5)
                y = erp_mean(agg)
                if np.isfinite(y).any(): ax.plot(tms, y, color="#111", lw=1.8)
                ax.set_title(f"AGGREGATE {'PRE' if ci == 0 else 'POST'} (n={len(agg)})", fontsize=8, color=NAVY); ax.tick_params(labelsize=6); ax.set_xlabel("ms", fontsize=7)
            caption(fig, "erp"); pdf.savefig(fig); plt.close()
        print(f"wrote {os.path.basename(path)}")
        return path

    # ==================================================================== classification + shared pages
    def classify(traj):  # 95% CI of the mean acute effect -- kept only for the CI display values
        vals = np.array([v for v in traj.values() if np.isfinite(v)], float)
        if len(vals) < 3: return np.nan, np.nan, np.nan
        mean = vals.mean(); sem = vals.std(ddof=1) / np.sqrt(len(vals)); t = stats.t.ppf(.975, len(vals) - 1)
        return mean, mean - t * sem, mean + t * sem

    def hotelling_klass(w):  # SAME responder test as the graph pipeline: Hotelling T^2 across the 4 bands (frontal No-Go)
        piv = w.pivot_table(index="visit", columns="band", values="l2").reindex(columns=BANDS).dropna()
        X = piv.values
        if X.shape[0] == 0: return "nonresponder"
        n, p = X.shape; mean = X.mean(0); net = float(mean.sum())
        if n <= p: return "nonresponder"                     # too few visits for the joint test
        S = np.cov(X, rowvar=False)
        try: Si = np.linalg.inv(S)
        except np.linalg.LinAlgError: return "nonresponder"
        T2 = n * mean @ Si @ mean; F = (n - p) / (p * (n - 1)) * T2; pv = float(stats.f.sf(F, p, n - p))
        if pv < 0.05: return "rapid responder" if net > 0 else "decreaser"
        return "nonresponder"

    def slope_class(d):  # -> slope, pval(of slope), traj{visit:acute}, klass (Hotelling T^2), (mean,lo,hi of acute)
        gg = d[(d.condition == "NoGo") & (d.region == "frontal") & (d.n_epochs >= THR)]
        w = gg.pivot_table(index=["visit", "band"], columns="phase", values="abs_power_uv2").reindex(columns=["pre", "post"]).dropna(subset=["pre", "post"]) if len(gg) else pd.DataFrame()
        if not len(w): return np.nan, np.nan, {}, "nonresponder", (np.nan, np.nan, np.nan)
        w = w.assign(l2=np.log2(w.post) - np.log2(w.pre)).reset_index()
        m = w.groupby("visit", as_index=False).l2.mean(); traj = dict(zip(m.visit.astype(int), m.l2))
        slope, pval = (np.nan, np.nan)
        if len(m) >= 3:
            lr = stats.linregress(m.visit, m.l2); slope, pval = lr.slope, lr.pvalue  # separate 'adaptiveness' metric
        klass = hotelling_klass(w)                            # responder type = Hotelling T^2 (matches the graphs)
        amean, alo, ahi = classify(traj)
        return slope, pval, traj, klass, (amean, alo, ahi)

    def erp_ci(eps):  # CI band across epochs
        if not eps: return None
        M = np.stack(eps); n = len(M); mean = M.mean(0)
        sem = M.std(0, ddof=1) / np.sqrt(n) if n > 1 else np.zeros(EP); t = stats.t.ppf(.975, max(1, n - 1))
        return mean, mean - t * sem, mean + t * sem, n

    N2WIN_R, P3WIN_R = (200.0, 350.0), (300.0, 500.0)

    def erp_components(mean):  # peak N200 (neg, 200-350 ms) and P300 (pos, 300-500 ms) from a mean ERP
        if mean is None or not np.isfinite(mean).any(): return np.nan, np.nan
        n2 = mean[(tms >= N2WIN_R[0]) & (tms <= N2WIN_R[1])]; p3 = mean[(tms >= P3WIN_R[0]) & (tms <= P3WIN_R[1])]
        return (float(np.nanmin(n2)) if n2.size else np.nan), (float(np.nanmax(p3)) if p3.size else np.nan)

    def draw_acute(pdf, d, vs, title):  # No-Go + Go acute-effect trajectory over the supplied visits
        fig = styled(title, "Study Trajectories"); g2 = gsp(fig, 2, 3, hspace=.5, wspace=.36, bottom=0.13)
        base = d[(d.region == "frontal") & (d.n_epochs >= THR)]

        def acute(cond):
            s = base[base.condition == cond]
            w = s.pivot_table(index=["visit", "band"], columns="phase", values="abs_power_uv2").reindex(columns=["pre", "post"]).dropna(subset=["pre", "post"]) if len(s) else pd.DataFrame()
            return w.assign(l2=np.log2(w.post) - np.log2(w.pre)).reset_index() if len(w) else pd.DataFrame()
        wN, wG = acute("NoGo"), acute("Go")
        if len(wN) or len(wG):
            for i, b in enumerate(BANDS):
                ax = fig.add_subplot(g2[i // 3, i % 3]); ax.axhline(0, color="grey", ls=":", lw=.7)
                for wdat, col, mk, lab in [(wN, BCOL[b], "o", "No-Go"), (wG, "#9aa0a6", "s", "Go")]:
                    sb = wdat[wdat.band == b].sort_values("visit") if len(wdat) else pd.DataFrame()
                    if len(sb): ax.scatter(sb.visit, sb.l2, color=col, s=15, marker=mk, zorder=3, label=lab); fitln(ax, sb.visit.values, sb.l2.values, col)
                ax.set_title(b); ax.set_xlabel("visit"); ax.set_ylabel("log2(POST/PRE)"); ax.grid(alpha=.3)
                if i == 0: ax.legend(fontsize=7)
            axm = fig.add_subplot(g2[1, 1]); axm.axhline(0, color="grey", ls=":", lw=.7)
            for wdat, col, mk, lab in [(wN, "#333", "o", "No-Go"), (wG, "#9aa0a6", "s", "Go")]:
                if len(wdat):
                    gm = wdat.groupby("visit", as_index=False).l2.mean(); axm.scatter(gm.visit, gm.l2, color=col, s=18, marker=mk, zorder=3, label=lab); fitln(axm, gm.visit.values, gm.l2.values, col)
            axm.set_title("4-band mean"); axm.set_xlabel("visit"); axm.set_ylabel("log2(POST/PRE)"); axm.grid(alpha=.3); axm.legend(fontsize=7)
        else:
            fig.text(0.5, 0.5, "Not enough paired PRE+POST sessions (>= 20 clean epochs each) to plot adaptiveness.", ha="center", color=GREY, fontsize=11)
        caption(fig, "acute"); pdf.savefig(fig); plt.close()

    def draw_levels(pdf, d, title, vs=None):  # overall + each band on its own graph (PRE vs POST)
        fig = styled(title, "Study Trajectories"); g2 = gsp(fig, 3, 2, hspace=.55, wspace=.3, bottom=0.09)
        sub = d[(d.condition == "NoGo") & (d.region == "frontal")]
        axA = fig.add_subplot(g2[0, :])
        for ph, col, ls in [("pre", "#7f8c8d", "--"), ("post", "#22406b", "-")]:
            s = sub[sub.phase == ph].groupby("visit").abs_power_uv2.sum()
            if len(s): axA.plot(s.index, s.values, ls, marker="o", ms=4, color=col, label=ph.upper())
        axA.set_title("Overall power (all bands combined, frontal No-Go)"); axA.set_ylabel("power (uV^2)"); axA.set_xlabel("visit"); axA.legend(fontsize=8); axA.grid(alpha=.3)
        if vs is not None: axA.set_xticks(vs)
        for i, b in enumerate(BANDS):
            ax = fig.add_subplot(g2[1 + i // 2, i % 2])
            for ph, ls, a in [("post", "-", 1.0), ("pre", "--", 0.6)]:
                s = sub[(sub.phase == ph) & (sub.band == b)].groupby("visit").abs_power_uv2.mean()
                if len(s): ax.plot(s.index, s.values, ls, marker="o", ms=3, lw=1.3, color=BCOL[b], alpha=a, label=ph.upper())
            ax.set_title(f"{b}  (solid = POST, dashed = PRE)", fontsize=9); ax.set_xlabel("visit"); ax.set_ylabel("power (uV^2)"); ax.grid(alpha=.3)
            if vs is not None: ax.set_xticks(vs)
            if i == 0: ax.legend(fontsize=7)
        caption(fig, "eeg_levels"); pdf.savefig(fig); plt.close()

    def band_bars_pages(pdf, beP, beQ, labels):  # two pages: absolute + relative (% of total) band power, PRE vs POST
        def rel_epochs(be):  # convert per-epoch band powers -> % of total (per epoch)
            bands = [b for b in BANDS if b in be and be[b] is not None and len(np.asarray(be[b]))]
            if len(bands) < len(BANDS): return {}
            M = np.vstack([np.asarray(be[b], float) for b in BANDS]); tot = M.sum(0); tot = np.where(tot == 0, np.nan, tot)
            return {b: np.asarray(be[b], float) / tot * 100 for b in BANDS}
        for rel in (False, True):
            dP = rel_epochs(beP) if rel else beP; dQ = rel_epochs(beQ) if rel else beQ
            ttl = "Band power (relative, % of total): PRE vs POST (frontal No-Go)" if rel else "Band power (absolute): PRE vs POST (frontal No-Go)"
            fig = styled(ttl, "Task Analysis"); g = gsp(fig, 2, 2, hspace=.4, wspace=.3, bottom=0.14)
            for i, b in enumerate(BANDS):
                ax = fig.add_subplot(g[i // 2, i % 2]); ax.set_title(b)
                for j, (lab, be) in enumerate([("PRE", dP), ("POST", dQ)]):
                    arr = be.get(b) if be else None
                    arr = np.asarray(arr, float)[np.isfinite(np.asarray(arr, float))] if arr is not None else None
                    if arr is not None and len(arr):
                        m = arr.mean(); sem = arr.std(ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else 0.0
                        ax.bar(j, m, yerr=1.96 * sem, capsize=4, color=["#7f8c8d", BCOL[b]][j])
                    else: ax.text(j, 0, labels.get(lab, "missing"), ha="center", va="bottom", color=GREY, fontsize=7)
                ax.set_xticks([0, 1]); ax.set_xticklabels(["PRE", "POST"]); ax.set_ylabel("% of total power" if rel else "power (uV^2)"); ax.grid(alpha=.3, axis="y")
            caption(fig, "bandbars_rel" if rel else "bandbars"); pdf.savefig(fig); plt.close()

    # ==================================================================== COMPLETE (per participant)
    def build_complete_participant(pid):
        print(f"Building complete report for P{pid} (processing {len(eeg_files(pid))} recordings; may take ~30s)...", flush=True)
        nogo_pre, nogo_post = [], []; rows = []; bep = {b: [] for b in BANDS}; beq = {b: [] for b in BANDS}
        for f in eeg_files(pid):
            r = process(f)
            if not r: continue
            rows += r["rows"]
            (nogo_pre if r["phase"] == "pre" else nogo_post).extend(r["nogo"])
            for b, arr in (r.get("band_epochs") or {}).items():
                (bep if r["phase"] == "pre" else beq)[b].extend(list(arr))
        if not rows: print(f"P{pid}: no usable EEG for complete report"); return None
        d = pd.DataFrame(rows); vs = sorted(int(v) for v in d.visit.unique())
        slope, pval, traj, klass, acute_ci = slope_class(d)
        path = os.path.join(REP, f"P{pid}_complete.pdf")
        with PdfPages(path) as pdf:
            pdf.savefig(page_intro(f"Participant {pid} - complete report", f"All {len(vs)} visits combined     -     classified as: {klass}")); plt.close()
            # ERP PRE vs POST: line graph (waveform, CI across epochs) + N200/P300 peak bars
            fig = styled("ERP: PRE vs POST (frontal No-Go)", "Task Analysis")
            ax = fig.add_axes([0.10, 0.16, 0.55, 0.68])          # line graph (kept - most useful)
            ax.axhline(0, color="grey", lw=.6); ax.axvline(0, color="grey", ls=":", lw=.8)
            ax.axvspan(*N2WIN_R, color="#3498db", alpha=.06); ax.axvspan(*P3WIN_R, color="#e67e22", alpha=.06)
            means = {}
            for eps, col, lab in [(nogo_pre, "#7f8c8d", "PRE"), (nogo_post, "#c0392b", "POST")]:
                c = erp_ci(eps)
                if c: mean, lo, hi, n = c; means[lab] = mean; ax.plot(tms, mean, color=col, lw=2, label=f"{lab} (n={n})"); ax.fill_between(tms, lo, hi, color=col, alpha=.2)
            ax.set_xlabel("time from stimulus (ms)"); ax.set_ylabel("frontal amp (uV)"); ax.legend(); ax.grid(alpha=.3)
            axb = fig.add_axes([0.74, 0.16, 0.20, 0.68])         # N200/P300 peak-amplitude bars
            n2p, p3p = erp_components(means.get("PRE")); n2q, p3q = erp_components(means.get("POST"))
            axb.bar([0, 0.8], [n2p, n2q], width=.75, color=["#7f8c8d", "#c0392b"])
            axb.bar([2.2, 3.0], [p3p, p3q], width=.75, color=["#7f8c8d", "#c0392b"])
            axb.axhline(0, color="k", lw=.7); axb.set_xticks([0.4, 2.6]); axb.set_xticklabels(["N200", "P300"])
            axb.set_ylabel("peak amplitude (uV)"); axb.set_title("Components", fontsize=9); axb.grid(alpha=.3, axis="y")
            axb.legend(handles=[Line2D([0], [0], color="#7f8c8d", lw=6, label="PRE"), Line2D([0], [0], color="#c0392b", lw=6, label="POST")], fontsize=7, loc="best")
            caption(fig, "erp"); pdf.savefig(fig); plt.close()
            # band bars PRE vs POST (absolute + relative) with CI across epochs
            band_bars_pages(pdf, {b: np.array(bep[b], float) for b in BANDS}, {b: np.array(beq[b], float) for b in BANDS}, {"PRE": "missing", "POST": "missing"})
            draw_acute(pdf, d, vs, "Acute-effect trajectory (all visits)")
            draw_levels(pdf, d, "EEG power over visits (PRE vs POST)")
        print(f"wrote {os.path.basename(path)}  (class: {klass})")
        return path

    # ==================================================================== COMPLETE (group / all participants)
    def build_complete(pids):
        cdir = os.path.join(os.path.dirname(REP), "_complete_cache"); os.makedirs(cdir, exist_ok=True)

        def agg_pid(pid):  # aggregate one participant (cached to disk -> resumable across runs)
            fs = eeg_files(pid)
            key = hashlib.md5(("v3|" + "|".join(f"{os.path.basename(f)}:{os.path.getsize(f)}" for f in fs)).encode()).hexdigest()[:10]
            cf = os.path.join(cdir, f"P{pid}_{key}.pkl")
            if os.path.exists(cf):
                try:
                    with open(cf, "rb") as fh: return pickle.load(fh)
                except Exception: pass
            erpc = {"frontal": {"pre": [], "post": []}, "posterior": {"pre": [], "post": []}}; rows = []
            for f in fs:
                r = process(f)
                if not r: continue
                rows += r["rows"]; ph = r["phase"]
                if ph in ("pre", "post"):
                    erpc["frontal"][ph].extend(r.get("nogo", [])); erpc["posterior"][ph].extend(r.get("nogo_par", []))
            if not rows: return None
            d = pd.DataFrame(rows)
            slope, pval, traj, klass, acute_ci = slope_class(d)     # responder class stays frontal (matches the graphs)
            rec = dict(slope=slope, pval=pval, klass=klass, traj=traj)
            for reg, sfx in [("frontal", ""), ("posterior", "_par")]:  # store BOTH regions (P3/P4 = parietal); ERP + spectral
                gg = d[(d.condition == "NoGo") & (d.region == reg) & (d.n_epochs >= THR)]
                sub = d[(d.condition == "NoGo") & (d.region == reg)]
                rec[f"erp_pre{sfx}"] = erp_mean(erpc[reg]["pre"]); rec[f"erp_post{sfx}"] = erp_mean(erpc[reg]["post"])
                rec[f"pre_over{sfx}"] = {int(v): x for v, x in sub[sub.phase == "pre"].groupby("visit").abs_power_uv2.sum().items()}
                rec[f"post_over{sfx}"] = {int(v): x for v, x in sub[sub.phase == "post"].groupby("visit").abs_power_uv2.sum().items()}
                bp = {}
                for ph in ["pre", "post"]:
                    for b in BANDS: bp[(ph, b)] = gg[(gg.phase == ph) & (gg.band == b)].abs_power_uv2.mean()
                rec[f"bp{sfx}"] = bp
            try:
                with open(cf, "wb") as fh: pickle.dump(rec, fh)
            except Exception: pass
            return rec

        perp = {}  # pid -> aggregates
        print(f"Building group complete report for {len(pids)} participant(s). "
              f"The first run processes every recording and can take a few minutes...", flush=True)
        for i, pid in enumerate(pids, 1):
            print(f"  [{i}/{len(pids)}] processing participant P{pid}...", flush=True)
            rec = agg_pid(pid)
            if rec is not None: perp[pid] = rec
        if not perp: print("no data for complete report"); return None
        print("  rendering pages...", flush=True)
        path = os.path.join(REP, "complete_report_all.pdf")

        def ci(mat):  # participants x time -> mean, lo, hi
            mat = np.array([m for m in mat if np.isfinite(m).any()]); n = len(mat)
            mean = np.nanmean(mat, 0); sem = np.nanstd(mat, 0, ddof=1) / np.sqrt(n); t = stats.t.ppf(.975, max(1, n - 1))
            return mean, mean - t * sem, mean + t * sem, n

        def group(k): return [p for p in perp if perp[p]["klass"] == k]
        counts = {k: len(group(k)) for k in KORD}

        def mean_ci_by_visit(grp, field):  # across-participant mean + 95% CI at each visit
            acc = {}
            for p in grp:
                for v, val in perp[p][field].items():
                    if np.isfinite(val): acc.setdefault(int(v), []).append(val)
            xs = sorted(acc); mean = []; lo = []; hi = []
            for v in xs:
                a = np.array(acc[v], float); mu = a.mean()
                h = (stats.t.ppf(.975, len(a) - 1) * a.std(ddof=1) / np.sqrt(len(a))) if len(a) > 1 else 0.0
                mean.append(mu); lo.append(mu - h); hi.append(mu + h)
            return xs, np.array(mean), np.array(lo), np.array(hi)
        with PdfPages(path) as pdf:
            sub = f"Group-level results  -  n={len(perp)} participants\n" + "     ".join(f"{counts[k]} {k}s" for k in KORD)
            pdf.savefig(page_intro("Complete report - all participants", sub)); plt.close()
            # 1) group ERP PRE vs POST with CI across participants
            fig = styled("Group ERP: PRE vs POST (frontal No-Go)", "Overall"); ax = fig.add_axes([0.12, 0.16, 0.8, 0.68])
            for key, col, lab in [("erp_pre", "#7f8c8d", "PRE"), ("erp_post", "#c0392b", "POST")]:
                mean, lo, hi, n = ci([perp[p][key] for p in perp]); ax.plot(tms, mean, color=col, lw=2, label=f"{lab} (n={n})"); ax.fill_between(tms, lo, hi, color=col, alpha=.2)
            ax.axhline(0, color="grey", lw=.6); ax.axvline(0, color="grey", ls=":", lw=.8); ax.set_xlabel("time from stimulus (ms)"); ax.set_ylabel("frontal amp (uV)"); ax.legend(); ax.grid(alpha=.3)
            caption(fig, "group_prepost"); pdf.savefig(fig); plt.close()
            # 2) group band power PRE vs POST per band with CI -- absolute + relative
            def bp_val(p, ph, b, rel, sfx=""):
                v = perp[p][f"bp{sfx}"].get((ph, b), np.nan)
                if rel:
                    tot = np.nansum([perp[p][f"bp{sfx}"].get((ph, bb), np.nan) for bb in BANDS])
                    v = (v / tot * 100) if (np.isfinite(v) and tot > 0) else np.nan
                return v
            for rel in (False, True):
                ttl = "Group band power (relative, % of total): PRE vs POST" if rel else "Group band power (absolute): PRE vs POST"
                fig = styled(ttl, "Overall"); g = gsp(fig, 2, 2, hspace=.4, wspace=.3, bottom=0.13)
                for i, b in enumerate(BANDS):
                    ax = fig.add_subplot(g[i // 2, i % 2])
                    for j, ph in enumerate(["pre", "post"]):
                        vals = np.array([bp_val(p, ph, b, rel) for p in perp], float); vals = vals[np.isfinite(vals)]
                        mean = vals.mean() if len(vals) else np.nan; sem = vals.std(ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0
                        ax.bar(j, mean, yerr=1.96 * sem, color=["#7f8c8d", BCOL[b]][j], capsize=4)
                    ax.set_xticks([0, 1]); ax.set_xticklabels(["PRE", "POST"]); ax.set_title(b); ax.set_ylabel("% of total power" if rel else "power (uV^2)"); ax.grid(alpha=.3, axis="y")
                caption(fig, "group_bandrel" if rel else "group_bandabs"); pdf.savefig(fig); plt.close()
            # 3) adaptiveness per participant (single neutral color - NOT by responder type)
            fig = styled("Adaptiveness per participant", "Overall"); ax = fig.add_axes([0.12, 0.18, 0.8, 0.62])
            ids = sorted(p for p in perp if np.isfinite(perp[p]["slope"]))
            if not ids:
                ax.text(0.5, 0.5, "Not enough paired data to estimate adaptiveness (needs >= 3 visits\nwith both PRE and POST usable per participant).", ha="center", va="center", color=GREY, fontsize=11, transform=ax.transAxes)
            for i, p in enumerate(ids):
                ax.bar(i, perp[p]["slope"], color="#22406b")
            ax.axhline(0, color="k", lw=.8); ax.set_xticks(range(len(ids))); ax.set_xticklabels([f"P{p}" for p in ids])
            ax.set_ylabel("adaptiveness slope (trend of acute effect over visits)"); ax.grid(alpha=.3, axis="y")
            caption(fig, "adaptiveness_bar"); pdf.savefig(fig); plt.close()
            # 4) ERP PRE vs POST, separated by responder type
            fig = styled("ERP: PRE vs POST by responder type", "Overall"); g = gsp(fig, 2, 2, hspace=.5, wspace=.28, bottom=0.12); a0 = None
            panels = [("All participants", list(perp))] + [(k.capitalize() + "s", group(k)) for k in KORD]
            for idx, (title, grp) in enumerate(panels):
                ax = fig.add_subplot(g[idx // 2, idx % 2], sharex=a0, sharey=a0) if a0 else fig.add_subplot(g[idx // 2, idx % 2]); a0 = a0 or ax
                ax.axhline(0, color="grey", lw=.5); ax.axvline(0, color="grey", ls=":", lw=.7)
                for key, col, lab in [("erp_pre", "#7f8c8d", "PRE"), ("erp_post", "#c0392b", "POST")]:
                    arrs = [perp[p][key] for p in grp]
                    if not arrs: continue
                    mean, lo, hi, n = ci(arrs)
                    if np.isfinite(mean).any():
                        ax.plot(tms, mean, color=col, lw=1.6, label=f"{lab} (n={n})"); ax.fill_between(tms, lo, hi, color=col, alpha=.18)
                ax.set_title(f"{title} (n={len(grp)})", fontsize=9); ax.set_xlabel("time from stimulus (ms)"); ax.set_ylabel("frontal amp (uV)"); ax.grid(alpha=.3); ax.legend(fontsize=6)
            caption(fig, "erp_byclass"); pdf.savefig(fig); plt.close()
            # 5) per-participant trajectories: faded individuals + thick group-mean trendline with CI
            fig = styled("Per-participant acute-effect trajectories", "Overall"); g = gsp(fig, 2, 2, hspace=.45, wspace=.28, bottom=0.12); a0 = None
            panels = [("All participants", list(perp), "#333")] + [(k.capitalize() + "s", group(k), KCOL[k]) for k in KORD]
            for idx, (title, grp, tcol) in enumerate(panels):
                ax = fig.add_subplot(g[idx // 2, idx % 2], sharex=a0, sharey=a0) if a0 else fig.add_subplot(g[idx // 2, idx % 2]); a0 = a0 or ax
                ax.axhline(0, color="grey", ls=":", lw=.7); nlines = 0
                for p in grp:
                    t = perp[p]["traj"]
                    if not t: continue
                    nlines += 1; xs = sorted(t); ax.plot(xs, [t[v] for v in xs], "-", lw=0.8, color=KCOL[perp[p]["klass"]], alpha=.22)
                xs, mean, lo, hi = mean_ci_by_visit(grp, "traj")
                if len(xs):
                    ax.fill_between(xs, lo, hi, color=tcol, alpha=.18)
                    ax.plot(xs, mean, "-", lw=2.6, color=tcol, zorder=5, label="group mean +/- 95% CI")
                    ax.legend(fontsize=6, loc="best")
                ax.set_title(f"{title} (n={nlines})", fontsize=9); ax.set_xlabel("visit"); ax.set_ylabel("log2(POST/PRE)"); ax.grid(alpha=.3)
            allv = [v for p in perp for v in perp[p]["traj"].values() if np.isfinite(v)]
            rr = max(1.5, np.nanpercentile(np.abs(allv), 90) * 2.0) if allv else 2.0
            a0.set_ylim(-rr, rr)  # readable range; extreme small-n CI bands may clip
            caption(fig, "classify"); pdf.savefig(fig); plt.close()
            # 6) EEG PRE & POST over visits: mean + 95% CI, all + by class
            fig = styled("EEG PRE & POST over visits", "Overall"); g = gsp(fig, 2, 2, hspace=.45, wspace=.28, bottom=0.12); a0 = None
            panels = [("All participants", list(perp))] + [(k.capitalize() + "s", group(k)) for k in KORD]
            for idx, (title, grp) in enumerate(panels):
                ax = fig.add_subplot(g[idx // 2, idx % 2], sharex=a0, sharey=a0) if a0 else fig.add_subplot(g[idx // 2, idx % 2]); a0 = a0 or ax
                for which, col, lab in [("pre_over", "#7f8c8d", "PRE"), ("post_over", "#22406b", "POST")]:
                    xx, mean, lo, hi = mean_ci_by_visit(grp, which)
                    if len(xx): ax.plot(xx, mean, "-o", ms=3, color=col, label=lab); ax.fill_between(xx, lo, hi, color=col, alpha=.18)
                ax.set_title(f"{title} (n={len(grp)})", fontsize=9); ax.set_xlabel("visit"); ax.set_ylabel("power (uV^2)"); ax.grid(alpha=.3)
                if idx == 0: ax.legend(fontsize=8)
            allpw = [v for p in perp for fld in ("pre_over", "post_over") for v in perp[p][fld].values() if np.isfinite(v)]
            top = max(50.0, np.nanpercentile(allpw, 92) * 1.8) if allpw else 300.0
            a0.set_ylim(0, top)  # power >= 0; readable range, extreme small-n CI bands may clip
            caption(fig, "over_visits"); pdf.savefig(fig); plt.close()

            # ================= additional pages: band bars by responder type + PARIETAL analyses =================
            def _bandbars_byclass(sfx, rlab):   # band power PRE vs POST, All + each responder type
                fig = styled(f"Band power by responder type ({rlab} No-Go)", "Overall"); g = gsp(fig, 2, 2, hspace=.5, wspace=.3, bottom=0.12)
                panels = [("All participants", list(perp))] + [(k.capitalize() + "s", group(k)) for k in KORD]
                x = np.arange(len(BANDS)); wbar = 0.38
                for idx, (title, grp) in enumerate(panels):
                    ax = fig.add_subplot(g[idx // 2, idx % 2])
                    for j, ph in enumerate(["pre", "post"]):
                        ms, es = [], []
                        for b in BANDS:
                            vals = np.array([bp_val(p, ph, b, False, sfx) for p in grp], float); vals = vals[np.isfinite(vals)]
                            ms.append(vals.mean() if len(vals) else np.nan); es.append(1.96 * vals.std(ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0)
                        ax.bar(x + (j - 0.5) * wbar, ms, wbar, yerr=es, capsize=3, color=["#7f8c8d", "#c0392b"][j], label=["PRE", "POST"][j])
                    ax.set_xticks(x); ax.set_xticklabels(BANDS, fontsize=8); ax.set_title(f"{title} (n={len(grp)})", fontsize=9); ax.set_ylabel("power (uV^2)"); ax.grid(alpha=.3, axis="y")
                    if idx == 0: ax.legend(fontsize=7)
                caption(fig, "group_bandabs"); pdf.savefig(fig); plt.close()

            def _erp_all(sfx, rlab):            # group ERP PRE vs POST, all participants
                fig = styled(f"Group ERP: PRE vs POST ({rlab} No-Go)", "Overall"); ax = fig.add_axes([0.12, 0.16, 0.8, 0.68])
                for key, col, lab in [(f"erp_pre{sfx}", "#7f8c8d", "PRE"), (f"erp_post{sfx}", "#c0392b", "POST")]:
                    mean, lo, hi, n = ci([perp[p][key] for p in perp]); ax.plot(tms, mean, color=col, lw=2, label=f"{lab} (n={n})"); ax.fill_between(tms, lo, hi, color=col, alpha=.2)
                ax.axhline(0, color="grey", lw=.6); ax.axvline(0, color="grey", ls=":", lw=.8); ax.set_xlabel("time from stimulus (ms)"); ax.set_ylabel(f"{rlab} amp (uV)"); ax.legend(); ax.grid(alpha=.3)
                caption(fig, "group_prepost"); pdf.savefig(fig); plt.close()

            def _erp_byclass(sfx, rlab):        # group ERP PRE vs POST, All + each responder type
                fig = styled(f"ERP: PRE vs POST by responder type ({rlab} No-Go)", "Overall"); g = gsp(fig, 2, 2, hspace=.5, wspace=.28, bottom=0.12); a0 = None
                panels = [("All participants", list(perp))] + [(k.capitalize() + "s", group(k)) for k in KORD]
                for idx, (title, grp) in enumerate(panels):
                    ax = fig.add_subplot(g[idx // 2, idx % 2], sharex=a0, sharey=a0) if a0 else fig.add_subplot(g[idx // 2, idx % 2]); a0 = a0 or ax
                    ax.axhline(0, color="grey", lw=.5); ax.axvline(0, color="grey", ls=":", lw=.7)
                    for key, col, lab in [(f"erp_pre{sfx}", "#7f8c8d", "PRE"), (f"erp_post{sfx}", "#c0392b", "POST")]:
                        arrs = [perp[p][key] for p in grp]
                        if not arrs: continue
                        mean, lo, hi, n = ci(arrs)
                        if np.isfinite(mean).any(): ax.plot(tms, mean, color=col, lw=1.6, label=f"{lab} (n={n})"); ax.fill_between(tms, lo, hi, color=col, alpha=.18)
                    ax.set_title(f"{title} (n={len(grp)})", fontsize=9); ax.set_xlabel("time from stimulus (ms)"); ax.set_ylabel(f"{rlab} amp (uV)"); ax.grid(alpha=.3); ax.legend(fontsize=6)
                caption(fig, "erp_byclass"); pdf.savefig(fig); plt.close()

            def _bandbars_all(sfx, rlab):       # group band power PRE vs POST, all participants
                fig = styled(f"Group band power (absolute): PRE vs POST ({rlab} No-Go)", "Overall"); g = gsp(fig, 2, 2, hspace=.4, wspace=.3, bottom=0.13)
                for i, b in enumerate(BANDS):
                    ax = fig.add_subplot(g[i // 2, i % 2])
                    for j, ph in enumerate(["pre", "post"]):
                        vals = np.array([bp_val(p, ph, b, False, sfx) for p in perp], float); vals = vals[np.isfinite(vals)]
                        mean = vals.mean() if len(vals) else np.nan; sem = vals.std(ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0
                        ax.bar(j, mean, yerr=1.96 * sem, color=["#7f8c8d", BCOL[b]][j], capsize=4)
                    ax.set_xticks([0, 1]); ax.set_xticklabels(["PRE", "POST"]); ax.set_title(b); ax.set_ylabel("power (uV^2)"); ax.grid(alpha=.3, axis="y")
                caption(fig, "group_bandabs"); pdf.savefig(fig); plt.close()

            def _over_visits(sfx, rlab):        # PRE & POST power over visits, All + each responder type
                fig = styled(f"EEG PRE & POST over visits ({rlab})", "Overall"); g = gsp(fig, 2, 2, hspace=.45, wspace=.28, bottom=0.12); a0 = None
                panels = [("All participants", list(perp))] + [(k.capitalize() + "s", group(k)) for k in KORD]
                for idx, (title, grp) in enumerate(panels):
                    ax = fig.add_subplot(g[idx // 2, idx % 2], sharex=a0, sharey=a0) if a0 else fig.add_subplot(g[idx // 2, idx % 2]); a0 = a0 or ax
                    for which, col, lab in [(f"pre_over{sfx}", "#7f8c8d", "PRE"), (f"post_over{sfx}", "#22406b", "POST")]:
                        xx, mean, lo, hi = mean_ci_by_visit(grp, which)
                        if len(xx): ax.plot(xx, mean, "-o", ms=3, color=col, label=lab); ax.fill_between(xx, lo, hi, color=col, alpha=.18)
                    ax.set_title(f"{title} (n={len(grp)})", fontsize=9); ax.set_xlabel("visit"); ax.set_ylabel("power (uV^2)"); ax.grid(alpha=.3)
                    if idx == 0: ax.legend(fontsize=8)
                caption(fig, "over_visits"); pdf.savefig(fig); plt.close()

            _bandbars_byclass("", "Frontal")     # NEW: frontal band bars split by responder type
            _erp_all("_par", "Parietal")         # NEW: parietal analyses (ERP + spectral, all + by responder type)
            _erp_byclass("_par", "Parietal")
            _bandbars_all("_par", "Parietal")
            _bandbars_byclass("_par", "Parietal")
            _over_visits("_par", "Parietal")
        print(f"wrote {os.path.basename(path)}  ({', '.join(f'{counts[k]} {k}' for k in KORD)})")
        return path

    # ---------- router ----------
    def _isall(x): return isinstance(x, str) and x.strip().lower() in ("all", "a", "*")

    def _visits_of(pid):
        return sorted({meta(os.path.basename(f))[1] for f in eeg_files(pid) if meta(os.path.basename(f))[1]})

    def _daily_batch(pids, vsel):
        for pid in pids:
            vs = _visits_of(pid) if (_isall(vsel) or vsel is None) else (list(vsel) if isinstance(vsel, (list, tuple, set)) else [vsel])
            for v in vs:
                p = build_daily(pid, v)
                if p: written.append(p)

    def _weekly_batch(pids, widxs):
        for pid in pids:
            for wi in widxs:
                p = build_weekly(pid, wi)
                if p: written.append(p)

    if kind is None and ask:
        try:
            kind = input("Report type - daily / weekly / complete / all: ").strip().lower()
        except Exception:
            kind = "complete"
    kind = (kind or "complete").lower()
    written = []
    ALLW = list(range(len(WEEKS)))

    # ==================================================================== EXPORT (added: dump numbers to CSV)
    # Runs the SAME process() used by every report, but instead of drawing PDFs it writes the
    # underlying numbers to two tidy CSV files so they can be analysed independently.
    if kind == "export":
        sess_recs, epoch_recs = [], []
        pids = all_pids if (participant is None or _isall(participant)) else [int(participant)]
        todo = [(pid, f) for pid in pids for f in eeg_files(pid)]
        print(f"Exporting band power from {len(todo)} recording(s)...", flush=True)
        for i, (pid, f) in enumerate(todo, 1):
            print(f"  [{i}/{len(todo)}] P{pid}  {os.path.basename(f)}", flush=True)
            r = process(f)
            if not r:
                print("        -> no usable EEG, skipped", flush=True); continue
            v = r["visit"]
            wk = next((n for n, lo, hi in WEEKS if v is not None and lo <= v <= hi), "unassigned")
            src = os.path.basename(f)
            for row in r["rows"]:
                sess_recs.append(dict(participant=pid, week=wk, source_file=src, **row))
            for b, arr in (r.get("band_epochs") or {}).items():
                for k, val in enumerate(np.asarray(arr, float)):
                    epoch_recs.append(dict(participant=pid, visit=v, phase=r["phase"], week=wk,
                                           condition="NoGo", region="frontal", band=b,
                                           epoch=k, abs_power_uv2=float(val), source_file=src))
        outdir = os.path.abspath(out_dir); os.makedirs(outdir, exist_ok=True)
        p1 = os.path.join(outdir, "band_power_session.csv")
        p2 = os.path.join(outdir, "band_power_epochs.csv")
        pd.DataFrame(sess_recs).to_csv(p1, index=False)
        pd.DataFrame(epoch_recs).to_csv(p2, index=False)
        print(f"\nwrote {p1}   ({len(sess_recs)} rows)")
        print(f"wrote {p2}   ({len(epoch_recs)} rows)")
        return [p1, p2]

    if kind == "all":  # every report type, every participant, every visit / week
        _daily_batch(all_pids, "all")
        _weekly_batch(all_pids, ALLW)
        for pid in all_pids:
            p = build_complete_participant(pid)
            if p: written.append(p)
        p = build_complete(all_pids)
        if p: written.append(p)
    elif kind.startswith("d"):
        if participant is None and ask:
            print("participants available:", all_pids, " (or type 'all')")
            try:
                r = input("Which participant # (or 'all'): ").strip().lower()
                participant = "all" if _isall(r) else (int(r) if r.isdigit() else all_pids[0])
            except Exception: participant = all_pids[0]
        pids = all_pids if (_isall(participant) or participant is None) else [participant]
        vsel = visit
        if _isall(participant): vsel = "all"
        elif visit is None and ask:
            byv = {}
            for f in eeg_files(participant):
                _, v, ph = meta(os.path.basename(f))
                if v: byv.setdefault(v, set()).add(ph)
            vlist = sorted(byv); paired = sorted(v for v in vlist if "pre" in byv[v] and "post" in byv[v])
            print(f"visits for P{participant}: {vlist}")
            print(f"visits with BOTH pre & post recorded (recommended): {paired}   (or type 'all')")
            try:
                r = input("Which visit # (or 'all'): ").strip().lower()
                vsel = "all" if _isall(r) else (int(r) if r.isdigit() else ((paired or vlist)[len(paired or vlist) // 2] if (paired or vlist) else None))
            except Exception: vsel = (paired or vlist)[len(paired or vlist) // 2] if (paired or vlist) else None
        _daily_batch(pids, vsel)
    elif kind.startswith("w"):
        if participant is None and ask:
            print("participants available:", all_pids, " (or type 'all')")
            try:
                r = input("Which participant # (or 'all'): ").strip().lower()
                participant = "all" if _isall(r) else (int(r) if r.isdigit() else all_pids[0])
            except Exception: participant = all_pids[0]
        pids = all_pids if (_isall(participant) or participant is None) else [participant]
        if _isall(participant) or _isall(week):
            widxs = ALLW
        elif week is None and ask:
            print("weeks:  " + "   ".join(f"{i+1}={n} (V{lo}-V{hi})" for i, (n, lo, hi) in enumerate(WEEKS)) + "   (or type 'all')")
            try:
                r = input("Which week # (1-4 or 'all'): ").strip().lower()
                widxs = ALLW if _isall(r) else ([min(max(int(r) - 1, 0), len(WEEKS) - 1)] if r.isdigit() else [0])
            except Exception: widxs = [0]
        else:
            widxs = [min(max((week or 1) - 1, 0), len(WEEKS) - 1)]
        _weekly_batch(pids, widxs)
    else:  # complete: individual (per participant) or group (all participants)
        if participant is None and ask:
            print("participants available:", all_pids)
            print("Enter a participant # for an INDIVIDUAL complete report, or 'all' for the GROUP (all-participants) report.")
            try:
                r = input("Participant # or 'all': ").strip().lower()
                participant = int(r) if r.isdigit() else "all"
            except Exception: participant = "all"
        if _isall(participant) or participant is None:
            p = build_complete(all_pids)
            if p: written.append(p)
        else:
            p = build_complete_participant(participant)
            if p: written.append(p)

    # ---------- bundle + links ----------
    if not written: return []
    zip_path = os.path.join(REP, zip_name)
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in written: zf.write(p, os.path.basename(p))
    print(f"bundled {len(written)} report(s) -> {zip_path}")
    if show_links:
        try:
            from IPython import get_ipython
            from IPython.display import display, HTML
            if get_ipython() is not None:
                import base64
                def _dl(pth, lab):
                    with open(pth, "rb") as fh: b64 = base64.b64encode(fh.read()).decode()
                    return f'<a download="{os.path.basename(pth)}" href="data:application/octet-stream;base64,{b64}">{lab}</a>'
                rows = [_dl(zip_path, "Download all (zip)")] + [_dl(p, os.path.basename(p)) for p in written]
                display(HTML("<div style='line-height:1.9;font-family:sans-serif'>" + "<br>".join(rows) + "</div>"))
        except Exception:
            pass
    return written


if __name__ == "__main__":
    generate_reports()