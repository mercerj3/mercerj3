"""Diagnose why recordings were reported as "no usable EEG, skipped".

Paste everything below the closing quotes into a notebook cell (DATA_FOLDER must
already be set by Step 1), list the recordings you want checked, and run.

It answers three questions per file, and the FIRST one that comes back zero is
the cause:
  1. does the Marker column actually use codes 5 (Go) and 6 (No-Go)?
  2. do those markers land on rows that carry EEG (not blank rows)?
  3. of the epochs that form, how many survive the 100 uV artifact test?
"""
import numpy as np, pandas as pd, zipfile, glob, os
from scipy import signal

# ---- list every recording you want to check (as it appears in the skip message)
FILES_TO_CHECK = ["*V8 POST*", "*V9 POST*", "*V11 PRE*"]

CH = ["F3", "F4", "C3", "Cz", "C4", "P3", "P4"]
REGIONS = [("frontal", ["F3", "F4"]), ("central", ["C3", "Cz", "C4"]), ("posterior", ["P3", "P4"])]
FS, NPRE, NPOST, ART = 250, 50, 200, 100.0
b_bp, a_bp = signal.butter(4, [0.1 / (FS / 2), 30 / (FS / 2)], btype="band")
b_n, a_n = signal.iirnotch(60.0, 30.0, FS)


def diagnose(pattern):
    hits = sorted(glob.glob(os.path.join(DATA_FOLDER, "**", pattern), recursive=True))
    if not hits:
        print(f"!! no file matching {pattern}\n"); return
    path = hits[0]
    print("=" * 70)
    print(os.path.basename(path), "\n")

    if path.endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            df = pd.read_csv(z.open(z.namelist()[0]))
    else:
        df = pd.read_csv(path)

    tcol = [c for c in df.columns if c.strip().lower().startswith("time")][0]
    mk = pd.to_numeric(df["Marker"], errors="coerce").fillna(0).to_numpy()
    t_all = pd.to_numeric(df[tcol], errors="coerce").to_numpy(float)
    n_go, n_nogo = int((mk == 5).sum()), int((mk == 6).sum())

    print("1) Marker values present in this file:")
    print(pd.Series(mk).value_counts().sort_index().to_string())
    print(f"\n   pipeline needs 5 (Go) and 6 (No-Go)  ->  found {n_go} Go, {n_nogo} No-Go")
    if n_go + n_nogo == 0:
        print("\n   *** CAUSE 1: this session did not use marker codes 5/6 ***\n"); return

    print("\n2) Do those markers land on rows that have EEG?")
    surviving = 0
    for region, chans in REGIONS:
        chans = [c for c in chans if c in df]
        if not chans:
            print(f"   {region:10s}: channels missing from file"); continue
        good = (df[chans].apply(pd.to_numeric, errors="coerce").notna().all(axis=1).to_numpy()
                & np.isfinite(t_all))
        kept = int(((mk == 5) | (mk == 6))[good].sum())
        surviving = max(surviving, kept)
        print(f"   {region:10s}: {int(good.sum()):>7,} usable samples | "
              f"{kept:>4} of {n_go + n_nogo} markers survive")
    if surviving == 0:
        print("\n   *** CAUSE 2: markers sit on rows where the EEG cells are blank ***\n"); return

    print("\n3) Of the epochs that form, how many are clean enough?  (every region)")
    print(f"   {'region':10s} {'epochs':>7} {'passed':>7} {'median uV':>10} {'best uV':>9} {'worst uV':>9}")
    total_kept = 0
    for region, chans in REGIONS:
        chans = [c for c in chans if c in df]
        if not chans:
            print(f"   {region:10s}  channels missing from file"); continue
        good = (df[chans].apply(pd.to_numeric, errors="coerce").notna().all(axis=1).to_numpy()
                & np.isfinite(t_all))
        if good.sum() < 500:
            print(f"   {region:10s}  too few usable samples"); continue
        ts = t_all[good]
        X = df[chans].apply(pd.to_numeric, errors="coerce").to_numpy(float)[good]
        m = mk[good]
        o = np.argsort(ts, kind="stable"); ts, X, m = ts[o], X[o], m[o]
        filt = np.empty_like(X)
        for k in range(X.shape[1]):
            filt[:, k] = signal.filtfilt(b_n, a_n, signal.filtfilt(b_bp, a_bp, X[:, k]))
        total = kept_ep = 0; ptps = []
        for code, t in zip(m, ts):
            if code not in (5, 6): continue
            s = int(np.searchsorted(ts, t)); a, bb = s - NPRE, s + NPOST
            if a < 0 or bb > len(filt): continue
            total += 1
            # the pipeline rejects on the LARGEST swing across the region's channels
            seg = filt[a:bb, :] - filt[a:bb, :][:NPRE].mean(0, keepdims=True)
            p = float(np.ptp(seg, 0).max()); ptps.append(p)
            if p <= ART: kept_ep += 1
        total_kept += kept_ep
        if ptps:
            print(f"   {region:10s} {total:>7} {kept_ep:>7} {np.median(ptps):>10.0f} "
                  f"{np.min(ptps):>9.0f} {np.max(ptps):>9.0f}")
        else:
            print(f"   {region:10s} {total:>7} {kept_ep:>7}        no epochs formed")
    print(f"\n   threshold is {ART:.0f} uV; a region needs >= 1 epoch to yield any number,"
          f"\n   and >= 20 for the analysis to use it")
    if total_kept == 0:
        print("\n   *** CAUSE 3: every epoch on every region exceeded the threshold ***\n")
    else:
        print(f"\n   {total_kept} epoch(s) survived somewhere - this file is not empty\n")


for pattern in FILES_TO_CHECK:
    diagnose(pattern)
