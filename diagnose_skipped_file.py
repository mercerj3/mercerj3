"""Diagnose why recordings were reported as "no usable EEG, skipped".

Paste everything below the closing quotes into a notebook cell (DATA_FOLDER must
already be set by Step 1), list the recordings you want checked, and run.

For each file it answers three questions, and the FIRST one that comes back zero
is the cause:
  1. does the Marker column actually use codes 5 (Go) and 6 (No-Go)?
  2. do those markers land on rows that carry EEG (not blank rows)?
  3. of the epochs that form, how many survive the 100 uV artifact test?

A summary table of every file checked is printed at the end.
"""
import numpy as np, pandas as pd, zipfile, glob, os, re
from scipy import signal

# ---- the recordings to check, as (participant, visit, phase) ----------------
# Matching is done on the parsed filename, exactly the way the pipeline does it,
# so capitalisation ("POST" vs "Post") and the date/time suffix do not matter.
FILES_TO_CHECK = [
    (1, 11, "PRE"),
    (1, 8, "POST"),
    (1, 9, "POST"),
    (5, 23, "POST"),
    (5, 7, "POST"),
    (7, 16, "POST"),
    (7, 9, "PRE"),
    (9, 19, "POST"),
    (9, 19, "PRE"),
]

CH = ["F3", "F4", "C3", "Cz", "C4", "P3", "P4"]
REGIONS = [("frontal", ["F3", "F4"]), ("central", ["C3", "Cz", "C4"]), ("posterior", ["P3", "P4"])]
FS, NPRE, NPOST, ART = 250, 50, 200, 100.0
b_bp, a_bp = signal.butter(4, [0.1 / (FS / 2), 30 / (FS / 2)], btype="band")
b_n, a_n = signal.iirnotch(60.0, 30.0, FS)

# ---- build the file list once, the same way the pipeline scans -------------
ALL_FILES = sorted({f for pat in ("**/*Go No-Go*.csv", "**/*Go No-Go*.csv.zip")
                    for f in glob.glob(os.path.join(DATA_FOLDER, pat), recursive=True)
                    if re.search(r"PIN\s*\d+", os.path.basename(f), re.I)})


def find_file(pin, visit, phase):
    """Locate one recording by participant / visit / phase, ignoring case."""
    hits = []
    for f in ALL_FILES:
        up = os.path.basename(f).upper()
        p, v = re.search(r"PIN\s*(\d+)", up), re.search(r"\bV(\d+)", up)
        ph = "POST" if "POST" in up else ("PRE" if "PRE" in up else "NA")
        if p and v and int(p.group(1)) == pin and int(v.group(1)) == visit and ph == phase.upper():
            hits.append(f)
    # the pipeline prefers a plain .csv over a .csv.zip when both exist
    hits.sort(key=lambda x: x.lower().endswith(".zip"))
    return hits[0] if hits else None


def diagnose(pin, visit, phase):
    label = f"P{pin} V{visit} {phase.upper()}"
    path = find_file(pin, visit, phase)
    if path is None:
        print(f"!! no file found for {label}\n")
        return dict(file=label, cause="file not found")

    print("=" * 72)
    print(f"{label}   {os.path.basename(path)}\n")
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
        print("\n   *** CAUSE 1: this session did not use marker codes 5/6 ***\n")
        return dict(file=label, cause="1: wrong marker codes", markers=0)

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
        print("\n   *** CAUSE 2: markers sit on rows where the EEG cells are blank ***\n")
        return dict(file=label, cause="2: markers on blank rows", markers=n_go + n_nogo)

    print("\n3) Of the epochs that form, how many are clean enough?  (every region)")
    print(f"   {'region':10s} {'epochs':>7} {'passed':>7} {'median uV':>10} {'best uV':>9} {'worst uV':>9}")
    total_kept, med_front, best_front = 0, np.nan, np.nan
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
            if region == "frontal":
                med_front, best_front = float(np.median(ptps)), float(np.min(ptps))
        else:
            print(f"   {region:10s} {total:>7} {kept_ep:>7}        no epochs formed")
    print(f"\n   threshold is {ART:.0f} uV; a region needs >= 1 epoch to yield any number,"
          f"\n   and >= 20 for the analysis to use it")
    if total_kept == 0:
        print("\n   *** CAUSE 3: every epoch on every region exceeded the threshold ***\n")
        cause = "3: all epochs over threshold"
    else:
        print(f"\n   {total_kept} epoch(s) survived somewhere - this file is not empty\n")
        cause = f"not empty ({total_kept} epochs survived)"
    return dict(file=label, cause=cause, markers=n_go + n_nogo,
                frontal_median_uV=round(med_front) if np.isfinite(med_front) else np.nan,
                frontal_best_uV=round(best_front) if np.isfinite(best_front) else np.nan)


summary = pd.DataFrame([diagnose(*spec) for spec in FILES_TO_CHECK])
print("=" * 72)
print("SUMMARY\n")
print(summary.to_string(index=False))
