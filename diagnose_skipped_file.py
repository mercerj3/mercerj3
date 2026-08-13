"""Diagnose why a recording was reported as "no usable EEG, skipped".

Paste the body of this file into a notebook cell (DATA_FOLDER must already be
set by Step 1), change FILE_PATTERN to the recording you want to check, and run.

It answers three questions in order, and the FIRST one that comes back zero is
your cause:
  1. does the Marker column actually use codes 5 (Go) and 6 (No-Go)?
  2. do those markers land on rows that carry EEG (not blank rows)?
  3. of the epochs that form, how many survive the 100 uV artifact test?
"""
# ---- Why was this recording skipped? ----
import numpy as np, pandas as pd, zipfile, glob, os
from scipy import signal

FILE_PATTERN = "*V8 POST*"          # <-- change to the file you want to check

CH   = ["F3", "F4", "C3", "Cz", "C4", "P3", "P4"]
FS, NPRE, NPOST, ART = 250, 50, 200, 100.0
b_bp, a_bp = signal.butter(4, [0.1 / (FS / 2), 30 / (FS / 2)], btype="band")
b_n,  a_n  = signal.iirnotch(60.0, 30.0, FS)

path = sorted(glob.glob(os.path.join(DATA_FOLDER, "**", FILE_PATTERN), recursive=True))[0]
print("checking:", os.path.basename(path), "\n")

if path.endswith(".zip"):
    with zipfile.ZipFile(path) as z:
        df = pd.read_csv(z.open(z.namelist()[0]))
else:
    df = pd.read_csv(path)

tcol = [c for c in df.columns if c.strip().lower().startswith("time")][0]
mk   = pd.to_numeric(df["Marker"], errors="coerce").fillna(0).to_numpy()

print("1) What values are in the Marker column?")
vc = pd.Series(mk).value_counts().sort_index()
print(vc.to_string(), "\n")
print("   the pipeline only uses Marker == 5 (Go) and Marker == 6 (No-Go)")
print(f"   Go markers   (5): {int((mk == 5).sum())}")
print(f"   NoGo markers (6): {int((mk == 6).sum())}\n")

print("2) Do those markers land on rows that HAVE eeg?")
t_all = pd.to_numeric(df[tcol], errors="coerce").to_numpy(float)
for region, chans in [("frontal", ["F3", "F4"]), ("central", ["C3", "Cz", "C4"]), ("posterior", ["P3", "P4"])]:
    chans = [c for c in chans if c in df]
    if not chans:
        print(f"   {region:10s}: channels missing from file"); continue
    good = df[chans].apply(pd.to_numeric, errors="coerce").notna().all(axis=1).to_numpy() & np.isfinite(t_all)
    kept = int(((mk == 5) | (mk == 6))[good].sum())
    print(f"   {region:10s}: {int(good.sum()):>7,} usable samples | "
          f"{kept:>4} of {int(((mk==5)|(mk==6)).sum())} markers survive")

print("\n3) Of the surviving epochs, how many are clean enough to keep?")
chans = [c for c in ["F3", "F4"] if c in df]
good  = df[chans].apply(pd.to_numeric, errors="coerce").notna().all(axis=1).to_numpy() & np.isfinite(t_all)
if good.sum() > 500 and chans:
    ts = t_all[good]; X = df[chans].apply(pd.to_numeric, errors="coerce").to_numpy(float)[good]; m = mk[good]
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
        seg = filt[a:bb, :] - filt[a:bb, :][:NPRE].mean(0, keepdims=True)
        p = float(np.ptp(seg, 0).max()); ptps.append(p)
        if p <= ART: kept_ep += 1
    print(f"   frontal epochs formed : {total}")
    print(f"   passed the 100 uV test: {kept_ep}   (need >= 1 to keep the file, >= 20 for analysis)")
    if ptps:
        print(f"   typical swing per epoch: median {np.median(ptps):.0f} uV, "
              f"worst {np.max(ptps):.0f} uV")
else:
    print("   not enough usable frontal samples to test")
