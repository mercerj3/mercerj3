# Workday plan: your own spectral analysis on top of the lab pipeline

**Your goal:** for each participant, trend frontal delta / theta / alpha / beta power across
visits, PRE vs POST exercise.

**Why you can't do that today with what you have:** your colleague's notebook is one
1,000-line function that goes *raw files → PDF*. Every number it computes lives inside that
function and is thrown away once the figure is drawn. There is no CSV, no table, nothing you
can plot yourself. That is the whole gap — not the signal processing, which is already done
and already correct.

**The fix is small.** The pipeline's internal `process()` already builds exactly the table you
want, one row per *participant × visit × phase × condition × region × band*. It just never
saves it. So one addition — `kind="export"` — writes it to CSV, and from there everything you
want is ordinary spreadsheet-style work on 8 columns.

**Ground rule for the day:** don't reimplement any signal processing. Every number you
report should come out of your colleague's code, so your results and her reports can never
disagree.

---

## What is already prepared for you in this repo

| file | what it is |
|---|---|
| `eeg_report.py` | Your colleague's pipeline, **unchanged**, plus one added `kind="export"` branch that writes the numbers to CSV. All three of her report types (`daily` / `weekly` / `complete`) were re-run after the change and produce identical output. |
| `eeg_trends.py` | New analysis helpers — loading, quality checks, trend figures, statistics. Short, commented, meant to be read. |
| `spectral_trends.ipynb` | The notebook you actually run, 12 steps, top to bottom. |

All three were tested end to end on synthetic recordings built to match your file format.

---

## Block 0 — Orientation (20 min)

1. Put `eeg_report.py`, `eeg_trends.py` and `spectral_trends.ipynb` in the same folder.
2. Open the notebook, run **Step 1** and **Step 2**.
3. Stop and read Step 2's output carefully.

Step 2 lists every participant, every visit found, and which visits have **both** PRE and
POST. Do not move on until that list matches what you believe you collected. A visit missing
here is a filename problem — the pipeline needs `PIN <n>`, `V<n>`, `PRE`/`POST` and the
literal text `Go No-Go` in the filename. Finding this now costs a minute; finding it after a
40-minute processing run costs 40 minutes.

**Done when:** you can say out loud how many usable PRE+POST visits each participant has.

---

## Block 1 — Understand what the pipeline actually did (45 min)

This is the block that makes the rest of the day defensible. Read the cheat sheet at the
bottom of this document ("What the pipeline does, precisely") against your own description of
the pipeline. Your summary was close — it is worth knowing the four places it wasn't:

- **The stats you described aren't in this notebook.** You said Friedman test. There is no
  Friedman test anywhere in her code. Her reports use a **Hotelling T² test** across the four
  bands to classify responders, plus a **linear regression** of the acute effect on visit
  number. If someone told you Friedman, it came from a different script. (A real Friedman test
  is provided in Step 10 of your notebook.)
- **Epochs are 1.0 s, and they start *before* the stimulus** — 200 ms before to 800 ms after.
  So "post-stimulus band power" actually includes 200 ms of pre-stimulus signal.
- **"Average power across epochs" happens in a specific order**: per-channel spectra first,
  then average across epochs, then average across the region's channels. Pooling channels
  before the FFT would keep only phase-locked activity and would be wrong for band power —
  her code is careful about this, and the comment at lines 167–172 says so.
- **0.48 Hz bin spacing is right** (250/512 = 0.4883 Hz), but bin spacing is not resolution.
  See Block 5.

**Done when:** you can point at the line that sets each parameter you'd have to report in a
methods section.

---

## Block 2 — Get the numbers out (45 min, mostly waiting)

Run **Step 3** with `participant=1` first. Inspect the CSV. Then change it to
`participant="all"` and re-run.

Expect roughly a second or two per recording — a few minutes for the full study. If it is
dramatically slower, the bottleneck is reading the raw CSVs, not the processing.

Then run **Step 4** and read the `inventory()` table.

**`visits_usable` is the most important number of your day.** It counts visits where *both*
PRE and POST survived with ≥20 clean epochs. It is the ceiling on every statistic that
follows. If a participant has 4 usable visits, no test in Step 9 can return a meaningful
p-value for them, and you should decide now to describe them descriptively rather than
statistically.

**Done when:** `band_power_session.csv` exists and you have written down each participant's
`visits_usable`.

---

## Block 3 — Your actual goal: the trend figures (60 min)

Run **Step 5** (prepare) then **Step 6** (per-participant figures).

Step 5 does four things — read its description in the notebook, especially the fourth:
it adds **relative power**, each band's share of total 0.5–30 Hz power.

**Run both the absolute and the relative version, and compare them.** Absolute power depends
on electrode contact, hair and gel, which change between visits. A rise in absolute power
across visits can be a better-seated electrode rather than a change in the brain. Relative
power divides that out. An effect that appears in both is credible; one that appears only in
absolute is a lead, not a finding.

Then run **Step 8** for the one-band-all-participants view, which is where between-person
consistency becomes visible.

**Done when:** you have looked at all four bands for every participant, in both absolute and
relative, and can name which participants show a consistent alpha trend.

---

## Block 4 — The acute effect (45 min)

Run **Step 7**. This pairs each visit's PRE with its POST and computes **log2(POST/PRE)**:
`0` = no change, `+1` = doubled, `−1` = halved.

The log matters. A rise 10→20 and a fall 20→10 are the same size change, but as raw
differences they're `+10` and `−10` against different baselines, so averaging them across
visits is meaningless. In log2 they are exactly `+1` and `−1`.

This is the same quantity your colleague's "adaptiveness" figures use, so your numbers should
line up with her weekly reports — a good cross-check. If they don't, it's because her acute
figures additionally require ≥20 clean epochs on *both* phases; `prepare()` applies the same
filter, so they should agree.

**Done when:** you can say, per participant, whether the acute effect grows or shrinks across
the study.

---

## Block 5 — Statistics (45 min)

Run **Step 9** (per participant) and **Step 10** (group).

Both tests are non-parametric, which is the right call: band power is skewed and you have few
visits, so a t-test's normality assumption isn't safe.

Three things to hold onto:

- **Read `n_visits` before any p-value.** With fewer than 6 visits the Wilcoxon test *cannot*
  return p < 0.03 even if every visit moves the same direction. A non-significant result there
  is a statement about your sample size, not about the brain.
- **Friedman needs complete blocks.** Every participant must have every visit. The helper
  keeps the largest usable run and reports what it kept — read `n_participants` and `visits`
  before the p-value.
- **PRE vs POST is not a Friedman question.** With two related conditions Friedman degenerates
  to a sign test; Wilcoxon signed-rank is correct, which is what Step 10's second cell uses.

**Done when:** `per_participant_stats.csv` exists and you've flagged which rows have enough
visits to interpret.

---

## Block 6 — Write down the limitations while they're fresh (30 min)

Read **Step 11** in the notebook and copy the relevant parts into your own notes. The one that
most affects your stated goal:

**Frontal delta is the weakest of your four bands.** With 1.0 s epochs and `NW=2`, the
multitaper estimate smooths the spectrum by ±2 Hz — an effective resolution of about **4 Hz**,
which is *wider than the entire delta band* (0.5–4 Hz). In practice delta is 7 frequency bins
spanning 0.98–3.91 Hz, heavily blended with theta. You can still trend it, but never
interpret a delta change without checking whether theta moved the same way, and expect a
reviewer to ask about this.

Alpha (10 bins) and theta (8 bins) are the well-resolved bands. If you need one headline
result, take it from alpha.

Then run **Step 12** to bundle everything into a zip.

---

## If you finish early

In rough order of value:

1. **Re-run Step 5 with `region="posterior"`.** Alpha is largest over parietal sites, so
   posterior alpha is often the cleaner signal — one line change, whole new result.
2. **Compare Go against No-Go** (`condition="Go"` in Step 5). Her reports centre on No-Go;
   a change appearing in both is general arousal, one specific to No-Go is inhibitory control.
3. **Use `band_power_epochs.csv`** for within-session error bars — it holds every individual
   clean epoch, so you can put a confidence interval on a single session rather than a point.
4. **Add week as a factor** — the `week` column already labels Control / Heat 1–3, which is
   probably the comparison your study is actually built around.

---

# What the pipeline does, precisely

Line numbers refer to `eeg_report.py`. FS = 250 Hz throughout.

| # | Step | Exactly what happens | Settings | Line |
|---|---|---|---|---|
| 1 | Find files | Recursive search for `*Go No-Go*.csv`/`.csv.zip` with `PIN <n>` in the name; participant, visit and PRE/POST parsed from the filename | — | 108, 132 |
| 2 | Read | Keeps `Time`, `Marker`, and F3 F4 C3 Cz C4 P3 P4 | — | 148 |
| 3 | Bandpass | 4th-order Butterworth, **0.1–30 Hz**, applied forwards and backwards (`filtfilt`) so no phase distortion | `butter(4, [0.1, 30])` | 98, 190 |
| 4 | Notch | 60 Hz mains notch, Q = 30 | `iirnotch(60, 30)` | 99 |
| 5 | Regions | frontal = F3,F4 · central = C3,Cz,C4 · posterior = P3,P4. Each region handled independently, so a dead parietal channel never kills frontal | — | 91, 182 |
| 6 | Epoching | Around each `Marker` 5 (Go) / 6 (No-Go): **50 samples before + 200 after = 250 samples = −200 to +800 ms** | `NPRE=50, NPOST=200` | 89, 199 |
| 7 | Baseline | Subtract the mean of the 50 pre-stimulus samples, per channel | — | 202 |
| 8 | Artifact rejection | Drop the epoch if peak-to-peak **> 100 µV** on *any* channel in that region | `ART=100` | 89, 203 |
| 9 | Spectrum | Multitaper (DPSS), **NW = 2.0, 3 tapers**, zero-padded 250 → **512** points → bin spacing **250/512 = 0.4883 Hz** | `NFFT=512, NW=2, KT=3` | 89, 100, 160 |
| 10 | Averaging | Mean across clean epochs **first**, then average the region's channels. Channels are pooled *after* the transform — pooling before it would keep only phase-locked activity and understate band power | — | 211 |
| 11 | PSD scaling | Divide by FS, then double every bin except DC and Nyquist (one-sided spectrum) | — | 211 |
| 12 | Band power | Sum the PSD bins inside each band × bin width → `abs_power_uv2` | delta 0.5–4 · theta 4–8 · alpha 8–13 · beta 13–30 Hz | 93, 214 |
| 13 | Her statistics | **Hotelling T²** across the four bands to classify responders; **linear regression** of acute effect on visit number for "adaptiveness". No Friedman test | `THR=20` min epochs | 512, 525 |

### Things worth knowing that the table doesn't show

- **No frequency bin lands exactly on a band edge**, so bands don't overlap and nothing is
  double-counted. Actual coverage: delta 0.98–3.91 Hz (7 bins) · theta 4.39–7.81 (8) ·
  alpha 8.30–12.70 (10) · beta 13.18–29.79 (35).
- **Beta sits on the filter edge.** The bandpass ends at 30 Hz, where the filter has already
  halved the power. Absolute beta is therefore underestimated — equally at every visit, so
  trends are fine, but don't quote an absolute beta value as unattenuated.
- **Bin spacing ≠ resolution.** 0.4883 Hz is the spacing after zero-padding; the true
  resolution is set by the taper — ±2 Hz, i.e. about 4 Hz wide. Zero-padding interpolates the
  spectrum, it does not add information.
- **Duplicate recordings** are de-duplicated per (visit, phase), preferring the plain `.csv`
  over a `.csv.zip`.
