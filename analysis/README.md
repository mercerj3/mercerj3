# Exercise Vitals — quick analyses

`exercise_vitals_analyses.R` assumes your data frame is called `exercise`:

```r
exercise <- Exercise_Vitals_1_9v2
```

Packages: `tidyverse`, `lme4`, `lmerTest`, `broom`, `broom.mixed`.

Three files, each runnable on its own:

| file | what it answers |
|---|---|
| `exercise_vitals_analyses.R` | the full set — QC, descriptives, pre→post, trends, ICC, EKG, FeNO, hydration, anthropometrics |
| `weight_check.R` | how much weight each participant loses per session, and how resting weight moves across the study |
| `weight_unit_discrepancy.R` | how far apart the lbs and kg columns are, and whether that disagreement is systematic |

## What's in the file

| § | Analysis |
|---|----------|
| 0 | Cleaning: drop the dictionary row, coerce numerics, trim the ` 009` PIN, derive MAP / pulse pressure / day-in-study / weight change, build a tidy `paired` frame |
| 1 | QC: visit inventory, missingness, lbs↔kg agreement (descriptive), physiologic range flags, within-person outliers |
| 2 | Descriptives overall and per participant |
| 3 | Pre→post response to a bout (naive paired *t* vs random-intercept model) |
| 4 | Training effect: resting vitals across `Visit_Num` |
| 5 | ICC — how much variance is the person vs the day |
| 6 | EKG: rhythm/result distributions, QTc formula verification, Bazett vs Fridericia, safety screen |
| 7 | FeNO pre→post, attempt-to-attempt reliability (Bland–Altman), ATS bands |
| 8 | Hydration: urine specific gravity, sweat loss, effect on HR response |
| 9 | Anthropometrics (baseline only, n = 9) — L/R symmetry |
| 10 | Figures: paired slopes, change-score boxplots, timelines, correlation heatmap |

## Data shape

- **214 visit rows, 9 participants (001–009), 23–25 visits each.** Repeated
  measures — cluster on `PIN` rather than treating rows as independent.
- Pre/post pairs exist for weight, RR, BP, HR, temp, SpO₂, FeNO (~166 complete
  pairs for the vitals, ~69 for FeNO).
- EKG on 124 of 214 visits; anthropometrics once per participant; ionized
  calcium on <10 visits; `Heart_Rate_Repeat_Pre` is entirely empty.

## Things to fix or confirm before publishing

1. **`PIN = " 009"`** on one row — leading space, would become a 10th
   participant in any `group_by()`. §0 trims it.
2. **`Respirations_Pre = 98`** (PIN 008, visit 11) — not a respiratory rate.
3. **`Temperature_Post` is labelled Celsius** in the dictionary row, but the
   values run 96.8–99.6 °F. The label looks wrong, not the data — confirm.
4. **`Weight_conv_*`** are QC ratio formulas (lbs/kg ≈ 2.2046), not
   measurements. They spill `#DIV/0!` down ~360 empty rows below the data;
   import with `na = c("", "NA", "#DIV/0!")` or §0 will coerce them to `NA`.

## How weight is handled

Every weight is used exactly as recorded. Nothing is corrected, excluded, or
back-converted on the basis of the lbs/kg ratio.

Weight appears twice at each timepoint (lbs and kg), and the two columns don't
always imply the same weight — median disagreement 0.07–0.09 lb, max 1.34 lb.
Rather than judge one column right, §0 computes the pre→post change down each
unit separately and averages the two:

```r
Weight_delta_from_lbs = (Weight_lbs_post - Weight_lbs_pre) / 2.20462
Weight_delta_from_kg  =  Weight_kg_post  - Weight_kg_pre
Weight_delta_kg       = mean of the two   # na.rm, so one unit alone still works
```

The disagreement largely cancels, because whatever offset separates the two
columns is present in both the pre and the post value of the same column. On
the 176 visits with all four weights:

| route | mean change | SD |
|---|---|---|
| via lbs | −1.2420 kg | 0.5223 |
| via kg | −1.2438 kg | 0.5206 |
| **averaged** | **−1.2429 kg** | **0.5188** |

The two routes correlate at r = 0.98 and their means differ by 0.002 kg, so
averaging is a small variance reduction, not a correction. §1c prints this
table so the choice is defensible in a methods section.

The same averaging is applied to the anthropometrics in §9, which are likewise
recorded in both inches and centimetres.
