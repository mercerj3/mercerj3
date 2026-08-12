# Exercise Vitals — quick analyses

`exercise_vitals_analyses.R` assumes your data frame is called `exercise`:

```r
exercise <- Exercise_Vitals_1_9v2
```

Packages: `tidyverse`, `lme4`, `lmerTest`, `broom`, `broom.mixed`.

## What's in the file

| § | Analysis |
|---|----------|
| 0 | Cleaning: drop the dictionary row, coerce numerics, trim the ` 009` PIN, derive MAP / pulse pressure / day-in-study, build a tidy `paired` frame |
| 1 | QC: visit inventory, missingness, lbs↔kg conversion check, physiologic range flags, within-person outliers |
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
   All 187 real rows pass the conversion check.
