# Ionized calcium report

`ionized_calcium_report.R` builds the multi-page landscape PDF from a
blood-gas export.

## Run it

```r
# from the folder that contains the workbook
source("R/ionized_calcium_report.R")
```

or from a shell:

```sh
Rscript R/ionized_calcium_report.R
```

The script reads `DATA_FILE` relative to the **current working directory**.
Either `setwd()` to the folder holding the workbook first, or put an absolute
path in the config block:

```r
DATA_FILE <- "~/Downloads/icalpin14.xlsx"        # or a full path
```

`file.exists(DATA_FILE)` should return `TRUE` before you run anything else.

## Packages

`readxl`, `dplyr`, `tidyr`, `ggplot2`, `scales`, `grid`. `lme4` only if
`FIT_MODEL <- TRUE`.

```r
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "scales"))
```

## Configuration

Everything that changes between runs is in the config block at the top:

| Setting | What it does |
| --- | --- |
| `DATA_FILE`, `SHEET`, `OUT_PDF` | input workbook, sheet, output PDF |
| `WEEK_MAP` | which visit belongs to which week block (the shaded bands) |
| `PH_MODE` | `"exclude"` (default, matches the reference report) or `"sensitivity"`, which adds a pH page and a pH-normalised iCa page |
| `PH_REF`, `PH_SLOPE` | reference pH and the ~5%-per-0.1-unit slope used for normalisation |
| `FIT_MODEL` | fit `iCa ~ phase + (1 | PIN)`; off by default because it is not estimable at this sample size |
| `MEASURES_USED` | analytes carried through the report |
| `FCOHB_FLAG`, `DPV_EXTREME` | thresholds used on the checks page |

## Expected columns

`PIN`, `Visit_Number`, `Pre_Post`, `Date (yyyy-mm-dd)`, `Ionized_Calcium`,
`pH`, `ctHb`, `Hct`, `cNa+`, `cCl-`, `FCOHb`, `FMetHb`, `Notes`.

Column matching ignores case, spaces and punctuation, so `cCl-` / `cCl−` /
`cCl` all resolve. Everything is read as text and coerced, because the export
stores the same column as a mix of numbers and strings.

`Hct` is detected as percent or fraction from its median and converted to a
fraction before it enters Dill-Costill.

## Method

```
%dPV (Dill-Costill) = 100*[(Hb_pre/Hb_post)*((1-Hct_post)/(1-Hct_pre)) - 1]
%dPV (Hb only)      = 100*[(Hb_pre/Hb_post) - 1]
iCa_PVadj,post      = iCa_post * (1 + %dPV/100)
```

Raw iCa is the concentration relevant to calcium-sensing receptor signalling.
PV-adjusted iCa estimates whether vascular calcium *content* changed
independently of haemoconcentration. Both are reported because they answer
different questions.
