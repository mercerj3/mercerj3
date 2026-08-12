# =============================================================================
# Exercise Vitals (PIN 001-009) -- quick analyses
#
# Assumes your data frame is called `exercise`
#   exercise <- Exercise_Vitals_1_9v2
#
# Structure of the data:
#   214 visit rows, 9 participants, 23-25 visits each  -> REPEATED MEASURES.
#   Almost every column is measured pre- and post-exercise, so the two
#   natural questions are (a) what does a bout do (pre -> post), and
#   (b) what does the program do (drift across Visit_Num).
#
# Section 0 is cleaning; sections 1-10 are the analyses. Run in order.
# =============================================================================

library(tidyverse)
library(lme4)
library(lmerTest)   # adds p-values to lmer summaries
library(broom)
library(broom.mixed)

theme_set(theme_minimal(base_size = 12))


# =============================================================================
# 0. CLEAN
# =============================================================================

# Row 2 of the spreadsheet is a data dictionary ("Numeric (mmHg)", etc.), not
# data. If your import kept it, this drops it. Also drops the ~360 trailing
# rows that only contain the Weight_conv formula spilling #DIV/0!.
exercise <- exercise %>%
  filter(!is.na(PIN), PIN != "Participant ID (001-009)")

# Anything that arrived as text because of the dictionary row / #DIV/0! gets
# forced back to numeric. `suppressWarnings` swallows the "NAs introduced"
# note, which is exactly what we want for "#DIV/0!".
num_vars <- c(
  "Visit_Num", "Weight_lbs_pre", "Weight_kg_pre", "Weight_conv_pre",
  "Respirations_Pre", "Systolic_Pre", "Diastolic_Pre", "Heart_Rate_Pre",
  "Temperature_Pre", "Oxygen_Saturation_Pre", "Urinalysis_USG",
  "EKG_HR", "PR_Interval", "QRS_Duration", "QT_Interval", "QTc_Interval",
  "P_Axis", "R_Axis", "T_Axis",
  "FeNO_Attempt1_Pre", "FeNO_Attempt1_Post",
  "FeNO_Attempt2_Pre", "FeNO_Attempt2_Post",
  "Weight_lbs_post", "Weight_kg_post", "Weight_conv_post",
  "Respirations_Post", "Systolic_Post", "Diastolic_Post", "Heart_Rate_Post",
  "Temperature_Post", "Oxygen_Saturation_Post"
)

exercise <- exercise %>%
  mutate(
    across(any_of(num_vars), ~ suppressWarnings(as.numeric(as.character(.x)))),
    # one PIN is stored as " 009" with a leading space -> would be a 10th
    # participant in every group_by() if left alone
    PIN        = factor(str_trim(as.character(PIN))),
    Visit_Date = as.Date(Visit_Date),
    EKG_Rhythm = factor(EKG_Rhythm, levels = 1:6, labels = c(
      "NSR", "Atrial arr.", "Sinus tachy", "Vent. arr.",
      "Sinus brady", "Undetermined")),
    EKG_Results = factor(EKG_Results, levels = c(1, 2, 3, 5), labels = c(
      "Normal", "Abnormal NCS", "Abnormal CS", "Unable to evaluate")),
    # derived
    MAP_Pre    = Diastolic_Pre  + (Systolic_Pre  - Diastolic_Pre)  / 3,
    MAP_Post   = Diastolic_Post + (Systolic_Post - Diastolic_Post) / 3,
    Pulse_Pressure_Pre  = Systolic_Pre  - Diastolic_Pre,
    Pulse_Pressure_Post = Systolic_Post - Diastolic_Post
  ) %>%
  group_by(PIN) %>%
  arrange(Visit_Date, .by_group = TRUE) %>%
  mutate(Day_In_Study = as.numeric(Visit_Date - min(Visit_Date, na.rm = TRUE))) %>%
  ungroup()

# The pre/post columns, paired up once so every later section can reuse them.
pair_spec <- tribble(
  ~measure,        ~pre,                     ~post,                     ~unit,
  "Heart rate",    "Heart_Rate_Pre",         "Heart_Rate_Post",         "bpm",
  "Systolic BP",   "Systolic_Pre",           "Systolic_Post",           "mmHg",
  "Diastolic BP",  "Diastolic_Pre",          "Diastolic_Post",          "mmHg",
  "MAP",           "MAP_Pre",                "MAP_Post",                "mmHg",
  "SpO2",          "Oxygen_Saturation_Pre",  "Oxygen_Saturation_Post",  "%",
  "Respirations",  "Respirations_Pre",       "Respirations_Post",       "rpm",
  "Temperature",   "Temperature_Pre",        "Temperature_Post",        "degF",
  "Weight",        "Weight_kg_pre",          "Weight_kg_post",          "kg",
  "FeNO",          "FeNO_Attempt1_Pre",      "FeNO_Attempt1_Post",      "ppb"
)

paired <- pair_spec %>%
  pmap(function(measure, pre, post, unit) {
    exercise %>%
      transmute(PIN, Visit_Num, Visit_Date, Day_In_Study,
                measure = measure, unit = unit,
                pre = .data[[pre]], post = .data[[post]])
  }) %>%
  list_rbind() %>%
  filter(!is.na(pre), !is.na(post)) %>%
  mutate(delta = post - pre,
         measure = factor(measure, levels = pair_spec$measure))


# =============================================================================
# 1. QC / DATA-INTEGRITY PASS  (do this before you believe any p-value)
# =============================================================================

# 1a. Visit inventory -- who has how many visits, over what window, and are
#     Visit_Num values contiguous?
exercise %>%
  group_by(PIN) %>%
  summarise(
    n_visits    = n(),
    first_visit = min(Visit_Date, na.rm = TRUE),
    last_visit  = max(Visit_Date, na.rm = TRUE),
    span_days   = as.numeric(max(Visit_Date, na.rm = TRUE) - min(Visit_Date, na.rm = TRUE)),
    visit_max   = max(Visit_Num, na.rm = TRUE),
    gaps        = paste(setdiff(seq_len(max(Visit_Num, na.rm = TRUE)), Visit_Num),
                        collapse = ", ")
  ) %>%
  print(n = Inf)

# 1b. Missingness by column -- tells you which analyses are actually powered.
#     Expect: anthropometrics ~11 rows (baseline only), ionized calcium <10,
#     Heart_Rate_Repeat_Pre completely empty, EKG ~124 of 214.
exercise %>%
  summarise(across(everything(), ~ sum(!is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "n_present") %>%
  mutate(pct_present = round(100 * n_present / nrow(exercise), 1)) %>%
  arrange(n_present) %>%
  print(n = Inf)

# 1c. Unit-conversion check. Weight_conv_* is lbs/kg and should sit at 2.2046.
#     Anything off means one of the two weights was mistyped.
exercise %>%
  mutate(ratio_pre  = Weight_lbs_pre  / Weight_kg_pre,
         ratio_post = Weight_lbs_post / Weight_kg_post) %>%
  filter(abs(ratio_pre - 2.20462) > 0.01 | abs(ratio_post - 2.20462) > 0.01) %>%
  select(PIN, Visit_Num, Weight_lbs_pre, Weight_kg_pre, ratio_pre,
         Weight_lbs_post, Weight_kg_post, ratio_post)

# 1d. Physiologic range flags. Known hit: PIN 008 visit 11 has
#     Respirations_Pre = 98, which is not a respiratory rate.
range_rules <- tribble(
  ~column,                   ~lo,  ~hi,
  "Heart_Rate_Pre",           30,  120,
  "Heart_Rate_Post",          40,  200,
  "Systolic_Pre",             80,  180,
  "Diastolic_Pre",            40,  110,
  "Respirations_Pre",          8,   30,
  "Respirations_Post",         8,   40,
  "Oxygen_Saturation_Pre",    90,  100,
  "Oxygen_Saturation_Post",   88,  100,
  "Temperature_Pre",          95,  101,
  "Temperature_Post",         95,  101,
  "Urinalysis_USG",        1.000, 1.035
)

range_rules %>%
  pmap(function(column, lo, hi) {
    exercise %>%
      filter(.data[[column]] < lo | .data[[column]] > hi) %>%
      transmute(PIN, Visit_Num, column = column, value = .data[[column]],
                expected = paste0(lo, "-", hi))
  }) %>%
  list_rbind()

# 1e. Within-person outliers: >4 SD from that participant's own mean. Catches
#     values that are physiologically plausible but wrong for this person.
paired %>%
  group_by(PIN, measure) %>%
  mutate(z = abs(pre - mean(pre)) / sd(pre)) %>%
  filter(z > 4) %>%
  select(PIN, Visit_Num, measure, pre, z) %>%
  arrange(desc(z))

# NOTE, not a code problem: the data dictionary labels Temperature_Post as
# Celsius, but the values run 96.8-99.6, i.e. Fahrenheit like Temperature_Pre.
# The dictionary is wrong, not the data -- but confirm before publishing.


# =============================================================================
# 2. DESCRIPTIVES
# =============================================================================

# 2a. Overall, per measure.
paired %>%
  group_by(measure, unit) %>%
  summarise(
    n_visits = n(), n_pin = n_distinct(PIN),
    pre_mean  = mean(pre),  pre_sd  = sd(pre),
    post_mean = mean(post), post_sd = sd(post),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# 2b. Per participant -- resting (pre) values. This is the table that shows you
#     how much of the total variance is just "which person is it".
paired %>%
  group_by(measure, PIN) %>%
  summarise(n = n(), mean = mean(pre), sd = sd(pre),
            cv_pct = 100 * sd(pre) / mean(pre), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  pivot_wider(id_cols = PIN, names_from = measure, values_from = mean)


# =============================================================================
# 3. THE MAIN EFFECT: pre -> post response to a bout
#
# Do NOT just run t.test() on all 214 rows -- the rows are ~24 per person, so
# a paired t-test treats one participant's 24 visits as 24 independent facts
# and understates the SE. Run both and compare; the mixed model is the one to
# report.
# =============================================================================

# 3a. Naive paired t-test (for contrast only).
paired %>%
  group_by(measure, unit) %>%
  summarise(broom::tidy(t.test(post, pre, paired = TRUE)), .groups = "drop") %>%
  select(measure, unit, mean_delta = estimate, conf.low, conf.high, p.value) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3)))

# 3b. Random-intercept model on the change score -- correct clustering.
delta_models <- paired %>%
  nest(data = -c(measure, unit)) %>%
  mutate(
    fit = map(data, ~ lmer(delta ~ 1 + (1 | PIN), data = .x)),
    out = map(fit, ~ broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE))
  ) %>%
  unnest(out) %>%
  select(measure, unit, mean_delta = estimate, conf.low, conf.high, p.value) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3)))

delta_models

# Expected shape of the result: HR ~ +21 bpm, SpO2 ~ -1.7%, weight ~ -1.2 kg
# (fluid loss), BP roughly flat / slightly down. If BP comes out flat, that is
# a real finding for a post-exercise measurement taken after recovery, not a
# null result to bury.

# 3c. Is the HR response stable across the study, or does it shrink? LRT of
#     adding a visit-number term. Significant => the same bout is producing a
#     different HR rise later in the program.
hr <- filter(paired, measure == "Heart rate")
m0 <- lmer(delta ~ 1 + (1 | PIN), data = hr, REML = FALSE)
m1 <- lmer(delta ~ 1 + Visit_Num + (1 | PIN), data = hr, REML = FALSE)
anova(m0, m1)

# Per-participant mean HR response, ranked.
hr %>%
  group_by(PIN) %>%
  summarise(n = n(), mean_delta = mean(delta), sd = sd(delta)) %>%
  arrange(desc(mean_delta))


# =============================================================================
# 4. TRAINING EFFECT: does resting physiology drift across the program?
#
# If this is a training study, the headline is resting HR / BP / weight
# trending over Visit_Num. Random slope lets each person have their own trend.
# =============================================================================

trend_fit <- function(var) {
  d <- exercise %>% select(PIN, Visit_Num, y = all_of(var)) %>% drop_na()
  m <- lmer(y ~ Visit_Num + (Visit_Num | PIN), data = d)
  broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE) %>%
    filter(term == "Visit_Num") %>%
    mutate(outcome = var, .before = 1)
}

# If a random-slope model fails to converge (only 9 clusters), fall back to
# (1 | PIN) -- swap the formula above.
c("Heart_Rate_Pre", "Systolic_Pre", "Diastolic_Pre", "Weight_kg_pre",
  "Oxygen_Saturation_Pre", "Respirations_Pre") %>%
  map(trend_fit) %>%
  list_rbind() %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3)))

# Same question for the *response* rather than the resting value: does the
# same bout produce a smaller HR rise as people get fitter?
lmer(delta ~ Visit_Num + (Visit_Num | PIN), data = hr) %>% summary()

ggplot(exercise, aes(Visit_Num, Heart_Rate_Pre, colour = PIN)) +
  geom_point(alpha = .6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .6) +
  labs(title = "Resting HR across visits", x = "Visit", y = "HR (bpm)")


# =============================================================================
# 5. HOW MUCH IS THE PERSON vs THE DAY?  (ICC)
#
# With 24 visits per person you can actually estimate this well. A high ICC
# means a single visit characterises the participant; a low one means you need
# repeats and single-timepoint screening would be noise.
# =============================================================================

icc_of <- function(var) {
  d <- exercise %>% select(PIN, y = all_of(var)) %>% drop_na()
  vc <- as.data.frame(VarCorr(lmer(y ~ 1 + (1 | PIN), data = d)))
  between <- vc$vcov[vc$grp == "PIN"]
  within  <- vc$vcov[vc$grp == "Residual"]
  tibble(outcome = var,
         sd_between = sqrt(between), sd_within = sqrt(within),
         ICC = between / (between + within))
}

c("Heart_Rate_Pre", "Systolic_Pre", "Diastolic_Pre", "Weight_kg_pre",
  "Temperature_Pre", "Oxygen_Saturation_Pre", "Urinalysis_USG",
  "QTc_Interval", "FeNO_Attempt1_Pre") %>%
  map(icc_of) %>%
  list_rbind() %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))


# =============================================================================
# 6. EKG
#
# 124 tracings. Heavily skewed: ~88% sinus bradycardia, ~84% read as
# "Abnormal NCS", resting EKG HR mean ~51 bpm. That is the athletic-heart
# pattern, and it is the most interesting descriptive finding in the file.
# =============================================================================

exercise %>% count(EKG_Rhythm) %>% mutate(pct = round(100 * n / sum(n), 1))
exercise %>% count(EKG_Results) %>% mutate(pct = round(100 * n / sum(n), 1))
exercise %>% count(PIN, EKG_Results) %>% pivot_wider(names_from = EKG_Results,
                                                     values_from = n, values_fill = 0)

# 6a. Verify the reported QTc. The vendor used Bazett (QT / sqrt(RR)); this
#     reproduces QTc_Interval to <1 ms on average, so any row that disagrees by
#     more than a few ms is a transcription error worth chasing.
ekg <- exercise %>%
  filter(!is.na(QT_Interval), !is.na(EKG_HR)) %>%
  mutate(
    RR         = 60 / EKG_HR,
    QTc_Bazett = QT_Interval / sqrt(RR),
    QTc_Frid   = QT_Interval / RR^(1/3),     # Fridericia -- better at low HR
    resid      = QTc_Bazett - QTc_Interval
  )

summary(ekg$resid)
ekg %>% filter(abs(resid) > 5) %>%
  select(PIN, Visit_Num, EKG_HR, QT_Interval, QTc_Interval, QTc_Bazett, resid)

# 6b. Bazett over-corrects at high HR and under-corrects at low HR, and this
#     cohort is bradycardic -- so check whether QTc still tracks HR. A non-zero
#     slope means the correction is not doing its job here and Fridericia is
#     the better summary.
summary(lmer(QTc_Interval ~ EKG_HR + (1 | PIN), data = ekg))
summary(lmer(QTc_Frid     ~ EKG_HR + (1 | PIN), data = ekg))

ekg %>%
  pivot_longer(c(QTc_Interval, QTc_Frid), names_to = "formula", values_to = "qtc") %>%
  ggplot(aes(EKG_HR, qtc, colour = formula)) +
  geom_point(alpha = .6) + geom_smooth(method = "lm") +
  labs(title = "QTc vs heart rate: is the correction working?",
       x = "EKG HR (bpm)", y = "QTc (ms)")

# 6c. Safety screen -- prolongation and conduction intervals.
ekg %>%
  summarise(
    n              = n(),
    qtc_over_450   = sum(QTc_Interval > 450, na.rm = TRUE),
    qtc_over_500   = sum(QTc_Interval > 500, na.rm = TRUE),
    pr_over_200    = sum(PR_Interval  > 200, na.rm = TRUE),
    qrs_over_120   = sum(QRS_Duration > 120, na.rm = TRUE),
    qtc_max        = max(QTc_Interval, na.rm = TRUE)
  )

# 6d. Within-person QTc stability -- day-to-day repeatability of the interval.
ekg %>%
  group_by(PIN) %>%
  summarise(n = n(), mean_qtc = mean(QTc_Interval), sd_qtc = sd(QTc_Interval),
            range = diff(range(QTc_Interval)))


# =============================================================================
# 7. FeNO (airway inflammation)
# =============================================================================

# 7a. Pre -> post, clustered.
feno <- filter(paired, measure == "FeNO")
summary(lmer(delta ~ 1 + (1 | PIN), data = feno))

# FeNO is right-skewed; log-transform and report a ratio instead of a
# difference if the residuals look bad.
feno %>%
  mutate(log_ratio = log(post / pre)) %>%
  { summary(lmer(log_ratio ~ 1 + (1 | PIN), data = .)) }

# 7b. Measurement reliability: attempt 1 vs attempt 2, Bland-Altman.
#     Only ~26 visits have a second attempt, so treat this as exploratory.
ba <- exercise %>%
  filter(!is.na(FeNO_Attempt1_Pre), !is.na(FeNO_Attempt2_Pre)) %>%
  mutate(avg = (FeNO_Attempt1_Pre + FeNO_Attempt2_Pre) / 2,
         dif = FeNO_Attempt1_Pre - FeNO_Attempt2_Pre)

ba %>% summarise(n = n(), bias = mean(dif), sd = sd(dif),
                 loa_lo = mean(dif) - 1.96 * sd(dif),
                 loa_hi = mean(dif) + 1.96 * sd(dif),
                 r = cor(FeNO_Attempt1_Pre, FeNO_Attempt2_Pre))

ggplot(ba, aes(avg, dif)) +
  geom_point() +
  geom_hline(yintercept = c(mean(ba$dif),
                            mean(ba$dif) + 1.96 * sd(ba$dif),
                            mean(ba$dif) - 1.96 * sd(ba$dif)),
             linetype = c("solid", "dashed", "dashed")) +
  labs(title = "FeNO attempt 1 vs 2 (Bland-Altman)",
       x = "Mean of attempts (ppb)", y = "Attempt 1 - Attempt 2 (ppb)")

# 7c. FeNO clinical bands (ATS): <25 low, 25-50 intermediate, >50 high.
exercise %>%
  filter(!is.na(FeNO_Attempt1_Pre)) %>%
  mutate(band = cut(FeNO_Attempt1_Pre, c(-Inf, 25, 50, Inf),
                    labels = c("<25", "25-50", ">50"))) %>%
  count(PIN, band) %>%
  pivot_wider(names_from = band, values_from = n, values_fill = 0)


# =============================================================================
# 8. HYDRATION
#
# Urine specific gravity is collected pre only, and weight loss across the bout
# is a sweat-loss proxy. Are people showing up to sessions dehydrated, and does
# arriving dehydrated change the HR response?
# =============================================================================

exercise %>%
  filter(!is.na(Urinalysis_USG)) %>%
  mutate(status = cut(Urinalysis_USG, c(-Inf, 1.010, 1.020, Inf),
                      labels = c("Euhydrated", "Mild", "Dehydrated (>1.020)"))) %>%
  count(PIN, status) %>%
  pivot_wider(names_from = status, values_from = n, values_fill = 0)

hydration <- exercise %>%
  mutate(pct_weight_loss = 100 * (Weight_kg_pre - Weight_kg_post) / Weight_kg_pre,
         hr_delta        = Heart_Rate_Post - Heart_Rate_Pre) %>%
  filter(!is.na(Urinalysis_USG), !is.na(hr_delta))

summary(lmer(hr_delta ~ scale(Urinalysis_USG) + (1 | PIN), data = hydration))
summary(lmer(pct_weight_loss ~ scale(Urinalysis_USG) + (1 | PIN),
             data = filter(hydration, !is.na(pct_weight_loss))))


# =============================================================================
# 9. ANTHROPOMETRICS
#
# Measured once per participant (baseline), so this is n = 9 -- describe, do
# not test. Left/right symmetry is the one thing worth checking.
# =============================================================================

anthro <- exercise %>%
  filter(!is.na(Anthro_Waist_cm) | !is.na(Anthro_Waist_in)) %>%
  transmute(
    PIN, Visit_Num,
    waist_cm = coalesce(Anthro_Waist_cm, Anthro_Waist_in * 2.54),
    hip_cm   = coalesce(Anthro_Hip_cm,   Anthro_Hip_in   * 2.54),
    whr      = waist_cm / hip_cm,
    rcalf    = coalesce(Anthro_Rcalf_cm, Anthro_Rcalf_in * 2.54),
    lcalf    = coalesce(Anthro_Lcalf_cm, Anthro_Lcalf_in * 2.54),
    rarm     = coalesce(Anthro_Rarm_cm,  Anthro_Rarm_in  * 2.54),
    larm     = coalesce(Anthro_Larm_cm,  Anthro_Larm_in  * 2.54),
    calf_asym = rcalf - lcalf,
    arm_asym  = rarm  - larm
  )

anthro
t.test(anthro$rcalf, anthro$lcalf, paired = TRUE)
t.test(anthro$rarm,  anthro$larm,  paired = TRUE)

# Baseline body composition vs mean HR response -- n = 9, so this is a
# hypothesis generator, nothing more.
anthro %>%
  left_join(hr %>% group_by(PIN) %>% summarise(mean_hr_delta = mean(delta)),
            by = "PIN") %>%
  ggplot(aes(whr, mean_hr_delta)) +
  geom_point(size = 3) + geom_text(aes(label = PIN), nudge_y = 1) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Waist-hip ratio vs mean HR response (n = 9)",
       x = "WHR", y = "Mean post - pre HR (bpm)")


# =============================================================================
# 10. FIGURES
# =============================================================================

# 10a. Paired slopes, faceted -- the single most informative plot in the set.
paired %>%
  pivot_longer(c(pre, post), names_to = "timepoint", values_to = "value") %>%
  mutate(timepoint = factor(timepoint, c("pre", "post"))) %>%
  ggplot(aes(timepoint, value, group = interaction(PIN, Visit_Num))) +
  geom_line(alpha = .15) +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               colour = "firebrick", linewidth = 1.2) +
  facet_wrap(~ measure, scales = "free_y") +
  labs(title = "Pre vs post by measure (grey = visits, red = mean)",
       x = NULL, y = NULL)

# 10b. Change scores by participant.
ggplot(paired, aes(PIN, delta)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.alpha = .4) +
  facet_wrap(~ measure, scales = "free_y") +
  labs(title = "Post - pre change by participant", x = NULL, y = "Delta")

# 10c. Resting-vitals timeline per participant.
exercise %>%
  select(PIN, Visit_Date, Heart_Rate_Pre, Systolic_Pre, Diastolic_Pre,
         Weight_kg_pre) %>%
  pivot_longer(-c(PIN, Visit_Date)) %>%
  ggplot(aes(Visit_Date, value, colour = PIN)) +
  geom_line(alpha = .8) +
  facet_wrap(~ name, scales = "free_y") +
  labs(title = "Resting vitals over the study", x = NULL, y = NULL)

# 10d. Correlation heatmap of the resting measures.
exercise %>%
  select(Heart_Rate_Pre, Systolic_Pre, Diastolic_Pre, Weight_kg_pre,
         Temperature_Pre, Oxygen_Saturation_Pre, Respirations_Pre,
         Urinalysis_USG, QTc_Interval) %>%
  cor(use = "pairwise.complete.obs") %>%
  as_tibble(rownames = "a") %>%
  pivot_longer(-a, names_to = "b", values_to = "r") %>%
  ggplot(aes(a, b, fill = r)) +
  geom_tile() +
  geom_text(aes(label = round(r, 2)), size = 3) +
  scale_fill_gradient2(limits = c(-1, 1)) +
  labs(title = "Resting measures, pairwise correlations", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
