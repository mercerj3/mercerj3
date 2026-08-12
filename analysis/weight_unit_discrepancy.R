# =============================================================================
# kg vs lbs DISCREPANCY -- standalone
#
# Runs on `exercise` alone; does not need the main analysis script.
#
# Weight is recorded twice at every timepoint, once in lbs and once in kg.
# The two do not always imply the same weight. This file characterises that
# disagreement. It changes nothing: no value is corrected, excluded, or
# back-converted anywhere below.
#
# The quantity throughout is
#     gap = recorded lbs  -  (recorded kg x 2.20462)
# i.e. how many pounds the two columns disagree by. Positive means the lbs
# column reads heavier than the kg column implies.
# =============================================================================

library(tidyverse)
theme_set(theme_minimal(base_size = 12))

K <- 2.20462  # lbs per kg


# ---- Build one row per weighing (2 per visit: pre and post) -----------------

gaps <- exercise %>%
  filter(!is.na(PIN), PIN != "Participant ID (001-009)") %>%
  mutate(PIN = factor(str_trim(as.character(PIN))),
         across(c(Weight_lbs_pre, Weight_kg_pre, Weight_lbs_post, Weight_kg_post),
                ~ suppressWarnings(as.numeric(as.character(.x))))) %>%
  select(PIN, Visit_Num, Visit_Date,
         Weight_lbs_pre, Weight_kg_pre, Weight_lbs_post, Weight_kg_post) %>%
  pivot_longer(
    -c(PIN, Visit_Num, Visit_Date),
    names_to = c("unit", "timepoint"),
    names_pattern = "Weight_(lbs|kg)_(pre|post)",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = unit, values_from = value) %>%
  rename(lbs_recorded = lbs, kg_recorded = kg) %>%
  drop_na(lbs_recorded, kg_recorded) %>%
  mutate(
    timepoint    = factor(str_to_title(timepoint), c("Pre", "Post")),
    lbs_implied  = kg_recorded * K,          # what the kg column implies
    gap_lbs      = lbs_recorded - lbs_implied,
    gap_pct      = 100 * gap_lbs / lbs_implied,
    mean_lbs     = (lbs_recorded + lbs_implied) / 2
  )


# ---- 1. Per participant -----------------------------------------------------

gaps %>%
  group_by(PIN) %>%
  summarise(
    n            = n(),
    mean_gap     = mean(gap_lbs),          # bias: is one column consistently high?
    sd_gap       = sd(gap_lbs),            # scatter: how repeatable is the gap?
    median_abs   = median(abs(gap_lbs)),
    max_abs      = max(abs(gap_lbs)),
    within_0.5lb = sprintf("%.0f%%", 100 * mean(abs(gap_lbs) <= 0.5))
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))


# ---- 2. All participants pooled ---------------------------------------------

gaps %>%
  summarise(
    n_participants = n_distinct(PIN),
    n_weighings    = n(),
    mean_gap       = mean(gap_lbs),
    sd_gap         = sd(gap_lbs),
    min_gap        = min(gap_lbs),
    max_gap        = max(gap_lbs),
    within_0.1lb   = sprintf("%.1f%%", 100 * mean(abs(gap_lbs) <= 0.1)),
    within_0.5lb   = sprintf("%.1f%%", 100 * mean(abs(gap_lbs) <= 0.5)),
    over_1lb       = sum(abs(gap_lbs) > 1)
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# split by pre vs post -- does the disagreement grow after exercise, when
# people are weighed in a hurry?
gaps %>%
  group_by(timepoint) %>%
  summarise(n = n(), mean_gap = mean(gap_lbs), sd_gap = sd(gap_lbs)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))


# ---- 3. Agreement plots -----------------------------------------------------
# Dashed line = perfect agreement (the two columns describe the same weight).

# per participant
ggplot(gaps, aes(lbs_implied, lbs_recorded)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(colour = timepoint), alpha = .8, size = 1.8) +
  facet_wrap(~ PIN, scales = "free") +
  labs(title = "Recorded lbs vs lbs implied by the kg column, by participant",
       subtitle = "Dashed line = perfect agreement. One point per weighing.",
       x = "kg column x 2.20462 (lbs)", y = "lbs column as recorded",
       colour = NULL)

# all participants
ggplot(gaps, aes(lbs_implied, lbs_recorded, colour = PIN)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = .8, size = 2) +
  coord_equal() +
  labs(title = "Recorded lbs vs lbs implied by the kg column, all participants",
       subtitle = "Dashed line = perfect agreement. One point per weighing.",
       x = "kg column x 2.20462 (lbs)", y = "lbs column as recorded")

# NOTE: at this scale the points sit on the line and the disagreement is
# invisible. Sections 4-5 zoom in on the gap itself, which is where it shows.


# ---- 4. The gap on its own --------------------------------------------------

# distribution per participant -- the plot to hand over
ggplot(gaps, aes(PIN, gap_lbs)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_jitter(aes(colour = timepoint), width = .18, alpha = .6, size = 1.4) +
  labs(title = "kg/lbs disagreement by participant",
       subtitle = "Dashed line = the two columns agree exactly",
       x = NULL, y = "Recorded lbs - implied lbs", colour = NULL)

# Bland-Altman: does the gap widen with body weight? A funnel shape would mean
# proportional error (a scale calibration issue); a flat band means it does not
# depend on how heavy the person is.
ggplot(gaps, aes(mean_lbs, gap_lbs)) +
  geom_hline(yintercept = mean(gaps$gap_lbs)) +
  geom_hline(yintercept = mean(gaps$gap_lbs) + c(-1.96, 1.96) * sd(gaps$gap_lbs),
             linetype = "dashed", colour = "firebrick") +
  geom_point(alpha = .6) +
  labs(title = "Bland-Altman: kg vs lbs columns",
       subtitle = "Solid = mean gap, red dashed = 95% limits of agreement",
       x = "Mean of the two (lbs)", y = "Recorded lbs - implied lbs")


# ---- 5. Is it systematic or is it noise? ------------------------------------
# Three checks. If all three come back flat, the disagreement is random
# recording noise and cannot bias any result.

# (a) does the gap scale with body weight?
cor.test(gaps$gap_lbs, gaps$kg_recorded)

# (b) does it drift across the study?
cor.test(gaps$gap_lbs, gaps$Visit_Num)

ggplot(gaps, aes(Visit_Num, gap_lbs, colour = PIN)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(alpha = .7) +
  facet_wrap(~ PIN) +
  guides(colour = "none") +
  labs(title = "Does the disagreement drift over the study?",
       x = "Visit", y = "Recorded lbs - implied lbs")

# (c) does any participant carry a consistent offset? A confidence interval
#     that excludes zero means that participant's two columns disagree
#     systematically, not randomly.
gaps %>%
  group_by(PIN) %>%
  summarise(broom::tidy(t.test(gap_lbs))) %>%
  select(PIN, mean_gap = estimate, conf.low, conf.high, p.value) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))


# ---- 6. How much of the gap is just rounding? -------------------------------
# Both columns are recorded to 1 decimal place. Rounding to 0.1 kg alone puts
# up to +/-0.05 kg (= +/-0.11 lbs) of slack between the columns before anyone
# has made a mistake. This compares the observed scatter against the scatter
# rounding could produce by itself.

rounding_sd <- sqrt((0.1 * K)^2 / 12 + 0.1^2 / 12)  # kg rounding + lbs rounding

tibble(
  observed_sd    = sd(gaps$gap_lbs),
  rounding_sd    = rounding_sd,
  ratio          = sd(gaps$gap_lbs) / rounding_sd
) %>%
  mutate(across(everything(), ~ round(.x, 3)))

# ratio ~1 => the columns disagree only as much as rounding forces them to.
# ratio >>1 => the two weights were read or transcribed independently, and the
# extra scatter is recording noise rather than arithmetic.
