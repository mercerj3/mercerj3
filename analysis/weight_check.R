# =============================================================================
# WEIGHT -- standalone check
#
# Runs on `exercise` alone; does not need the main analysis script.
#
# Every weight is used exactly as recorded. Nothing is corrected, excluded, or
# back-converted. Weight is captured twice at each timepoint (lbs and kg), so
# each quantity below is computed down each unit column separately and the two
# results are averaged -- neither column is treated as the authoritative one.
# =============================================================================

library(tidyverse)
theme_set(theme_minimal(base_size = 12))

K <- 2.20462  # lbs per kg


# ---- Derive both routes, then average ---------------------------------------

wt <- exercise %>%
  filter(!is.na(PIN), PIN != "Participant ID (001-009)") %>%
  mutate(
    PIN = factor(str_trim(as.character(PIN))),
    across(c(Weight_lbs_pre, Weight_kg_pre, Weight_lbs_post, Weight_kg_post),
           ~ suppressWarnings(as.numeric(as.character(.x)))),

    # resting (pre) weight, each unit -> kg
    rest_via_lbs = Weight_lbs_pre / K,
    rest_via_kg  = Weight_kg_pre,
    rest_kg      = rowMeans(cbind(rest_via_lbs, rest_via_kg), na.rm = TRUE),

    # change across the bout, each unit -> kg
    bout_via_lbs = (Weight_lbs_post - Weight_lbs_pre) / K,
    bout_via_kg  =  Weight_kg_post  - Weight_kg_pre,
    bout_kg      = rowMeans(cbind(bout_via_lbs, bout_via_kg), na.rm = TRUE),

    across(c(rest_kg, bout_kg), ~ if_else(is.nan(.x), NA_real_, .x)),
    bout_pct = 100 * bout_kg / rest_kg   # negative = weight lost
  )


# ---- 1. Bout weight change, per participant ---------------------------------
# What one session does. Negative = lost weight (fluid).
# %% of body mass is the number to look at: >2% is where hydration guidance
# starts to flag a session.

bout_by_pin <- wt %>%
  filter(!is.na(bout_kg)) %>%
  group_by(PIN) %>%
  summarise(
    n_visits  = n(),
    mean_kg   = mean(bout_kg),
    sd_kg     = sd(bout_kg),
    min_kg    = min(bout_kg),
    max_kg    = max(bout_kg),
    mean_pct  = mean(bout_pct),
    worst_pct = min(bout_pct),
    visits_over_2pct = sum(bout_pct <= -2)
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

bout_by_pin

# same thing pooled across everyone
wt %>%
  filter(!is.na(bout_kg)) %>%
  summarise(n_participants = n_distinct(PIN), n_visits = n(),
            mean_kg = mean(bout_kg), sd_kg = sd(bout_kg),
            mean_pct = mean(bout_pct),
            visits_over_2pct = sum(bout_pct <= -2)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))


# ---- 2. Resting weight across the study, per participant --------------------
# Separate question: is body weight drifting over the program?

study_by_pin <- wt %>%
  filter(!is.na(rest_kg)) %>%
  arrange(PIN, Visit_Date) %>%
  group_by(PIN) %>%
  summarise(
    n_visits  = n(),
    first_kg  = first(rest_kg),
    last_kg   = last(rest_kg),
    change_kg = last(rest_kg) - first(rest_kg),
    change_pct = 100 * (last(rest_kg) - first(rest_kg)) / first(rest_kg),
    mean_kg   = mean(rest_kg),
    sd_kg     = sd(rest_kg)
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

study_by_pin


# ---- 3. Do the two unit columns agree? --------------------------------------
# Dashed line = perfect agreement. Points on the line mean the lbs column and
# the kg column tell the same story about the bout, so the averaging in the
# tables above is a formality rather than a judgement call.

# per participant
ggplot(filter(wt, !is.na(bout_via_kg), !is.na(bout_via_lbs)),
       aes(bout_via_kg, bout_via_lbs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = .7, size = 1.8) +
  facet_wrap(~ PIN) +
  coord_equal() +
  labs(title = "Bout weight change: kg column vs lbs column, by participant",
       subtitle = "Dashed line = perfect agreement. Each point is one visit.",
       x = "Change computed from kg (kg)",
       y = "Change computed from lbs (kg)")

# all participants together
ggplot(filter(wt, !is.na(bout_via_kg), !is.na(bout_via_lbs)),
       aes(bout_via_kg, bout_via_lbs, colour = PIN)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = .8, size = 2) +
  coord_equal() +
  labs(title = "Bout weight change: kg column vs lbs column, all participants",
       subtitle = "Dashed line = perfect agreement. Each point is one visit.",
       x = "Change computed from kg (kg)",
       y = "Change computed from lbs (kg)")

# the disagreement as a single number, per participant
wt %>%
  filter(!is.na(bout_via_kg), !is.na(bout_via_lbs)) %>%
  group_by(PIN) %>%
  summarise(n = n(),
            mean_gap_kg = mean(bout_via_lbs - bout_via_kg),
            max_gap_kg  = max(abs(bout_via_lbs - bout_via_kg)),
            r = cor(bout_via_lbs, bout_via_kg)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))


# ---- 4. The two plots the boss actually wants --------------------------------

# Bout change per participant, every visit shown.
ggplot(filter(wt, !is.na(bout_pct)), aes(PIN, bout_pct)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_hline(yintercept = -2, linetype = "dashed", colour = "firebrick") +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_jitter(width = .15, alpha = .5, size = 1.5) +
  labs(title = "Weight lost per session, by participant",
       subtitle = "Red dashed line = 2% of body mass",
       x = NULL, y = "Change across bout (% body mass)")

# Resting weight over the study.
ggplot(filter(wt, !is.na(rest_kg)), aes(Visit_Date, rest_kg, colour = PIN)) +
  geom_line() +
  geom_point(size = 1) +
  labs(title = "Resting weight across the study",
       x = NULL, y = "Weight (kg)")
