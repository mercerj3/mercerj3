###############################################################################
# Predicted probabilities of ever having had a mammogram
# across intersectional strata - NULL MODEL
#
# Run after model_null in mammograms_key2.R:
#   model_null <- glmer(hadmam ~ (1|stratum) + year,
#                       data = mammograms_sp1, family = binomial)
#
# ORDER OF OPERATIONS (the part that makes the CI correct)
# -------------------------------------------------------
# For each of N_SIMS draws:
#   1. draw beta*  ~ MVN(beta_hat, vcov(model))   intercept + all year dummies
#   2. draw u*_j   ~ N(u_hat_j, condSD_j)         one per stratum
#   3. eta_jy = beta*_0 + beta*_y + u*_j          every stratum x year cell
#   4. p_jy   = plogis(eta_jy)                    probability scale
#   5. p_bar_j = sum_y ( w_y * p_jy )             AVERAGE OVER YEARS, INSIDE THE DRAW
# Then, across the N_SIMS values of p_bar_j:
#   6. point estimate = median
#   7. 95% CI = 2.5th and 97.5th percentiles
#
# Averaging happens at step 5, BEFORE any percentile is taken. Averaging the
# lwr/upr bounds afterwards instead gives intervals roughly twice too wide,
# because it keeps each year's coefficient uncertainty rather than averaging
# it away. Years are not processed in chronological order - order is irrelevant
# to a weighted mean; w_y is simply each year's share of respondents.
###############################################################################

library(dplyr)
library(lme4)
library(ggplot2)
library(MASS)          # mvrnorm
# NOTE: MASS masks dplyr::select, so select() is namespaced throughout.

## ---------------------------------------------------------------------------
## OPTIONS
## ---------------------------------------------------------------------------
N_SIMS    <- 2000
SEED      <- 456

# "average"   - marginalise over survey years (recommended for the study period)
# "reference" - hold year at its first level; simpler, condition on one year
YEAR_MODE <- "average"

# Only used when YEAR_MODE == "average":
# "common"  - every stratum weighted by the OVERALL year distribution, so strata
#             are compared at an identical year mix (recommended: differences
#             between strata are not contaminated by differences in when each
#             stratum happened to be sampled)
# "stratum" - each stratum weighted by its OWN year distribution, i.e. that
#             stratum's actual observed experience
YEAR_WEIGHTS <- "common"

REF_YEAR  <- levels(mammograms_sp1$year)[1]
OUT_PNG   <- "null_strata_predicted.png"
OUT_CSV   <- "null_strata_predicted.csv"

## ---------------------------------------------------------------------------
## 1. Stratum x year grid, with the averaging weights
## ---------------------------------------------------------------------------
strata_info <- mammograms_sp1 %>%
  ungroup() %>%
  distinct(stratum, stratum_label, race, sxorient, education, strataN) %>%
  arrange(stratum)

stopifnot(nrow(strata_info) == nlevels(droplevels(mammograms_sp1$stratum)))

if (YEAR_MODE == "reference") {
  grid <- strata_info %>%
    mutate(year = factor(REF_YEAR, levels = levels(mammograms_sp1$year)), w = 1)
} else {
  year_share <- mammograms_sp1 %>% ungroup() %>%
    count(year, name = "n_year") %>% mutate(w_common = n_year / sum(n_year))

  grid <- mammograms_sp1 %>% ungroup() %>%
    count(stratum, year, name = "n_cell") %>%
    group_by(stratum) %>% mutate(w_stratum = n_cell / sum(n_cell)) %>% ungroup() %>%
    left_join(year_share, by = "year") %>%
    mutate(w = if (YEAR_WEIGHTS == "common") w_common else w_stratum) %>%
    # renormalise within stratum so weights sum to exactly 1
    group_by(stratum) %>% mutate(w = w / sum(w)) %>% ungroup() %>%
    left_join(strata_info, by = "stratum")
}

stopifnot(all(abs(tapply(grid$w, grid$stratum, sum) - 1) < 1e-10))

## ---------------------------------------------------------------------------
## 2. Simulate  (steps 1-5 above)
## ---------------------------------------------------------------------------
set.seed(SEED)

beta_hat <- fixef(model_null)
V        <- as.matrix(vcov(model_null))
X        <- model.matrix(~ year, data = grid)          # design for every cell

u        <- ranef(model_null, condVar = TRUE)$stratum
u_sd     <- sqrt(attr(u, "postVar")[1, 1, ])
idx      <- match(as.character(grid$stratum), rownames(u))
stopifnot(!anyNA(idx))

beta_draws <- MASS::mvrnorm(N_SIMS, beta_hat, V)                             # 1
u_draws    <- sapply(idx, function(j) rnorm(N_SIMS, u[j, 1], u_sd[j]))       # 2
probs      <- plogis(beta_draws %*% t(X) + u_draws)                          # 3-4

strata_ids <- unique(as.character(grid$stratum))
p_bar <- sapply(strata_ids, function(sj) {                                   # 5
  k <- which(as.character(grid$stratum) == sj)
  as.vector(probs[, k, drop = FALSE] %*% grid$w[k])
})   # N_SIMS x n_strata: one averaged probability per draw per stratum

## ---------------------------------------------------------------------------
## 3. Summarise across draws  (steps 6-7)
## ---------------------------------------------------------------------------
strata_pred <- data.frame(
  stratum = strata_ids,
  pred    = apply(p_bar, 2, median),
  lwr     = apply(p_bar, 2, quantile, 0.025),
  upr     = apply(p_bar, 2, quantile, 0.975),
  row.names = NULL
) %>%
  left_join(strata_info %>% mutate(stratum = as.character(stratum)), by = "stratum") %>%
  mutate(
    # n printed in the label: exact, and avoids a size scale that cannot span
    # a 33 -> 146,742 range legibly
    lab_n = paste0(stratum_label, "  (n = ",
                   format(strataN, big.mark = ",", trim = TRUE), ")"),
    lab_n = reorder(lab_n, pred)
  ) %>%
  arrange(desc(pred))

## ---------------------------------------------------------------------------
## 4. Overall predicted probability (fixed effects only) for the reference line
## ---------------------------------------------------------------------------
if (YEAR_MODE == "reference") {
  grand_mean <- predict(model_null,
                        newdata = data.frame(year = factor(REF_YEAR,
                                             levels = levels(mammograms_sp1$year))),
                        re.form = NA, type = "response")
} else {
  yr_levels  <- levels(mammograms_sp1$year)
  yr_p       <- predict(model_null,
                        newdata = data.frame(year = factor(yr_levels, levels = yr_levels)),
                        re.form = NA, type = "response")
  yr_w       <- if (YEAR_WEIGHTS == "common") {
                  as.vector(table(mammograms_sp1$year))
                } else as.vector(table(mammograms_sp1$year))
  grand_mean <- weighted.mean(yr_p, w = yr_w)
}

## ---------------------------------------------------------------------------
## 5. Plot
## ---------------------------------------------------------------------------
subtitle_txt <- if (YEAR_MODE == "reference") {
  paste0("Intersectional strata, null model (random intercept only), ", REF_YEAR)
} else {
  paste0("Intersectional strata, null model (random intercept only)\n",
         "Averaged across survey years (",
         if (YEAR_WEIGHTS == "common") "common year mix" else "each stratum's own year mix", ")")
}

p_null <- ggplot(strata_pred, aes(x = pred, y = lab_n)) +
  geom_vline(xintercept = grand_mean, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0,
                 colour = "grey55", linewidth = 0.5) +
  geom_point(colour = "#2c7fb8", size = 2.3) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Predicted probability of ever having had a mammogram",
    subtitle = paste0(subtitle_txt, "\nDashed line = overall predicted probability"),
    x = "Predicted probability (95% CI)",
    y = NULL,
    caption = "Stratum = Race | Sexual orientation | Education"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    # theme_minimal leaves plot.background$fill NULL, so a saved PNG is
    # TRANSPARENT. Composited on a dark background the black text vanishes.
    # Set both backgrounds white explicitly (ggsave(bg=) below reinforces it).
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    text          = element_text(colour = "black"),
    axis.text     = element_text(colour = "black"),   # grey30 by default
    axis.title    = element_text(colour = "black"),
    plot.title    = element_text(colour = "black", face = "bold"),
    plot.subtitle = element_text(colour = "black"),
    plot.caption  = element_text(colour = "black"),
    panel.grid.major.y = element_line(colour = "grey92"),
    panel.grid.minor   = element_blank()
  )

print(p_null)
ggsave(OUT_PNG, p_null, width = 9.5, height = 7.5, dpi = 300, bg = "white")

## ---------------------------------------------------------------------------
## 6. The numbers behind the plot
## ---------------------------------------------------------------------------
strata_table <- strata_pred %>%
  dplyr::select(stratum, stratum_label, race, sxorient, education,
                strataN, pred, lwr, upr) %>%
  mutate(across(c(pred, lwr, upr), ~ round(.x, 4)))

print(dplyr::as_tibble(strata_table), n = 30)
write.csv(strata_table, OUT_CSV, row.names = FALSE)

cat("\nStrata:", nrow(strata_table),
    "| predicted range:", round(min(strata_table$pred), 3), "-",
    round(max(strata_table$pred), 3),
    "| overall:", round(grand_mean, 3), "\n")
