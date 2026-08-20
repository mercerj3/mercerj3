###############################################################################
# Predicted probabilities of ever having had a mammogram
# across intersectional strata - NULL MODEL
#
# Runs after model_null in mammograms_key2.R:
#   model_null <- glmer(hadmam ~ (1|stratum) + year, data = mammograms_sp1,
#                       family = binomial)
###############################################################################

library(dplyr)
library(merTools)
library(ggplot2)

# NOTE: merTools loads arm/MASS, which mask dplyr::select. Namespace it below.

## ---- 1. One row per stratum -------------------------------------------------
# The null model contains `year` as a fixed effect, so a predicted probability
# is only defined once year is fixed. Held at the reference year; the strata
# ORDER and SPACING are unaffected by this choice, only the overall level.
ref_year <- levels(mammograms_sp1$year)[1]   # "2018"

strata_nd <- mammograms_sp1 %>%
  ungroup() %>%
  distinct(stratum, stratum_label, race, sxorient, education, strataN) %>%
  mutate(year = factor(ref_year, levels = levels(mammograms_sp1$year))) %>%
  arrange(stratum)

stopifnot(nrow(strata_nd) == nlevels(mammograms_sp1$stratum))   # expect 30

## ---- 2. Predicted probability + 95% CI, including the random effect ---------
# include.resid.var = FALSE matches the rest of your script: uncertainty in the
# stratum mean, not the spread of individuals within it.
set.seed(123)
pi_null <- predictInterval(model_null,
                           newdata = as.data.frame(strata_nd),
                           level = 0.95,
                           type = "probability",
                           include.resid.var = FALSE,
                           n.sims = 2000)

strata_pred <- bind_cols(strata_nd, pi_null) %>%
  dplyr::rename(pred = fit) %>%
  mutate(stratum_label = reorder(stratum_label, pred))   # caterpillar ordering

## ---- 3. Grand mean (fixed part only) for the reference line ----------------
grand_mean <- predict(model_null, newdata = strata_nd[1, ],
                      re.form = NA, type = "response")

## ---- 4. Plot ---------------------------------------------------------------
p_null <- ggplot(strata_pred, aes(x = pred, y = stratum_label)) +
  geom_vline(xintercept = grand_mean, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0,
                 colour = "grey55", linewidth = 0.5) +
  geom_point(aes(size = strataN), colour = "#2c7fb8") +
  scale_size_continuous(range = c(1.5, 5), name = "Stratum n",
                        labels = scales::comma) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Predicted probability of ever having had a mammogram",
    subtitle = paste0("Intersectional strata, null model (random intercept only), ",
                      ref_year, "\nDashed line = overall predicted probability"),
    x = "Predicted probability (95% CI)",
    y = NULL,
    caption = "Stratum = Race | Sexual orientation | Education"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_line(colour = "grey92"),
        panel.grid.minor   = element_blank(),
        plot.title         = element_text(face = "bold"))

p_null
ggsave("null_strata_predicted.png", p_null, width = 9, height = 8, dpi = 300)

## ---- 5. The numbers behind the plot ----------------------------------------
strata_table <- strata_pred %>%
  dplyr::select(stratum, stratum_label, strataN, pred, lwr, upr) %>%
  arrange(desc(pred))
print(strata_table, n = 30)
# write.csv(strata_table, "null_strata_predicted.csv", row.names = FALSE)
