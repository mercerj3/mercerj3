###############################################################################
# Predicted probabilities of ever having had a mammogram across intersectional
# strata - NULL MODEL, averaged across survey years.
#
# Variant of plot_null_strata.R: instead of holding year at the reference
# level, this marginalises over the observed year distribution. The averaging
# happens INSIDE each simulation draw, so the CI is correct - averaging the
# lwr/upr bounds after the fact would not be.
###############################################################################

library(dplyr); library(lme4); library(merTools); library(ggplot2)
## ---- year-AVERAGED predicted probability per stratum ------------------------
# all stratum x year combinations, with the observed number of respondents
grid <- mammograms_sp1 %>% ungroup() %>%
  count(stratum, stratum_label, year, name = "n_cell") %>%
  group_by(stratum) %>% mutate(w = n_cell / sum(n_cell)) %>% ungroup()

set.seed(456); nsim <- 2000
beta_hat <- fixef(model_null); V <- as.matrix(vcov(model_null))
X <- model.matrix(~ year, data = grid)
u_hat <- ranef(model_null, condVar = TRUE)$stratum
u_sd  <- sqrt(attr(u_hat, "postVar")[1, 1, ])
idx   <- match(as.character(grid$stratum), rownames(u_hat))

beta_draws <- MASS::mvrnorm(nsim, beta_hat, V)
u_draws    <- sapply(idx, function(j) rnorm(nsim, u_hat[j, 1], u_sd[j]))
probs      <- plogis(beta_draws %*% t(X) + u_draws)     # nsim x ncell

# average across years WITHIN each simulation draw, weighted by observed n
strata_ids <- unique(grid$stratum)
avg <- sapply(strata_ids, function(sj) {
  k <- which(grid$stratum == sj)
  as.vector(probs[, k, drop = FALSE] %*% grid$w[k])     # weighted mean per draw
})

strata_pred <- mammograms_sp1 %>% ungroup() %>%
  distinct(stratum, stratum_label, strataN) %>%
  right_join(data.frame(stratum = strata_ids,
                        pred = apply(avg, 2, median),
                        lwr  = apply(avg, 2, quantile, .025),
                        upr  = apply(avg, 2, quantile, .975)), by = "stratum") %>%
  mutate(stratum_label = reorder(stratum_label, pred))

grand_mean <- weighted.mean(
  predict(model_null, newdata = data.frame(year = levels(mammograms_sp1$year)),
          re.form = NA, type = "response"),
  w = as.vector(table(mammograms_sp1$year)))

cat("year-averaged range:", round(min(strata_pred$pred),3), "-",
    round(max(strata_pred$pred),3), "| grand mean:", round(grand_mean,3), "\n")

## ---- plot, all text black ---------------------------------------------------
p <- ggplot(strata_pred, aes(x = pred, y = stratum_label)) +
  geom_vline(xintercept = grand_mean, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, colour = "grey55", linewidth = .5) +
  geom_point(aes(size = strataN), colour = "#2c7fb8") +
  scale_size_continuous(range = c(1.5, 5), name = "Stratum n", labels = scales::comma) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Predicted probability of ever having had a mammogram",
       subtitle = "Intersectional strata, null model, averaged across survey years\nDashed line = overall predicted probability",
       x = "Predicted probability (95% CI)", y = NULL,
       caption = "Stratum = Race | Sexual orientation | Education") +
  theme_minimal(base_size = 11) +
  theme(
    text          = element_text(colour = "black"),   # titles, legend, caption
    axis.text     = element_text(colour = "black"),   # tick labels (grey30 by default)
    axis.title    = element_text(colour = "black"),
    plot.title    = element_text(colour = "black", face = "bold"),
    plot.subtitle = element_text(colour = "black"),
    plot.caption  = element_text(colour = "black"),
    legend.text   = element_text(colour = "black"),
    legend.title  = element_text(colour = "black"),
    strip.text    = element_text(colour = "black"),
    panel.grid.major.y = element_line(colour = "grey92"),
    panel.grid.minor   = element_blank()
  )
p
ggsave("null_strata_yearavg.png", p, width = 9, height = 8, dpi = 300)
