# ============================================================
# Visit 1 Demographics — summary table + figures for presentation
# ============================================================

pkgs <- c("dplyr", "tidyr", "forcats", "ggplot2", "scales", "patchwork")
new_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new_pkgs)) install.packages(new_pkgs)
# gt is optional -- used only to render a polished Table 1 image.
# install.packages("gt") if you don't already have it.
has_gt <- requireNamespace("gt", quietly = TRUE)

library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(scales)
library(patchwork)

# ------------------------------------------------------------
# 1. Data entry
# ------------------------------------------------------------

demo_raw <- tibble::tribble(
  ~PIN, ~Visit_Date,   ~Height_in, ~Height_cm, ~Veteran, ~Gender, ~Age_at_Inclusion, ~Race,   ~Ethnicity,
  "001", "2025-07-21",       71.5,      181.6,        0,       0,                39,   "1",            0,
  "002", "2025-07-25",       71.2,      180.0,        0,       0,                34,   "1",            0,
  "003", "2025-09-08",       72.5,      184.0,        0,       0,                30,   "1,3",          1,
  "004", "2025-09-09",       72.5,      184.0,        0,       0,                32,   "1",            0,
  "005", "2025-12-02",       67.5,      171.0,        0,       0,                33,   "1",            0,
  "006", "2025-12-01",       68.0,      172.7,        0,       0,                36,   "1",            0,
  "007", "2026-03-04",       71.0,      181.0,        1,       0,                46,   "1",            0,
  "008", "2026-03-02",       63.5,      161.3,        0,       1,                29,   "1",            1,
  "009", "2026-04-30",       63.5,      161.0,        0,       1,                26,   "1",            1
)

race_labels <- c(`1` = "White", `2` = "Black/AA", `3` = "Asian",
                  `4` = "NH/PI", `5` = "AI/AN", `6` = "Other")

demo <- demo_raw %>%
  mutate(
    Visit_Date = as.Date(Visit_Date),
    Gender     = factor(Gender, levels = c(0, 1), labels = c("Male", "Female")),
    Veteran    = factor(Veteran, levels = c(0, 1), labels = c("No", "Yes")),
    Ethnicity  = factor(Ethnicity, levels = c(0, 1),
                         labels = c("Not Hispanic/Latino", "Hispanic/Latino")),
    n_races    = lengths(strsplit(Race, ",")),
    Race_Group = if_else(n_races > 1, "Multiracial",
                          unname(race_labels[Race])),
    Race_Group = fct_infreq(Race_Group)
  )

# One row per race code selected, for a multi-response race chart
# (participant 003 selected both White and Asian).
demo_race_long <- demo %>%
  separate_rows(Race, sep = ",") %>%
  mutate(Race_Label = unname(race_labels[Race]))

# ------------------------------------------------------------
# 2. Table 1 — demographic summary
# ------------------------------------------------------------

n_total <- nrow(demo)

fmt_pct <- function(n, denom = n_total) sprintf("%d (%.0f%%)", n, 100 * n / denom)

cat_summary <- function(data, var, label) {
  data %>%
    count(.data[[var]]) %>%
    transmute(Characteristic = label,
              Level = as.character(.data[[var]]),
              Summary = fmt_pct(n))
}

table1 <- bind_rows(
  tibble(Characteristic = "N", Level = "", Summary = as.character(n_total)),
  tibble(Characteristic = "Age at inclusion, years",
         Level = "Mean (SD)",
         Summary = sprintf("%.1f (%.1f)", mean(demo$Age_at_Inclusion), sd(demo$Age_at_Inclusion))),
  tibble(Characteristic = "Age at inclusion, years",
         Level = "Median (range)",
         Summary = sprintf("%.0f (%d–%d)", median(demo$Age_at_Inclusion),
                            min(demo$Age_at_Inclusion), max(demo$Age_at_Inclusion))),
  tibble(Characteristic = "Height, cm",
         Level = "Mean (SD)",
         Summary = sprintf("%.1f (%.1f)", mean(demo$Height_cm), sd(demo$Height_cm))),
  cat_summary(demo, "Gender", "Gender"),
  cat_summary(demo, "Veteran", "Veteran status"),
  cat_summary(demo, "Ethnicity", "Ethnicity"),
  demo %>% count(Race_Group) %>%
    transmute(Characteristic = "Race", Level = as.character(Race_Group), Summary = fmt_pct(n))
)

print(table1, n = Inf)
write.csv(table1, "table1_demographics.csv", row.names = FALSE)

if (has_gt) {
  library(gt)
  table1_gt <- table1 %>%
    gt(groupname_col = "Characteristic") %>%
    tab_header(title = "Table 1. Visit 1 Demographic Characteristics",
               subtitle = paste0("N = ", n_total)) %>%
    cols_label(Level = "", Summary = "n (%) / value") %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_row_groups()) %>%
    tab_options(table.font.size = 14, heading.align = "left")
  gtsave(table1_gt, "table1_demographics.html")
  # PNG export needs webshot2 + Chrome (install.packages("webshot2")); falls
  # back silently to the HTML file above if that isn't set up.
  tryCatch(
    gtsave(table1_gt, "table1_demographics.png", vwidth = 900, vheight = 900),
    error = function(e) message("PNG export skipped (install 'webshot2' for image export): ", conditionMessage(e))
  )
}

# ------------------------------------------------------------
# 3. Palette + shared theme
# ------------------------------------------------------------

pal <- c("#2a78d6", "#1baf7a", "#eda100", "#008300", "#4a3aa7", "#e34948")

theme_presentation <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "#0b0b0b"),
    plot.subtitle = element_text(size = 12, color = "#52514e"),
    axis.title = element_text(color = "#52514e"),
    axis.text = element_text(color = "#52514e"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#e1e0d9"),
    legend.position = "none"
  )
theme_set(theme_presentation)

bar_chart <- function(data, var, title, subtitle = NULL, denom = n_total) {
  d <- data %>% count(.data[[var]]) %>%
    mutate(pct = n / denom, lab = sprintf("%d (%.0f%%)", n, 100 * pct))
  ggplot(d, aes(x = .data[[var]], y = n, fill = .data[[var]])) +
    geom_col(width = 0.6) +
    geom_text(aes(label = lab), vjust = -0.6, size = 4.2, color = "#0b0b0b") +
    scale_fill_manual(values = pal) +
    scale_y_continuous(limits = c(0, denom), expand = expansion(mult = c(0, 0.15))) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Participants (n)")
}

# ------------------------------------------------------------
# 4. Categorical charts
# ------------------------------------------------------------

p_gender    <- bar_chart(demo, "Gender", "Gender")
p_veteran   <- bar_chart(demo, "Veteran", "Veteran Status")
p_ethnicity <- bar_chart(demo, "Ethnicity", "Ethnicity")

p_race <- demo_race_long %>%
  count(Race_Label) %>%
  mutate(lab = sprintf("%d (%.0f%%)", n, 100 * n / n_total)) %>%
  ggplot(aes(x = fct_reorder(Race_Label, n, .desc = TRUE), y = n, fill = Race_Label)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = lab), vjust = -0.6, size = 4.2, color = "#0b0b0b") +
  scale_fill_manual(values = pal) +
  scale_y_continuous(limits = c(0, n_total), expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Race", subtitle = "One participant reported more than one race",
       x = NULL, y = "Participants (n)")

# ------------------------------------------------------------
# 5. Continuous variables — dot plots (n = 9, so show every point)
# ------------------------------------------------------------

dot_plot <- function(data, var, title, ylab, unit = "") {
  m <- mean(data[[var]])
  rng <- diff(range(data[[var]]))
  ggplot(data, aes(x = "", y = .data[[var]])) +
    geom_hline(yintercept = m, linetype = "dashed", color = "#898781") +
    geom_jitter(width = 0.05, size = 4, color = pal[1], alpha = 0.85) +
    annotate("label", x = 1, y = m + 0.08 * rng,
             label = sprintf("mean = %.1f%s", m, unit),
             size = 4, color = "#52514e", label.size = 0, fill = "#fcfcfb") +
    scale_x_discrete(expand = expansion(add = c(0.7, 0.7))) +
    labs(title = title, x = NULL, y = ylab)
}

p_age    <- dot_plot(demo, "Age_at_Inclusion", "Age at Inclusion", "Years", " yrs")
p_height <- dot_plot(demo, "Height_cm", "Height", "Centimeters", " cm")

# ------------------------------------------------------------
# 6. Combined overview figure + individual exports
# ------------------------------------------------------------

overview <- (p_age | p_height) / (p_gender | p_veteran) / (p_ethnicity | p_race) +
  plot_annotation(
    title = "Visit 1 Demographics Overview",
    subtitle = paste0("N = ", n_total, " participants"),
    theme = theme(plot.title = element_text(face = "bold", size = 18),
                  plot.subtitle = element_text(size = 13, color = "#52514e"))
  )

ggsave("demographics_overview.png", overview, width = 10, height = 12, dpi = 300, bg = "white")

for (nm in c("p_gender", "p_veteran", "p_ethnicity", "p_race", "p_age", "p_height")) {
  ggsave(paste0(sub("^p_", "", nm), ".png"), get(nm), width = 6, height = 5, dpi = 300, bg = "white")
}

if (interactive()) print(overview)
