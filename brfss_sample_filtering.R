library(haven)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)

# ---------------------------------------------------------------------------
# Confirm the age variable before filtering
# AGE80 = age in years, top-coded at 80 (BRFSS _AGE80); no 7/9 missing codes
# ---------------------------------------------------------------------------
summary(mydata$AGE80)
table(mydata$AGE80, useNA = "ifany")

# ---------------------------------------------------------------------------
# Filtering pipeline
# 1) 2018-2023, 2) female sex only, 3) ages 40-75, 4) mammogram data,
# 5) education data, 6) race data, 7) sexual orientation data
# ---------------------------------------------------------------------------

# Step 1: filter out 2017
mydata_f <- mydata %>%
  filter(SURVYEAR != 2017)

table(mydata_f$SURVYEAR)

# Step 2: women only
mydata_f <- mydata_f %>%
  filter(SEX == 2)

table(mydata_f$SEX)

# Step 3: ages 40 through 75 (screening-eligible age range)
mydata_f <- mydata_f %>%
  filter(!is.na(AGE80), AGE80 >= 40, AGE80 <= 75)

summary(mydata_f$AGE80)

# Step 4: has mammogram data (non-missing)
mydata_f <- mydata_f %>%
  filter(!is.na(HADMAM99))

table(mydata_f$HADMAM99)

# Step 5: has education data (non-missing)
mydata_f <- mydata_f %>%
  filter(!is.na(EDUCAG99))

table(mydata_f$EDUCAG99)

# Step 6: has race data (non-missing)
mydata_f <- mydata_f %>%
  filter(!is.na(RACE99))

table(mydata_f$RACE99)

# Step 7: has LGBO (derived from SXORIENT) non-missing
mydata_f <- mydata_f %>%
  filter(!is.na(LGBO))

table(mydata_f$LGBO)


# ---------------------------------------------------------------------------
# 1) Re-run the filtering pipeline, capturing sample size at each step
# ---------------------------------------------------------------------------
mydata_f <- mydata
steps <- tibble::tibble(step = "Full sample", n = nrow(mydata))

mydata_f <- mydata_f %>% dplyr::filter(SURVYEAR != 2017)
steps <- bind_rows(steps, tibble::tibble(step = "Exclude 2017", n = nrow(mydata_f)))

mydata_f <- mydata_f %>% dplyr::filter(SEX == 2)
steps <- bind_rows(steps, tibble::tibble(step = "Women only", n = nrow(mydata_f)))

mydata_f <- mydata_f %>% dplyr::filter(!is.na(AGE80), AGE80 >= 40, AGE80 <= 75)
steps <- bind_rows(steps, tibble::tibble(step = "Ages 40-75", n = nrow(mydata_f)))

mydata_f <- mydata_f %>% dplyr::filter(!is.na(HADMAM99))
steps <- bind_rows(steps, tibble::tibble(step = "Non-missing HADMAM99", n = nrow(mydata_f)))

mydata_f <- mydata_f %>% dplyr::filter(!is.na(EDUCAG99))
steps <- bind_rows(steps, tibble::tibble(step = "Non-missing EDUCAG99", n = nrow(mydata_f)))

mydata_f <- mydata_f %>% dplyr::filter(!is.na(RACE99))
steps <- bind_rows(steps, tibble::tibble(step = "Non-missing RACE99", n = nrow(mydata_f)))

mydata_f <- mydata_f %>% dplyr::filter(!is.na(LGBO))
steps <- bind_rows(steps, tibble::tibble(step = "Non-missing LGBO", n = nrow(mydata_f)))

# 2) Excluded N at each step
steps <- steps %>%
  mutate(
    step = factor(step, levels = step),
    n_excluded = replace_na(lag(n) - n, 0)
  )

# 3) Single-hue ordinal ramp (funnel stages = magnitude, not separate categories)
bar_fill <- colorRampPalette(c("#86b6ef", "#184f95"))(nrow(steps))

# 4) Attrition / funnel chart
ggplot(steps, aes(x = step, y = n)) +
  geom_col(fill = bar_fill, width = 0.62) +
  geom_text(aes(label = comma(n)),
            vjust = -0.6, size = 3.6, color = "#0b0b0b", fontface = "bold") +
  geom_text(data = filter(steps, n_excluded > 0),
            aes(y = n, label = paste0("−", comma(n_excluded), " excluded")),
            vjust = 2.0, size = 3, color = "#898781") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Sample Attrition Across Sequential Exclusion Criteria",
    subtitle = "BRFSS 2018–2023 — women ages 40–75, sequential non-missing filters",
    x = NULL, y = "Remaining sample size (N)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "#e1e0d9", linewidth = 0.3),
    axis.text.x  = element_text(angle = 30, hjust = 1, color = "#52514e"),
    axis.text.y  = element_text(color = "#898781"),
    axis.title.y = element_text(color = "#52514e"),
    plot.title    = element_text(face = "bold", color = "#0b0b0b"),
    plot.subtitle = element_text(color = "#52514e", size = 10)
  )


# ---------------------------------------------------------------------------
# Save the filtered analytic sample as an SPSS (.sav) file
# ---------------------------------------------------------------------------
write_sav(mydata_f, "~/Desktop/mydata_f.sav")

check <- read_sav("~/Desktop/mydata_f.sav")
dim(check)
