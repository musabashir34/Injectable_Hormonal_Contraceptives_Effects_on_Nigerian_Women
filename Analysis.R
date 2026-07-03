## =============================================================================
## Effects of DMPA-SC and NET-EN on Anthropometric Adiposity Indices
## Full analysis pipeline: descriptive stats, visualization, robust within-group
## tests (Wilcoxon + paired t-test, FDR-corrected, with effect sizes), and the
## between-formulation comparison (ANCOVA, adjusted for baseline value and age)
## =============================================================================

## ---- 0. Packages -----------------------------------------------------------
required_pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "broom", "effectsize",
                   "rstatix", "knitr", "patchwork")
to_install <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(effectsize)   # for Cohen's dz / rank_biserial if preferred over manual calc
library(rstatix)      # tidy wilcox_test / t_test with effect sizes
library(knitr)
library(patchwork)    # combine ggplots side by side

theme_set(theme_minimal(base_size = 12))

## ---- 1. Read raw data and build long/wide structures -----------------------
## Update the path below to point at your raw data file
raw <- read_excel("raw_data.xlsx", sheet = "Sheet1")

## The sheet stacks two groups (ID resets to 1 between them): first block =
## DMPA-SC, second block = NET-EN. Split on the ID reset.
id_reset <- which(diff(raw$ID) < 0)[1]          # last row of group 1
df1 <- raw[1:id_reset, ]
df2 <- raw[(id_reset + 1):nrow(raw), ]

df1$group <- "DMPA-SC"
df2$group <- "NET-EN"
dat <- bind_rows(df1, df2) %>% filter(!is.na(ID))

## ---- 2. Derive anthropometric indices --------------------------------------
dat <- dat %>%
  mutate(
    height_m = Heightcm / 100,
    WHR_Bf   = `WCcm/Bf` / `HCcm/Bf`,
    WHR_Af   = `WCcm/Af` / `HCcm/Af`,
    WHtR_Bf  = `WCcm/Bf` / Heightcm,
    WHtR_Af  = `WCcm/Af` / Heightcm,
    BAI_Bf   = `HCcm/Bf` / (height_m^1.5) - 18,
    BAI_Af   = `HCcm/Af` / (height_m^1.5) - 18,
    group    = factor(group, levels = c("DMPA-SC", "NET-EN"))
  )

## Variable map: short label -> (before column, after column)
var_map <- list(
  NC   = c("NCcm/Bf",  "NCcm/Af"),
  WC   = c("WCcm/Bf",  "WCcm/Af"),
  HC   = c("HCcm/Bf",  "HCcm/Af"),
  WHR  = c("WHR_Bf",   "WHR_Af"),
  WHtR = c("WHtR_Bf",  "WHtR_Af"),
  BAI  = c("BAI_Bf",   "BAI_Af")
)

## Long format: one row per subject x variable x timepoint (used for plotting
## and for the paired tests / ANCOVA below)
long_dat <- bind_rows(lapply(names(var_map), function(v) {
  cols <- var_map[[v]]
  tibble(
    ID      = dat$ID,
    group   = dat$group,
    age     = dat$`age(yr)`,
    variable = v,
    before  = dat[[cols[1]]],
    after   = dat[[cols[2]]]
  )
}))

## =============================================================================
## PART A: DESCRIPTIVE STATISTICS
## =============================================================================

## ---- A1. Summary table (mean, SD, n) per variable, group, and timepoint ----
desc_tbl <- long_dat %>%
  pivot_longer(cols = c(before, after), names_to = "timepoint", values_to = "value") %>%
  group_by(group, variable, timepoint) %>%
  summarise(
    n    = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(timepoint = factor(timepoint, levels = c("before", "after"),
                            labels = c("Before", "After")))

print(kable(desc_tbl, digits = 2, caption = "Descriptive statistics by group, variable and timepoint"))
write.csv(desc_tbl, "descriptive_statistics.csv", row.names = FALSE)

## ---- A2. Visualization 1: Before/after paired boxplots, faceted by variable
plot_dat <- long_dat %>%
  pivot_longer(cols = c(before, after), names_to = "timepoint", values_to = "value") %>%
  mutate(timepoint = factor(timepoint, levels = c("before", "after"),
                            labels = c("Before", "After")))

p_box <- ggplot(plot_dat, aes(x = timepoint, y = value, fill = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = 21, position = position_dodge(0.75), width = 0.6) +
  facet_wrap(~ variable, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("DMPA-SC" = "#2C7FB8", "NET-EN" = "#D95F02")) +
  labs(
    x = NULL, y = "Value", fill = "Formulation") +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))
print(p_box)
ggsave("Figure-1.png", p_box, width = 10, height = 7, dpi = 300)


## ---- A4. Visualization 2: Mean change (95% CI) forest-style plot ----------
change_summary <- long_dat %>%
  mutate(delta = after - before) %>%
  group_by(group, variable) %>%
  summarise(
    n        = n(),
    mean_d   = mean(delta, na.rm = TRUE),
    se_d     = sd(delta, na.rm = TRUE) / sqrt(n),
    ci_lo    = mean_d - 1.96 * se_d,
    ci_hi    = mean_d + 1.96 * se_d,
    .groups  = "drop"
  )

p_forest <- ggplot(change_summary, aes(x = mean_d, y = variable, color = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi),
                  position = position_dodge(width = 0.5), size = 0.7) +
  scale_color_manual(values = c("DMPA-SC" = "#2C7FB8", "NET-EN" = "#D95F02")) +
  labs(
       x = "Mean change (After \u2212 Before)", y = NULL, color = "Formulation") +
  theme(legend.position = "top")
print(p_forest)
ggsave("Figure-2.png", p_forest, width = 8, height = 5, dpi = 300)

## =============================================================================
## PART B: WITHIN-GROUP INFERENTIAL TESTS
## Wilcoxon signed-rank (primary) + paired t-test (secondary), with effect
## sizes (matched-pairs rank-biserial r and Cohen's dz), FDR-corrected across
## the 6 outcomes within each group.
## =============================================================================

run_within_group_tests <- function(data, grp) {
  sub <- data %>% filter(group == grp)
  
  results <- lapply(names(var_map), function(v) {
    d <- sub %>% filter(variable == v)
    before <- d$before
    after  <- d$after
    delta  <- after - before
    
    ## Normality check on the difference scores (informs which test to trust)
    sw <- shapiro.test(delta)
    
    ## Paired t-test
    tt <- t.test(after, before, paired = TRUE)
    dz <- mean(delta) / sd(delta)                 # Cohen's dz for paired data
    
    ## Wilcoxon signed-rank test (handles ties via normal approximation)
    wt <- wilcox.test(after, before, paired = TRUE, conf.int = TRUE, exact = FALSE)
    
    ## Matched-pairs rank-biserial correlation (effect size for Wilcoxon)
    nz       <- delta[delta != 0]
    ranks    <- rank(abs(nz))
    Wpos     <- sum(ranks[nz > 0])
    Wneg     <- sum(ranks[nz < 0])
    r_rb     <- (Wpos - Wneg) / (Wpos + Wneg)
    
    tibble(
      group         = grp,
      variable      = v,
      n             = length(before),
      before_mean   = mean(before), before_sd = sd(before),
      after_mean    = mean(after),  after_sd  = sd(after),
      mean_diff     = mean(delta),
      ci_lo         = tt$conf.int[1], ci_hi = tt$conf.int[2],
      shapiro_p     = sw$p.value,
      t_stat        = unname(tt$statistic), p_ttest = tt$p.value, cohens_dz = dz,
      W_stat        = unname(wt$statistic), p_wilcoxon = wt$p.value, rank_biserial_r = r_rb
    )
  })
  
  out <- bind_rows(results)
  ## Benjamini-Hochberg FDR correction across the 6 outcomes, within this group
  out$p_ttest_FDR    <- p.adjust(out$p_ttest,    method = "BH")
  out$p_wilcoxon_FDR <- p.adjust(out$p_wilcoxon, method = "BH")
  out
}

table1_dmpa  <- run_within_group_tests(long_dat, "DMPA-SC")
table2_neten <- run_within_group_tests(long_dat, "NET-EN")

print(kable(table1_dmpa  %>% select(variable, n, before_mean, after_mean, mean_diff,
                                    p_wilcoxon_FDR, rank_biserial_r, cohens_dz),
            digits = 4, caption = "Table 1: DMPA-SC \u2014 Wilcoxon signed-rank, FDR-adjusted"))
print(kable(table2_neten %>% select(variable, n, before_mean, after_mean, mean_diff,
                                    p_wilcoxon_FDR, rank_biserial_r, cohens_dz),
            digits = 4, caption = "Table 2: NET-EN \u2014 Wilcoxon signed-rank, FDR-adjusted"))

write.csv(table1_dmpa,  "table-1.csv",  row.names = FALSE)
write.csv(table2_neten, "table-2.csv", row.names = FALSE)

## =============================================================================
## PART C: BETWEEN-FORMULATION COMPARISON
## ANCOVA: post-value ~ baseline + age + group, per outcome.
## This is the test of whether DMPA-SC and NET-EN differ.
## =============================================================================

run_ancova <- function(v) {
  cols <- var_map[[v]]
  d <- dat %>%
    transmute(group, age = `age(yr)`,
              baseline = .data[[cols[1]]],
              post     = .data[[cols[2]]])
  
  fit <- lm(post ~ baseline + age + group, data = d)
  tt  <- broom::tidy(fit, conf.int = TRUE)
  grp_row <- tt %>% filter(grepl("^group", term))
  
  ## Standardised beta for the group effect (effect size)
  std_beta <- coef(fit)["groupNET-EN"] * sd(as.numeric(d$group) - 1) / sd(d$post)
  
  tibble(
    variable        = v,
    group_estimate  = grp_row$estimate,
    group_se        = grp_row$std.error,
    t_value         = grp_row$statistic,
    p_value         = grp_row$p.value,
    std_beta        = std_beta,
    r_squared       = summary(fit)$r.squared
  )
}

ancova_results <- bind_rows(lapply(names(var_map), run_ancova))
ancova_results$p_value_FDR <- p.adjust(ancova_results$p_value, method = "BH")

print(kable(ancova_results, digits = 4,
            caption = "Table 3: ANCOVA \u2014 DMPA-SC vs NET-EN, adjusted for baseline value and age"))
write.csv(ancova_results, "table-3.csv", row.names = FALSE)

## =============================================================================
## PART D: VISUALIZATION OF THE COMPARATIVE (ANCOVA) RESULT
## =============================================================================

p_ancova <- ggplot(ancova_results, aes(x = group_estimate, y = variable)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = group_estimate - 1.96 * group_se,
                      xmax = group_estimate + 1.96 * group_se),
                  color = "#7570B3", size = 0.7) +
  labs(
       x = "Adjusted mean difference (NET-EN \u2212 DMPA-SC)", y = NULL)
print(p_ancova)
ggsave("Figure-3.png", p_ancova, width = 8, height = 5, dpi = 300)

## =============================================================================
## END OF SCRIPT
## Outputs written to working directory:
##   descriptive_statistics.csv, table1_dmpa_results.csv, table2_neten_results.csv,
##   table3_ancova_results.csv, mixed_model_sensitivity.csv
##   fig1_boxplots_before_after.png, fig2_spaghetti_trajectories.png,
##   fig3_mean_change_forest.png, fig4_ancova_group_effect.png
## =============================================================================
