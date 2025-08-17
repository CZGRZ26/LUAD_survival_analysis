# lung-adenocarcinoma-RWE-GEO/
#   ├── data/
#   │   ├── raw/                   # Original GEO and clinical files
#   │   └── processed/             # Cleaned expression + clinical data
#   ├── scripts/
#   │   ├── 01_download_data.R
# │   ├── 02_clean_merge_data.R
# │   ├── 03_survival_analysis.R
# │   ├── 04_visualizations.R
# ├── results/
#   │   ├── figures/
#   │   └── tables/
#   ├── report/
#   │   └── RWE_LUAD_summary.Rmd   # Final analysis report (R Markdown)
# ├── README.md
# └── renv/                     # Optional: reproducible R environment
#   



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

# Load packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GEOquery")

# ==========================================
# LUAD Survival Analysis using GEO Data
# ==========================================
# --- 1. Load packages ---
library(GEOquery)
library(survival)
library(survminer)
library(dplyr)

# --- 2. Download GEO dataset ---
gse <- getGEO("GSE68465", GSEMatrix = TRUE)[[1]]
expr <- exprs(gse)    # expression matrix
clin <- pData(gse)    # clinical data

# --- 3. Prepare survival variables ---
clin$time <- as.numeric(clin$`months_to_last_contact_or_death:ch1`) * 30.44
vstatus <- tolower(trimws(clin$`vital_status:ch1`))
clin$status <- ifelse(vstatus %in% c("dead", "deceased"), 1,
                      ifelse(vstatus %in% c("alive", "living"), 0, NA))

clin_clean <- clin %>%
  filter(!is.na(status) & !is.na(time)) %>%
  # drop rare disease stage levels
  filter(`disease_stage:ch1` != "pNXpT1")

# --- 4. Get platform annotation ---
gpl <- getGEO(annotation(gse), AnnotGPL = TRUE)
gpl_table <- Table(gpl)

# --- 5. Find EGFR probes ---
egfr_probes <- gpl_table$ID[grep("^EGFR$", gpl_table$`Gene symbol`, ignore.case = TRUE)]
if (length(egfr_probes) == 0) stop("No probes found for EGFR!")

# --- 6. Extract EGFR expression (log2 if desired) ---
egfr_expr <- colMeans(expr[egfr_probes, , drop = FALSE], na.rm = TRUE)
egfr_expr <- egfr_expr[match(rownames(clin_clean), names(egfr_expr))]
clin_clean$EGFR <- log2(egfr_expr + 1)  # log2 transformation

# --- 7. Create high/low EGFR groups ---
clin_clean$EGFR_group <- factor(ifelse(clin_clean$EGFR >= median(clin_clean$EGFR, na.rm = TRUE),
                                       "High", "Low"))

# --- 8. Clean clinical covariates for Cox ---
clin_clean$age <- as.numeric(clin_clean$`age:ch1`)
clin_clean$Sex <- factor(clin_clean$`Sex:ch1`)
clin_clean$smoking_history <- factor(clin_clean$`smoking_history:ch1`)
clin_clean$disease_stage <- factor(clin_clean$`disease_stage:ch1`)

# --- 9. Univariable Cox model ---
cox_fit <- coxph(Surv(time, status) ~ EGFR_group, data = clin_clean)

# --- 10. Kaplan–Meier plot ---
fit_km <- survfit(Surv(time, status) ~ EGFR_group, data = clin_clean)

hr <- round(exp(coef(cox_fit)), 2)
ci <- round(exp(confint(cox_fit)), 2)
hr_label <- paste0("HR = ", hr, " (", ci[1], "–", ci[2], ")")

km_plot <- ggsurvplot(
  fit_km,
  data = clin_clean,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  legend.title = "EGFR Expression",
  legend.labs = c("High", "Low"),
  palette = c("#1f77b4", "#ff7f0e"),
  surv.median.line = "hv",
  ggtheme = theme_minimal(),
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE
)

km_plot$plot <- km_plot$plot +
  annotate(
    "text",
    x = max(clin_clean$time, na.rm = TRUE) * 0.7,
    y = 0.1,
    label = hr_label,
    size = 4
  )

print(km_plot)

# --- 11. Multivariable Cox model ---
cox_multi <- coxph(
  Surv(time, status) ~ EGFR_group + age + Sex + smoking_history + disease_stage,
  data = clin_clean
)
summary(cox_multi)

# --- 12. Forest plot for multivariable Cox ---
forest_plot <- ggforest(
  cox_multi,
  data = clin_clean,
  main = "Multivariable Cox",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1
)
print(forest_plot)

library(cowplot)  # or patchwork

# --- 1. KM plot with risk table ---
km_plot_full <- ggsurvplot(
  fit_km,
  data = clin_clean,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  legend.title = "EGFR Expression",
  legend.labs = c("High", "Low"),
  palette = c("#1f77b4", "#ff7f0e"),
  surv.median.line = "hv",
  ggtheme = theme_minimal(),
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE
)

# Convert ggplot objects from ggsurvplot
km_plot_main <- km_plot_full$plot
km_risk_table <- km_plot_full$table

# --- 2. Forest plot ---
forest_plot <- ggforest(
  cox_multi,
  data = clin_clean,
  main = "Multivariable Cox",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1
)

# --- 3. Combine KM + risk table + forest plot ---
combined_panel <- plot_grid(
  km_plot_main,
  km_risk_table,
  forest_plot,
  ncol = 1,
  rel_heights = c(3, 1, 3)  # adjust to balance space
)

# --- 4. Save for publication ---
ggsave("Survival_Figure.png", combined_panel, width = 8, height = 12, dpi = 300)