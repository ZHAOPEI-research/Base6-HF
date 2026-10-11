rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Public-repository path configuration
# Restricted external-cohort data and derived individual-level outputs
# must remain outside the public repository.
.base6_external_output_dir <- Sys.getenv("BASE6_EXTERNAL_OUTPUT_DIR")

if (.base6_external_output_dir == "") {
  stop(
    paste0(
      "Set BASE6_EXTERNAL_OUTPUT_DIR to an authorized local directory ",
      "for restricted external-cohort analysis outputs."
    )
  )
}

if (!dir.exists(.base6_external_output_dir)) {
  dir.create(
    .base6_external_output_dir,
    recursive = TRUE
  )
}

setwd(.base6_external_output_dir)

# ============================================================
# 1. Load STEP 52 results
# ============================================================

load(
  "Chinese_GWTG_STEP11P_SUBGROUP_CI_2000BOOT.RData"
)

stopifnot(
  exists("final_subgroup_CI")
)

dat <- final_subgroup_CI


# ============================================================
# 2. Standardize subgroup labels and order
# ============================================================

dat <- dat %>%
  mutate(

    subgroup_label =
      case_when(

        subgroup == "Sex" &
          level == "Female" ~
          paste0(
            "Sex: Female (n=",
            n,
            ", deaths=",
            deaths,
            ")"
          ),

        subgroup == "Sex" &
          level == "Male" ~
          paste0(
            "Sex: Male (n=",
            n,
            ", deaths=",
            deaths,
            ")"
          ),

        subgroup == "Age" &
          level == ">=70" ~
          paste0(
            "Age: ≥70 years (n=",
            n,
            ", deaths=",
            deaths,
            ")"
          ),

        subgroup == "LVEF" &
          level == "<50%" ~
          paste0(
            "LVEF: <50% (n=",
            n,
            ", deaths=",
            deaths,
            ")"
          ),

        subgroup == "LVEF" &
          level == ">=50%" ~
          paste0(
            "LVEF: ≥50% (n=",
            n,
            ", deaths=",
            deaths,
            ")"
          ),

        subgroup == "CKD" &
          level == "No CKD" ~
          paste0(
            "CKD: No CKD (n=",
            n,
            ", deaths=",
            deaths,
            ")"
          ),

        subgroup == "CKD" &
          level == "CKD" ~
          paste0(
            "CKD: CKD (n=",
            n,
            ", deaths=",
            deaths,
            ")"
          ),

        subgroup == "DM" &
          level == "No DM" ~
          paste0(
            "DM: No DM (n=",
            n,
            ", deaths=",
            deaths,
            ")"
          ),

        subgroup == "AF" &
          level == "AF" ~
          paste0(
            "AF: AF (n=",
            n,
            ", deaths=",
            deaths,
            ")"
          ),

        subgroup == "AF" &
          level == "No AF" ~
          paste0(
            "AF: No AF (n=",
            n,
            ", deaths=",
            deaths,
            ")"
          ),

        TRUE ~
          paste(
            subgroup,
            level,
            sep = ": "
          )
      )
  )


# Desired top-to-bottom order
subgroup_order <- c(
  dat$subgroup_label[
    dat$subgroup == "Sex" &
      dat$level == "Female"
  ],

  dat$subgroup_label[
    dat$subgroup == "Sex" &
      dat$level == "Male"
  ],

  dat$subgroup_label[
    dat$subgroup == "Age" &
      dat$level == ">=70"
  ],

  dat$subgroup_label[
    dat$subgroup == "LVEF" &
      dat$level == "<50%"
  ],

  dat$subgroup_label[
    dat$subgroup == "LVEF" &
      dat$level == ">=50%"
  ],

  dat$subgroup_label[
    dat$subgroup == "CKD" &
      dat$level == "No CKD"
  ],

  dat$subgroup_label[
    dat$subgroup == "CKD" &
      dat$level == "CKD"
  ],

  dat$subgroup_label[
    dat$subgroup == "DM" &
      dat$level == "No DM"
  ],

  dat$subgroup_label[
    dat$subgroup == "AF" &
      dat$level == "AF"
  ],

  dat$subgroup_label[
    dat$subgroup == "AF" &
      dat$level == "No AF"
  ]
)

dat$subgroup_label <- factor(
  dat$subgroup_label,
  levels = rev(subgroup_order)
)


# ============================================================
# 3. Build long AUC dataset
# ============================================================

auc_long <- bind_rows(

  dat %>%
    transmute(
      subgroup_label,
      Model = "Base6",
      estimate = AUC_Base6,
      lower = AUC_Base6_low,
      upper = AUC_Base6_high
    ),

  dat %>%
    transmute(
      subgroup_label,
      Model = "MAGGIC",
      estimate = AUC_MAGGIC,
      lower = AUC_MAGGIC_low,
      upper = AUC_MAGGIC_high
    ),

  dat %>%
    transmute(
      subgroup_label,
      Model = "GWTG-HF",
      estimate = AUC_GWTG,
      lower = AUC_GWTG_low,
      upper = AUC_GWTG_high
    )
)


# ============================================================
# 4. Build long C-index dataset
# ============================================================

c_long <- bind_rows(

  dat %>%
    transmute(
      subgroup_label,
      Model = "Base6",
      estimate = C_Base6,
      lower = C_Base6_low,
      upper = C_Base6_high
    ),

  dat %>%
    transmute(
      subgroup_label,
      Model = "MAGGIC",
      estimate = C_MAGGIC,
      lower = C_MAGGIC_low,
      upper = C_MAGGIC_high
    ),

  dat %>%
    transmute(
      subgroup_label,
      Model = "GWTG-HF",
      estimate = C_GWTG,
      lower = C_GWTG_low,
      upper = C_GWTG_high
    )
)


# ============================================================
# 5.Model order + colors + shapes
# ============================================================

model_order <- c(
  "Base6",
  "GWTG-HF",
  "MAGGIC"
)

auc_long$Model <- factor(
  auc_long$Model,
  levels = model_order
)

c_long$Model <- factor(
  c_long$Model,
  levels = model_order
)

model_colors <- c(
  "Base6"   = "#0072B2",
  "GWTG-HF" = "#E69F00",
  "MAGGIC"  = "#D55E00"
)

model_shapes <- c(
  "Base6"   = 16,   # circle
  "GWTG-HF" = 15,   # square
  "MAGGIC"  = 17    # triangle
)

pd <- position_dodge(
  width = 0.62
)

# ============================================================
# 6. Common plotting settings
# ============================================================

pd <- position_dodge(
  width = 0.62
)

common_theme <- theme_classic(
  base_size = 10
) +
  theme(
    axis.title.y =
      element_blank(),

    axis.text.y =
      element_text(
        size = 8.5
      ),

    axis.text.x =
      element_text(
        size = 9
      ),

    axis.title.x =
      element_text(
        size = 10
      ),

    legend.title =
      element_blank(),

    legend.position =
      "bottom",

    legend.text =
      element_text(
        size = 9
      ),

    plot.title =
      element_text(
        face = "bold",
        size = 11,
        hjust = 0
      ),

    plot.margin =
      margin(
        5,
        8,
        5,
        5
      )
  )


# ============================================================
# 7. Panel A: 364-day AUC
# ============================================================

# ============================================================
# Panel A: 364-day AUC
# ============================================================

p_auc <- ggplot(
  auc_long,
  aes(
    x = estimate,
    y = subgroup_label,
    color = Model,
    shape = Model
  )
) +

  geom_vline(
    xintercept = 0.5,
    linewidth = 0.35,
    linetype = "dashed",
    color = "grey40"
  ) +

  geom_errorbar(
    aes(
      xmin = lower,
      xmax = upper
    ),
    orientation = "y",
    position = pd,
    width = 0,
    linewidth = 0.55
  ) +

  geom_point(
    position = pd,
    size = 2.5,
    stroke = 0.8
  ) +

  scale_color_manual(
    values = model_colors,
    breaks = model_order
  ) +

  scale_shape_manual(
    values = model_shapes,
    breaks = model_order
  ) +

  scale_x_continuous(
    limits = c(0.30, 0.90),
    breaks = seq(0.3, 0.9, 0.1),
    expand = expansion(
      mult = c(0, 0.01)
    )
  ) +

  labs(
    x = "364-day AUC",
    title = "A"
  ) +

  common_theme


# ============================================================
# 8. Panel B: Harrell C-index
# ============================================================

# ============================================================
# Panel B: Harrell C-index
# ============================================================

p_c <- ggplot(
  c_long,
  aes(
    x = estimate,
    y = subgroup_label,
    color = Model,
    shape = Model
  )
) +

  geom_vline(
    xintercept = 0.5,
    linewidth = 0.35,
    linetype = "dashed",
    color = "grey40"
  ) +

  geom_errorbar(
    aes(
      xmin = lower,
      xmax = upper
    ),
    orientation = "y",
    position = pd,
    width = 0,
    linewidth = 0.55
  ) +

  geom_point(
    position = pd,
    size = 2.5,
    stroke = 0.8
  ) +

  scale_color_manual(
    values = model_colors,
    breaks = model_order
  ) +

  scale_shape_manual(
    values = model_shapes,
    breaks = model_order
  ) +

  scale_x_continuous(
    limits = c(0.30, 0.90),
    breaks = seq(0.3, 0.9, 0.1),
    expand = expansion(
      mult = c(0, 0.01)
    )
  ) +

  labs(
    x = "Harrell C-index",
    title = "B"
  ) +

  common_theme


# ============================================================
# 9. Combine
# ============================================================

figS4 <-
  p_auc +
  p_c +
  plot_layout(
    ncol = 2,
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )


# ============================================================
# 10. Display
# ============================================================

print(figS4)


# ============================================================
# 11. Save publication versions
# ============================================================

ggsave(
  filename =
    "FINAL_FigureS4_SubgroupDiscrimination_CI.pdf",
  plot = figS4,
  width = 12,
  height = 6.6,
  units = "in"
)

ggsave(
  filename =
    "FINAL_FigureS4_SubgroupDiscrimination_CI.tiff",
  plot = figS4,
  width = 12,
  height = 6.6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename =
    "FINAL_FigureS4_SubgroupDiscrimination_CI.png",
  plot = figS4,
  width = 12,
  height = 6.6,
  units = "in",
  dpi = 600
)


# ============================================================
# 12. Save source object
# ============================================================

save(
  dat,
  auc_long,
  c_long,
  p_auc,
  p_c,
  figS4,
  file =
    "FINAL_FigureS4_SubgroupDiscrimination_CI.RData"
)


cat(
  "\n============================================\n"
)

cat(
  "Figure S4 with 95% CI completed\n"
)

cat(
  "============================================\n"
)
