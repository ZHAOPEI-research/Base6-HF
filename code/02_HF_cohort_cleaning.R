# ============================================================
# Step 1
# NHANES 2011-2018 HF研究
# 原始变量类型、编码、缺失与异常值审计
# ============================================================

rm(list = ls())

library(dplyr)

# ------------------------------------------------------------
# 1. 读取上一阶段生成的数据
# ------------------------------------------------------------

load("NHANES_2011_2018_RAW_HFstudy.RData")

# 确认对象存在
if (!exists("final_all")) {
  stop("未找到 final_all，请检查 RData 文件。")
}

dat0 <- final_all

cat("============================================\n")
cat("Step 1: 原始数据审计\n")
cat("============================================\n")

cat("总样本量:", nrow(dat0), "\n")
cat("唯一 SEQN:", n_distinct(dat0$seqn), "\n\n")


# ============================================================
# 2. 研究核心变量清单
# ============================================================

audit_vars <- c(

  # ID / 周期 / survey
  "seqn",
  "cycle",
  "sddsrvyr",
  "wtint2yr",
  "wtmec2yr",
  "sdmvstra",
  "sdmvpsu",

  # 人口学
  "riagendr",
  "ridageyr",
  "ridreth3",
  "dmdeduc2",
  "dmdmartl",
  "indfmpir",

  # HF / cardiovascular
  "mcq160b",
  "mcq160c",
  "mcq160d",
  "mcq160e",
  "mcq160f",

  # BP
  "bpxsy1",
  "bpxsy2",
  "bpxsy3",
  "bpxdi1",
  "bpxdi2",
  "bpxdi3",
  "bpxpls",

  # BMI
  "bmxbmi",

  # CBC
  "lbxhgb",
  "lbxrdw",
  "lbxwbcsi",
  "lbdneno",
  "lbdlymno",

  # glucose / HbA1c
  "lbxglu",
  "lbdglusi",
  "lbxgh",

  # biochemistry
  "lbxscr",
  "lbdscrsi",
  "lbxsal",
  "lbdsalsi",
  "lbxsbu",
  "lbdsbusi",
  "lbxsnasi",

  # kidney
  "urdact",

  # inflammation
  "lbxhscrp",

  # lipid
  "lbdhdd",
  "lbdhddsi",
  "lbxtr",
  "lbdtrsi",
  "lbdldl",
  "lbdldnsi",
  "lbxtc",
  "lbdtcsi",

  # hypertension / diabetes / kidney
  "bpq020",
  "bpq040a",
  "bpq050a",
  "diq010",
  "diq050",
  "diq070",
  "kiq022",
  "kiq025",

  # smoking / alcohol
  "smq020",
  "smq040",
  "alq111",
  "alq121",

  # physical activity / function
  "paq605",
  "paq620",
  "paq650",
  "paq665",
  "pad680",
  "pfq061b",
  "pfq061c",

  # healthcare
  "hud080",

  # medication
  "med_cnt",
  "aceiarb_use",
  "bb_use",
  "mra_use",
  "sglt2i_use",
  "vericiguat_use",

  # mortality
  "eligstat",
  "mortstat",
  "ucod_leading",
  "permth_int",
  "permth_exm"
)

audit_vars <- intersect(
  audit_vars,
  names(dat0)
)


# ============================================================
# 3. 变量类型
# ============================================================

cat("\n============================================\n")
cat("A. 变量类型\n")
cat("============================================\n")

type_table <- data.frame(
  variable = audit_vars,
  class = sapply(
    dat0[audit_vars],
    function(x) paste(class(x), collapse = "/")
  ),
  stringsAsFactors = FALSE
)

print(
  type_table,
  row.names = FALSE
)


# ============================================================
# 4. 缺失情况
# ============================================================

cat("\n============================================\n")
cat("B. 缺失情况\n")
cat("============================================\n")

missing_table <- data.frame(
  variable = audit_vars,

  n_missing = sapply(
    dat0[audit_vars],
    function(x) sum(is.na(x))
  ),

  pct_missing = round(
    sapply(
      dat0[audit_vars],
      function(x) mean(is.na(x)) * 100
    ),
    2
  ),

  stringsAsFactors = FALSE
)

missing_table <- missing_table %>%
  arrange(
    desc(pct_missing)
  )

print(
  missing_table,
  row.names = FALSE
)


# ============================================================
# 5. 各周期缺失情况
# ============================================================

cat("\n============================================\n")
cat("C. 各周期核心变量缺失率\n")
cat("============================================\n")

core_missing_vars <- c(
  "mcq160b",
  "ridageyr",
  "riagendr",

  "bpxsy1",
  "bpxsy2",
  "bpxsy3",

  "bmxbmi",

  "lbxhgb",
  "lbxrdw",
  "lbxwbcsi",
  "lbdneno",

  "lbxscr",
  "lbxsal",
  "lbxsbu",
  "lbxsnasi",

  "pfq061b",

  "med_cnt",

  "mortstat",
  "permth_exm"
)

core_missing_vars <- intersect(
  core_missing_vars,
  names(dat0)
)

cycle_missing <- dat0 %>%
  group_by(cycle) %>%
  summarise(
    across(
      all_of(core_missing_vars),
      ~ round(mean(is.na(.x)) * 100, 2),
      .names = "{.col}_missing_pct"
    ),
    .groups = "drop"
  )

print(cycle_missing)


# ============================================================
# 6. 离散/问卷变量编码审计
# ============================================================

cat("\n============================================\n")
cat("D. 问卷变量编码审计\n")
cat("============================================\n")

categorical_vars <- c(

  "riagendr",
  "ridreth3",
  "dmdeduc2",
  "dmdmartl",

  "mcq160b",
  "mcq160c",
  "mcq160d",
  "mcq160e",
  "mcq160f",

  "bpq020",
  "bpq040a",
  "bpq050a",

  "diq010",
  "diq050",
  "diq070",

  "kiq022",
  "kiq025",

  "smq020",
  "smq040",

  "alq111",

  "paq605",
  "paq620",
  "paq650",
  "paq665",

  "pfq061b",
  "pfq061c",

  "aceiarb_use",
  "bb_use",
  "mra_use",
  "sglt2i_use",

  "mortstat",
  "eligstat"
)

categorical_vars <- intersect(
  categorical_vars,
  names(dat0)
)


for (v in categorical_vars) {

  cat("\n--------------------------------------------\n")
  cat(v, "\n")

  print(
    table(
      dat0[[v]],
      useNA = "ifany"
    )
  )
}


# ============================================================
# 7. 按周期查看关键问卷变量编码
# ============================================================

cat("\n============================================\n")
cat("E. 按周期检查关键变量编码\n")
cat("============================================\n")

key_code_vars <- c(
  "mcq160b",
  "pfq061b",
  "bpq020",
  "diq010",
  "smq020"
)

key_code_vars <- intersect(
  key_code_vars,
  names(dat0)
)

for (v in key_code_vars) {

  cat("\n##############################\n")
  cat("变量:", v, "\n")

  print(
    dat0 %>%
      count(
        cycle,
        .data[[v]],
        name = "n"
      ) %>%
      arrange(
        cycle,
        .data[[v]]
      )
  )
}


# ============================================================
# 8. 连续变量范围审计
# ============================================================

cat("\n============================================\n")
cat("F. 连续变量范围\n")
cat("============================================\n")

continuous_vars <- c(

  "ridageyr",

  "bpxsy1",
  "bpxsy2",
  "bpxsy3",
  "bpxdi1",
  "bpxdi2",
  "bpxdi3",
  "bpxpls",

  "bmxbmi",

  "lbxhgb",
  "lbxrdw",
  "lbxwbcsi",
  "lbdneno",
  "lbdlymno",

  "lbxglu",
  "lbxgh",

  "lbxscr",
  "lbxsal",
  "lbxsbu",
  "lbxsnasi",

  "urdact",

  "lbxhscrp",

  "lbdhdd",
  "lbxtr",
  "lbdldl",
  "lbxtc",

  "pad680",

  "med_cnt",

  "permth_int",
  "permth_exm"
)

continuous_vars <- intersect(
  continuous_vars,
  names(dat0)
)


range_table <- do.call(
  rbind,
  lapply(
    continuous_vars,
    function(v) {

      x <- suppressWarnings(
        as.numeric(
          as.character(
            dat0[[v]]
          )
        )
      )

      data.frame(
        variable = v,
        n_nonmissing = sum(!is.na(x)),
        min = ifelse(
          all(is.na(x)),
          NA,
          min(x, na.rm = TRUE)
        ),
        p01 = ifelse(
          all(is.na(x)),
          NA,
          quantile(
            x,
            0.01,
            na.rm = TRUE
          )
        ),
        median = ifelse(
          all(is.na(x)),
          NA,
          median(
            x,
            na.rm = TRUE
          )
        ),
        p99 = ifelse(
          all(is.na(x)),
          NA,
          quantile(
            x,
            0.99,
            na.rm = TRUE
          )
        ),
        max = ifelse(
          all(is.na(x)),
          NA,
          max(x, na.rm = TRUE)
        )
      )
    }
  )
)

print(
  range_table,
  row.names = FALSE
)


# ============================================================
# 9. 重点检查 HF 定义
# ============================================================

cat("\n============================================\n")
cat("G. HF候选人群初步检查\n")
cat("============================================\n")

hf_check <- dat0 %>%
  mutate(
    mcq160b_num =
      suppressWarnings(
        as.numeric(
          as.character(
            mcq160b
          )
        )
      )
  ) %>%
  count(
    cycle,
    mcq160b_num,
    name = "n"
  )

print(hf_check)


# ============================================================
# 10. PFQ061B 特别检查
# ============================================================

cat("\n============================================\n")
cat("H. PFQ061B 行走困难编码\n")
cat("============================================\n")

pfq_check <- dat0 %>%
  count(
    cycle,
    pfq061b,
    name = "n"
  ) %>%
  arrange(
    cycle,
    pfq061b
  )

print(pfq_check)


# ============================================================
# 11. 中性粒/WBC 原始逻辑检查
# ============================================================

cat("\n============================================\n")
cat("I. NWR 原始变量检查\n")
cat("============================================\n")

if (
  all(
    c(
      "lbdneno",
      "lbxwbcsi"
    ) %in%
    names(dat0)
  )
) {

  nwr_check <- dat0 %>%
    mutate(

      neut =
        suppressWarnings(
          as.numeric(
            as.character(
              lbdneno
            )
          )
        ),

      wbc =
        suppressWarnings(
          as.numeric(
            as.character(
              lbxwbcsi
            )
          )
        ),

      nwr_raw =
        neut / wbc
    ) %>%
    summarise(

      neut_missing =
        sum(
          is.na(neut)
        ),

      wbc_missing =
        sum(
          is.na(wbc)
        ),

      wbc_zero =
        sum(
          wbc == 0,
          na.rm = TRUE
        ),

      nwr_gt1 =
        sum(
          nwr_raw > 1,
          na.rm = TRUE
        ),

      nwr_lt0 =
        sum(
          nwr_raw < 0,
          na.rm = TRUE
        ),

      nwr_min =
        min(
          nwr_raw,
          na.rm = TRUE
        ),

      nwr_max =
        max(
          nwr_raw,
          na.rm = TRUE
        )
    )

  print(nwr_check)
}


# ============================================================
# 12. Creatinine 单位一致性检查
# ============================================================

cat("\n============================================\n")
cat("J. Creatinine 单位一致性检查\n")
cat("============================================\n")

if (
  all(
    c(
      "lbxscr",
      "lbdscrsi"
    ) %in%
    names(dat0)
  )
) {

  cr_check <- dat0 %>%
    transmute(

      scr_mgdl =
        suppressWarnings(
          as.numeric(
            as.character(
              lbxscr
            )
          )
        ),

      scr_umol =
        suppressWarnings(
          as.numeric(
            as.character(
              lbdscrsi
            )
          )
        ),

      ratio =
        scr_umol / scr_mgdl
    ) %>%
    summarise(

      n_pair =
        sum(
          !is.na(scr_mgdl) &
            !is.na(scr_umol)
        ),

      mgdl_min =
        min(
          scr_mgdl,
          na.rm = TRUE
        ),

      mgdl_max =
        max(
          scr_mgdl,
          na.rm = TRUE
        ),

      umol_min =
        min(
          scr_umol,
          na.rm = TRUE
        ),

      umol_max =
        max(
          scr_umol,
          na.rm = TRUE
        ),

      conversion_median =
        median(
          ratio,
          na.rm = TRUE
        )
    )

  print(cr_check)
}


# ============================================================
# 13. Survey weight 基础审计
# ============================================================

cat("\n============================================\n")
cat("K. Survey design变量检查\n")
cat("============================================\n")

survey_check <- dat0 %>%
  group_by(cycle) %>%
  summarise(

    n = n(),

    wtint_missing =
      sum(
        is.na(
          wtint2yr
        )
      ),

    wtmec_missing =
      sum(
        is.na(
          wtmec2yr
        )
      ),

    wtmec_zero =
      sum(
        suppressWarnings(
          as.numeric(
            as.character(
              wtmec2yr
            )
          )
        ) <= 0,
        na.rm = TRUE
      ),

    strata_missing =
      sum(
        is.na(
          sdmvstra
        )
      ),

    psu_missing =
      sum(
        is.na(
          sdmvpsu
        )
      ),

    .groups = "drop"
  )

print(survey_check)


# ============================================================
# 14. 保存审计对象
# ============================================================

save(
  dat0,
  type_table,
  missing_table,
  cycle_missing,
  range_table,
  file =
    "NHANES_2011_2018_STEP1_AUDIT.RData"
)

cat("\n============================================\n")
cat("Step 1 审计完成\n")
cat("============================================\n")
cat(
  "已保存：NHANES_2011_2018_STEP1_AUDIT.RData\n"
)
