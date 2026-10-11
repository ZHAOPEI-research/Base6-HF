# ============================================================
# STEP 4C
# MICE convergence and distribution diagnostics
# ============================================================

rm(list = ls())

library(mice)
library(lattice)

analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
setwd(analysis_dir)

load("NHANES_2011_2018_STEP4B_MICE_REVISED.RData")

cat("============================================\n")
cat("STEP 4C: MICE diagnostics\n")
cat("============================================\n")


# ============================================================
# 1. 基本信息
# ============================================================

cat("\n============================================\n")
cat("A. MICE基本信息\n")
cat("============================================\n")

print(imp)

cat(
  "\nImputations:",
  imp$m,
  "\n"
)

cat(
  "Iterations:",
  imp$iteration,
  "\n"
)


# ============================================================
# 2. 自动 convergence statistics
# ============================================================

cat("\n============================================\n")
cat("B. Convergence statistics\n")
cat("============================================\n")

conv <- convergence(imp)

print(conv)

# 保存
write.csv(
  conv,
  "MICE_convergence_statistics.csv",
  row.names = FALSE
)


# ============================================================
# 3. 查看每个变量的 chain mean / SD
# ============================================================

vars_diag <- c(
  "sbp_mean",
  "bpxpls",
  "bmxbmi",
  "lbxhgb",
  "lbxrdw",
  "lbdneno",
  "lbxwbcsi",
  "scr_mgdl",
  "lbxsal",
  "lbxsbu",
  "lbxsnasi"
)

vars_diag <- intersect(
  vars_diag,
  names(imp$imp)
)


# ============================================================
# 4. Trace plots
#
# 会依次在绘图窗口生成
# ============================================================

cat("\n============================================\n")
cat("C. Trace plots\n")
cat("============================================\n")

for (v in vars_diag) {

  cat(
    "绘制:",
    v,
    "\n"
  )

  print(
    plot(
      imp,
      v
    )
  )
}


# ============================================================
# 5. Density plots
#
# 蓝/红等颜色由 lattice 自动控制
# 观察 observed vs imputed distribution
# ============================================================

cat("\n============================================\n")
cat("D. Density plots\n")
cat("============================================\n")

for (v in vars_diag) {

  cat(
    "Density:",
    v,
    "\n"
  )

  print(
    densityplot(
      imp,
      as.formula(
        paste0(
          "~",
          v
        )
      )
    )
  )
}


# ============================================================
# 6. Strip plots
# ============================================================

cat("\n============================================\n")
cat("E. Strip plots\n")
cat("============================================\n")

for (v in vars_diag) {

  print(
    stripplot(
      imp,
      as.formula(
        paste0(
          v,
          " ~ .imp"
        )
      ),
      pch = 20,
      cex = 0.5
    )
  )
}


# ============================================================
# 7. Chain mean最终5次迭代稳定性
#
# 对每个变量比较最后5次迭代均值
# ============================================================

cat("\n============================================\n")
cat("F. 最后5次迭代 chain mean 稳定性\n")
cat("============================================\n")

chain_mean <- imp$chainMean

stability_results <- list()

for (v in vars_diag) {

  if (!v %in% dimnames(chain_mean)[[1]]) {
    next
  }

  x <- chain_mean[v, , , drop = FALSE]

  # dimensions:
  # variable x iteration x chain

  iter_n <- dim(x)[2]

  start_iter <- max(
    1,
    iter_n - 4
  )

  last5 <- x[
    1,
    start_iter:iter_n,
    ,
    drop = FALSE
  ]

  chain_means_final <- apply(
    last5,
    3,
    mean,
    na.rm = TRUE
  )

  stability_results[[v]] <- data.frame(
    variable = v,
    chain = seq_along(chain_means_final),
    final5_mean = chain_means_final
  )
}


stability_table <- do.call(
  rbind,
  stability_results
)

print(
  stability_table
)

write.csv(
  stability_table,
  "MICE_chain_final5_means.csv",
  row.names = FALSE
)


# ============================================================
# 8. 插补数据集之间的总体均值波动
# ============================================================

cat("\n============================================\n")
cat("G. 30个插补数据集之间的均值波动\n")
cat("============================================\n")

completed_means <- list()

for (m in 1:imp$m) {

  d <-
    complete(
      imp,
      m
    )

  one <- data.frame(
    imputation = m
  )

  for (v in vars_diag) {

    one[[v]] <-
      mean(
        d[[v]],
        na.rm = TRUE
      )
  }

  completed_means[[m]] <-
    one
}

completed_means <-
  do.call(
    rbind,
    completed_means
  )


mean_variation <- data.frame(
  variable = vars_diag,

  mean_of_means =
    sapply(
      vars_diag,
      function(v)
        mean(
          completed_means[[v]]
        )
    ),

  sd_between_imputations =
    sapply(
      vars_diag,
      function(v)
        sd(
          completed_means[[v]]
        )
    ),

  min_mean =
    sapply(
      vars_diag,
      function(v)
        min(
          completed_means[[v]]
        )
    ),

  max_mean =
    sapply(
      vars_diag,
      function(v)
        max(
          completed_means[[v]]
        )
    )
)


print(
  mean_variation,
  row.names = FALSE
)

write.csv(
  mean_variation,
  "MICE_between_imputation_means.csv",
  row.names = FALSE
)


# ============================================================
# 9. 保存诊断对象
# ============================================================

save(
  conv,
  stability_table,
  completed_means,
  mean_variation,
  file =
    "NHANES_2011_2018_STEP4C_MICE_DIAGNOSTICS.RData"
)


cat("\n============================================\n")
cat("STEP 4C 完成\n")
cat("============================================\n")
