# Lab 1：练习 1--3
# 运行：Rscript --vanilla lab1_exercises_1_3.R
# 或在 RStudio 中打开本文件并点击 Source（文件编码为 UTF-8）。
# 只使用 R 自带的 base / stats / utils，不需要额外安装程序包。

# 将结果保存到本脚本旁边的 lab1_results 文件夹。
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_file <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[1])
} else if (!is.null(source_file)) {
  source_file
} else {
  NULL
}
output_dir <- file.path(if (is.null(script_file)) getwd() else
                         dirname(normalizePath(script_file)), "lab1_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# 练习 1：相互独立的 Bernoulli(0.3)、Bernoulli(0.6) 数据 -----------
# 题目未指定 n，选取 20、100、1000，比较小样本和大样本。
# 先产生全部数据，再做置换，避免置换次数影响后续数据的生成。
set.seed(111)
sample_sizes <- c(20L, 100L, 1000L)
bernoulli_data <- lapply(sample_sizes, function(n) {
  data.frame(x = rbinom(n, size = 1, prob = 0.3),
             y = rbinom(n, size = 1, prob = 0.6))
})

set.seed(112)
N1 <- 100000L
exercise1 <- do.call(rbind, lapply(seq_along(sample_sizes), function(j) {
  x <- bernoulli_data[[j]]$x
  y <- bernoulli_data[[j]]$y
  n <- length(x)
  if (sd(x) == 0 || sd(y) == 0) {
    stop("A sample is constant; its correlation is undefined.")
  }
  r <- cor(x, y)
  t_test <- cor.test(x, y, method = "pearson", alternative = "two.sided")
  z <- sqrt(n - 2) * r             # 按讲义使用 sqrt(n - 2)，不是 sqrt(n)
  p_z <- 2 * pnorm(-abs(z))

  # 置换时均值、方差不变，使用与相关系数成正比的整数统计量：
  # T = n * sum(x * y) - sum(x) * sum(y)。整数比较避免边界舍入误差。
  sx <- sum(x)
  sy <- sum(y)
  T0 <- n * sum(x * y) - sx * sy
  T_perm <- replicate(N1, {
    x_perm <- x[sample.int(n)]
    n * sum(x_perm * y) - sx * sy
  })
  extreme <- sum(abs(T_perm) >= abs(T0))
  p_perm <- extreme / N1           # 严格使用讲义公式 (1)

  # 二元数据可用超几何分布给出同一统计量的精确置换 p 值，作核对。
  # 固定两个边际和时，K = sum(x_perm * y) ~ Hypergeom(sx, n-sx, sy)。
  k <- seq.int(max(0, sx + sy - n), min(sx, sy))
  p_exact <- sum(dhyper(k, sx, n - sx, sy)[
    abs(n * k - sx * sy) >= abs(T0)
  ])
  write.csv(bernoulli_data[[j]],
            file.path(output_dir, paste0("exercise1_data_n", n, ".csv")),
            row.names = FALSE)
  data.frame(n = n, r = r, t = unname(t_test$statistic), z = z,
             p_t = t_test$p.value, p_z = p_z, N = N1,
             extreme = extreme, p_permutation = p_perm,
             MCSE = sqrt(p_perm * (1 - p_perm) / N1),
             p_exact_reference = p_exact)
}))

# 练习 2：Spearman 双侧置换检验 ------------------------------
x <- c(-11, -2, 5, 6, 8, 9)
y <- c(-3, 2, 4, 1, 12, 5)
n <- length(x)
rank_x <- rank(x)
rank_y <- rank(y)
rho0 <- cor(rank_x, rank_y)
spearman_exact <- cor.test(x, y, method = "spearman", exact = TRUE)

# 对固定秩序列，用 n 倍的中心化乘积和表示 Spearman 系数的分子。
# 与直接比较 abs(cor(rank_x_perm, rank_y)) >= abs(rho0) 等价。
rank_offset <- sum(rank_x) * sum(rank_y)
S0 <- n * sum(rank_x * rank_y) - rank_offset
set.seed(222)
permutation_counts <- c(1000L, 10000L, 100000L)
exercise2 <- do.call(rbind, lapply(permutation_counts, function(N) {
  # 三个 N 各自执行指定次数的随机置换；一次置换内部不放回，
  # 不同次置换之间允许重复，这是随机置换模拟的正常方式。
  S_perm <- replicate(N, {
    n * sum(rank_x[sample.int(n)] * rank_y) - rank_offset
  })
  extreme <- sum(abs(S_perm) >= abs(S0))
  p <- extreme / N
  data.frame(N = N, extreme = extreme, p_permutation = p,
             MCSE = sqrt(p * (1 - p) / N),
             p_exact = spearman_exact$p.value,
             absolute_error = abs(p - spearman_exact$p.value))
}))

# 穷举 6! = 720 个排列，用于核验练习 2、3 的精确结果。
all_permutations <- function(v) {
  if (length(v) == 1L) return(matrix(v, nrow = 1L))
  do.call(rbind, lapply(seq_along(v), function(i) {
    cbind(v[i], all_permutations(v[-i]))
  }))
}
perms <- all_permutations(seq_len(n))
S_all <- apply(perms, 1L, function(idx) {
  n * sum(rank_x[idx] * rank_y) - rank_offset
})
spearman_extreme <- sum(abs(S_all) >= abs(S0))
spearman_enumerated <- spearman_extreme / nrow(perms)

# 练习 3：逆序数 D 的置换检验 --------------------------------
# 先按 x 从小到大排序，D = sum_{i<j} I(y_i > y_j)。
y_sorted <- y[order(x)]
pairs <- combn(seq_len(n), 2L)
count_inversions <- function(v) {
  sum(v[pairs[1L, ]] > v[pairs[2L, ]])
}
D0 <- count_inversions(y_sorted)
M <- choose(n, 2)                 # 总对数 15
tau0 <- 1 - 2 * D0 / M

set.seed(333)
N3 <- 100000L
D_perm <- replicate(N3, count_inversions(y_sorted[sample.int(n)]))
D_all <- apply(perms, 1L, function(idx) count_inversions(y_sorted[idx]))

# 正相关单侧备择：D 越小越极端，p = P(D_perm <= D0)。
# 双侧备择：以 M/2 为中心，等价于比较 |Kendall tau|。
one_tail <- function(D) D <= D0
two_tail <- function(D) abs(M - 2 * D) >= abs(M - 2 * D0)
exercise3 <- do.call(rbind, lapply(c("positive_one_sided", "two_sided"), function(test) {
  is_extreme <- if (test == "positive_one_sided") one_tail else two_tail
  extreme <- sum(is_extreme(D_perm))
  exact_extreme <- sum(is_extreme(D_all))
  p <- extreme / N3
  data.frame(test = test, D_observed = D0, tau_observed = tau0,
             N = N3, extreme = extreme, p_permutation = p,
             MCSE = sqrt(p * (1 - p) / N3),
             exact_extreme = exact_extreme, total_permutations = nrow(perms),
             p_exact = exact_extreme / nrow(perms))
}))

# 正确性检查：穷举结果必须与 R 内置精确检验一致。
kendall_one <- cor.test(x, y, method = "kendall", exact = TRUE,
                        alternative = "greater")
kendall_two <- cor.test(x, y, method = "kendall", exact = TRUE,
                        alternative = "two.sided")
stopifnot(
  nrow(unique(perms)) == factorial(n),
  isTRUE(all.equal(spearman_enumerated, spearman_exact$p.value)),
  D0 == 3, isTRUE(all.equal(tau0, unname(kendall_two$estimate))),
  isTRUE(all.equal(exercise3$p_exact[1], kendall_one$p.value)),
  isTRUE(all.equal(exercise3$p_exact[2], kendall_two$p.value))
)

# 导出原始结果，便于复核和修改参数后重跑。
write.csv(exercise1, file.path(output_dir, "exercise1.csv"), row.names = FALSE)
write.csv(exercise2, file.path(output_dir, "exercise2.csv"), row.names = FALSE)
write.csv(exercise3, file.path(output_dir, "exercise3.csv"), row.names = FALSE)
transcript <- capture.output({
  cat("Lab 1 - Exercises 1-3\n", R.version.string, "\n", sep = "")
  cat("Seeds: exercise 1 data=111, permutation=112; exercise 2=222; exercise 3=333.\n")
  cat("\nExercise 1: independent Bernoulli data; two-sided tests\n")
  print(exercise1, digits = 9, row.names = FALSE)
  cat("\nExercise 2: Spearman rho =", rho0, "\n")
  print(exercise2, digits = 9, row.names = FALSE)
  cat("Exact enumeration:", spearman_extreme, "/", nrow(perms),
      "=", spearman_enumerated, "\n")
  cat("\nExercise 3: inversion count D =", D0, "; Kendall tau =", tau0, "\n")
  print(exercise3, digits = 9, row.names = FALSE)
  cat("\nAll exact-enumeration checks passed.\n\n")
  print(sessionInfo())
})
writeLines(transcript, file.path(output_dir, "results.txt"), useBytes = TRUE)
cat(paste(transcript, collapse = "\n"), "\n")
