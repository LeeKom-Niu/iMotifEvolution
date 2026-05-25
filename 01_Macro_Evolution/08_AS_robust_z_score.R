# ==================================================================
# 祖先状态重建 · 最终版 (PNG+PDF + 三大演化节点净增量分析)
# 策略：用外群节点代替父节点
# ==================================================================

library(ape)
library(phytools)

# ---------- 1. 读取数据 ----------
tree <- read.tree("species_tree.nwk")
dat  <- read.delim("species_taxid_sampled_final_400.tsv", header = TRUE, 
                   sep = "\t", stringsAsFactors = FALSE)

# ---------- 2. 名称清洗与匹配 ----------
tree$tip.label <- gsub(":.*", "", tree$tip.label)
rownames(dat) <- dat$organism_name_asm_name

common <- intersect(tree$tip.label, rownames(dat))
cat("树中叶节点:", length(tree$tip.label), " 表格中物种:", nrow(dat),
    " 共有:", length(common), "\n")

missing <- setdiff(tree$tip.label, rownames(dat))
if (length(missing) > 0) {
  cat("移除树中缺失物种:", length(missing), "\n")
  tree <- drop.tip(tree, missing)
}
dat <- dat[tree$tip.label, ]

# ---------- 3. 定根 ----------
tree <- midpoint.root(tree)
cat("树已定根，总节点数:", tree$Nnode + Ntip(tree), "\n")

# ---------- 4. 构建性状向量 ----------
trait <- setNames(as.numeric(dat[tree$tip.label, "robust_z_score"]), tree$tip.label)
if (any(is.na(trait))) {
  na_tips <- names(trait)[is.na(trait)]
  cat("移除", length(na_tips), "个 NA 值物种\n")
  tree <- drop.tip(tree, na_tips)
  trait <- trait[!is.na(trait)]
}

# ---------- 5. 祖先状态重建 ----------
fit <- fastAnc(tree, trait, vars = TRUE, CI = TRUE)
ace     <- fit$ace
ace_CI95 <- fit$CI95
cat("祖先状态重建完成\n")

# ---------- 6. 可视化 ----------
obj <- contMap(tree, trait, plot = FALSE)
obj <- setMap(obj, colors = colorRampPalette(c("blue", "white", "red"))(100))

pdf("ancestral_state_reconstruction.pdf", width = 18, height = 30)
plot(obj, fsize = 0.2, lwd = 1.5, outline = FALSE)
dev.off()

png("ancestral_state_reconstruction.png", width = 18, height = 30, units = "in", res = 300)
plot(obj, fsize = 0.2, lwd = 1.5, outline = FALSE)
dev.off()
cat("图形已保存 PDF 和 PNG\n")

# ---------- 7. 辅助函数 ----------
# 获取指定类群的 MRCA 节点号和祖先值
get_node_val <- function(groups, label) {
  tips <- tree$tip.label[dat$group %in% groups]
  if (length(tips) < 2) return(NULL)
  node <- getMRCA(tree, tips)
  node_idx <- as.character(node)
  list(node = node, label = label,
       z = ace[node_idx],
       ci_low = ace_CI95[node_idx, 1],
       ci_up  = ace_CI95[node_idx, 2])
}

# 计算两个节点的净增量和置信区间
calc_delta <- function(target, reference) {
  if (is.null(target) || is.null(reference)) return(NULL)
  dz <- target$z - reference$z
  d_low <- target$ci_low - reference$ci_up
  d_up  <- target$ci_up - reference$ci_low
  c(dz, d_low, d_up)
}

# ---------- 8. 定义节点 ----------
# 1) 真核冠群（所有真核生物）
euk_node <- get_node_val(c("Fungi","Plant","Invertebrate","Vertebrate Other",
                           "Mammalian","Protozoa"), "Eukaryota_crown")
# 原核冠群
prok_node <- get_node_val(c("Archaea","Bacteria"), "Prokaryota_crown")

# 2) 动物冠群（后生动物）
metazoa_node <- get_node_val(c("Invertebrate","Vertebrate Other","Mammalian"),
                             "Metazoa_crown")
# 非动物真核（除动物外的所有真核：真菌+植物+原生动物，这里简化，如有其他需调整）
non_metazoa_euk_node <- get_node_val(c("Fungi","Plant","Protozoa"),
                                     "NonMetazoa_Eukaryota_crown")

# 3) 脊椎动物冠群
vert_node <- get_node_val(c("Vertebrate Other","Mammalian"),
                          "Vertebrata_crown")
# 非脊椎动物后口动物/无脊椎动物（作为外群）
# 树中有 Invertebrate 类群，我们选 Invertebrate + maybe 其他非脊椎动物
non_vert_node <- get_node_val(c("Invertebrate"), "Invertebrate_crown")

# 如果没有独立的 Invertebrate 节点，也可以用非脊椎动物+脊椎动物的父节点，但这里 Invertebrate 存在

# ---------- 9. 计算净增量 ----------
res <- list()

calc_and_add <- function(target, ref, target_label, ref_label) {
  delta <- calc_delta(target, ref)
  if (!is.null(delta)) {
    res[[length(res) + 1]] <<- data.frame(
      Transition = paste0(target_label, " - ", ref_label),
      Target_Z = target$z,
      Reference_Z = ref$z,
      Delta_Z = delta[1],
      Delta_Lower = delta[2],
      Delta_Upper = delta[3],
      stringsAsFactors = FALSE
    )
  }
}

calc_and_add(euk_node, prok_node, "Eukaryota", "Prokaryota")
calc_and_add(metazoa_node, non_metazoa_euk_node, "Metazoa", "NonMetazoa_Euk")
calc_and_add(vert_node, non_vert_node, "Vertebrata", "Invertebrate")

key_df <- do.call(rbind, res)
write.csv(key_df, "key_nodes_ancestral_Z.csv", row.names = FALSE)

# ---------- 10. 输出结论 ----------
cat("\n===== 三大演化过渡净增量分析 =====\n")
print(key_df, row.names = FALSE)

cat("\n===== 演化结论 =====\n")
for (i in 1:nrow(key_df)) {
  trans <- key_df$Transition[i]
  dz <- key_df$Delta_Z[i]
  dl <- key_df$Delta_Lower[i]
  du <- key_df$Delta_Upper[i]
  if (is.na(dz)) {
    cat(trans, ": 数据不足\n")
  } else if (dl > 0) {
    cat(sprintf("%s: ΔZ = %.2f [%.2f, %.2f] → 净增加显著，支持独立增强\n", trans, dz, dl, du))
  } else if (du < 0) {
    cat(sprintf("%s: ΔZ = %.2f [%.2f, %.2f] → 净减少显著\n", trans, dz, dl, du))
  } else {
    cat(sprintf("%s: ΔZ = %.2f [%.2f, %.2f] → 置信区间包含0，未显著偏离继承\n", trans, dz, dl, du))
  }
}
cat("若ΔZ的95%置信区间完全大于0，表明该演化过渡是一次独立的i-Motif富集事件，而非单纯继承。\n")
