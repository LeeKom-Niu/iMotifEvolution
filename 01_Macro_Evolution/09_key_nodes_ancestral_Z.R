library(ggplot2)

# 数据（基于 400 物种分析，全数据跑完后替换数值即可）
df <- data.frame(
  Transition = c("Eukaryota - Prokaryota", "Metazoa - NonMetazoa_Euk", "Vertebrata - Invertebrate"),
  Delta_Z = c(10.36, 2.05, 4.43),
  Lower = c(6.99, -1.45, 0.88),
  Upper = c(13.73, 5.55, 7.99)
)
df$Transition <- factor(df$Transition, levels = rev(df$Transition))

p <- ggplot(df, aes(x = Delta_Z, y = Transition)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 4, color = "#2c7bb6") +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.3, linewidth = 1.5, 
                color = "#2c7bb6", orientation = "y") +
  labs(x = expression(Delta ~ "robust Z-score (net increase)"),
       y = "",
       title = "Ancestral state reconstruction of i-Motif enrichment") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

# 同时输出 PDF 和 PNG
ggsave("delta_z_forest.pdf", plot = p, width = 8, height = 4)
ggsave("delta_z_forest.png", plot = p, width = 8, height = 4, dpi = 300)
