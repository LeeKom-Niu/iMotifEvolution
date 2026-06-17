# ==================================================================
# Ancestral state reconstruction · Final (PNG+PDF + net increase analysis for three major evolutionary nodes)
# Strategy: use outgroup nodes instead of ancestral nodes
# ==================================================================

library(ape)
library(phytools)

# ---------- 1. Read data ----------
tree <- read.tree("species_tree.nwk")
dat  <- read.delim("species_taxid_sampled_final_400.tsv", header = TRUE, 
                   sep = "\t", stringsAsFactors = FALSE)

# ---------- 2. Name cleaning and matching ----------
tree$tip.label <- gsub(":.*", "", tree$tip.label)
rownames(dat) <- dat$organism_name_asm_name

common <- intersect(tree$tip.label, rownames(dat))
cat("Tips in tree:", length(tree$tip.label), " Species in table:", nrow(dat),
    " Common:", length(common), "\n")

missing <- setdiff(tree$tip.label, rownames(dat))
if (length(missing) > 0) {
  cat("Removing missing species from tree:", length(missing), "\n")
  tree <- drop.tip(tree, missing)
}
dat <- dat[tree$tip.label, ]

# ---------- 3. Rooting ----------
tree <- midpoint.root(tree)
cat("Tree rooted, total nodes:", tree$Nnode + Ntip(tree), "\n")

# ---------- 4. Build trait vector ----------
trait <- setNames(as.numeric(dat[tree$tip.label, "robust_z_score"]), tree$tip.label)
if (any(is.na(trait))) {
  na_tips <- names(trait)[is.na(trait)]
  cat("Removing", length(na_tips), "species with NA values\n")
  tree <- drop.tip(tree, na_tips)
  trait <- trait[!is.na(trait)]
}

# ---------- 5. Ancestral state reconstruction ----------
fit <- fastAnc(tree, trait, vars = TRUE, CI = TRUE)
ace     <- fit$ace
ace_CI95 <- fit$CI95
cat("Ancestral state reconstruction completed\n")

# ---------- 6. Visualization ----------
obj <- contMap(tree, trait, plot = FALSE)
obj <- setMap(obj, colors = colorRampPalette(c("blue", "white", "red"))(100))

pdf("ancestral_state_reconstruction.pdf", width = 18, height = 30)
plot(obj, fsize = 0.2, lwd = 1.5, outline = FALSE)
dev.off()

png("ancestral_state_reconstruction.png", width = 18, height = 30, units = "in", res = 300)
plot(obj, fsize = 0.2, lwd = 1.5, outline = FALSE)
dev.off()
cat("Figures saved as PDF and PNG\n")

# ---------- 7. Helper functions ----------
# Get MRCA node number and ancestral value for a specified group
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

# Calculate net increase and confidence interval between two nodes
calc_delta <- function(target, reference) {
  if (is.null(target) || is.null(reference)) return(NULL)
  dz <- target$z - reference$z
  d_low <- target$ci_low - reference$ci_up
  d_up  <- target$ci_up - reference$ci_low
  c(dz, d_low, d_up)
}

# ---------- 8. Define nodes ----------
# 1) Eukaryotic crown group (all eukaryotes)
euk_node <- get_node_val(c("Fungi","Plant","Invertebrate","Vertebrate Other",
                           "Mammalian","Protozoa"), "Eukaryota_crown")
# Prokaryotic crown group
prok_node <- get_node_val(c("Archaea","Bacteria"), "Prokaryota_crown")

# 2) Animal crown group (Metazoa)
metazoa_node <- get_node_val(c("Invertebrate","Vertebrate Other","Mammalian"),
                             "Metazoa_crown")
# Non-animal eukaryotes (all eukaryotes except animals: Fungi + Plant + Protozoa, simplified; adjust if needed)
non_metazoa_euk_node <- get_node_val(c("Fungi","Plant","Protozoa"),
                                     "NonMetazoa_Eukaryota_crown")

# 3) Vertebrate crown group
vert_node <- get_node_val(c("Vertebrate Other","Mammalian"),
                          "Vertebrata_crown")
# Non-vertebrate deuterostomes/invertebrates (as outgroup)
# The tree has Invertebrate group; select Invertebrate + other non-vertebrates
non_vert_node <- get_node_val(c("Invertebrate"), "Invertebrate_crown")

# If no independent Invertebrate node exists, use the ancestral node of non-vertebrate + vertebrate. Here Invertebrate is available.

# ---------- 9. Calculate net increase ----------
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

# ---------- 10. Output conclusions ----------
cat("\n===== Net increase analysis for three major evolutionary transitions =====\n")
print(key_df, row.names = FALSE)

cat("\n===== Evolutionary Conclusions =====\n")
for (i in 1:nrow(key_df)) {
  trans <- key_df$Transition[i]
  dz <- key_df$Delta_Z[i]
  dl <- key_df$Delta_Lower[i]
  du <- key_df$Delta_Upper[i]
  if (is.na(dz)) {
    cat(trans, ": Insufficient data\n")
  } else if (dl > 0) {
    cat(sprintf("%s: ΔZ = %.2f [%.2f, %.2f] → Significant net increase, supports independent enrichment\n", trans, dz, dl, du))
  } else if (du < 0) {
    cat(sprintf("%s: ΔZ = %.2f [%.2f, %.2f] → Significant net decrease\n", trans, dz, dl, du))
  } else {
    cat(sprintf("%s: ΔZ = %.2f [%.2f, %.2f] → Confidence interval contains 0, no significant deviation from inheritance\n", trans, dz, dl, du))
  }
}
cat("If the 95% CI of ΔZ is entirely above 0, the evolutionary transition represents an independent i-Motif enrichment event, not mere inheritance.\n")
