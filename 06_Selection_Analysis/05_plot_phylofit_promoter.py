#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# set font
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

# read branch length data
df = pd.read_csv("branch_lengths.tsv", sep="\t")

# pivot table
pivot = df.pivot_table(index=["group", "region"], columns="part", values="hs1_length").reset_index()

# extract only the promoter region
promoter = pivot[pivot["region"] == "promoter"].set_index("group")

# compute promoter stem/loop ratio
promoter_ratio = promoter["stems"] / promoter["loops"]

# organize into a DataFrame
merged = pd.DataFrame({
    "ratio": promoter_ratio
}).reset_index()
merged["group"] = pd.Categorical(merged["group"], categories=["hominid", "homininae", "hominini"], ordered=True)
merged = merged.sort_values("group")

print("Promoter stem/loop ratios:")
print(merged)

# define colors (consistent with the HKA figure)
colors = {
    "hominid": "#E6A8B6",
    "homininae": "#7A4F2C",
    "hominini": "#66C2A5"
}

# plot
fig, ax = plt.subplots(figsize=(6, 5))

# bars (with different colors)
bars = []
for i, row in merged.iterrows():
    bar = ax.bar(row["group"], row["ratio"], color=colors[row["group"]], edgecolor="black", linewidth=0.8)
    bars.append(bar)

# add value labels
for bar, val in zip(bars, merged["ratio"]):
    # bar is a BarContainer, need to take the first rectangle
    rect = bar[0]
    ax.text(rect.get_x() + rect.get_width()/2, rect.get_height() + 0.01, f"{val:.3f}",
            ha="center", va="bottom", fontsize=10)

# reference line y=1
ax.axhline(y=1, color="gray", linestyle="dashed", linewidth=1, alpha=0.7)

# axis labels
ax.set_ylabel("Stem / Loop branch length ratio", fontsize=13)
ax.set_xlabel("Evolutionary group", fontsize=13)
ax.set_title("PhyloFit stem-loop ratios in promoters", fontsize=14, pad=15)

# polish
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="y", linestyle="--", alpha=0.5)

plt.tight_layout()
plt.savefig("Figure_phylofit_promoter.pdf", format="pdf")
plt.savefig("Figure_phylofit_promoter.png", dpi=300)
print("Figure saved as Figure_phylofit_promoter.pdf and .png")
plt.show()
