#!/usr/bin/env python3
"""
venn_combined.py
Draw two Venn diagrams side by side in one figure, consistent style.
Add bold uppercase letters A / B above the outside of the subpanels,
Keep the subpanel titles, placed one line below the labels (on a different line from the labels).
Data are based on actual run results.
"""

import matplotlib.pyplot as plt

# ========== Data 1: iMab CUT&Tag vs ATAC (HEK293T) ==========
A1_total = 8046       # iM promoter peak total
B1_total = 35053      # ATAC peak total
both1 = 7825          # intersection

only_A1 = A1_total - both1   # 221
only_B1 = B1_total - both1   # 27228

# ========== Data 2: Predicted iM vs ATAC (T2T, score ≥ 0.5) ==========
A2_total = 22660      # iM promoter peak total
B2_total = 50883      # ATAC peak total
both2 = 14937         # intersection

only_A2 = A2_total - both2   # 7723
only_B2 = B2_total - both2   # 35946

# ========== Plotting function (draw on specified axes) ==========
def draw_venn_on_ax(ax, only_A, only_B, both, panel_label, title):
    r = 1.0
    cx1, cy1 = -0.5, 0.0
    cx2, cy2 = 0.5, 0.0

    # Draw the two circles
    circle1 = plt.Circle((cx1, cy1), r, color='#E69F00', alpha=0.4)
    circle2 = plt.Circle((cx2, cy2), r, color='#56B4E9', alpha=0.4)
    ax.add_patch(circle1)
    ax.add_patch(circle2)

    # Set labels
    ax.text(cx1, cy1 - r - 0.2, 'iM near TSS (±1kb)', ha='center', va='top', fontsize=11)
    ax.text(cx2, cy2 - r - 0.2, 'ATAC', ha='center', va='top', fontsize=11)

    # Number labels
    ax.text(cx1 - 0.55, cy1, f'{only_A}', ha='center', va='center', fontsize=13, fontweight='bold')
    ax.text(cx2 + 0.55, cy2, f'{only_B}', ha='center', va='center', fontsize=13, fontweight='bold')
    ax.text((cx1 + cx2) / 2, cy1 - 0.05, f'{both}', ha='center', va='center', fontsize=13, fontweight='bold')

    ax.set_xlim(-2, 2)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.axis('off')

    # Subpanel title (placed above the outside of the panel, on a different line from the A/B labels)
    # Note: set_title still works after axis('off'), so we use ax.set_title here
    #       and increase pad so the title sits one line below the A/B labels
    ax.set_title(title, fontsize=13, pad=18)

    # Add bold letters A / B at the top outside the subpanel (above the title)
    ax.text(-0.05, 1.14, panel_label,
            transform=ax.transAxes,
            fontsize=18, fontweight='bold',
            ha='left', va='bottom', color='black')

# ========== Create figure with two subpanels side by side ==========
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Draw the two Venn diagrams, passing panel_label and title
draw_venn_on_ax(ax1, only_A1, only_B1, both1, 'A', 'iMab CUT&Tag vs ATAC (HEK293T)')
draw_venn_on_ax(ax2, only_A2, only_B2, both2, 'B', 'Predicted iM vs ATAC (T2T, score≥0.5)')

# Adjust layout and save
plt.tight_layout()
plt.subplots_adjust(top=0.88)   # top margin, space for the A / B letters and titles
plt.savefig("venn_combined.pdf", dpi=300, bbox_inches='tight')
plt.savefig("venn_combined.png", dpi=300, bbox_inches='tight')
print("Generated: venn_combined.pdf and venn_combined.png")
