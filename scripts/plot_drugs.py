"""Scatter grid: chosen drugs (GDSC AUC) vs per-line signalling features.
usage: python plot_drugs.py W features.csv out.png feat1 feat2 ...
"""
import sys, re, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

W, FEAT, OUT, cols = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
BLUE, GRID, AXIS, INK, INK2, SURF = "#2a78d6", "#e1e0d9", "#c3c2b7", "#0b0b0b", "#52514e", "#fcfcfb"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "font.size": 8, "axes.titlesize": 8.5})
norm = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper())

DRUGS = [("Doramapimod", "p38 / JNK2"), ("TAK-715", "p38α/β"), ("VX-702", "p38"),
         ("CHIR-99021", "GSK3α/β"), ("SB216763", "GSK3α/β")]
feat = pd.read_csv(FEAT, index_col=0); feat["key"] = feat.index.map(norm)
g = pd.read_parquet(f"{W}/gdsc/gdsc.parquet"); g["key"] = g.CELL_LINE_NAME.map(norm)
g = g[g.key.isin(feat.key)]

rows = []
for dn, tgt in DRUGS:
    sub = g[g.DRUG_NAME == dn]
    # use the GDSC release that covers more of our lines
    ds = sub.groupby("DATASET").key.nunique().idxmax()
    auc = sub[sub.DATASET == ds].groupby("key").AUC.mean()
    rows.append((dn, tgt, ds, auc))

fig, axs = plt.subplots(len(rows), len(cols), figsize=(2.6 * len(cols), 2.2 * len(rows)), squeeze=False)
stats = []
for i, (dn, tgt, ds, auc) in enumerate(rows):
    d = feat.join(auc.rename("AUC"), on="key").dropna(subset=["AUC"])
    for j, c in enumerate(cols):
        ax = axs[i, j]
        x, y = d[c], d["AUC"]
        rho, p = spearmanr(x, y, nan_policy="omit")
        stats.append({"drug": dn, "target": tgt, "dataset": ds, "feature": c, "n": int(x.notna().sum()), "rho": rho, "p": p})
        ax.scatter(x, y, s=18, color=BLUE, edgecolor=SURF, lw=0.8, zorder=3)
        ax.set_title(f"ρ = {rho:.2f}, p = {p:.2g}, n = {x.notna().sum()}", color=INK)
        if j == 0: ax.set_ylabel(f"{dn} ({tgt})\n{ds} AUC", color=INK, fontsize=8)
        if i == len(rows) - 1: ax.set_xlabel(c, color=INK)
fig.suptitle("p38 and GSK3 inhibitor sensitivity vs EGF-response features (lower AUC = more sensitive)", fontsize=10, fontweight="bold")
fig.tight_layout(); fig.savefig(OUT, dpi=170)
pd.DataFrame(stats).to_csv(OUT.replace(".png", ".csv"), index=False)
print(pd.DataFrame(stats).round(3).to_string())
