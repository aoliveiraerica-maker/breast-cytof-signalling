import sys, glob, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

S = sys.argv[1]; R = f"{S}/results"
BLUE, ORANGE, RED, GRID, AXIS, INK, INK2, SURF = "#2a78d6", "#eb6834", "#e34948", "#e1e0d9", "#c3c2b7", "#0b0b0b", "#52514e", "#fcfcfb"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold"})

feat = pd.read_csv(f"{R}/per_line_features.csv", index_col=0)
tc = pd.read_csv(f"{R}/egf_timecourse_medians.csv")
corr = pd.read_csv(f"{R}/correlations.csv")

# 1. time courses
fig, ax = plt.subplots(figsize=(7, 4.2))
for cl, g in tc.groupby("cell_line"):
    ax.plot(g.time, g["p.ERK"], color=AXIS, lw=1)
mid = tc.groupby("time")["p.ERK"].median()
ax.plot(mid.index, mid.values, color=BLUE, lw=2.5, label="median across lines")
ax.set(xlabel="minutes after EGF", ylabel="median pERK per line (arcsinh)", title="EGF-induced pERK: every line (grey) peaks early, at different heights")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(f"{R}/fig1_pERK_timecourses.png", dpi=200); plt.close()

# 2. two lines, similar mean, different tail
f = feat.copy()
best = None
for a in f.index:
    for b in f.index:
        if a >= b: continue
        dm = abs(f.loc[a, "pERK_mean"] - f.loc[b, "pERK_mean"])
        if dm < 0.1:
            dt = abs(f.loc[a, "pERK_tail"] - f.loc[b, "pERK_tail"])
            if best is None or dt > best[0]: best = (dt, a, b)
_, a, b = best
peaks = {}
for cl in (a, b):
    d = pd.read_parquet(f"{S}/reduced/{cl}.parquet")
    peaks[cl] = d[(d.treatment == "EGF") & (d.time == f.loc[cl, "peak_time"])]["p.ERK"].values
thr = pd.read_csv(f"{R}/thresholds.csv", index_col=0).loc["pERK", "value"]
fig, ax = plt.subplots(figsize=(7, 4.2))
bins = np.linspace(0, 8, 80)
for cl, col in ((a, BLUE), (b, ORANGE)):
    ax.hist(peaks[cl], bins=bins, density=True, histtype="step", lw=2, color=col,
            label=f"{cl}: mean {f.loc[cl,'pERK_mean']:.2f}, {100*f.loc[cl,'pERK_tail']:.1f}% of cells above dashed line")
ax.axvline(thr, color=INK2, ls="--", lw=1); ax.text(thr + 0.05, ax.get_ylim()[1] * 0.92, "pooled 95th pct", color=INK2, fontsize=9)
ax.set(xlabel="pERK at peak (arcsinh)", ylabel="density", title="Same mean, different high-pERK tail")
ax.legend(frameon=False, fontsize=9, loc="upper left")
fig.tight_layout(); fig.savefig(f"{R}/fig2_same_mean_different_tail.png", dpi=200); plt.close()

# 3. mean vs tail
fig, ax = plt.subplots(figsize=(6, 4.6))
ax.scatter(f.pERK_mean, 100 * f.pERK_tail, s=40, color=BLUE, edgecolor=SURF, lw=1.5, zorder=3)
for cl in f.pERK_tail_resid.abs().nlargest(5).index.tolist() + [a, b]:
    ax.annotate(cl, (f.loc[cl, "pERK_mean"], 100 * f.loc[cl, "pERK_tail"]), fontsize=8, color=INK2, xytext=(4, 3), textcoords="offset points")
rho = f[["pERK_mean", "pERK_tail"]].corr("spearman").iloc[0, 1]
ax.set(xlabel="mean pERK at peak", ylabel="% cells above pooled 95th pct", title=f"Mean and tail move together (Spearman ρ = {rho:.2f})")
fig.tight_layout(); fig.savefig(f"{R}/fig3_mean_vs_tail.png", dpi=200); plt.close()

# 4. heatmap of rho
corr["row"] = corr.marker + " · " + corr.metric
corr["col"] = corr.pathway + ": " + corr.drug.str.replace(r" \(GDSC\d\)", "", regex=True)
order_rows = [f"{m} · {k}" for m in ("pERK", "pAKT", "pS6") for k in ("mean", "induction", "tail", "upper_spread", "tail_resid")]
cols = corr.drop_duplicates("col").col.tolist()
H = corr.pivot(index="row", columns="col", values="rho").loc[order_rows, cols]
P = corr.pivot(index="row", columns="col", values="q_BH").loc[order_rows, cols]
cmap = LinearSegmentedColormap.from_list("div", [BLUE, "#f0efec", RED])
fig, ax = plt.subplots(figsize=(13, 6.5))
im = ax.imshow(H.values, cmap=cmap, vmin=-0.8, vmax=0.8, aspect="auto")
ax.grid(False)
ax.set_xticks(range(len(cols)), cols, rotation=60, ha="right", fontsize=8); ax.set_yticks(range(len(order_rows)), order_rows, fontsize=8)
for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        if P.values[i, j] < 0.1: ax.text(j, i, "*", ha="center", va="center", color=INK, fontsize=12)
cb = fig.colorbar(im, ax=ax, shrink=0.7); cb.set_label("Spearman ρ with GDSC AUC\n(negative = higher signal, more sensitive)")
ax.set_title("Signalling feature vs drug AUC across lines  (* = FDR q < 0.1)")
fig.tight_layout(); fig.savefig(f"{R}/fig4_correlation_heatmap.png", dpi=200); plt.close()
print("pair:", a, b)
