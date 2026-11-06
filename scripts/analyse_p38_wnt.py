"""p38 / GSK3-beta-catenin readouts in the EGF time course: dynamics, heterogeneity, cross-talk.
usage: python analyse_p38_wnt.py W   (W/reduced2/*.parquet -> W/results/)
"""
import sys, glob, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

W = sys.argv[1]; OUT = f"{W}/results"; os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, AQUA, GRID, AXIS, INK, INK2, SURF = "#2a78d6", "#eb6834", "#1baf7a", "#e1e0d9", "#c3c2b7", "#0b0b0b", "#52514e", "#fcfcfb"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold"})

M = {"pERK": "p.ERK", "p-p38": "p.p38", "p-MAPKAPK2": "p.MAPKAPK2", "p-MKK3/6": "p.MKK3.MKK6",
     "p-GSK3β": "p.GSK3b", "β-catenin": "b.CATENIN", "p-JNK": "p.JNK", "p-AKT": "p.Akt.Ser473.", "p-S6": "p.S6"}

egf = pd.concat([pd.read_parquet(f).query("treatment == 'EGF'") for f in sorted(glob.glob(f"{W}/reduced2/*.parquet"))], ignore_index=True)
egf["cell_line"] = egf.cell_line.astype(str).str.strip()
lines = sorted(egf.cell_line.unique())
# keep only time points measured in most lines (12, 14, 35 min exist for 1-2 lines only)
cov = egf.groupby("time").cell_line.nunique()
egf = egf[egf.time.isin(cov[cov >= 0.75 * len(lines)].index)]
BCAT_FLOOR = egf["b.CATENIN"].min() + 1e-3   # ~19% of cells have no detectable beta-catenin

# ---- 1. time courses, normalised to each line's own time 0 ----
tc = egf.groupby(["cell_line", "time"])[[*M.values(), "GAPDH"]].median().reset_index()
show = ["pERK", "p-p38", "p-MAPKAPK2", "p-MKK3/6", "p-JNK", "p-GSK3β", "β-catenin", "p-AKT"]
fig, axs = plt.subplots(2, 4, figsize=(13, 5.8), sharex=True)
for ax, k in zip(axs.flat, show):
    c = M[k]
    for cl, g in tc.groupby("cell_line"):
        ax.plot(g.time, g[c] - g[c].iloc[0], color=AXIS, lw=0.8)
    d = tc.assign(delta=tc.groupby("cell_line")[c].transform(lambda s: s - s.iloc[0])).groupby("time").delta.median()
    ax.plot(d.index, d.values, color=BLUE, lw=2.2)
    ax.axhline(0, color=INK2, lw=0.7)
    ax.set_title(k)
for ax in axs[1]: ax.set_xlabel("minutes after EGF")
for ax in axs[:, 0]: ax.set_ylabel("change in median vs t=0")
fig.suptitle("Does EGF move p38 and GSK3/β-catenin?  grey = each line, blue = median across lines", fontweight="bold")
fig.tight_layout(); fig.savefig(f"{OUT}/fig6_timecourses_p38_wnt.png", dpi=180); plt.close()

# ---- 2. per-line features at that line's pERK peak ----
post = tc[tc.time > 0]
peak = post.loc[post.groupby("cell_line")["p.ERK"].idxmax(), ["cell_line", "time"]].set_index("cell_line").time
pk = egf[egf.time == egf.cell_line.map(peak)]
base = egf[egf.time == 0]

def pooled_q(df, col, q, n=3000):
    return df.groupby("cell_line", group_keys=False).apply(lambda g: g.sample(min(len(g), n), random_state=0))[col].quantile(q)
thr = {k: pooled_q(pk, c, 0.95) for k, c in M.items()}

def resid(y, x):  # remove the part explained by GAPDH (cell size / staining)
    return y - np.polyval(np.polyfit(x, y, 1), x)

rows = []
for cl in lines:
    p, b = pk[pk.cell_line == cl], base[base.cell_line == cl]
    r = {"cell_line": cl, "peak_time": peak[cl], "n": len(p)}
    for k, c in M.items():
        r[f"{k}_mean"] = p[c].mean()
        r[f"{k}_induction"] = p[c].mean() - b[c].mean()
        r[f"{k}_tail"] = (p[c] > thr[k]).mean()
    hi = p["p.ERK"] > thr["pERK"]
    for k in ("p-p38", "p-MAPKAPK2", "p-JNK", "p-AKT", "p-S6", "β-catenin"):
        c = M[k]
        r[f"corr_pERK_{k}"] = spearmanr(p["p.ERK"], p[c])[0]
        r[f"corrGAPDHadj_pERK_{k}"] = spearmanr(resid(p["p.ERK"].values, p.GAPDH.values), resid(p[c].values, p.GAPDH.values))[0]
        allhi = (p[c] > thr[k]).mean()
        r[f"enrich_{k}_in_pERKhigh"] = (p.loc[hi, c] > thr[k]).mean() / allhi if hi.sum() >= 20 and allhi > 0 else np.nan
    q = p[p["b.CATENIN"] > BCAT_FLOOR]          # cells with detectable beta-catenin only
    r["bcat_frac_detectable"] = len(q) / len(p)
    r["corr_pGSK3b_bcat"] = spearmanr(q["p.GSK3b"], q["b.CATENIN"])[0]
    r["corrGAPDHadj_pGSK3b_bcat"] = spearmanr(resid(q["p.GSK3b"].values, q.GAPDH.values), resid(q["b.CATENIN"].values, q.GAPDH.values))[0]
    r["corrGAPDHadj_pERK_GAPDH_control"] = spearmanr(p["p.ERK"], p["GAPDH"])[0]
    rows.append(r)
F = pd.DataFrame(rows).set_index("cell_line")
F.to_csv(f"{OUT}/p38_wnt_per_line_features.csv")

# ---- 3. cross-talk: within-cell correlation of pERK with other nodes, after removing GAPDH ----
ks = ["p-AKT", "p-S6", "p-p38", "p-MAPKAPK2", "p-JNK", "β-catenin"]
fig, ax = plt.subplots(figsize=(8, 4.4))
for i, k in enumerate(ks):
    v = F[f"corrGAPDHadj_pERK_{k}"].values
    ax.scatter(np.full(len(v), i) + np.random.default_rng(i).uniform(-0.18, 0.18, len(v)), v, s=14, color=BLUE, alpha=0.75, edgecolor="none")
    ax.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color=INK, lw=2)
ax.axhline(0, color=INK2, lw=0.7)
ax.set_xticks(range(len(ks)), ks); ax.set_ylabel("within-cell Spearman ρ with pERK\n(GAPDH-adjusted)")
ax.set_title("At the pERK peak, which nodes rise in the same cells as pERK?  (each dot = one line)")
fig.tight_layout(); fig.savefig(f"{OUT}/fig7_crosstalk_with_pERK.png", dpi=180); plt.close()

# ---- 4. heterogeneity: mean vs tail across lines for p-p38, p-MAPKAPK2, β-catenin ----
fig, axs = plt.subplots(1, 3, figsize=(13, 4.2))
for ax, k in zip(axs, ["p-p38", "p-MAPKAPK2", "β-catenin"]):
    x, y = F[f"{k}_mean"], 100 * F[f"{k}_tail"]
    ax.scatter(x, y, s=30, color=BLUE, edgecolor=SURF, lw=1, zorder=3)
    res = y - np.polyval(np.polyfit(x, y, 2), x)
    for cl in res.abs().nlargest(4).index:
        ax.annotate(cl, (x[cl], y[cl]), fontsize=7.5, color=INK2, xytext=(4, 3), textcoords="offset points")
    ax.set(xlabel=f"mean {k} at pERK peak", ylabel="% cells above pooled 95th pct",
           title=f"{k}: ρ(mean, tail) = {spearmanr(x, y)[0]:.2f}")
fig.suptitle("Line-to-line heterogeneity: is the high tail just the mean?", fontweight="bold")
fig.tight_layout(); fig.savefig(f"{OUT}/fig8_mean_vs_tail_p38_bcat.png", dpi=180); plt.close()

# ---- 5. GSK3 -> beta-catenin, two example lines with the most different coupling ----
ex = [F.corrGAPDHadj_pGSK3b_bcat.idxmax(), F.corrGAPDHadj_pGSK3b_bcat.idxmin()]
fig, axs = plt.subplots(1, 3, figsize=(13, 4.2))
v = F.corrGAPDHadj_pGSK3b_bcat.sort_values()
axs[0].barh(range(len(v)), v.values, color=BLUE, height=0.7); axs[0].set_yticks(range(len(v)), v.index, fontsize=6)
axs[0].axvline(0, color=INK2, lw=0.7); axs[0].set_xlabel("within-cell ρ, p-GSK3β vs β-catenin (GAPDH-adj.)")
axs[0].set_title("Coupling per line (cells with detectable β-catenin)")
for ax, cl in zip(axs[1:], ex):
    p = pk[(pk.cell_line == cl) & (pk["b.CATENIN"] > BCAT_FLOOR)]
    ax.hexbin(p["p.GSK3b"], p["b.CATENIN"], gridsize=45, cmap="Blues", mincnt=1)
    ax.set(xlabel="p-GSK3β", ylabel="β-catenin", title=f"{cl}: ρ_adj = {F.loc[cl,'corrGAPDHadj_pGSK3b_bcat']:.2f}")
fig.tight_layout(); fig.savefig(f"{OUT}/fig9_GSK3_bcatenin.png", dpi=180); plt.close()

pd.set_option("display.width", 220)
print("peak-time induction (median across lines):")
print(F[[f"{k}_induction" for k in M]].median().round(3))
print("\nwithin-cell corr with pERK, raw vs GAPDH-adjusted (median across lines):")
print(pd.DataFrame({k: [F[f"corr_pERK_{k}"].median(), F[f"corrGAPDHadj_pERK_{k}"].median(), F[f"enrich_{k}_in_pERKhigh"].median()] for k in ks},
                   index=["raw", "GAPDH-adj", "enrichment"]).round(2))
print("control: pERK vs GAPDH raw corr, median", F.corrGAPDHadj_pERK_GAPDH_control.median().round(2))
print("\nGSK3b-bcat coupling raw/adj median:", F.corr_pGSK3b_bcat.median().round(2), F.corrGAPDHadj_pGSK3b_bcat.median().round(2), "range adj", F.corrGAPDHadj_pGSK3b_bcat.min().round(2), F.corrGAPDHadj_pGSK3b_bcat.max().round(2))
for k in ["p-p38", "p-MAPKAPK2", "β-catenin", "pERK"]:
    print(f"{k}: rho(mean,tail)={spearmanr(F[f'{k}_mean'], F[f'{k}_tail'])[0]:.2f}  tail range {F[f'{k}_tail'].min():.3f}-{F[f'{k}_tail'].max():.3f}")
print("examples:", ex, " beta-cat detectable frac median/min:", F.bcat_frac_detectable.median().round(2), F.bcat_frac_detectable.min().round(2))
print(F.peak_time.value_counts().sort_index().to_dict())
