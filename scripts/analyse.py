"""Mean vs tail of EGF-induced signalling in 44 breast lines, vs GDSC drug sensitivity."""
import sys, glob, re, os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata

S = sys.argv[1]
OUT = f"{S}/results"; os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(0)
norm = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper())

MARK = {"pERK": "p.ERK", "pS6": "p.S6", "pAKT": "p.Akt.Ser473."}

# ---------- load EGF time courses ----------
frames = []
for f in sorted(glob.glob(f"{S}/reduced/*.parquet")):
    d = pd.read_parquet(f)
    frames.append(d[d.treatment == "EGF"])
egf = pd.concat(frames, ignore_index=True)
egf["cell_line"] = egf.cell_line.astype(str).str.strip()
lines = sorted(egf.cell_line.unique())
print("lines:", len(lines), "EGF cells:", len(egf))

# time course of median pERK
tc = egf.groupby(["cell_line", "time"])[list(MARK.values())].median().reset_index()
tc.to_csv(f"{OUT}/egf_timecourse_medians.csv", index=False)

# peak time per line = EGF time point (>0) with highest median pERK
post = tc[tc.time > 0]
peak = post.loc[post.groupby("cell_line")["p.ERK"].idxmax(), ["cell_line", "time"]].rename(columns={"time": "peak_time"})
egf = egf.merge(peak, on="cell_line")
pk = egf[egf.time == egf.peak_time]
base = egf[egf.time == 0]

# pooled thresholds: equal number of cells per line so big lines don't dominate
def pooled_q(df, col, q, n=3000):
    s = df.groupby("cell_line", group_keys=False).apply(lambda g: g.sample(min(len(g), n), random_state=0))
    return s[col].quantile(q)

thr = {k: pooled_q(pk, c, 0.95) for k, c in MARK.items()}
pd.Series(thr, name="value").to_csv(f"{OUT}/thresholds.csv")
print("pooled 95th pct thresholds at peak:", {k: round(v, 3) for k, v in thr.items()})

rows = []
for cl in lines:
    p, b = pk[pk.cell_line == cl], base[base.cell_line == cl]
    r = {"cell_line": cl, "peak_time": p.peak_time.iloc[0], "n_peak_cells": len(p)}
    for k, c in MARK.items():
        x = p[c].values
        r[f"{k}_mean"] = x.mean()
        r[f"{k}_induction"] = x.mean() - b[c].mean()
        r[f"{k}_tail"] = (x > thr[k]).mean()
        r[f"{k}_upper_spread"] = np.quantile(x, 0.95) - np.median(x)   # location-free width of upper tail
    # cross-talk: among pERK-high cells, how enriched are pS6-high / pAKT-high cells?
    hi = p["p.ERK"].values > thr["pERK"]
    for k in ("pS6", "pAKT"):
        allhi = (p[MARK[k]].values > thr[k]).mean()
        r[f"coact_{k}"] = ((p[MARK[k]].values[hi] > thr[k]).mean() / allhi) if hi.sum() >= 20 and allhi > 0 else np.nan
        r[f"cellcorr_pERK_{k}"] = spearmanr(p["p.ERK"], p[MARK[k]])[0]
    r["cleavedCas_base_mean"] = b["cleavedCas"].mean()
    rows.append(r)
feat = pd.DataFrame(rows).set_index("cell_line")

# tail residual: how much tail is left once the mean is accounted for (location-free)
for k in MARK:
    m, t = feat[f"{k}_mean"], feat[f"{k}_tail"]
    coef = np.polyfit(m, t, 2)
    feat[f"{k}_tail_resid"] = t - np.polyval(coef, m)
feat.to_csv(f"{OUT}/per_line_features.csv")

# ---------- replicate check at baseline (EGF t0 has 2 files per line) ----------
thr0 = pooled_q(base, "p.ERK", 0.95)
rep = []
for cl, g in base.groupby("cell_line"):
    fids = g.fileID.unique()
    if len(fids) < 2: continue
    a, b = [g[g.fileID == f]["p.ERK"].values for f in fids[:2]]
    rep.append({"cell_line": cl, "mean_a": a.mean(), "mean_b": b.mean(),
                "tail_a": (a > thr0).mean(), "tail_b": (b > thr0).mean()})
rep = pd.DataFrame(rep)
snr = {}
for m in ("mean", "tail"):
    between = rep[[f"{m}_a", f"{m}_b"]].mean(axis=1).std()
    within = (rep[f"{m}_a"] - rep[f"{m}_b"]).abs().mean() / np.sqrt(2)
    snr[m] = between / within
    print(f"baseline pERK {m}: between-line SD / within-replicate SD = {snr[m]:.2f}  (n={len(rep)})")
rep.to_csv(f"{OUT}/baseline_replicates.csv", index=False)

# ---------- GDSC ----------
g = pd.read_parquet(f"{S}/gdsc/gdsc.parquet")
g["key"] = g.CELL_LINE_NAME.map(norm)
DRUGS = {  # pre-specified, by pathway; GDSC2 preferred, GDSC1 fallback
    "MEK": ["Trametinib", "Selumetinib", "PD0325901", "Refametinib"],
    "ERK": ["SCH772984", "Ulixertinib", "VX-11e"],
    "EGFR/HER2": ["Gefitinib", "Erlotinib", "Afatinib", "Lapatinib", "Sapitinib"],
    "PI3K": ["Pictilisib", "Alpelisib", "Taselisib", "AZD8186"],
    "AKT": ["MK-2206", "Capivasertib", "Afuresertib", "Uprosertib"],
    "mTOR": ["AZD8055", "AZD2014", "Rapamycin", "OSI-027"],
}
resp = {}
for path, names in DRUGS.items():
    for dn in names:
        sub = g[g.DRUG_NAME == dn]
        if sub.empty: continue
        ds = "GDSC2" if (sub.DATASET == "GDSC2").any() else "GDSC1"
        sub = sub[sub.DATASET == ds]
        resp[(path, f"{dn} ({ds})")] = sub.groupby("key").AUC.mean()
resp = pd.DataFrame(resp)
resp.columns = [f"{a}|{b}" for a, b in resp.columns]
feat["key"] = feat.index.map(norm)
merged = feat.join(resp, on="key")
merged.to_csv(f"{OUT}/features_with_gdsc_auc.csv")
print("lines with GDSC data:", merged[resp.columns].notna().any(axis=1).sum())

# ---------- correlations ----------
def partial_spearman(x, y, z):
    ok = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    x, y, z = rankdata(x[ok]), rankdata(y[ok]), rankdata(z[ok])
    rx = x - np.polyval(np.polyfit(z, x, 1), z); ry = y - np.polyval(np.polyfit(z, y, 1), z)
    return np.corrcoef(rx, ry)[0, 1]

METRICS = ["mean", "induction", "tail", "upper_spread", "tail_resid"]
res = []
for col in resp.columns:
    path, drug = col.split("|")
    y = merged[col].values.astype(float)
    for k in MARK:
        for m in METRICS:
            x = merged[f"{k}_{m}"].values.astype(float)
            ok = ~np.isnan(y)
            rho, p = spearmanr(x[ok], y[ok])
            r = {"pathway": path, "drug": drug, "marker": k, "metric": m, "n": ok.sum(), "rho": rho, "p": p}
            if m == "tail":
                r["partial_rho_tail_given_mean"] = partial_spearman(x, y, merged[f"{k}_mean"].values.astype(float))
                # bootstrap: |rho(tail)| - |rho(mean)|
                xm = merged[f"{k}_mean"].values.astype(float)[ok]; xt = x[ok]; yy = y[ok]; n = len(yy)
                diffs = []
                for _ in range(2000):
                    i = rng.integers(0, n, n)
                    diffs.append(abs(spearmanr(xt[i], yy[i])[0]) - abs(spearmanr(xm[i], yy[i])[0]))
                r["boot_diff_abs_rho_tail_minus_mean_lo"], r["boot_diff_hi"] = np.nanpercentile(diffs, [2.5, 97.5])
            res.append(r)
res = pd.DataFrame(res)
# Benjamini-Hochberg across all tests
o = np.argsort(res.p.values); q = np.empty(len(o)); m_ = len(o)
q[o] = np.minimum.accumulate((res.p.values[o] * m_ / np.arange(1, m_ + 1))[::-1])[::-1]
res["q_BH"] = np.minimum(q, 1)
res.to_csv(f"{OUT}/correlations.csv", index=False)

pd.set_option("display.width", 200, "display.max_rows", 200)
print("\nmetric inter-correlation (Spearman, across lines):")
for k in MARK:
    print(k); print(feat[[f"{k}_{m}" for m in METRICS]].corr("spearman").round(2))
print("\npERK metrics vs MEK/ERK drugs:")
print(res[(res.marker == "pERK") & res.pathway.isin(["MEK", "ERK"])].round(3).to_string())
print("\nall tests with p<0.05:")
print(res[res.p < 0.05].sort_values("p").round(3).to_string())
