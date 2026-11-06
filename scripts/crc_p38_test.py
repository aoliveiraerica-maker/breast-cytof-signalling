"""Hypothesis: p38alpha (MAPK14) is a brake on MAPK/ERK; removing it in MAPK-high (KRAS/BRAF-mutant)
colorectal cancer pushes ERK too high -> toxicity.

Test 1 (CRISPR): MAPK14 knockout fitness effect by RAS/RAF status, Sanger+Broad combined screens
               (Pacini et al. 2021, batch-corrected log fold change; more negative = more essential).
Test 2 (drugs):  GDSC p38 inhibitor response by RAS/RAF status.
usage: python crc_p38_test.py W     (needs W/cmp, W/depmap, W/gdsc)  -> W/results/
"""
import sys, os, re, glob, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, kruskal, spearmanr
import statsmodels.formula.api as smf

W = sys.argv[1]; OUT = f"{W}/results"; os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, AQUA, GRID, AXIS, INK, INK2, SURF = "#2a78d6", "#eb6834", "#1baf7a", "#e1e0d9", "#c3c2b7", "#0b0b0b", "#52514e", "#fcfcfb"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold"})
GROUPS = ["RAS/RAF WT", "KRAS/NRAS mut", "BRAF V600E"]
GCOL = dict(zip(GROUPS, [BLUE, ORANGE, AQUA]))

# ---------- annotation: tissue, MSI, RAS/RAF hotspot status ----------
ml = pd.read_csv(f"{W}/cmp/model_list_latest.csv.gz", low_memory=False)
mut = pd.read_csv(f"{W}/cmp/mutations_summary_latest.csv.gz", low_memory=False)
pm = mut.protein_mutation.astype(str)
ras = (mut.gene_symbol.isin(["KRAS", "NRAS"]) & pm.str.match(r"p\.(G12|G13|Q61|K117|A146)"))
braf = (mut.gene_symbol == "BRAF") & pm.str.startswith("p.V600")
ann = ml[["model_id", "model_name", "tissue", "cancer_type", "msi_status"]].drop_duplicates("model_id").set_index("model_id")
ann["has_mut_data"] = ann.index.isin(mut.model_id.unique())
ann["RAS_mut"] = ann.index.isin(mut.loc[ras, "model_id"])
ann["BRAF_mut"] = ann.index.isin(mut.loc[braf, "model_id"])
ann["group"] = np.select([ann.BRAF_mut, ann.RAS_mut], ["BRAF V600E", "KRAS/NRAS mut"], "RAS/RAF WT")
ann.loc[ann.BRAF_mut & ann.RAS_mut, "group"] = "both"
ann.loc[~ann.has_mut_data, "group"] = np.nan
ann["CRC"] = ann.cancer_type.astype(str).str.contains("Colorectal", case=False)

# ---------- Test 1: CRISPR ----------
D = glob.glob(f"{W}/depmap/Project_score_combined*/*fold_change_values*.tsv")[0]
raw = pd.read_csv(D, sep="\t", header=None, low_memory=False)
hdr = raw.iloc[:4, 1:].T; hdr.columns = ["model_name", "model_id", "source", "qc_pass"]
fc = raw.iloc[4:].set_index(0); fc.columns = range(fc.shape[1]); fc = fc.apply(pd.to_numeric, errors="coerce")
GENES = ["MAPK14", "MAPK11", "MAPKAPK2", "MAP2K3", "MAP2K6", "DUSP6", "DUSP4", "MAPK1", "KRAS", "BRAF"]
keep = hdr.qc_pass.astype(str).str.upper().eq("TRUE").values
sub = fc.loc[GENES, keep].T
sub["model_id"] = hdr.model_id.values[keep]
ko = sub.groupby("model_id").mean()            # average duplicate screens of the same line
ko = ko.join(ann, how="inner")
ko.to_csv(f"{OUT}/crispr_ko_selected_genes.csv")

crc = ko[ko.CRC & ko.group.isin(GROUPS)]
print("CRC lines with CRISPR + mutation data:", len(crc), crc.group.value_counts().to_dict())
rows = []
for gname in GENES:
    x = {g: crc.loc[crc.group == g, gname].dropna() for g in GROUPS}
    r = {"gene": gname, **{f"median_{g}": x[g].median() for g in GROUPS}, **{f"n_{g}": len(x[g]) for g in GROUPS}}
    r["p_RASmut_vs_WT"] = mannwhitneyu(x["KRAS/NRAS mut"], x["RAS/RAF WT"]).pvalue
    r["p_BRAF_vs_WT"] = mannwhitneyu(x["BRAF V600E"], x["RAS/RAF WT"]).pvalue
    r["p_anyMAPKmut_vs_WT"] = mannwhitneyu(pd.concat([x["KRAS/NRAS mut"], x["BRAF V600E"]]), x["RAS/RAF WT"]).pvalue
    rows.append(r)
crc_stats = pd.DataFrame(rows); crc_stats.to_csv(f"{OUT}/crispr_crc_group_tests.csv", index=False)

# MSI-adjusted model for MAPK14 in CRC
m = smf.ols("MAPK14 ~ C(group, Treatment('RAS/RAF WT')) + C(msi_status)", data=crc.dropna(subset=["MAPK14", "msi_status"])).fit()
print("\nMAPK14 ~ group + MSI (CRC):"); print(m.summary().tables[1])

# pan-cancer, adjusting for tissue
pan = ko[ko.group.isin(GROUPS)].dropna(subset=["MAPK14"])
pan = pan[pan.groupby("tissue").tissue.transform("size") >= 10]
mp = smf.ols("MAPK14 ~ C(group, Treatment('RAS/RAF WT')) + C(tissue)", data=pan).fit()
print(f"\nPan-cancer MAPK14 ~ group + tissue (n={len(pan)}):")
print(mp.summary().tables[1].data[:3])
# in CRC: does needing p38a go with needing the MAPK pathway?
print("\nCRC: Spearman of MAPK14 KO effect with ...")
for gname in ["MAPK1", "KRAS", "BRAF", "DUSP6", "MAPK11", "MAPKAPK2"]:
    r, p = spearmanr(crc.MAPK14, crc[gname], nan_policy="omit"); print(f"  {gname}: rho={r:.2f} p={p:.3f}")

fig, axs = plt.subplots(2, 5, figsize=(14, 6.4), sharey=False)
for ax, gname in zip(axs.flat, GENES):
    for i, g in enumerate(GROUPS):
        v = crc.loc[crc.group == g, gname].dropna().values
        ax.scatter(i + np.random.default_rng(i).uniform(-0.17, 0.17, len(v)), v, s=16, color=GCOL[g], alpha=0.85, edgecolor="none")
        ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=INK, lw=2)
    ax.axhline(0, color=INK2, lw=0.7)
    p = crc_stats.set_index("gene").loc[gname, "p_anyMAPKmut_vs_WT"]
    ax.set_title(f"{gname}  (mut vs WT p = {p:.2g})", fontsize=9)
    ax.set_xticks(range(3), ["WT", "KRAS/\nNRAS", "BRAF\nV600E"], fontsize=8)
for ax in axs[:, 0]: ax.set_ylabel("knockout effect (log fold change)\nmore negative = more essential")
fig.suptitle(f"CRC cell lines, CRISPR knockout (Sanger+Broad, n = {len(crc)}): is p38α (MAPK14) a liability in RAS/RAF-mutant lines?", fontweight="bold")
fig.tight_layout(); fig.savefig(f"{OUT}/fig10_crc_crispr_MAPK14.png", dpi=180); plt.close()

# ---------- Test 2: GDSC drugs ----------
g = pd.read_parquet(f"{W}/gdsc/gdsc.parquet")
g = g.join(ann[["group", "msi_status", "CRC"]], on="SANGER_MODEL_ID")
g = g[(g.TCGA_DESC == "COREAD") & g.group.isin(GROUPS)]
DRUGS = [("Doramapimod", "GDSC2"), ("TAK-715", "GDSC1"), ("VX-702", "GDSC1"), ("Trametinib", "GDSC2"), ("Dabrafenib", "GDSC2"), ("SB216763", "GDSC2")]
rows = []
fig, axs = plt.subplots(1, len(DRUGS), figsize=(15, 4))
for ax, (dn, ds) in zip(axs, DRUGS):
    s = g[(g.DRUG_NAME == dn) & (g.DATASET == ds)].groupby(["SANGER_MODEL_ID", "group"]).agg(AUC=("AUC", "mean"), LN_IC50=("LN_IC50", "mean"), MAX=("MAX_CONC", "max")).reset_index()
    x = {gg: s.loc[s.group == gg, "AUC"] for gg in GROUPS}
    mutv = pd.concat([x["KRAS/NRAS mut"], x["BRAF V600E"]])
    rows.append({"drug": dn, "dataset": ds, **{f"medianAUC_{gg}": x[gg].median() for gg in GROUPS}, **{f"n_{gg}": len(x[gg]) for gg in GROUPS},
                 "p_anyMAPKmut_vs_WT": mannwhitneyu(mutv, x["RAS/RAF WT"]).pvalue,
                 "frac_IC50_below_max": (s.LN_IC50 < np.log(s.MAX)).mean()})
    for i, gg in enumerate(GROUPS):
        v = x[gg].values
        ax.scatter(i + np.random.default_rng(i).uniform(-0.17, 0.17, len(v)), v, s=16, color=GCOL[gg], alpha=0.85, edgecolor="none")
        ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=INK, lw=2)
    ax.set_title(f"{dn}\n(mut vs WT p = {rows[-1]['p_anyMAPKmut_vs_WT']:.2g})", fontsize=9)
    ax.set_xticks(range(3), ["WT", "KRAS/\nNRAS", "BRAF\nV600E"], fontsize=8)
axs[0].set_ylabel("GDSC AUC (lower = more sensitive)")
fig.suptitle("CRC cell lines, GDSC: p38 inhibitors vs controls (trametinib = MEK, dabrafenib = BRAF, SB216763 = GSK3)", fontweight="bold")
fig.tight_layout(); fig.savefig(f"{OUT}/fig11_crc_gdsc_p38i.png", dpi=180); plt.close()
gd = pd.DataFrame(rows); gd.to_csv(f"{OUT}/gdsc_crc_group_tests.csv", index=False)

pd.set_option("display.width", 250, "display.max_columns", 30)
print("\nCRISPR CRC group medians / tests:"); print(crc_stats.round(3).to_string(index=False))
print("\nGDSC CRC:"); print(gd.round(3).to_string(index=False))
