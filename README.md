# Breast cancer single-cell signalling: mean vs tail

Does the **fraction of cells pushed into a very high signalling state** (the "tail")
explain drug sensitivity across breast cancer cell lines better than the
**population mean**?

Data: the Tognetti/Bodenmiller mass cytometry dataset behind the
*Single Cell Signaling in Breast Cancer DREAM Challenge*
([Synapse syn20366914](https://www.synapse.org/Synapse:syn20366914);
Tognetti et al., *Cell Systems* 2021; Gabor et al., *Mol Syst Biol* 2021, e10402).
Drug sensitivity: [GDSC](https://www.cancerrxgene.org) release 8.5 (GDSC2, GDSC1 as fallback), AUC.

## What was done

1. **Download** (`scripts/fetch.py`): the 44 cell lines in
   `single_cell_phospo/complete_cell_lines` (syn20613593). Each file is reduced to
   the columns used here (pERK, pMEK, p90RSK, pAKT S473/T308, pS6, pS6K, p4EBP1,
   cleaved caspase-3, Ki-67, pRB, GAPDH) and saved as Parquet. Values are already
   arcsinh-transformed in the challenge files.
2. **Features per line** (`scripts/analyse.py`), EGF-only condition:
   - peak time = EGF time point with the highest median pERK
   - at the peak, for pERK, pAKT (S473) and pS6:
     - `mean`, `induction` (peak mean minus time-0 mean)
     - `tail` = fraction of cells above the pooled 95th percentile (pooled with equal cells per line)
     - `upper_spread` = 95th percentile minus median (location-free tail width)
     - `tail_resid` = tail left over after regressing it on the mean (quadratic)
   - cross-talk: co-activation of pS6/pAKT within high-pERK cells, and within-cell correlations
   - replicate check: the two time-0 EGF files per line, between-line SD / within-replicate SD
3. **Drug link**: Spearman correlation of each feature with GDSC AUC for 24
   pre-specified MEK, ERK, EGFR/HER2, PI3K, AKT and mTOR drugs (34 lines matched).
   Benjamini–Hochberg FDR across all 360 tests; bootstrap CI for
   |ρ(tail)| − |ρ(mean)|; partial correlation of tail given mean.
4. **Figures** (`scripts/plots.py`).

## First results (Nov 2026)

- Tail and mean are almost the same measurement across lines (pERK ρ = 0.87, fig3).
- But lines with equal means can have very different tails, e.g. HCC1187 vs MDA-MB-157:
  1.5 % vs 16 % of cells above the threshold (fig2).
- The tail is a reproducible line property: between-line / within-replicate ratio 4.8 (tail) vs 3.4 (mean) at baseline.
- **pERK vs MEK/ERK inhibitors: the tail never beat the mean** (bootstrap CI for the
  difference includes 0 for all drugs). Best: PD0325901, ρ ≈ −0.43 for both.
- **Only FDR-surviving hit:** pAKT `tail_resid` vs MK-2206 (AKT inhibitor) AUC,
  ρ = −0.60, q = 0.07. Lines with more high-pAKT cells than their mean predicts are
  more sensitive; the same direction holds across other AKT/PI3K/mTOR drugs (fig4).
  Exploratory: 23/360 tests had p < 0.05 (~18 expected by chance).

### Caveats
- n = 34 lines with drug data: limited power.
- Not yet adjusted for subtype (luminal / HER2+ / basal).
- Some lines barely respond to EGF, so their "peak" time is noisy.
- Within-cell correlations between phospho markers are partly driven by cell size/staining.

## Next steps
- Add subtype and test the pAKT finding within subtypes.
- Use inhibitor conditions (does the high-pERK tail persist under MEK inhibition?).
- Join DepMap (CRISPR / PRISM) by cell line name.

## Rerunning

Needs Python with `pandas numpy scipy matplotlib pyarrow openpyxl`, and a Synapse
account with access to the challenge data. Scripts take one argument: a working
folder `W` that holds `reduced/`, `gdsc/` and gets `results/`.

```bash
python scripts/fetch.py W          # downloads to W/data, writes W/reduced/*.parquet
# GDSC: download GDSC1/GDSC2_fitted_dose_response_27Oct23.xlsx from
# https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/ into W/gdsc/ and combine
# them into W/gdsc/gdsc.parquet (pandas: read_excel both, concat, to_parquet)
python scripts/analyse.py W
python scripts/plots.py W
```

`fetch.py` uses plain `curl` against the Synapse REST API and expects the
Synapse access token to be supplied as an `Authorization: Bearer` header
(in the Claude cloud environment this is injected by the network proxy).
Raw data are not stored in this repository.

## p38 / GSK3 / β-catenin (added Nov 2026)

Scripts: `fetch2.py` (re-download keeping p-p38, p-MAPKAPK2, p-MKK3/6, p-GSK3β,
β-catenin, p-JNK, p-CREB → `W/reduced2/`), `plot_drugs.py`, `analyse_p38_wnt.py`.

- **Drugs (fig5):** p38 inhibitors (doramapimod, TAK-715, VX-702) and GSK3 inhibitors
  (CHIR-99021, SB216763) barely affect breast lines in GDSC (median AUC 0.95–0.99;
  IC50 reached in ≤ 8 % of lines), so there is no response to correlate. No feature
  correlated (all p > 0.05). SB202190 is not in GDSC.
- **Dynamics (fig6):** EGF raises p-p38 (sustained, peaks ~10–20 min) and p-MAPKAPK2,
  modestly p-JNK and p-GSK3β; **β-catenin does not change** within 60 min.
- **Cross-talk (fig7):** within the same cell at the pERK peak, pERK co-varies with p-AKT
  (median ρ 0.54) and p-p38 / p-MAPKAPK2 (~0.45), less with p-JNK (0.20) and hardly with
  β-catenin (0.09). GAPDH-adjusted, but caution: every phospho marker co-varies with
  cell size/staining (pERK vs GAPDH ρ ≈ 0.5), so absolute values are inflated; the
  ranking is more informative.
- **Heterogeneity (fig8):** as for pERK, the high tail of p-p38, p-MAPKAPK2 and
  β-catenin mostly tracks the mean across lines (ρ 0.77–0.95).
- **GSK3β → β-catenin (fig9):** weak positive within-cell coupling (median ρ 0.12), only
  HCC1599 negative. ~16 % of cells have no detectable β-catenin (excluded).
- Time points 12, 14 and 35 min exist for only 1–2 lines and are dropped.

## CRC test of the p38α "inhibitor of the inhibitor" hypothesis (added Nov 2026)

Hypothesis: p38α (MAPK14) restrains MAPK/ERK; removing it in MAPK-high (KRAS/BRAF-mutant)
colorectal cancer pushes signalling too high → toxicity. Script: `crc_p38_test.py`.

Data: Sanger Cell Model Passports (mutations, MSI), combined Sanger+Broad CRISPR screens
(`Project_score_combined_Sanger_v1_Broad_20Q2_20210311`, Pacini et al. 2021; batch-corrected
log fold change, QC-passed screens, duplicates averaged), GDSC 8.5. 39 CRC lines with CRISPR
(12 RAS/RAF WT, 20 KRAS/NRAS hotspot, 7 BRAF V600E). DepMap's own portal blocks automated
downloads, so the Sanger-hosted combined release was used.

- **CRISPR (fig10):** MAPK14 knockout is mildly more deleterious in **BRAF V600E** CRC
  (median LFC −0.11 vs +0.08 in WT; p = 0.01; p = 0.043 adjusting for MSI; holds when any one
  line is dropped). No difference for KRAS/NRAS-mutant lines. Across CRC lines, MAPK14
  dependency tracks **BRAF dependency** (ρ = 0.51, p = 0.001; also within non-BRAF lines, ρ = 0.42).
  MAP2K3 (MKK3, upstream of p38) shows the same direction (p = 0.014).
  **Effect size is small** (KRAS knockout in KRAS-mutant lines: −2.3).
- Positive controls behave as expected (KRAS in KRAS-mutant, BRAF and MAPK1 in BRAF-mutant).
  The ERK phosphatase **DUSP6** (a direct "ERK brake") is not essential in any group.
- **GDSC (fig11):** p38 inhibitors barely kill CRC lines (doramapimod, VX-702: IC50 almost
  never reached; TAK-715 in 16 %). TAK-715 trends the *opposite* way (mutants less sensitive,
  p = 0.03). Trametinib and dabrafenib controls work as expected.
- Interpretation: consistent with p38α mattering more where cells depend on BRAF→ERK, but
  CRISPR fitness cannot tell "hyperactivation toxicity" from "p38α supports growth". Needs a
  pERK/DUSP readout after MAPK14 loss and a MEK-inhibitor rescue.
