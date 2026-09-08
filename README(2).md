# Virtual Screening Tutorial: Drug-Likeness & Toxicity Risk Assessment

A hands-on pipeline for screening a compound library for **oral drug-likeness**
(Lipinski's Rule of Five) and **structural toxicity alerts** (mutagenicity,
tumorigenicity, reproductive/developmental toxicity, irritancy) using
**OpenBabel** and **DataWarrior**.

**Input:** `ranked-smiles-1500-1.csv` — 1,500 SMILES strings (single `Smiles` column)
**Output:** A 3D SDF file loadable in DataWarrior, RO5 descriptor table, and
DataWarrior toxicity-risk annotations.

---

## 1. Background — What are we actually screening for?

Before running commands, it helps to know what each parameter means and *why*
it matters in early-stage drug discovery.

### 1.1 Lipinski's Rule of Five (RO5)

Proposed by Christopher Lipinski (Pfizer, 1997), RO5 is a heuristic that flags
molecules likely to have **poor oral absorption or permeation** based on
simple physicochemical properties. A compound is considered "drug-like" if it
violates **no more than one** of the following:

| Property | Threshold | Why it matters |
|---|---|---|
| **Molecular Weight (MW)** | ≤ 500 Da | Larger molecules diffuse poorly across membranes |
| **cLogP** (calculated octanol-water partition coefficient) | ≤ 5 | Too lipophilic → poor solubility, high metabolism/toxicity risk |
| **H-Bond Donors (HBD)** | ≤ 5 (sum of OH + NH) | Too many donors hinder passive membrane diffusion |
| **H-Bond Acceptors (HBA)** | ≤ 10 (sum of N + O) | Same as above — excess polarity impairs permeability |

Some later extensions add **Rotatable Bonds ≤ 10** and **TPSA (Topological
Polar Surface Area) ≤ 140 Å²** as secondary flexibility/permeability filters —
both are included in this pipeline's output table for reference.

> RO5 is a *filter*, not a verdict. It flags candidates for further
> scrutiny — many approved oral drugs (and essentially all biologics, natural
> products, and PROTACs) violate it.

### 1.2 Toxicity Risk Parameters (DataWarrior)

DataWarrior's built-in **Toxicity Risk Predictor** matches a molecule's
substructures against a curated fragment library of known toxicophores
(originally derived from the Registry of Toxic Effects of Chemical Substances,
RTECS, and expanded in later literature). Each compound is flagged
**High, Low, or None** risk for four independent categories:

| Parameter | Definition |
|---|---|
| **Mutagenicity / Carcinogenicity** | Structural fragments statistically associated with the potential to cause genetic mutations that can lead to cancer (e.g., aromatic nitro groups, epoxides, certain aromatic amines, alkylating groups). |
| **Tumorigenicity** | Fragments associated with tumor formation *not necessarily via a mutagenic mechanism* (e.g., promoting uncontrolled cell proliferation) — a distinct risk class from carcinogenicity, since not all tumor-promoting agents damage DNA directly. |
| **Reproductive Effect** | Substructures linked to developmental/reproductive toxicity — effects on fertility, teratogenicity (birth defects), or embryo-fetal development. |
| **Irritant** | Fragments associated with irritation of skin, eyes, or mucous membranes (e.g., reactive electrophiles, strong acids/bases as substructural motifs, aldehydes). |

These are **fragment-based statistical alerts**, not experimental data or
mechanistic proof — they are a fast triage step to deprioritize compounds
before committing to costly *in vitro*/*in vivo* toxicology.

---

## 2. Pipeline Overview

```
ranked-smiles-1500-1.csv
        │
        ▼
  [Step 1] Extract SMILES → compounds.smi  (Python/awk)
        │
        ▼
  [Step 2] OpenBabel: SMILES → 3D SDF      (obabel --gen3d)
        │
        ▼
  [Step 3] Load compounds_3D.sdf in DataWarrior
        │
        ├──▶ [Step 4] Lipinski RO5 (built-in "Druglikeness" / calculated properties)
        │
        └──▶ [Step 5] Toxicity Risk Predictor
                (Mutagenic, Tumorigenic, Reproductive, Irritant)
        │
        ▼
  [Step 6] Filter / export annotated results
```

---

## 3. Installation

### 3.1 OpenBabel (command line, converts chemical file formats)

**Ubuntu / Debian / WSL:**
```bash
sudo apt-get update
sudo apt-get install -y openbabel
obabel -V        # verify: should print "Open Babel 3.x.x"
```

**macOS (Homebrew):**
```bash
brew install open-babel
```

**Windows:** Download the installer from
[openbabel.org/wiki/Category:Installation](http://openbabel.org/wiki/Category:Installation)
and ensure `obabel.exe` is added to your PATH.

**Conda (cross-platform, recommended if you already use conda):**
```bash
conda install -c conda-forge openbabel
```

### 3.2 DataWarrior (GUI, no license required — free for any use)

1. Go to **https://openmolecules.org/datawarrior/download.html**
2. Download the installer for your OS (Windows/macOS/Linux — it's a Java app, so it runs anywhere with Java installed).
3. Install and launch it. No account or license key needed.

### 3.3 (Optional) Python + RDKit — for scripting descriptor calculations yourself

```bash
pip install rdkit pandas --break-system-packages   # or: conda install -c conda-forge rdkit
```

---

## 4. Step-by-Step Walkthrough

### Step 1 — Prepare the `.smi` file from your CSV

DataWarrior and OpenBabel both read plain `.smi` files: one line per molecule,
`SMILES<TAB>ID`.

```python
import csv

with open('ranked-smiles-1500-1.csv') as f:
    reader = csv.reader(f)
    next(reader)  # skip header "Smiles"
    rows = [row[0].strip() for row in reader if row and row[0].strip()]

with open('compounds.smi', 'w') as out:
    for i, smi in enumerate(rows, start=1):
        out.write(f"{smi}\tCPD{i:04d}\n")

print(f"Wrote {len(rows)} molecules to compounds.smi")
```

Run it:
```bash
python3 make_smi.py
wc -l compounds.smi     # should read 1500
```

> **Tip:** Always give each compound a stable ID (`CPD0001`, `CPD0002`, ...)
> before converting formats — some tools drop or auto-renumber names
> otherwise, and you lose the ability to trace a result back to its original
> CSV row.

### Step 2 — Convert SMILES → 3D SDF with OpenBabel

DataWarrior can technically import `.smi` directly, but generating 3D
coordinates up front is good practice — it lets DataWarrior compute
3D-dependent descriptors later, and SDF is the more universal exchange format
for cheminformatics pipelines.

```bash
obabel compounds.smi -O compounds_3D.sdf --gen3d -h
```

Flag reference:
- `-O compounds_3D.sdf` — output file (format inferred from extension)
- `--gen3d` — generate 3D coordinates (OpenBabel builds a rough 3D
  conformer using its force-field-based coordinate generator)
- `-h` — add explicit hydrogens (needed for correct valence/descriptor
  calculations downstream)

Expected terminal output:
```
1500 molecules converted
```
(You may see benign warnings like `Conflicting single bond directions around
double bond` — these relate to ambiguous cis/trans stereo notation in a few
input SMILES and do not affect 2D connectivity or standard descriptor
calculations.)

**Sanity check:**
```bash
grep -c '\$\$\$\$' compounds_3D.sdf   # each molecule block ends in $$$$
# should print 1500
```

### Step 3 — Load into DataWarrior

1. Open DataWarrior.
2. `File → Open...` → select `compounds_3D.sdf`.
3. DataWarrior auto-detects the structure column and displays 2D depictions
   in the spreadsheet-like grid, with `CompoundID` as a row label.

### Step 4 — Compute Lipinski RO5 properties

1. Right-click the structure column header → **New Column with Chemical
   Descriptor... → "cLogP"**, repeat for `Molecular Weight`, and use
   `Chemistry → Calculate → Molecule Properties` for **H-Donors**,
   **H-Acceptors**, **Rotatable Bonds**, and **TPSA** (all are available as
   preset properties — check each box you need).
2. Optionally use `Chemistry → Filter → Druglikeness Filter` (or manually
   combine the four columns with a formula column) to flag violations.
3. Add a calculated column:
   `RO5 Violations = (MW>500) + (cLogP>5) + (HBD>5) + (HBA>10)`
   via `Data → New Column with Formula...`.

### Step 5 — Toxicity Risk screening

1. With the structure column selected: `Chemistry → Predict Toxicity Risks`.
2. DataWarrior adds four new columns, each valued **None / Low / High**:
   - `Mutagenic`
   - `Tumorigenic`
   - `Irritant`
   - `Reproductive Effective` (DataWarrior's exact label for reproductive/developmental toxicity)
3. Each cell is colored green (none), yellow (low/borderline), or red (high)
   for quick visual triage.
4. Hover/click a flagged cell to see which substructure fragment(s) triggered
   the alert.

### Step 6 — Filter and export your hit list

Typical "clean" filter for a screening funnel:
```
RO5 Violations ≤ 1
AND Mutagenic = None
AND Tumorigenic = None
AND Reproductive Effective = None
AND Irritant = None
```
Apply via `Data → New Row List... → filter by column values`, then
`File → Save Special... → Save Visible Rows As...` to export your
shortlisted, "clean" compounds as a new SDF or CSV.

---

## 5. Files in this Repository

| File | Description |
|---|---|
| `ranked-smiles-1500-1.csv` | Original input — 1,500 SMILES |
| `compounds.smi` | SMILES + ID, tab-separated (Step 1 output) |
| `compounds_3D.sdf` | 3D structures with explicit H, ready for DataWarrior (Step 2 output) |
| `ro5_results.csv` | RDKit-computed MW/cLogP/HBD/HBA/RotB/TPSA + RO5 pass/fail per compound (companion cross-check to DataWarrior's own calculation) |
| `README.md` | This file |

---

## 6. Notes, Caveats & Good Practice

- **Toxicity alerts are structural, not mechanistic.** A "High" flag means
  the molecule contains a fragment statistically overrepresented in known
  toxic compounds — it is *not* a prediction from an assay or a physics-based
  model. Always treat hits as candidates for confirmatory testing, not
  final verdicts.
- **RO5 is for oral small molecules only.** It doesn't apply to injectables,
  PROTACs, macrocycles, or biologics — don't over-filter a library that's
  intentionally exploring "beyond Rule of Five" (bRO5) chemical space.
- **3D generation is a single low-energy conformer**, not a full conformer
  ensemble — sufficient for 2D-topology-based descriptors (which is all RO5
  and DataWarrior's toxicity alerts use) but not for anything conformer-
  sensitive (e.g., docking) without further refinement (e.g., `--conformer`,
  or energy minimization with `obminimize`).
- **Version your reagent/library files.** Re-running `obabel --gen3d` is not
  deterministic run-to-run (coordinate generation has a randomized seed by
  default) — pin a `--seed` if reproducible 3D coordinates matter for your
  workflow: `obabel compounds.smi -O compounds_3D.sdf --gen3d -h --seed 42`.

---

## 7. Quick Reference — All Commands

```bash
# 1. Install
sudo apt-get install -y openbabel

# 2. Build .smi from CSV (see Step 1 script above)
python3 make_smi.py

# 3. Convert to 3D SDF for DataWarrior
obabel compounds.smi -O compounds_3D.sdf --gen3d -h

# 4. Verify
grep -c '\$\$\$\$' compounds_3D.sdf

# 5. Open DataWarrior → File → Open → compounds_3D.sdf
#    → Chemistry → Predict Toxicity Risks
#    → Chemistry → Calculate → Molecule Properties (MW, cLogP, HBD, HBA, RotB, TPSA)
```
