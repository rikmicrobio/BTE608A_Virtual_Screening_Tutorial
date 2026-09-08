# Virtual Screening Tutorial: Drug-Likeness & Toxicity Risk Assessment

This is a hands-on walkthrough for screening a compound library for oral
drug-likeness (Lipinski's Rule of Five) and structural toxicity alerts
(mutagenicity, tumorigenicity, reproductive toxicity, irritancy), using
OpenBabel and DataWarrior.

**Input:** `ranked-smiles-1500-1.csv`, 1,500 SMILES strings in a single `Smiles` column.
**Output:** A 3D SDF file you can load straight into DataWarrior, plus a
descriptor table and toxicity risk annotations.

---

## 1. What is SMILES notation?

Before anything else, let's cover the format all of this is built on, since
everything downstream (the .smi file, the SDF conversion, the descriptors)
starts from it.

SMILES stands for Simplified Molecular Input Line Entry System. It's a way
of writing a 3D molecular structure as a single line of plain text, using
letters, numbers, and a handful of symbols. Chemists needed something that
was compact enough to store in a spreadsheet or database, but still precise
enough that software could rebuild the exact 2D structure from it. That's
what SMILES does.

Here's how to read one, piece by piece:

**Atoms** are written as their element symbols. Carbon is so common in
organic molecules that it doesn't even need brackets: `C` is a carbon atom.
Other atoms usually don't need brackets either if they have a "normal"
number of bonds and no charge, for example `O` for oxygen or `N` for
nitrogen. Atoms that need extra information (a charge, an isotope, an exact
hydrogen count) get wrapped in square brackets, like `[NH2+]` or `[13C]`.

**Bonds** between atoms are shown by just writing the atoms next to each
other. `CC` is two carbons joined by a single bond, which is ethane. A
double bond is written with `=`, so `C=C` is ethylene. A triple bond uses
`#`, so `C#C` is acetylene. If you don't write a bond symbol at all, SMILES
assumes a single bond.

**Branches** are written in parentheses. `CC(C)C` means a carbon chain where
the second carbon has a branch coming off it, in this case that's
isobutane. You can nest branches inside branches too.

**Rings** are shown with matching numbers. Each atom that opens or closes a
ring gets a number right after it, and SMILES connects the two atoms sharing
that number with a bond. `C1CCCCC1` is a six-membered ring of carbons
(cyclohexane): the `1` after the first `C` and the `1` after the last `C`
tell the software to close the ring there.

**Aromaticity** (rings like benzene, where the bonding is delocalized) is
usually written with lowercase letters instead of uppercase. Benzene is
`c1ccccc1`, all lowercase, versus `C1CCCCC1` for the fully saturated
cyclohexane ring above. This distinction matters a lot for descriptor
calculations and toxicity fragment matching, so it's worth training your eye
to spot it.

**Stereochemistry** shows up as `/` and `\` around double bonds (to mark
cis/trans, also called E/Z), and `@` or `@@` on an atom in brackets (to mark
a stereocenter's 3D handedness, R or S). This is exactly what triggered the
CreateCisTrans warning discussed in Step 2, so it's worth knowing these
symbols exist even before you get there.

**Charges** go inside the brackets too, as `+` or `-`. `[O-]` is a
negatively charged oxygen, `[NH4+]` is a positively charged ammonium.

Let's put it together with a real molecule from this library. Row 1 of the
CSV is:
```
OB(c1ccccc1)O
```
Reading left to right: `O` is a hydroxyl oxygen, `B` is boron, `(c1ccccc1)`
is a branch containing an aromatic six membered ring (a benzene ring), and
the final `O` is a second hydroxyl oxygen. Put together, this is phenylboronic
acid, a benzene ring with a boronic acid group `B(OH)2` attached. This kind
of building block shows up constantly in Suzuki coupling chemistry, which is
a strong hint that this library is a boronic acid fragment set, matching
what we saw when every compound sailed through RO5 with room to spare (small
fragments, low molecular weight).

A couple of things worth knowing as you work with SMILES:

- **The same molecule can have more than one valid SMILES string.** Atom
  order and starting point aren't fixed, so `CCO` and `OCC` both describe
  ethanol. Software usually generates one "canonical" SMILES per molecule so
  you always get a consistent string back, but don't be surprised if a tool
  rewrites your input into a different-looking (but chemically identical)
  string.
- **SMILES describes 2D connectivity, not 3D shape.** That's exactly why
  Step 2 of this pipeline exists: OpenBabel's job is to take the flat,
  text-based connectivity in a SMILES string and build an actual 3D
  conformer out of it.
- **Not every string of letters is valid SMILES.** If a bracket, ring
  number, or branch parenthesis doesn't match up, the parser will fail on
  that line. This is one of the first things to check if a conversion tool
  reports fewer output molecules than you fed it.

If you want to explore this interactively, DataWarrior lets you type or
paste a SMILES string directly into a structure field and will draw it for
you immediately, which is a good way to build intuition before diving into
the full library.

---

## 2. Background: what are we actually screening for?

Before jumping into commands, it's worth knowing what each parameter means
and why chemists care about it.

### 2.1 Lipinski's Rule of Five (RO5)

Christopher Lipinski (Pfizer) proposed this back in 1997 as a quick heuristic
for flagging molecules that are likely to have poor oral absorption or
permeability, just from their basic physicochemical properties. A compound is
considered "drug-like" if it violates no more than one of the following:

| Property | Threshold | Why it matters |
|---|---|---|
| Molecular Weight (MW) | ≤ 500 Da | Larger molecules diffuse poorly across membranes |
| cLogP (calculated octanol-water partition coefficient) | ≤ 5 | Too lipophilic means poor solubility and higher metabolism/toxicity risk |
| H-Bond Donors (HBD) | ≤ 5 (sum of OH + NH) | Too many donors hinder passive membrane diffusion |
| H-Bond Acceptors (HBA) | ≤ 10 (sum of N + O) | Same idea, excess polarity impairs permeability |
| TPSA (Topological Polar Surface Area) | 50 to 140 Å² | Below ~50 tends to favor CNS/BBB penetration; above ~140 usually means poor oral absorption. This isn't one of Lipinski's original four, but it's the most common add-on filter |

Rotatable bonds (≤ 10) is another common add-on that captures molecular
flexibility, and it's included in the output table too.

RO5 is a filter, not a verdict. It's meant to flag candidates for a closer
look, not disqualify them outright. Plenty of approved oral drugs (and pretty
much all biologics, natural products, and PROTACs) break this rule and work
just fine.

### 2.2 Toxicity risk parameters (DataWarrior)

DataWarrior's built-in Toxicity Risk Predictor checks a molecule's
substructures against a curated library of known toxicophores (fragments
originally pulled from RTECS and expanded on since). Each compound gets
flagged High, Low, or None for four categories:

**Mutagenicity / Carcinogenicity**
Structural fragments statistically associated with the potential to cause
genetic mutations that can lead to cancer. Think aromatic nitro groups,
epoxides, certain aromatic amines, alkylating groups.

**Tumorigenicity**
Fragments linked to tumor formation that don't necessarily work through a
mutagenic mechanism, for example by promoting uncontrolled cell
proliferation. It's treated as a separate risk class from carcinogenicity
because not every tumor-promoting agent damages DNA directly.

**Reproductive effect**
Substructures linked to developmental or reproductive toxicity: effects on
fertility, teratogenicity (birth defects), or embryo-fetal development.

**Irritant**
Fragments associated with irritation of skin, eyes, or mucous membranes,
things like reactive electrophiles, strong acid/base motifs, or aldehydes.

Keep in mind these are fragment-based statistical alerts, not experimental
data and not a mechanistic prediction. They're a fast way to triage a library
before spending time and money on actual toxicology testing, nothing more.

---

## 3. Pipeline overview

```
<img width="1240" height="848" alt="Image" src="https://github.com/user-attachments/assets/4944006a-5b91-421c-9296-3d5fdbd39d45" />
```

---

## 4. Installation

### 4.1 OpenBabel (command line tool for converting chemical file formats)

Ubuntu / Debian / WSL:
```bash
sudo apt-get update
sudo apt-get install -y openbabel
obabel -V
```
This should print something like "Open Babel 3.x.x" if it worked.

macOS (Homebrew):
```bash
brew install open-babel
```

Windows: grab the installer from
openbabel.org/wiki/Category:Installation and make sure `obabel.exe` ends up
on your PATH.

Conda, if that's already your setup:
```bash
conda install -c conda-forge openbabel
```

### 4.2 DataWarrior (the GUI app, free to use, no license needed)

1. Go to openmolecules.org/datawarrior/download.html
2. Download the installer for your OS. It's a Java app, so it runs anywhere Java is installed.
3. Install it and launch it. No account, no license key.

### 4.3 Python + RDKit (optional, if you want to script the descriptors yourself)

```bash
pip install rdkit pandas --break-system-packages
```

---

## 5. Step by step

### Step 1: Build the .smi file from your CSV

Both DataWarrior and OpenBabel read plain .smi files, one line per molecule,
SMILES and ID separated by a tab.

```python
import csv

with open('ranked-smiles-1500-1.csv') as f:
    reader = csv.reader(f)
    next(reader)  # skip header row "Smiles"
    rows = [row[0].strip() for row in reader if row and row[0].strip()]

with open('compounds.smi', 'w') as out:
    for i, smi in enumerate(rows, start=1):
        out.write(f"{smi}\tCPD{i:04d}\n")

print(f"Wrote {len(rows)} molecules to compounds.smi")
```

Run it:
```bash
python3 make_smi.py
wc -l compounds.smi
```
That last line should read 1500.

Give every compound a stable ID (CPD0001, CPD0002, and so on) before you
convert formats. Some tools will drop or renumber names on their own, and
then you lose the ability to trace a result back to its row in the original
CSV.

### Step 2: Convert SMILES to 3D SDF with OpenBabel

DataWarrior can technically read a .smi file directly, but generating 3D
coordinates up front is good practice. It lets DataWarrior compute
3D-dependent descriptors later, and SDF is the more universal exchange format
in cheminformatics.

```bash
obabel compounds.smi -O compounds_3D.sdf --gen3d -h
```

What the flags do:
- `-O compounds_3D.sdf` sets the output file (OpenBabel figures out the format from the extension)
- `--gen3d` builds 3D coordinates for each molecule
- `-h` adds explicit hydrogens, which you need for correct valence and descriptor calculations later

You should see something like:
```
1500 molecules converted
```

**About the CreateCisTrans warning:**

If you're seeing a wall of this in your terminal:
```
*** Open Babel Warning  in CreateCisTrans
  Error in cis/trans stereochemistry specified for the double bond
```

This is OpenBabel telling you that a handful of SMILES in your input have a
double bond stereo notation (the `/` and `\` slashes) that it can't resolve
cleanly, usually because the two substituents needed to define "cis" or
"trans" aren't both present around that bond, or the notation is written
ambiguously. OpenBabel just drops the stereo assignment for that bond and
keeps going. It's a warning, not an error that stops the run, and it doesn't
affect the 2D connectivity or the standard descriptors used for RO5 and the
DataWarrior toxicity screens. Your molecule count should still match your
input count.

A few ways to handle it:

1. **Ignore it.** If you're only doing RO5 and toxicity alert screening
   (both are 2D/topology based), this warning has no effect on your results.
   It's just noise in the terminal.
2. **Find which molecules triggered it,** if you want to double check them:
   ```bash
   obabel compounds.smi -O compounds_3D.sdf --gen3d -h --log 2> obabel_log.txt
   ```
   Then look through `obabel_log.txt` alongside your `.smi` file to spot the
   problem entries (they'll usually be molecules with `/` or `\` in their
   SMILES).
3. **Quiet the terminal output** without losing the conversion, using `-q`:
   ```bash
   obabel compounds.smi -O compounds_3D.sdf --gen3d -h -q
   ```
4. **Fix the source SMILES** if it matters for your project, by re-drawing
   the affected double bond stereochemistry in a tool like DataWarrior or
   RDKit's SMILES sanitizer, then re-export a corrected SMILES string.

Quick sanity check once it's done:
```bash
grep -c '\$\$\$\$' compounds_3D.sdf
```
Each molecule block in an SDF ends with `$$$$`, so this count should also
read 1500.

### Step 3: Load it into DataWarrior

1. Open DataWarrior.
2. File, then Open, and pick `compounds_3D.sdf`.
3. DataWarrior will auto-detect the structure column and show 2D depictions
   in its spreadsheet grid, with CompoundID as the row label.

### Step 4: Calculate Lipinski RO5 properties

1. Right-click the structure column header, choose New Column with Chemical
   Descriptor, and add cLogP. Repeat for Molecular Weight, and use
   Chemistry, then Calculate, then Molecule Properties for H-Donors,
   H-Acceptors, Rotatable Bonds, and TPSA (these are all preset options,
   just check the boxes you need).
2. You can also use Chemistry, then Filter, then Druglikeness Filter, or
   just build your own formula column combining the four core RO5 values.
3. Add a calculated violation count via Data, then New Column with Formula:
   `RO5 Violations = (MW>500) + (cLogP>5) + (HBD>5) + (HBA>10)`

### Step 5: Run the toxicity risk screen

1. With your structure column selected, go to Chemistry, then Predict
   Toxicity Risks.
2. DataWarrior adds four columns, each valued None, Low, or High:
   - Mutagenic
   - Tumorigenic
   - Irritant
   - Reproductive Effective (that's DataWarrior's exact label for reproductive/developmental toxicity)
3. Cells are color coded, green for none, yellow for low or borderline, red
   for high, so you can scan a large table quickly.
4. Click on a flagged cell to see which fragment triggered the alert.

### Step 6: Filter and export your shortlist

A reasonable "clean" filter for a screening funnel looks like this:
```
RO5 Violations <= 1
AND Mutagenic = None
AND Tumorigenic = None
AND Reproductive Effective = None
AND Irritant = None
```
Apply it through Data, then New Row List, filtering by column values, and
then File, Save Special, Save Visible Rows As, to export your shortlisted
compounds as a new SDF or CSV.

---

## 6. Files in this repo

| File | Description |
|---|---|
| `ranked-smiles-1500-1.csv` | Original input, 1,500 SMILES |
| `compounds.smi` | SMILES plus ID, tab separated (Step 1 output) |
| `compounds_3D.sdf` | 3D structures with explicit hydrogens, ready for DataWarrior (Step 2 output) |
| `ro5_results.csv` | RDKit computed MW, cLogP, HBD, HBA, rotatable bonds, TPSA, and RO5 pass/fail per compound, useful as a cross check against DataWarrior's own numbers |
| `README.md` | This file |

---

## 7. Notes and good practice

Toxicity alerts are structural, not mechanistic. A "High" flag means the
molecule contains a fragment that shows up disproportionately often in known
toxic compounds. It's not a result from an assay or a physics based model.
Treat hits as candidates worth a closer look, not a final verdict.

RO5 is meant for oral small molecules. It doesn't really apply to
injectables, PROTACs, macrocycles, or biologics, so don't over-filter a
library that's intentionally exploring beyond Rule of Five chemical space.

The 3D structures generated here are a single low energy conformer, not a
full conformer ensemble. That's fine for 2D topology based descriptors
(which is all RO5 and the DataWarrior toxicity alerts actually use), but not
enough for anything conformer sensitive like docking, without further
refinement, for example `obminimize` or generating multiple conformers.

Coordinate generation in obabel --gen3d isn't deterministic by default. If
you need reproducible 3D coordinates run to run, pin a seed:
```bash
obabel compounds.smi -O compounds_3D.sdf --gen3d -h --seed 42
```

---

## 8. Quick reference, all commands

```bash
# 1. Install
sudo apt-get install -y openbabel

# 2. Build .smi from CSV (see Step 1 script above)
python3 make_smi.py

# 3. Convert to 3D SDF for DataWarrior
obabel compounds.smi -O compounds_3D.sdf --gen3d -h

# (optional) quiet the cis/trans warnings if they're cluttering your terminal
obabel compounds.smi -O compounds_3D.sdf --gen3d -h -q

# 4. Verify the molecule count
grep -c '\$\$\$\$' compounds_3D.sdf

# 5. In DataWarrior: File > Open > compounds_3D.sdf
#    then Chemistry > Predict Toxicity Risks
#    then Chemistry > Calculate > Molecule Properties (MW, cLogP, HBD, HBA, RotB, TPSA)
```
