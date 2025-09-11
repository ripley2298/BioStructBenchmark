# BioStructBenchmark 🧬

A comprehensive toolkit for benchmarking computational structure predictions against experimental reality, with specialized focus on **DNA-protein complexes**. Features advanced multi-frame alignment, **critical functional residue interaction analysis**, hydrogen bond network analysis with structural correspondence mapping, and X3DNA-DSSR integration for protein-DNA binding interface validation.

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Command Line Usage](#command-line-usage)
- [Critical Interaction Analysis](#critical-interaction-analysis)
- [Module Architecture](#module-architecture)
- [Analysis Capabilities](#analysis-capabilities)
- [Output Files](#output-files)
- [Best Practices](#best-practices)
- [API Usage](#api-usage)
- [Troubleshooting](#troubleshooting)
- [Testing and Development](#testing-and-development)
- [Contributing](#contributing)
- [License](#license)

## Overview

BioStructBenchmark compares experimental DNA-protein complex structures with their computationally predicted counterparts (AlphaFold3, RoseTTAFold, etc.) to identify systematic prediction errors. The toolkit provides comprehensive structural analysis through multiple reference frame alignments, **critical functional residue interaction verification**, detailed error decomposition, and extensive DNA geometry analysis.

### Why BioStructBenchmark?

- **🎯 Multi-frame alignment**: Three reference frames to distinguish positioning vs structural errors
- **🔬 Critical interaction analysis**: Verifies ≥3 essential DNA-binding residue interactions (Arg NH₃⁺→phosphate, Gln→base)
- **🧬 Advanced H-bond analysis**: Structural correspondence mapping fixes false negatives from residue numbering shifts
- **⚗️ X3DNA-DSSR integration**: 5 critical parameters for protein-DNA binding interface validation
- **🔗 Intelligent structure pairing**: Automatic matching of experimental/predicted pairs across naming conventions
- **📊 Dual RMSD reporting**: Both all-atom and backbone-only RMSD for comprehensive assessment
- **🎨 Publication-quality outputs**: Advanced visualizations, heatmaps, dashboards, and scientific reports

## Key Features

### 🎯 Multi-Frame Alignment Analysis
- **Full structure alignment**: Overall structural accuracy assessment
- **DNA-to-protein positioning**: Evaluates DNA placement relative to protein
- **DNA standalone**: Assesses DNA structure independent of protein context
- **Dual RMSD calculation**: All-atom (primary) + backbone-only (additional) for each frame

### 🔬 Critical Functional Residue Interaction Analysis ⭐
**Verifies ≥3 critical DNA-binding interactions essential for protein-DNA recognition:**

| Residue Type | Functional Atoms | Target | Threshold | Biological Role |
|--------------|------------------|--------|-----------|-----------------|
| **Arginine** | NH1, NH2 | DNA phosphate | 3.5Å | Electrostatic backbone contacts |
| **Lysine** | NZ | DNA phosphate | 3.5Å | Positive charge neutralization |
| **Glutamine** | NE2, OE1 | DNA bases | 3.2Å | Sequence-specific recognition |
| **Asparagine** | ND2, OD1 | DNA bases | 3.2Å | Hydrogen bond networks |
| **Serine** | OG | DNA phosphate | 3.2Å | Hydroxyl-phosphate contacts |
| **Threonine** | OG1 | DNA phosphate | 3.2Å | Backbone stabilization |


### 📊 Comprehensive Metrics
- **Dual RMSD calculation**: All-atom (backbone + side chains/bases) and backbone-only for all analyses
- Per-residue/nucleotide RMSD with gap-tolerant sequence alignment
- Error decomposition (translational vs rotational components)
- B-factor vs pLDDT confidence metric comparison
- Consensus error mapping across multiple structures
- Mutation impact analysis with local RMSD effects

### 🧬 Advanced Hydrogen Bond Analysis
- **Structural correspondence mapping** - Fixes critical false negatives from residue numbering differences
- Protein-DNA interface H-bond networks with geometric validation (distance ≤3.5Å, angle ≥120°)
- Conservation rate and prediction accuracy metrics relative to actual interactions
- Comprehensive H-bond classification (strong/medium/weak)
- Export detailed H-bond tables and network statistics

### ⚗️ X3DNA-DSSR Integration (5 Critical Parameters)
**Essential protein-DNA binding interface parameters:**
- **Base pairs**: Canonical base pairs at protein-binding interface
- **Helical twist**: Average twist per base pair step (affects groove dimensions)
- **Major groove width**: Primary protein-binding interface dimensions  
- **Minor groove width**: DNA bending and protein recognition effects
- **Stacking energy**: Interface stability indicator
- **Threshold alerts**: Flags deviations >5° (twist), >0.5Å (grooves), >2 kcal/mol (energy)

### 📈 Visualization
- Publication-quality plots and dashboards with dual RMSD displays
- Memory-efficient per-residue RMSD heatmaps (all-atom and backbone)
- Critical interaction conservation plots
- Correlation analysis between RMSD and secondary metrics
- DNA geometry comparison plots
- Hydrogen bond network visualization
- 3D structure alignment visualization

## Installation

### Requirements

- Python 3.8+
- BioPython ≥ 1.79
- NumPy ≥ 1.20
- Pandas ≥ 1.3
- Matplotlib ≥ 3.4
- SciPy ≥ 1.7

### Optional Dependencies

- **X3DNA-DSSR** (for critical protein-DNA binding interface analysis) - https://x3dna.org/
- **py3Dmol** (for interactive 3D visualization)
- **Seaborn** (for enhanced plotting)

### Install from Source

```bash
git clone https://github.com/yourusername/biostructbenchmark
cd biostructbenchmark
pip install -e .
```

### Verify Installation

```bash
biostructbenchmark --version
# BioStructBenchmark v0.1.0

# Run tests
pytest tests/
```

## Quick Start

### Basic Single Structure Comparison

```bash
# Simple RMSD calculation with dual reporting
biostructbenchmark observed.pdb predicted.pdb

# With output directory
biostructbenchmark -e experimental.pdb -p predicted.pdb -o results/
```

### Multi-Frame Alignment with Critical Interactions

```bash
# Comprehensive analysis with critical residue verification
biostructbenchmark -e exp.pdb -p pred.pdb --multi-frame --hbond --save-aligned

# This generates:
# 1. Full structure alignment (all-atom + backbone RMSD)
# 2. DNA positioning relative to protein
# 3. DNA standalone structural comparison
# 4. Critical functional residue interaction analysis
# 5. ≥3 critical interactions verification (PASS/FAIL)
```

### Batch Processing with Intelligent Pairing

```bash
# Process entire directories (automatic structure pairing)
biostructbenchmark -e experimental_dir/ -p predicted_dir/ -o results/ --all-benchmarks

# Intelligent pairing handles naming differences:
# experimental: p456_02_experimental.pdb ↔ predicted: p456_02_alphafold3.cif
# experimental: 1abc_exp.pdb ↔ predicted: 1abc_pred.cif

# With critical interaction analysis and visualization
biostructbenchmark -e exp/ -p pred/ -o out/ --hbond --dssr --visualize --parallel 4
```

## Command Line Usage

### Basic Syntax

```bash
biostructbenchmark [OPTIONS] [LEGACY_FILES]
```

### Primary Options

| Option | Description |
|--------|-------------|
| `-e, --experimental PATH` | Path to experimental structure(s) |
| `-p, --predicted PATH` | Path to predicted structure(s) |
| `-o, --output PATH` | Output directory (default: ./biostructbenchmark_results) |
| `-v, --version` | Show version and exit |

### Analysis Options

| Option | Description |
|--------|-------------|
| `--all-benchmarks` | Run all available analyses including critical interactions |
| `--multi-frame` | Perform multi-frame alignment (3 reference frames) |
| `--rmsd-only` | Basic RMSD analysis only (fastest) |
| `--hbond` | Hydrogen bond network analysis with structural correspondence |
| `--dssr` | X3DNA-DSSR analysis (5 critical protein-DNA binding parameters) |
| `--bfactor` | B-factor vs confidence metric analysis |
| `--consensus` | Consensus error mapping (requires multiple pairs) |
| `--mutations` | Detect and analyze mutations |
| `--visualize` | Generate publication-quality plots with dual RMSD |

### Advanced Options

| Option | Description |
|--------|-------------|
| `--reference-frame {full,protein,dna,multi}` | Alignment reference frame |
| `--output-format {csv,json,both}` | Output data format |
| `--rmsd-threshold FLOAT` | RMSD threshold for consensus (default: 3.0 Å) |
| `--save-aligned` | Save aligned PDB structures |
| `--export-all` | Export all intermediate files |
| `--parallel N` | Number of parallel processes |
| `--verbose` | Enable verbose output |
| `--quiet` | Suppress non-essential output |

### Example Commands

```bash
# Critical interaction analysis with interface verification
biostructbenchmark -e exp.pdb -p pred.pdb --hbond --dssr --multi-frame --visualize

# Batch analysis with functional residue verification
biostructbenchmark -e exp_dir/ -p pred_dir/ --all-benchmarks --visualize --parallel 4

# Quick screening with dual RMSD assessment
biostructbenchmark -e exp.pdb -p pred.pdb --rmsd-only --visualize

# Publication-ready comprehensive analysis
biostructbenchmark -e exp_dir/ -p pred_dir/ --all-benchmarks --save-aligned --export-all
```

## Critical Interaction Analysis

### Overview

The **Critical Functional Residue Interaction Analysis** addresses a key limitation in structural biology validation: ensuring that computationally predicted protein-DNA complexes maintain the essential molecular interactions required for biological function.

### Scientific Approach

**Problem Addressed**: Generic RMSD calculations fail to assess whether a predicted structure preserves the specific amino acid-DNA contacts critical for protein function.

**Solution**: Systematic verification of ≥3 critical DNA-binding residue interactions using:
1. **Specific functional atoms** (not generic residue centers)
2. **Biologically relevant distance thresholds** (not arbitrary cutoffs)
3. **Actual DNA-binding residues** (not all potential binders)

### Critical Interaction Types

#### Electrostatic Interactions
- **Arginine NH₃⁺ → DNA phosphate** (3.5Å threshold)
- **Lysine NH₃⁺ → DNA phosphate** (3.5Å threshold)

#### Hydrogen Bond Networks
- **Glutamine amide → DNA base** (3.2Å threshold)
- **Asparagine amide → DNA base** (3.2Å threshold)

#### Hydroxyl Contacts
- **Serine OH → DNA phosphate** (3.2Å threshold)
- **Threonine OH → DNA phosphate** (3.2Å threshold)

### Analysis Process

1. **Identify critical residues** in experimental structure
2. **Map corresponding residues** in predicted structure
3. **Verify functional atom contacts** within distance thresholds
4. **Calculate conservation rate** relative to actual DNA-binders
5. **Assess ≥3 requirement** for functional adequacy

### Output Interpretation

#### Conservation Rate Reporting
```
Total potential DNA-binding residues: 114
Actual DNA-binding residues (experimental): 8
Conservation rate (vs DNA-binders): 37.5% (3/8)  ← Biologically meaningful
Conservation rate (vs all residues): 2.6% (3/114) ← Misleading
```

#### Assessment Results
- **PASS**: ≥3 critical interactions conserved
- **FAIL**: <3 critical interactions conserved

### Example Results

```bash
# Real analysis output (p456_02 experimental vs AlphaFold3)
CONSERVED INTERACTIONS (3/8):
  A:SER:23 (Ser_OH_to_phosphate) - Error: 0.145 Å
  A:GLN:194 (Gln_amide_to_base) - Error: 0.071 Å  
  A:GLN:196 (Gln_amide_to_base) - Error: 0.150 Å

MISSING INTERACTIONS (5/8):
  A:ARG:82 (Arg_NH3_to_phosphate) - 3.260 Å (experimental)
  A:LYS:102 (Lys_NH3_to_phosphate) - 3.276 Å (experimental)
  ...

Assessment: PASS (≥3 requirement met)
Conservation rate: 37.5% (biologically reasonable)
```

## Module Architecture

### Core Modules (`biostructbenchmark/core/`)

#### `alignment.py`
- **Multi-frame alignment system**: Three reference frames for comprehensive assessment
- **Dual RMSD calculation**: All-atom (primary) + backbone-only (additional) for each frame
- **Gap-tolerant sequence alignment**: Handles missing residues in experimental structures
- **Intelligent structure correspondence**: Maps equivalent residues across numbering schemes

**Key Functions**:
- `perform_multi_frame_alignment()`: Execute all three reference frame alignments with dual RMSD
- `calculate_residue_rmsd()`: Calculate both all-atom and backbone RMSD for each residue/nucleotide
- `export_residue_rmsd_csv()`: Export dual RMSD data in separate ALL_ATOM_RMSD.csv and BACKBONE_RMSD.csv files

**Key Classes**:
- `ResidueRMSD`: Enhanced with `backbone_rmsd` and `backbone_atom_count` fields
- `AlignmentResult`: Multi-frame results with dual RMSD reporting
- `MultiFrameAlignmentResult`: Complete three-frame analysis

#### `metrics.py`
- **Error decomposition**: Separate translation/rotation components
- **Comprehensive statistics**: Per-chain, per-molecule type analysis
- **Dual RMSD metrics**: Statistics for both all-atom and backbone-only RMSD

#### `io.py`
- **Robust file handling**: PDB/CIF format detection and validation
- **Structure validation**: Missing residue/atom handling
- **Format flexibility**: Supports mixed file types in batch processing

### Analysis Modules (`biostructbenchmark/analysis/`)

#### `hbond.py` ⭐ **Enhanced with Critical Interaction Analysis**
- **Critical functional residue analysis**: Verifies ≥3 essential DNA-binding interactions
- **Structural correspondence mapping**: Fixes false negatives from residue numbering differences
- **Biologically meaningful reporting**: Conservation rate relative to actual DNA-binders

**New Classes**:
- `CriticalInteraction`: Container for functional residue interaction data
- **Enhanced `HBondAnalyzer`**: Now includes critical interaction methods

**Key Methods**:
- `analyze_critical_dna_binding_interactions()`: Main critical interaction analysis
- `verify_critical_interactions_requirement()`: ≥3 critical interactions verification
- `_find_critical_interaction()`: Specific functional atom contact detection

**Output Files**:
- `*_critical_interactions.csv`: Detailed critical interaction data
- `*_critical_verification.json`: Pass/fail assessment with metrics

#### `dssr.py` ⭐ **X3DNA-DSSR Integration**
- **5 critical parameters**: Base pairs, helical twist, groove widths, stacking energy
- **Threshold alerts**: Biological significance-based deviation warnings
- **Publication-quality reports**: Formatted comparison outputs

#### `bfactor.py`, `consensus.py`, `mutations.py`
- **B-factor analysis**: Experimental B-factors vs predicted confidence
- **Consensus error mapping**: Systematic error identification across multiple structures
- **Mutation impact analysis**: Local structural effects of sequence changes

### Visualization Modules (`biostructbenchmark/visualization/`)

#### `residue_plots.py` ⭐ **Enhanced for Dual RMSD**
- **Dual RMSD heatmaps**: Separate visualizations for all-atom and backbone RMSD
- **Critical interaction plots**: Conservation rate visualization
- **Publication dashboards**: Comprehensive analysis summaries
- **Memory-efficient processing**: Handles large datasets without memory issues

## Analysis Capabilities

### 1. Multi-Frame Alignment with Dual RMSD ⭐

```python
# Three reference frames, each with dual RMSD calculation
result = perform_multi_frame_alignment(exp_path, pred_path)

# All-atom RMSD (primary metric)
print(f"Full structure all-atom: {result.full_structure.overall_rmsd:.2f} Å")
print(f"DNA positioning all-atom: {result.dna_to_protein.overall_rmsd:.2f} Å")
print(f"DNA standalone all-atom: {result.dna_to_dna.overall_rmsd:.2f} Å")

# Backbone RMSD (additional metric)
for residue in result.full_structure.residue_rmsds:
    if residue.backbone_rmsd:
        print(f"{residue.residue_id}: All-atom={residue.rmsd:.3f}Å, Backbone={residue.backbone_rmsd:.3f}Å")
```

### 2. Critical Functional Residue Interaction Analysis ⭐

```python
# Analyze critical DNA-binding interactions
analyzer = HBondAnalyzer()
critical_interactions = analyzer.analyze_critical_dna_binding_interactions(
    exp_structure, pred_structure, correspondence_map)

# Verify ≥3 requirement
verification = analyzer.verify_critical_interactions_requirement(critical_interactions)

print(f"Meets ≥3 requirement: {verification['meets_requirement']}")
print(f"Conservation rate (vs DNA-binders): {verification['conservation_rate_dna_binders']:.1%}")

# Export results
analyzer.export_results(comparison, statistics, critical_interactions, 
                       verification, output_dir, pair_id)
```

### 3. X3DNA-DSSR Critical Parameters ⭐

**5 Most Critical Parameters for Protein-DNA Binding Interface**:

| Parameter | Threshold | Biological Significance |
|-----------|-----------|------------------------|
| **Base Pairs** | >2 pairs | Interface binding specificity |
| **Helical Twist** | >5° | Groove dimension alterations |
| **Major Groove Width** | >0.5Å | Primary protein recognition surface |
| **Minor Groove Width** | >0.5Å | DNA bending and flexibility |
| **Stacking Energy** | >2 kcal/mol | Interface stability |

### 4. Hydrogen Bond Network Analysis

- **Structural correspondence mapping**: Matches bonds across numbering schemes
- **Conservation assessment**: Fraction of experimental bonds preserved
- **Network topology**: Interface stability and connectivity analysis

## Output Files

### Directory Structure

```
biostructbenchmark_results/
├── structure1_vs_structure2/
│   ├── alignments/                           # Aligned structures (--save-aligned)
│   ├── visualizations/                       # Publication-quality plots
│   ├── hydrogen_bonds/                       # H-bond analysis
│   │   ├── *_critical_interactions.csv      # ⭐ Critical residue interactions
│   │   ├── *_critical_verification.json     # ⭐ ≥3 interactions assessment
│   │   ├── *_hbond_details.csv             # Detailed H-bond data
│   │   └── *_hbond_statistics.json         # Network statistics
│   ├── dssr_analysis/                       # X3DNA-DSSR critical parameters
│   ├── ALL_ATOM_RMSD_*.csv                 # ⭐ All-atom RMSD (primary)
│   ├── BACKBONE_RMSD_*.csv                 # ⭐ Backbone RMSD (additional)
│   ├── rmsd_full_structure.csv              # Combined dual RMSD
│   ├── rmsd_dna_to_protein.csv             # DNA positioning analysis
│   ├── rmsd_dna_standalone.csv             # DNA structure analysis
│   └── multi_frame_analysis.json           # Comprehensive summary
└── batch_analysis_report.json               # Batch summary
```

### Key Output Files ⭐

#### Dual RMSD Files
- **`ALL_ATOM_RMSD_*.csv`**: Primary metric including all atoms (backbone + side chains/bases)
- **`BACKBONE_RMSD_*.csv`**: Additional metric for backbone geometry assessment
- **Frame-specific files**: Separate files for each reference frame (full, DNA-to-protein, DNA-standalone)

#### Critical Interaction Files
- **`*_critical_interactions.csv`**: Detailed critical residue interaction analysis
- **`*_critical_verification.json`**: Pass/fail assessment with conservation rates

#### CSV Structure
```csv
# ALL_ATOM_RMSD_full_experimental.csv
unit_id,unit_type,all_atom_rmsd,atom_count,molecule_type,unit_class
A:ARG:156,ARG,1.707,11,protein,amino_acid
B:DT:4,DT,1.118,20,dna,nucleotide

# BACKBONE_RMSD_full_experimental.csv  
unit_id,unit_type,backbone_rmsd,backbone_atom_count,molecule_type,unit_class
A:ARG:156,ARG,1.233,4,protein,amino_acid
B:DT:4,DT,1.454,11,dna,nucleotide
```

## Best Practices

### 1. Structure Preparation

```bash
# Clean PDB files before analysis
grep -v "HOH\|HETATM" input.pdb > cleaned.pdb

# Ensure consistent chain naming between experimental and predicted
```

### 2. Choosing Analysis Modes

- **Quick screening**: `--rmsd-only --visualize` for rapid dual RMSD assessment
- **Critical interface analysis**: `--hbond --dssr --multi-frame` for functional verification
- **Publication figures**: `--all-benchmarks --visualize --export-all` for comprehensive outputs
- **Batch processing**: `biostructbenchmark -e exp_dir/ -p pred_dir/ --all-benchmarks --parallel 4`

### 3. Interpreting Results

#### RMSD Thresholds (both all-atom and backbone)
- < 2.0 Å: Excellent (near-experimental accuracy)
- 2.0-3.5 Å: Good (reliable for most analyses)
- 3.5-5.0 Å: Moderate (useful with caution)
- \> 5.0 Å: Poor (significant deviations)

#### Critical Interaction Assessment ⭐
- **Conservation Rate (vs DNA-binders)**: 
  - \>50%: Excellent functional preservation
  - 30-50%: Good functional preservation
  - <30%: Poor functional preservation
- **≥3 Requirement**: 
  - PASS: Adequate for functional interface
  - FAIL: Functional deficiencies likely

#### Dual RMSD Interpretation
- **All-atom vs Backbone**: Compare for insight into side chain vs backbone accuracy
- **All-atom higher**: Side chain positioning issues
- **Backbone higher**: Main chain geometry problems
- **Similar values**: Uniform prediction quality

## API Usage

### Python Integration

```python
from biostructbenchmark.core.alignment import perform_multi_frame_alignment
from biostructbenchmark.analysis.hbond import HBondAnalyzer
from biostructbenchmark.analysis.dssr import analyze_protein_dna_complexes

# Multi-frame alignment with dual RMSD
result = perform_multi_frame_alignment(
    experimental_path="exp.pdb",
    predicted_path="pred.pdb",
    output_dir="results/"
)

# Access dual RMSD results
for residue in result.full_structure.residue_rmsds:
    print(f"{residue.residue_id}: All-atom={residue.rmsd:.3f}Å, Backbone={residue.backbone_rmsd:.3f}Å")

# Critical interaction analysis
analyzer = HBondAnalyzer()
comparison, statistics, critical_interactions, verification = analyzer.analyze_structures_with_correspondence(
    "experimental.pdb", "predicted.pdb", correspondence_map
)

print(f"Critical interactions requirement: {'PASS' if verification['meets_requirement'] else 'FAIL'}")
print(f"Conservation rate (vs DNA-binders): {verification['conservation_rate_dna_binders']:.1%}")

# Export comprehensive results
analyzer.export_results(comparison, statistics, critical_interactions, 
                       verification, output_dir, pair_id)
```

## Troubleshooting

### Common Issues and Solutions

#### "No critical interactions found"
```bash
# Verify structure contains DNA
biostructbenchmark -e exp.pdb -p pred.pdb --verbose

# Check for proper protein-DNA complex
grep "^ATOM.*P   D" exp.pdb  # Look for DNA phosphate atoms
```

#### "Conservation rate appears low"
- **Check if using correct metric**: Use conservation rate vs DNA-binders, not vs all residues
- **Verify structure quality**: Very low rates (<10%) may indicate major structural problems
- **Check correspondence mapping**: Use `--verbose` to see residue matching details

#### Dual RMSD files not generated
```bash
# Ensure you're running recent version with dual RMSD support
biostructbenchmark --version

# Check output directory for ALL_ATOM_RMSD_*.csv and BACKBONE_RMSD_*.csv files
ls -la results/*/ALL_ATOM_RMSD_*.csv
```

## Testing and Development

### Test Organization

```
tests/
├── data/                              # Test structures
├── test_*_integration.py              # Integration tests
├── test_alignment.py                  # Core alignment tests
├── test_hbond_critical.py            # ⭐ Critical interaction tests
└── outputs/                           # Test results with dual RMSD
```

### Running Tests

```bash
# All tests including critical interaction analysis
pytest tests/

# Critical interaction specific tests
pytest tests/ -k "critical"

# Integration tests with dual RMSD verification
pytest tests/test_*integration*.py -v
```

## Contributing

We welcome contributions! Key areas for enhancement:

1. **Additional critical interactions**: Other functionally important residue-DNA contacts
2. **Improved thresholds**: Refinement of distance/angle criteria
3. **Extended analysis**: Integration with other structural validation tools
4. **Performance optimization**: Faster critical interaction detection algorithms

### Code Style
- Follow PEP 8
- Add type hints for function arguments
- Include comprehensive docstrings for critical interaction methods
- Add unit tests for new critical interaction features

## Citation

If you use BioStructBenchmark in your research, please cite:

```bibtex
@software{biostructbenchmark,
  title = {BioStructBenchmark: A toolkit for benchmarking DNA-protein structure predictions with critical functional residue interaction analysis},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/yourusername/biostructbenchmark}
}
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Acknowledgments

- BioPython team for excellent structure handling tools
- X3DNA-DSSR developers for DNA geometry analysis capabilities
- AlphaFold team for advancing structure prediction
- All contributors and users of BioStructBenchmark

---

*Built with coffee, BioPython, and an unhealthy obsession with functional structural accuracy* ☕🧬

For questions, issues, or suggestions, please open an issue on [GitHub](https://github.com/yourusername/biostructbenchmark).
