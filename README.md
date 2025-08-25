# BioStructBenchmark 🧬

A comprehensive toolkit for benchmarking computational structure predictions against experimental reality, with specialized focus on DNA-protein complexes. Features advanced multi-frame alignment, hydrogen bond network analysis with structural correspondence mapping, and X3DNA-DSSR integration for critical protein-DNA binding interface parameters.

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Command Line Usage](#command-line-usage)
- [Module Architecture](#module-architecture)
- [Analysis Capabilities](#analysis-capabilities)
- [Output Files](#output-files)
- [Best Practices](#best-practices)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

## Overview

BioStructBenchmark compares experimental DNA-protein complex structures with their computationally predicted counterparts (AlphaFold, RoseTTAFold, etc.) to identify systematic prediction errors. The toolkit provides comprehensive structural analysis through multiple reference frame alignments, detailed error decomposition, and extensive DNA geometry analysis.

### Why BioStructBenchmark?

- **🎯 Multi-frame alignment**: Three reference frames to distinguish positioning vs structural errors
- **🧬 Advanced H-bond analysis**: Structural correspondence mapping fixes false negatives from residue numbering shifts
- **⚗️ X3DNA-DSSR integration**: 5 critical parameters for protein-DNA binding interface validation
- **🔗 Intelligent structure pairing**: Automatic matching of experimental/predicted pairs across naming conventions
- **📊 Comprehensive metrics**: Error decomposition, per-residue analysis, DNA geometry, interface statistics
- **🎨 Publication-quality outputs**: Advanced visualizations, heatmaps, dashboards, and scientific reports

## Key Features

### 🎯 Multi-Frame Alignment Analysis
- **Full structure alignment**: Overall structural accuracy assessment
- **DNA-to-protein positioning**: Evaluates DNA placement relative to protein
- **DNA standalone**: Assesses DNA structure independent of protein context

### 📊 Comprehensive Metrics
- Per-residue/nucleotide RMSD calculations with gap handling
- Error decomposition (translational vs rotational components)
- B-factor vs pLDDT confidence metric comparison
- Consensus error mapping across multiple structures
- Mutation impact analysis with local RMSD effects

### 🧬 Advanced Hydrogen Bond Analysis
- **Structural correspondence mapping** - Fixes critical false negatives from residue numbering differences
- Protein-DNA interface H-bond networks with geometric validation (distance ≤3.5Å, angle ≥120°)
- Conservation rate and prediction accuracy metrics
- Comprehensive H-bond classification (strong/medium/weak)
- Export detailed H-bond tables and network statistics

### ⚗️ X3DNA-DSSR Integration (5 Critical Parameters)
- **Base pairs**: Canonical base pairs at protein-binding interface
- **Helical twist**: Average twist per base pair step (affects groove dimensions)
- **Major groove width**: Primary protein-binding interface dimensions  
- **Minor groove width**: DNA bending and protein recognition effects
- **Stacking energy**: Interface stability indicator
- **Threshold alerts**: Flags deviations >5° (twist), >0.5Å (grooves), >2 kcal/mol (energy)

### 📈 Visualization
- Publication-quality plots and dashboards  
- Memory-efficient per-residue RMSD heatmaps
- Correlation analysis between RMSD and secondary metrics
- Comprehensive analysis dashboards with 2x2 subplot layouts
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

### Install from PyPI

```bash
pip install biostructbenchmark
```

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
# Simple RMSD calculation
biostructbenchmark observed.pdb predicted.pdb

# With output directory
biostructbenchmark -e experimental.pdb -p predicted.pdb -o results/
```

### Multi-Frame Alignment Analysis

```bash
# Perform all three reference frame alignments
biostructbenchmark -e exp.pdb -p pred.pdb --multi-frame --save-aligned

# This generates:
# 1. Full structure alignment
# 2. DNA positioning relative to protein
# 3. Standalone DNA structure comparison
```

### Batch Processing with Intelligent Pairing

```bash
# Process entire directories (automatic structure pairing)
biostructbenchmark -e experimental_dir/ -p predicted_dir/ -o results/ --all-benchmarks

# Intelligent pairing handles naming differences:
# experimental: p456_02_experimental.pdb ↔ predicted: p456_02_alphafold3.cif
# experimental: 1abc_exp.pdb ↔ predicted: 1abc_pred.cif

# With parallel processing and hydrogen bond analysis
biostructbenchmark -e exp/ -p pred/ -o out/ --hbond --dssr --parallel 4
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
| `--all-benchmarks` | Run all available analyses |
| `--multi-frame` | Perform multi-frame alignment (3 reference frames) |
| `--rmsd-only` | Basic RMSD analysis only (fastest) |
| `--hbond` | Hydrogen bond network analysis with structural correspondence |
| `--dssr` | X3DNA-DSSR analysis (5 critical protein-DNA binding parameters) |
| `--bfactor` | B-factor vs confidence metric analysis |
| `--consensus` | Consensus error mapping (requires multiple pairs) |
| `--mutations` | Detect and analyze mutations |
| `--visualize` | Generate publication-quality plots |

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
# Full analysis with all features and visualizations
biostructbenchmark -e exp.pdb -p pred.pdb --all-benchmarks --visualize

# Protein-DNA interface analysis with critical parameters
biostructbenchmark -e exp.pdb -p pred.pdb --hbond --dssr --multi-frame --visualize

# Batch with hydrogen bond and DNA structural analysis
biostructbenchmark -e exp_dir/ -p pred_dir/ --hbond --dssr --bfactor --visualize

# High-throughput processing with comprehensive analysis
biostructbenchmark -e structures/ -p predictions/ --parallel 8 --all-benchmarks --visualize

# Quick assessment with intelligent structure pairing
biostructbenchmark -e ~/experimental/ -p ~/predicted_alphafold3/ --rmsd-only --visualize
```

## Module Architecture

### Core Modules (`biostructbenchmark/core/`)

#### `alignment.py`
- **Purpose**: Structure superposition and RMSD calculations
- **Key Functions**:
  - `perform_multi_frame_alignment()`: Execute all three reference frame alignments
  - `align_structures_by_reference_frame()`: Flexible alignment with different references
  - `calculate_rmsd()`: RMSD calculation for atoms or coordinates
  - `calculate_per_residue_rmsd_for_subset()`: Per-residue RMSD for specific molecule types
- **Classes**:
  - `AlignmentResult`: Container for alignment results
  - `MultiFrameAlignmentResult`: Results from all three alignments
  - `ResidueRMSD`: Per-residue RMSD data

#### `metrics.py`
- **Purpose**: Error decomposition and comprehensive metrics
- **Key Functions**:
  - `decompose_structural_error()`: Separate translation/rotation components
  - `calculate_per_chain_metrics()`: Chain-specific analysis
  - `generate_comprehensive_metrics()`: Full metrics report
- **Key Metrics**:
  - Translation error magnitude
  - Rotation error (Frobenius norm)
  - Per-molecule type statistics
  - Backbone vs side-chain RMSD

#### `io.py`
- **Purpose**: Robust structure file handling
- **Features**:
  - Automatic PDB/CIF format detection
  - Missing residue/atom handling
  - Chain mismatch resolution
  - Structure validation

### Analysis Modules (`biostructbenchmark/analysis/`)

#### `hbond.py` ⭐ **Advanced H-Bond Analysis**
- **Purpose**: Hydrogen bond network analysis with structural correspondence mapping
- **Features**:
  - Geometric H-bond detection (distance ≤3.5Å, angle ≥120°)
  - Structural correspondence integration for proper bond matching
  - Conservation rate and prediction accuracy metrics
  - Detailed H-bond classification and network statistics
- **Classes**:
  - `HBondAnalyzer`: Main analysis class with correspondence support
  - `HydrogenBond`: H-bond container with geometric parameters
  - `HBondComparison`: Network comparison results
  - `HBondStatistics`: Comprehensive statistics

#### `dssr.py` ⭐ **X3DNA-DSSR Integration**
- **Purpose**: Critical protein-DNA binding interface parameters (5 key metrics)
- **Parameters**:
  - Base pairs: Canonical base pairs at binding interface
  - Helical twist: Average twist per base pair step (affects groove dimensions)
  - Major groove width: Primary protein-binding interface
  - Minor groove width: DNA bending effects
  - Stacking energy: Interface stability indicator
- **Features**:
  - Threshold alerts: >5° twist, >0.5Å groove, >2 kcal/mol energy deviations
  - Publication-quality comparison reports
  - Seamless integration with existing pipeline
- **Classes**:
  - `DSSRAnalyzer`: Main DSSR wrapper
  - `DSSRParameters`: Critical parameter container

#### `bfactor.py`
- **Purpose**: B-factor and confidence metric analysis
- **Features**:
  - B-factor extraction from experimental structures
  - pLDDT extraction from AlphaFold predictions
  - Correlation analysis
  - Disorder region identification

#### `consensus.py`
- **Purpose**: Identify consistent prediction errors
- **Requirements**: Minimum 3 structure pairs
- **Outputs**:
  - Consensus error positions
  - Statistical significance
  - Error frequency maps

#### `mutations.py`
- **Purpose**: Mutation detection and impact analysis
- **Features**:
  - Automatic mutation detection
  - Local RMSD impact calculation
  - Mutation clustering analysis

### Visualization Modules (`biostructbenchmark/visualization/`)

#### `residue_plots.py`
- **Purpose**: Comprehensive residue-level visualization with heatmaps and correlations
- **Generates**:
  - Per-residue RMSD heatmaps (memory-efficient using imshow)
  - Correlation matrices between RMSD and secondary data
  - Interactive dashboards for comprehensive analysis
  - Chain-by-chain comparison plots
- **Features**:
  - Unified module combining all plotting capabilities
  - Legacy compatibility with `PublicationPlotter` interface
  - Memory-efficient processing for large datasets

#### `structure.py`
- **Purpose**: 3D structure visualization
- **Features**:
  - Matplotlib 3D backbone traces
  - Per-residue RMSD coloring
  - Alignment visualization
  - Summary statistics export

#### `curves_plots.py`
- **Purpose**: DNA geometry visualization
- **Plots**:
  - Base pair parameter comparisons
  - Groove geometry analysis
  - Hydrogen bond networks
  - Parameter correlation matrices

## Analysis Capabilities

### 1. Multi-Frame Alignment Analysis

The cornerstone feature providing three complementary perspectives:

```python
# Alignment 1: Full structure
# Answers: How accurate is the overall prediction?
full_rmsd = result.full_structure.overall_rmsd

# Alignment 2: DNA aligned using protein reference
# Answers: Is the DNA correctly positioned relative to the protein?
positioning_rmsd = result.dna_to_protein.overall_rmsd

# Alignment 3: DNA-only alignment
# Answers: Is the DNA structure itself accurate?
dna_rmsd = result.dna_to_dna.overall_rmsd
```

**Interpretation Guide**:
- High positioning RMSD + Low DNA RMSD = Interface prediction problem
- Low positioning RMSD + High DNA RMSD = DNA geometry problem
- Both high = Systematic prediction issues

### 2. DNA Geometry Analysis

Comprehensive DNA structure characterization:

- **Base Pair Parameters**: Quantify intra-base pair distortions
- **Step Parameters**: Inter-base pair relationships
- **Groove Analysis**: Major/minor groove dimensions
- **Global Metrics**: Bending, writhe, handedness

### 3. Advanced Hydrogen Bond Network Analysis ⭐


**Features**:
- Geometric H-bond detection (distance ≤3.5Å, D-H-A angle ≥120°)
- **Correspondence mapping**: Matches H-bonds across different residue numbering schemes
- Conservation rate: Fraction of experimental bonds preserved in predictions
- Prediction accuracy: Fraction of predicted bonds that are correct
- Network topology and interface stability assessment

### 4. X3DNA-DSSR Critical Parameters ⭐

**5 Most Critical Parameters for Protein-DNA Binding Interface**:

| Parameter | Threshold | Biological Significance |
|-----------|-----------|------------------------|
| **Base Pairs** | >2 pairs | Interface binding specificity |
| **Helical Twist** | >5° | Groove dimension alterations |
| **Major Groove Width** | >0.5Å | Primary protein recognition surface |
| **Minor Groove Width** | >0.5Å | DNA bending and flexibility |
| **Stacking Energy** | >2 kcal/mol | Interface stability |

## Output Files

### Directory Structure

```
biostructbenchmark_results/
├── structure1_vs_structure2/
│   ├── alignments/                    # Aligned PDB files (with --save-aligned)
│   │   ├── aligned_1_full_to_experimental.pdb
│   │   ├── aligned_2_dna_to_protein_reference.pdb
│   │   ├── aligned_3_dna_to_dna.pdb
│   │   └── experimental_reference.pdb
│   ├── visualizations/               # Publication-quality visualizations
│   │   ├── residue_heatmap.png      # Per-residue RMSD heatmap
│   │   ├── residue_correlation.png   # Correlation matrix
│   │   ├── residue_chains.png       # Chain comparison
│   │   ├── residue_dashboard.png    # Comprehensive dashboard
│   │   ├── residue_summary.csv      # Visualization data
│   │   └── analysis_summary.png     # Analysis overview
│   ├── hydrogen_bonds/              # H-bond analysis with correspondence
│   │   ├── structure_vs_structure_hbond_details.csv
│   │   ├── structure_vs_structure_hbond_summary.csv
│   │   └── structure_vs_structure_hbond_statistics.json
│   ├── dssr_analysis/               # X3DNA-DSSR critical parameters
│   │   ├── dssr_parameters.csv      # 5 critical parameters
│   │   └── structure_vs_structure_dssr_comparison.txt
│   ├── rmsd_full_structure.csv       # Per-residue RMSD (3 reference frames)
│   ├── rmsd_dna_to_protein.csv      
│   ├── rmsd_dna_standalone.csv      
│   ├── multi_frame_analysis.json     # Comprehensive summary
│   └── bfactor_comparison.csv       
└── batch_analysis_report.json        # Batch summary with intelligent pairing
```

### Key Output Files

#### `multi_frame_analysis.json`
Complete analysis summary including:
- Overall RMSD values for all three alignments
- Per-molecule type statistics
- Worst-performing residues
- Automated interpretation

#### `hydrogen_bonds/` Directory ⭐
**Advanced H-bond analysis with structural correspondence**:
- `*_hbond_details.csv`: Complete H-bond listing with donor/acceptor atoms, distances, angles
- `*_hbond_summary.csv`: Network comparison metrics (conservation rate, prediction accuracy)
- `*_hbond_statistics.json`: Comprehensive statistics with correspondence mapping results

#### `dssr_analysis/` Directory ⭐
**X3DNA-DSSR critical parameters for protein-DNA binding**:
- `dssr_parameters.csv`: 5 critical parameters (base pairs, twist, groove widths, stacking energy)
- `*_dssr_comparison.txt`: Formatted comparison report with threshold alerts

## Best Practices

### 1. Structure Preparation

```bash
# Clean PDB files before analysis
# Remove water molecules, heteroatoms if needed
grep -v "HOH\|HETATM" input.pdb > cleaned.pdb

# Ensure consistent chain naming
# Experimental: Chain A (protein), Chain B (DNA)
# Predicted: Same chain naming scheme
```

### 2. Choosing Analysis Modes

- **Quick assessment**: Use `--rmsd-only --visualize` for rapid screening with heatmaps
- **Publication figures**: Always use `--visualize --export-all` for comprehensive outputs  
- **Protein-DNA interface focus**: Use `--hbond --dssr --multi-frame --visualize` for critical binding parameters
- **Comprehensive analysis**: Use `--all-benchmarks --visualize` for complete analysis with dashboards
- **Batch processing**: Use intelligent pairing: `biostructbenchmark -e exp_dir/ -p pred_dir/ --all-benchmarks`

### 3. Interpreting Results

#### RMSD Thresholds
- < 2.0 Å: Excellent (near-experimental accuracy)
- 2.0-3.5 Å: Good (reliable for most analyses)
- 3.5-5.0 Å: Moderate (useful with caution)
- \> 5.0 Å: Poor (significant deviations)

#### Hydrogen Bond Analysis Results ⭐
- **Conservation Rate**: 40-70% = Good, <40% = Poor interface prediction
- **Prediction Accuracy**: >70% = Excellent, 50-70% = Good, <50% = Poor
- **0% Rates**: Likely residue numbering issue (correspondence mapping fixes this)

#### X3DNA-DSSR Critical Thresholds ⭐
- **Base Pairs**: >2 difference = Interface binding changes
- **Helical Twist**: >5° deviation = Groove dimension alterations
- **Groove Widths**: >0.5Å deviation = Altered protein recognition
- **Stacking Energy**: >2 kcal/mol = Interface stability issues

### 4. Batch Processing

```bash
# For large datasets, use parallel processing
biostructbenchmark -e exp_dir/ -p pred_dir/ --parallel 8

# Monitor progress with verbose mode
biostructbenchmark -e exp/ -p pred/ --verbose --parallel 4

# For HPC environments
export OMP_NUM_THREADS=1  # Prevent nested parallelism
biostructbenchmark -e exp/ -p pred/ --parallel $SLURM_CPUS_PER_TASK
```

### 5. Memory Management

- Process large structures individually rather than in batch
- Use `--rmsd-only` for initial screening of large datasets
- Consider splitting very large complexes into domains

## Troubleshooting

### Common Issues and Solutions

#### "X3DNA-DSSR not found"
```bash
# Install X3DNA-DSSR from https://x3dna.org/
# Add to PATH or specify explicitly
export PATH=$PATH:/path/to/dssr/bin
dssr --version  # Test installation
```

#### "No matching structure pairs found"
```bash
# Check file naming patterns - intelligent pairing handles:
# exp: p456_02_experimental.pdb ↔ pred: p456_02_alphafold3.cif
# Use --verbose to see pairing details
biostructbenchmark -e exp_dir/ -p pred_dir/ --verbose
```

#### "No common residues found"
- Check chain naming consistency
- Verify sequence alignment
- Use `--verbose` to see detailed matching

#### Memory errors with large structures
```bash
# Reduce parallel processes
biostructbenchmark -e exp/ -p pred/ --parallel 2

# Process individually
for file in exp/*.pdb; do
    biostructbenchmark -e $file -p pred/$(basename $file) -o results/
done
```

#### Visualization errors
```bash
# Check matplotlib backend
export MPLBACKEND=Agg  # For headless systems

# Install optional dependencies
pip install seaborn py3Dmol

# For detailed visualizations, ensure residue data is available
biostructbenchmark -e exp.pdb -p pred.pdb --visualize  # Auto-includes residue data
```

#### No visualization outputs generated
- Ensure you use the `--visualize` flag
- Check that analysis completed successfully (RMSD calculation)
- For detailed heatmaps, residue data must be available
- Use `--verbose` to see detailed visualization status

## API Usage

### Python Integration

```python
from biostructbenchmark.core.alignment import perform_multi_frame_alignment
from biostructbenchmark.analysis.hbond import HBondAnalyzer
from biostructbenchmark.analysis.dssr import analyze_protein_dna_complexes
from biostructbenchmark.visualization.residue_plots import create_residue_analysis

# Multi-frame alignment
result = perform_multi_frame_alignment(
    experimental_path="exp.pdb",
    predicted_path="pred.pdb",
    output_dir="results/"
)

# Access results
print(f"Full RMSD: {result.full_structure.overall_rmsd:.2f} Å")
print(f"DNA positioning: {result.dna_to_protein.overall_rmsd:.2f} Å")
print(f"DNA structure: {result.dna_to_dna.overall_rmsd:.2f} Å")

# Advanced H-bond analysis with correspondence mapping
analyzer = HBondAnalyzer()
comparison, statistics = analyzer.analyze_structures_with_correspondence(
    "experimental.pdb", "predicted.pdb", correspondence_map
)
print(f"H-bond conservation rate: {statistics.conservation_rate:.2%}")
print(f"H-bond prediction accuracy: {statistics.prediction_accuracy:.2%}")

# X3DNA-DSSR critical parameters
df = analyze_protein_dna_complexes(
    experimental_pdbs=["exp.pdb"],
    predicted_pdbs=["pred.pdb"],
    output_csv="dssr_results.csv"
)

# Publication-quality visualizations
viz_paths = create_residue_analysis(
    residue_data=result.full_structure.residue_rmsds,
    output_dir=Path("visualizations/"),
    analysis_data={"hbond": hbond_data, "dssr": dssr_data}
)
```

## Contributing

We welcome contributions! Please follow these guidelines:

### Code Style
- Follow PEP 8
- Use type hints for function arguments
- Add comprehensive docstrings
- Include unit tests for new features

### Pull Request Process
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Testing
```bash
# Run all tests
pytest tests/

# Run specific test module
pytest tests/test_alignment.py -v

# Check coverage
pytest --cov=biostructbenchmark tests/
```

## Citation

If you use BioStructBenchmark in your research, please cite:

```bibtex
@software{biostructbenchmark,
  title = {BioStructBenchmark: A toolkit for benchmarking DNA-protein structure predictions},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/yourusername/biostructbenchmark}
}
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Acknowledgments

- BioPython team for excellent structure handling tools
- CURVES+ developers for DNA geometry analysis
- AlphaFold team for advancing structure prediction
- All contributors and users of BioStructBenchmark

---

*Built with coffee, BioPython, and an unhealthy obsession with structural accuracy* 🐱

For questions, issues, or suggestions, please open an issue on [GitHub](https://github.com/yourusername/biostructbenchmark).

---

## Complete Command Line Interface Reference

### Synopsis
```bash
biostructbenchmark [OPTIONS] [experimental_file predicted_file]
```

### Input/Output Options
```bash
-e, --experimental PATH     Path to experimental structure file or directory
-p, --predicted PATH        Path to predicted structure file or directory  
-o, --output PATH          Output directory (default: ./biostructbenchmark_results)
-v, --version              Show version and exit
--help                     Show help message and exit
```

### Analysis Selection
```bash
# Core analyses
--all-benchmarks           Enable all available analyses
--multi-frame              Multi-frame alignment (3 reference frames)
--rmsd-only               Basic RMSD analysis only (fastest)

# Advanced analyses ⭐
--hbond                   Hydrogen bond analysis with structural correspondence
--dssr                    X3DNA-DSSR critical protein-DNA binding parameters
--bfactor                 B-factor vs confidence correlation analysis
--consensus               Consensus error mapping (requires ≥3 structure pairs)
--mutations               Mutation detection and impact analysis

# Visualization
--visualize               Generate publication-quality plots and dashboards
```

### Analysis Parameters
```bash
--reference-frame FRAME   Alignment reference: {full,protein,dna,multi}
--rmsd-threshold FLOAT    RMSD threshold for consensus analysis (default: 3.0)
--output-format FORMAT    Output format: {csv,json,both} (default: both)
```

### Output Control
```bash
--save-aligned            Save aligned PDB structures to alignments/ directory
--export-all              Export all intermediate analysis files  
--quiet                   Suppress non-essential output
--verbose                 Enable detailed progress information
```

### Performance Options
```bash
--parallel N              Number of parallel processes for batch analysis
```

### Usage Examples

#### Basic Analysis
```bash
# Single structure pair
biostructbenchmark experimental.pdb predicted.pdb

# With explicit paths and output
biostructbenchmark -e exp.pdb -p pred.pdb -o results/

# Quick RMSD screening
biostructbenchmark -e exp.pdb -p pred.pdb --rmsd-only --visualize
```

#### Advanced Protein-DNA Interface Analysis ⭐
```bash
# Critical binding interface parameters
biostructbenchmark -e exp.pdb -p pred.pdb --hbond --dssr --multi-frame --visualize

# H-bond analysis with structural correspondence (fixes false negatives)
biostructbenchmark -e exp.pdb -p pred.pdb --hbond --export-all

# X3DNA-DSSR critical parameters (5 key metrics)
biostructbenchmark -e exp.pdb -p pred.pdb --dssr --visualize
```

#### Batch Processing with Intelligent Pairing
```bash
# Automatic structure pairing across naming conventions
biostructbenchmark -e experimental_dir/ -p predicted_dir/ --all-benchmarks

# High-throughput with parallel processing
biostructbenchmark -e exp_dir/ -p pred_dir/ -o results/ --parallel 8 --all-benchmarks

# Interface-focused batch analysis
biostructbenchmark -e exp_dir/ -p pred_dir/ --hbond --dssr --visualize --parallel 4
```

#### Comprehensive Analysis Pipeline
```bash
# Complete analysis with all features
biostructbenchmark -e exp_dir/ -p pred_dir/ -o comprehensive_results/ \
    --all-benchmarks --visualize --save-aligned --export-all --parallel 6

# Research-grade analysis with detailed outputs
biostructbenchmark -e experimental/ -p alphafold_predictions/ \
    --multi-frame --hbond --dssr --bfactor --consensus --visualize \
    --verbose --export-all
```

### File Format Support
- **Input**: PDB (.pdb), mmCIF (.cif), compressed files (.gz)
- **Output**: CSV, JSON, PNG (visualizations), PDB (aligned structures)
- **Intelligent Pairing**: Handles naming differences automatically
  - `p456_02_experimental.pdb` ↔ `p456_02_alphafold3.cif`
  - `1abc_exp.pdb` ↔ `1abc_pred.pdb`
  - `structure_experimental.pdb` ↔ `structure_predicted.cif`

### Exit Codes
- **0**: Successful completion
- **1**: Analysis errors or failures
- **130**: Interrupted by user (Ctrl+C)

### Environment Variables
```bash
export DSSR_EXEC=/path/to/dssr          # X3DNA-DSSR executable path
export MPLBACKEND=Agg                   # Matplotlib backend for headless systems
export OMP_NUM_THREADS=1                # Prevent nested parallelism in HPC
```

### Advanced Features

#### Multi-Frame Alignment Interpretation
- **Full structure**: Overall prediction accuracy
- **DNA positioning**: DNA placement relative to protein  
- **DNA standalone**: Intrinsic DNA structure quality

#### Hydrogen Bond Analysis Innovation ⭐
- **Structural correspondence mapping** fixes false negatives from residue numbering shifts
- **Before**: 0% conservation rate (false negatives)
- **After**: 40-70% conservation rate (true accuracy)

#### X3DNA-DSSR Critical Parameters ⭐
- **Base Pairs**: Interface binding specificity
- **Helical Twist**: Groove dimension effects (>5° threshold)
- **Major Groove Width**: Primary protein recognition (>0.5Å threshold)  
- **Minor Groove Width**: DNA bending effects (>0.5Å threshold)
- **Stacking Energy**: Interface stability (>2 kcal/mol threshold)

### Notes
- Use `--verbose` to see detailed progress and pairing information
- Use `--visualize` for publication-quality plots and dashboards  
- Combine `--hbond --dssr` for comprehensive protein-DNA interface analysis
- All analyses are compatible and can be combined for comprehensive evaluation

## Testing and Development

### Test Organization

The repository tests are organized in the `tests/` directory:

```
tests/
├── README.md                          # Test documentation  
├── data/                              # Test data files
│   ├── proteins_pdb/                  # PDB format test structures
│   ├── proteins_cif/                  # CIF format test structures
│   ├── experimental/                  # Experimental structures
│   └── predicted_alphafold3/          # AlphaFold predicted structures
├── demos/                             # Demo scripts and examples
│   ├── alignment_demo.py              # Basic alignment demonstration
│   └── comprehensive_alignment_demo.py # Comprehensive alignment example
├── outputs/                           # Test outputs and results
│   ├── alignment_outputs*/            # Alignment test results
│   ├── test_*_outputs/               # Integration test outputs
│   └── htmlcov/                      # HTML coverage reports
└── test_*.py                         # Test scripts
```

### Test Categories

#### Unit Tests
- `test_io.py` - I/O functionality tests
- `test_metrics.py` - Metrics calculation tests
- `test_alignment.py` - Core alignment tests
- `test_visualization.py` - Visualization tests

#### Integration Tests
- `test_core_analysis_integration.py` - Core + Analysis integration
- `test_core_integration.py` - Core module integration
- `test_integration_comprehensive.py` - Comprehensive integration tests
- `test_pca_integration.py` - PCA analysis integration
- `test_dssr_integration.py` - DSSR nucleic acid analysis (requires DSSR)
- `test_curves_setup.py` - CURVES+ setup verification

#### Analysis Tests
- `test_core_analysis_basic.py` - Basic analysis pipeline
- `test_alignment_comprehensive.py` - Comprehensive alignment testing
- `test_dna_p_alignment.py` - DNA-protein alignment tests

### Running Tests

#### All Tests
```bash
# From repository root
python -m pytest tests/

# With coverage
python -m pytest tests/ --cov=biostructbenchmark --cov-report=html
```

#### Specific Test Categories
```bash
# Unit tests only
python -m pytest tests/test_io.py tests/test_metrics.py

# Integration tests
python -m pytest tests/test_*integration*.py

# PCA analysis tests
python tests/test_pca_integration.py
```

#### Individual Integration Tests
```bash
# Core-Analysis integration
python tests/test_core_analysis_integration.py

# PCA analysis
python tests/test_pca_integration.py

# DSSR setup (requires DSSR license)
python tests/test_dssr_integration.py
```

### Test Data

#### Structure Files
- **1bom.pdb/cif** - Small protein structure for basic tests
- **2r4g.pdb/cif** - Medium complexity structure
- **p456_02_experimental.pdb** - Experimental DNA-protein complex
- **p456_02_predicted.cif** - AlphaFold predicted structure

#### Invalid Files
- `empty.cif` - Empty file for error handling tests
- `invalid.pdb` - Malformed PDB for parser tests
- `no_extension` - File without extension for type detection tests

### Demo Scripts

#### alignment_demo.py
Basic demonstration of structure alignment:
```bash
python tests/demos/alignment_demo.py
```

#### comprehensive_alignment_demo.py
Full pipeline demonstration with analysis:
```bash
python tests/demos/comprehensive_alignment_demo.py
```

### External Dependencies

Some tests require external tools:

- **DSSR** (3DNA suite) - For nucleic acid analysis
  - Test: `test_dssr_integration.py`
  - Setup: Run `./setup_dssr.sh` after obtaining license from https://x3dna.org/

- **CURVES+** - Alternative nucleic acid analysis (deprecated in favor of DSSR)
  - Test: `test_curves_setup.py`
  - Setup: `./install_curves_plus.sh`

### New Analysis Modules

The toolkit includes several advanced analysis capabilities:

#### Principal Component Analysis (PCA)
- **Purpose**: Identify outlier structures and error patterns
- **Module**: `biostructbenchmark.analysis.pca`
- **Test**: `tests/test_pca_integration.py`
- **Features**:
  - Structure-level outlier detection
  - Residue-level error pattern analysis
  - Feature importance identification
  - Comprehensive visualization plots

#### DSSR Integration (Nucleic Acid Analysis)
- **Purpose**: Modern replacement for CURVES+ using 3DNA-DSSR
- **Module**: `biostructbenchmark.analysis.dssr`
- **Test**: `tests/test_dssr_integration.py`
- **Features**:
  - Base pair and base step parameter analysis
  - Groove geometry measurements
  - Protein-nucleic acid contact detection
  - JSON-based output parsing

#### Enhanced Analysis Pipeline
- **B-factor Analysis**: Compare experimental B-factors with predicted confidence
- **Consensus Analysis**: Identify systematically mispredicted regions
- **Mutation Analysis**: Detect and analyze structural impact of mutations
- **Secondary Structure**: Compare predicted vs experimental secondary structure

### Adding New Tests

1. **Unit Tests**: Add to existing `test_*.py` files or create new ones
2. **Integration Tests**: Follow naming pattern `test_*_integration.py`
3. **Test Data**: Add to appropriate `data/` subdirectory
4. **Outputs**: Will be automatically organized in `outputs/`

### Notes

- Tests are organized by functionality and complexity
- All test outputs are preserved for analysis and debugging
- Demo scripts serve as both examples and integration tests
- External tool tests gracefully handle missing dependencies