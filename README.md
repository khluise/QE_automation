## Structural and Electronic Properties Analysis by Quantum ESPRESSO

### Automation using Python

---

## Introduction

Quantum ESPRESSO (QE) is an open-source suite widely used for electronic-structure calculations and materials modeling at the nanoscale. Based on density functional theory (DFT), plane-wave basis sets, and pseudopotentials, QE enables simulations of total energy, structural optimization, electronic band structure, phonon dispersion, and more.

While QE supports high-performance simulations with strong flexibility, its command-line interface can be less user-friendly compared to GUI-based software like VASP or Materials Studio. To address this, a Python-based tool was developed to automate QE input generation, execution, and result processing with minimal user input, while preserving QE's customizability.

This tool automates the setup and execution of:

1. Geometry Optimization
2. Band Structure Calculation
3. Density of States (DOS)
4. Projected DOS (PDOS)
5. Elastic Properties (via `thermo_pw`)
6. SOC (Spin-Orbit Coupling) effects on Band, DOS, PDOS

Input files are generated from CIF (Crystallographic Information Files), and output post-processing is handled automatically.

---

## Prerequisites

- Quantum ESPRESSO installed
- `thermo_pw` installed for elastic calculations
- Python 3 with `numpy` and `matplotlib`

---

## Usage

1. Clone the repo and place all four Python files in your working directory.
2. Place your CIF and pseudopotential files in the same folder (for SOC use fully relativistic pseudopotentials).
3. Open `QE_automation.py` and scroll to the end of the script to select which calculation(s) to run by commenting/uncommenting lines.
4. Run the script with:

```bash
python3 QE_automation.py
```

**Important**:

- Before running any non-SOC calculation or `elastic_run()`, run `get_tuned_scf()` at least once.
- For SOC calculations, run `soc.get_tuned_scf()` first.

---

## Workflow Summary
![](./workflow.png){width="6.6929in"
height="8.7299in"}

### 1. CIF to QE Input

- Lattice type and ibrav determined from cell parameters and space group
- Symmetry operations applied to atomic positions, duplicates removed
- K-points estimated using Monkhorst-Pack grid
- QE input sections: `&CONTROL`, `&SYSTEM`, `&ELECTRONS`, `ATOMIC_SPECIES`, `ATOMIC_POSITIONS`, `K_POINTS`

### 2. Parameter Optimization

- **K-points**: Converged by monitoring total energy changes using increasingly dense grids
- **ecutwfc**: Converged using fixed k-grid and varying cutoffs; ecutrho = 8 × ecutwfc
- **degauss**: Converged by scanning values and checking energy variations
- **vc-relax**: Run with optimized parameters to refine cell and atomic positions

### 3. Band Structure

- SCF: Compute ground-state density
- NSCF: Use band path and number of bands based on valence electrons
- bands.x: Processes output to generate CSV and plot band structure

### 4. DOS

- SCF & NSCF: Run with dense k-grid and tetrahedra\_opt smearing
- dos.x: Computes DOS using fine energy steps; outputs CSV and plots

### 5. PDOS

- SCF & NSCF similar to DOS
- projwfc.x: Projects states onto atomic orbitals; output parsed to CSV and plotted

### 6. SOC

- Requires `noncolin = .TRUE.` and `lspinorb = .TRUE.` in `&SYSTEM`
- SOC SCF, NSCF, bands, DOS, PDOS all supported (with FR pseudopotentials)

### 7. Elastic Constants (via Thermo\_pw)

- Automatically parses elastic constants
- Computes advanced mechanical properties including:
  - Bulk & Shear modulus
  - Poisson's ratio
  - Debye temperature
  - Melting temperature
  - Thermal conductivity (Clarke, Cahill, and lattice methods)
  - Anisotropy indexes, etc.

---

## Discussion

This tool reduces the complexity of QE-based simulations by automating tedious setup tasks and offering a clear workflow. It lowers the entry barrier for new users while supporting high-throughput studies and expert-level customizations.

With modules for elastic and SOC-aware electronic structure calculations, the tool provides a comprehensive and streamlined platform for materials research. Internally consistent file structures, error prevention via validation, and automatic plotting and data export further enhance usability.

---

## Future Work

- Enhanced runtime feedback and parameter validation
- Logging system with error handling and job resumption
- Support for phonon dispersion, Tc prediction
- Graphical User Interface (GUI) for ease of use

---

**Note**: For mathematical equations and extended derivations used in the Elastic module, refer to the supplementary PDF documentation.
---



