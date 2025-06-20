# Pipeline

This repository contains scripts to process molecular structures and perform DFTB+ calculations. The pipeline consists of two main scripts: `get_xyz.py` and `run_dftb.py`. This README provides instructions to set up and execute these scripts.

## Prerequisites

- **Python 3.13** with the following packages:
  - `pandas`
  - `numpy`
  - `rdkit`
  - `tqdm`
  - `stk`
- **DFTB+**: Ensure DFTB+ is installed on your system. You can download it from [DFTB+ official website](https://dftbplus.org/).
- A CSV file containing SMILES strings (e.g., `pubchem_dataset.csv` with columns named `smiles` and `id`).

## Setup

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/gabrielogia/Pipeline.git
   cd Pipeline
   ```

2. **Install Python Dependencies**:
   Install the required Python packages using pip:
   ```bash
   pip install pandas numpy rdkit tqdm stk
   ```

3. **Install DFTB+**:
   - Follow the installation instructions for DFTB+ from its official documentation.
   - Ensure the `dftb+` executable is accessible in your system's PATH.

4. **Set Environment Variables**:
   Before running the scripts, configure the following environment variables:
   ```bash
   export OMP_NUM_THREADS=1
   export DFTBPLUS_PATH=$(dirname $(which dftb+))
   ```

## Usage

The pipeline involves two steps: generating XYZ coordinates from SMILES strings and running DFTB+ calculations.

### Step 1: Generate XYZ Coordinates
The `get_xyz.py` script converts SMILES strings from a CSV file into XYZ coordinate files.

**Command**:
```bash
python get_xyz.py
```

**Input**:
- A CSV file (e.g., `pubchem_dataset.csv`) with columns named `smiles` and `id`.

**Output**:
- A directory named `xyz_files` containing XYZ files for each molecule, named sequentially (e.g., `0.xyz`, `1.xyz`, etc.).

**Details**:
- The script uses RDKit to generate 3D conformations.
- It processes each SMILES string, optimizes the molecular geometry, and saves the coordinates in XYZ format.
- Progress is displayed using a tqdm progress bar.

### Step 2: Run DFTB+ Calculations
The `run_dftb.py` script executes DFTB+ calculations on the XYZ files generated in Step 1.

**Command**:
```bash
python run_dftb.py
```

**Input**:
- The `xyz_files` directory containing XYZ files from the previous step.

**Output**:
- A directory named `dftb_output` containing:
  - DFTB+ output files for each molecule (e.g., `0_dftb.out`, `1_dftb.out`, etc.).

**Details**:
- The script generates a DFTB+ input file (`dftb_in.hsd`) for each XYZ file.
- It runs DFTB+ with a predefined Hamiltonian (SCC-DFTB with specific parameters).
- The total energy is extracted from the DFTB+ output and compiled into `results.csv`.

## Example Workflow

1. Prepare an `input.csv` file with a `SMILES` column.
2. Run:
   ```bash
   export OMP_NUM_THREADS=1
   export DFTBPLUS_PATH=$(dirname $(which dftb+))
   python get_xyz.py
   python run_dftb.py
   ```
3. Check the `xyz_files` directory for XYZ files and the `dftb_output` directory for DFTB+ results.

## Notes

- Ensure the input .csv file is in the same directory as the scripts.
- The scripts assume the DFTB+ executable is correctly installed and accessible.
- The DFTB+ parameters in `run_dftb.py` are set for a specific use case; modify them as required for your calculations.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.