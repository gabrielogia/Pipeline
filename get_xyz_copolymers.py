import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
import logging
from joblib import Parallel, delayed
import stk
import numpy as np

# Set up logging to suppress UFFTYPER warnings and capture errors
logging.getLogger('rdkit').setLevel(logging.ERROR)  # Suppress RDKit warnings
logging.basicConfig(filename='xyz_generation_errors.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Load and merge data
df = pd.read_csv('pubchem_dataset.csv')

# Create output directory if it doesn't exist
output_dir = 'xyz_files'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Write XYZ file function unchanged
def write_xyz_file(mol, filename):
    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()
    with open(filename, 'w') as f:
        f.write(f"{num_atoms}\n")
        f.write(f"Molecule ID: {os.path.basename(filename).split('.')[0]}\n")
        for i in range(num_atoms):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            symbol = atom.GetSymbol()
            f.write(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

# Function to process a single molecule row (for monomer .xyz generation)
def process_row(row):
    mol_id = row['id']
    sml = row['smiles']
    
    try:
        m = Chem.MolFromSmiles(sml, sanitize=True)
        if m is None:
            logging.error(f"Failed to create molecule for ID {mol_id} with SMILES {sml}")
            return f"Skipping ID {mol_id}: Invalid SMILES"

        m_h = Chem.AddHs(m)
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        params.maxIterations = 1000
        params.numThreads = 1
        params.randomSeed = 42
        
        if AllChem.EmbedMolecule(m_h, params) == -1:
            logging.error(f"Embedding failed for ID {mol_id} with SMILES {sml}")
            return f"Skipping ID {mol_id}: Embedding failed"
        
        xyz_filename = os.path.join(output_dir, f"{mol_id}.xyz")
        write_xyz_file(m_h, xyz_filename)
        
    except Exception as e:
        logging.error(f"Error processing ID {mol_id} with SMILES {sml}: {str(e)}")
        return f"Skipping ID {mol_id}: Exception occurred - {str(e)}"

# Replace only the first occurrence of 'C=C' with '[Br]C(Br)'
def replace_first_cce(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    
    # Find the first C=C double bond
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'C':
                # Instead of SMILES replacement, modify the molecule
                mol = Chem.RWMol(mol)
                bond.SetBondType(Chem.BondType.SINGLE)
                # Add Br atoms
                br1 = mol.AddAtom(Chem.Atom('Br'))
                br2 = mol.AddAtom(Chem.Atom('Br'))
                mol.AddBond(atom1.GetIdx(), br1, Chem.BondType.SINGLE)
                mol.AddBond(atom2.GetIdx(), br2, Chem.BondType.SINGLE)
                return Chem.MolToSmiles(mol)
    
    return smiles  # No double bond found, return original

# Fixed monomers and their IDs (assuming these IDs exist in df)
fixed_monomer_1 = 'C=CC(=O)O'
fixed_monomer_2 = 'C=CC(=O)OCCO'

# Fix the fixed monomers smiles to bromo form
fixed1_bromo = replace_first_cce(fixed_monomer_1)
fixed2_bromo = replace_first_cce(fixed_monomer_2)

# Filter dataframe for candidates that contain exactly one 'C=C' (to avoid valence problems)
df_candidates = df[df['smiles'].str.count('C=C') == 1].copy()

# We will take 50 candidates containing one 'C=C' for each fixed monomer
# For reproducibility, sort and take top 50
df_candidates = df_candidates.reset_index(drop=True)

candidates_smls = df_candidates['smiles'].tolist()[:250]

# Replace first 'C=C' occurrence with bromo for candidates
candidates_bromo = [replace_first_cce(s) for s in candidates_smls]

# Function to build copolymer and save xyz file
def build_and_save_copolymer(sml1_bromo, sml2_bromo, out_name):
    try:
        bb1 = stk.BuildingBlock(sml1_bromo, [stk.BromoFactory()])
        bb2 = stk.BuildingBlock(sml2_bromo, [stk.BromoFactory()])
        
        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2),
                repeating_unit='AB',
                num_repeating_units=3,
                optimizer=stk.Collapser(scale_steps=False),
            ),
        )
        
        rdkit_polymer = polymer.to_rdkit_mol()
        rdkit_polymer = Chem.AddHs(rdkit_polymer)
        Chem.SanitizeMol(rdkit_polymer)
        
        # Embed if necessary (sometimes stk embedding is enough, but we do ETKDG for 3D coords)
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        params.maxIterations = 1000
        params.numThreads = 1
        params.randomSeed = 42
        if AllChem.EmbedMolecule(rdkit_polymer, params) == -1:
            logging.warning(f"Embedding failed for copolymer {out_name}")
        
        AllChem.MMFFOptimizeMolecule(rdkit_polymer)
        
        polymer = polymer.with_position_matrix(
            position_matrix=rdkit_polymer.GetConformer().GetPositions()
        )
        
        # Write .xyz file
        xyz_filename = os.path.join(output_dir, f"{out_name}.xyz")
        write_xyz_file(polymer.to_rdkit_mol(), xyz_filename)
        print(f"Saved copolymer: {xyz_filename}")
        
    except Exception as e:
        logging.error(f"Error building copolymer {out_name}: {str(e)}")

# First: generate xyz files for all molecules in original df (monomers)
print("Generating .xyz files for original monomers...")
results = Parallel(n_jobs=-1, backend='loky')(
    delayed(process_row)(row) for _, row in tqdm(df.iterrows(), total=len(df))
)

# Function to wrap build_and_save_copolymer for parallel execution
def process_copolymer(sml1_bromo, sml2_bromo, out_name):
    try:
        build_and_save_copolymer(sml1_bromo, sml2_bromo, out_name)
        return f"Completed copolymer: {out_name}"
    except Exception as e:
        logging.error(f"Failed to process copolymer {out_name}: {str(e)}")
        return f"Skipping copolymer {out_name}: Exception occurred - {str(e)}"

# Generate tasks for all copolymers
copolymer_tasks = []

# Tasks for fixed1 + candidates
for i, cand_bromo in enumerate(candidates_bromo, start=1):
    copolymer_tasks.append((fixed1_bromo, cand_bromo, f"copoly_fixed1_cand{i}"))

# Tasks for fixed2 + candidates
for i, cand_bromo in enumerate(candidates_bromo, start=1):
    copolymer_tasks.append((fixed2_bromo, cand_bromo, f"copoly_fixed2_cand{i}"))

# Task for fixed1 + fixed2
copolymer_tasks.append((fixed1_bromo, fixed2_bromo, "copoly_fixed1_fixed2"))

# Parallel execution of copolymer generation
print("Building all copolymers in parallel...")
results = Parallel(n_jobs=-1, backend='loky')(
    delayed(process_copolymer)(sml1_bromo, sml2_bromo, out_name)
    for sml1_bromo, sml2_bromo, out_name in tqdm(copolymer_tasks, total=len(copolymer_tasks))
)