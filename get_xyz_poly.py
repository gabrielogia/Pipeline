import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
import logging
from joblib import Parallel, delayed
import stk

# Set up logging to suppress UFFTYPER warnings and capture errors
logging.getLogger('rdkit').setLevel(logging.ERROR)  # Suppress RDKit warnings
logging.basicConfig(filename='xyz_generation_errors.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Load data
df = pd.read_csv('pubchem_dataset.csv')

# Create output directory if it doesn't exist
output_dir = 'xyz_files'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Write XYZ file function
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

# Function to check if a SMILES is an acrylate
def is_acrylate(smiles):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return False
    Chem.SanitizeMol(mol, catchErrors=True)
    return mol.HasSubstructMatch(Chem.MolFromSmarts('C=C-C(=O)O'))

# Improved function to replace C=C in acrylate group with single bond and add Br atoms
def replace_first_acrylate_cce(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    acrylate_pattern = Chem.MolFromSmarts('C=C-C(=O)O')
    if not mol.HasSubstructMatch(acrylate_pattern):
        raise ValueError("No acrylate group (C=C-C(=O)O) found in the molecule")

    matches = mol.GetSubstructMatches(acrylate_pattern)
    match = matches[0]
    c1_idx, c2_idx = match[0], match[1]

    bond = mol.GetBondBetweenAtoms(c1_idx, c2_idx)
    if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
        raise ValueError("Expected double bond not found in acrylate group")

    rw_mol = Chem.RWMol(mol)
    bond.SetBondType(Chem.BondType.SINGLE)

    br1 = rw_mol.AddAtom(Chem.Atom('Br'))
    br2 = rw_mol.AddAtom(Chem.Atom('Br'))
    rw_mol.AddBond(c1_idx, br1, Chem.BondType.SINGLE)
    rw_mol.AddBond(c2_idx, br2, Chem.BondType.SINGLE)

    Chem.SanitizeMol(rw_mol)
    return Chem.MolToSmiles(rw_mol, isomericSmiles=True)

# Function to process a single monomer row (for monomer .xyz generation)
def process_monomer_row(row):
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
        
        xyz_filename = os.path.join(output_dir, f"monomer_{mol_id}.xyz")
        write_xyz_file(m_h, xyz_filename)
        return f"Saved monomer: {xyz_filename}"
    
    except Exception as e:
        logging.error(f"Error processing ID {mol_id} with SMILES {sml}: {str(e)}")
        return f"Skipping ID {mol_id}: Exception occurred - {str(e)}"

# Function to build homopolymer and save xyz file
def build_and_save_homopolymer(sml, mol_id):
    try:
        # Check if the monomer is an acrylate
        if not is_acrylate(sml):
            logging.warning(f"ID {mol_id} with SMILES {sml} is not an acrylate")
            return f"Skipping ID {mol_id}: Not an acrylate"

        # Replace C=C with Br atoms
        sml_bromo = replace_first_acrylate_cce(sml)
        
        # Create building block
        bb1 = stk.BuildingBlock(sml_bromo, [stk.BromoFactory()])
        bb2 = stk.BuildingBlock(sml_bromo, [stk.BromoFactory()])
        
        # Build homopolymer
        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1,bb2),
                repeating_unit='AB',
                num_repeating_units=4,
                optimizer=stk.Collapser(scale_steps=False),
            ),
        )
        
        # Convert to RDKit molecule
        rdkit_polymer = polymer.to_rdkit_mol()
        rdkit_polymer = Chem.AddHs(rdkit_polymer)
        Chem.SanitizeMol(rdkit_polymer)
        
        # Remove bromine atoms
        rw_mol = Chem.RWMol(rdkit_polymer)
        bromine_atoms = [atom.GetIdx() for atom in rw_mol.GetAtoms() if atom.GetSymbol() == 'Br']
        for idx in sorted(bromine_atoms, reverse=True):
            rw_mol.RemoveAtom(idx)
        rdkit_polymer = rw_mol.GetMol()
        Chem.SanitizeMol(rdkit_polymer)
        
        # Embed and optimize
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        params.maxIterations = 1000
        params.numThreads = 1
        params.randomSeed = 42
        if AllChem.EmbedMolecule(rdkit_polymer, params) == -1:
            logging.warning(f"Embedding failed for homopolymer {mol_id}")
        
        AllChem.MMFFOptimizeMolecule(rdkit_polymer)
        
        # Update polymer with new coordinates
        polymer = polymer.with_position_matrix(
            position_matrix=rdkit_polymer.GetConformer().GetPositions()
        )
        
        # Write .xyz file
        xyz_filename = os.path.join(output_dir, f"homopoly_{mol_id}.xyz")
        write_xyz_file(rdkit_polymer, xyz_filename)
        return f"Saved homopolymer: {xyz_filename}"
    
    except Exception as e:
        logging.error(f"Error building homopolymer {mol_id} with SMILES {sml}: {str(e)}")
        return f"Skipping ID {mol_id}: Exception occurred - {str(e)}"

# Filter dataset for acrylate monomers
print("Filtering dataset for acrylate monomers...")
df_acrylates = df[df['smiles'].apply(is_acrylate)].copy()
print(f"Found {len(df_acrylates)} acrylate monomers")

# Generate homopolymers for all acrylate monomers
print("Building homopolymers in parallel...")
homopolymer_results = Parallel(n_jobs=-1, backend='loky')(
    delayed(build_and_save_homopolymer)(row['smiles'], row['id'])
    for _, row in tqdm(df_acrylates.iterrows(), total=len(df_acrylates))
)

# Print summary
print("Processing complete.")
print(f"Homopolymer XYZ files generated: {sum('Saved homopolymer' in r for r in homopolymer_results)}")
print("Check 'xyz_generation_errors.log' for any errors.")