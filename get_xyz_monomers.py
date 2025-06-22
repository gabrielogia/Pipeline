import pandas as pd
import os, logging
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
from joblib import Parallel, delayed

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

# First: generate xyz files for all molecules in original df (monomers)
print("Generating .xyz files for original monomers...")
results = Parallel(n_jobs=-1, backend='loky')(
    delayed(process_row)(row) for _, row in tqdm(df.iterrows(), total=len(df))
)