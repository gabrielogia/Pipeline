import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
import logging
from joblib import Parallel, delayed

# Set up logging to suppress UFFTYPER warnings and capture errors
logging.getLogger('rdkit').setLevel(logging.ERROR)  # Suppress RDKit warnings
logging.basicConfig(filename='xyz_generation_errors.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Load and merge data
df = pd.read_csv('PubChem_compound_text_Acrylate.csv')
df = df.rename(columns={' cid': 'id'})
pub_chem_db = pd.read_csv('pubchem_db.txt', sep='\t', names=['id', 'smiles'])
df = df.merge(pub_chem_db, on='id')

# Create output directory if it doesn't exist
output_dir = 'xyz_files'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Function to write XYZ file
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

# Function to process a single row
def process_row(row):
    mol_id = row['id']
    sml = row['smiles']
    
    try:
        # Create molecule from SMILES
        m = Chem.MolFromSmiles(sml, sanitize=True)
        if m is None:
            logging.error(f"Failed to create molecule for ID {mol_id} with SMILES {sml}")
            return f"Skipping ID {mol_id}: Invalid SMILES"
    
        # Add hydrogens
        m_h = Chem.AddHs(m)
    
        # Set up embedding parameters
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        params.maxIterations = 1000
        params.numThreads = 1  # Use 1 thread per molecule to avoid RDKit threading issues
        params.randomSeed = 42
    
        # Embed molecule
        if AllChem.EmbedMolecule(m_h, params) == -1:
            logging.error(f"Embedding failed for ID {mol_id} with SMILES {sml}")
            return f"Skipping ID {mol_id}: Embedding failed"
    
        # Write .xyz file
        xyz_filename = os.path.join(output_dir, f"{mol_id}.xyz")
        write_xyz_file(m_h, xyz_filename)
        
    except Exception as e:
        logging.error(f"Error processing ID {mol_id} with SMILES {sml}: {str(e)}")
        return f"Skipping ID {mol_id}: Exception occurred - {str(e)}"

# Parallel processing with joblib
results = Parallel(n_jobs=-1, backend='loky')(
    delayed(process_row)(row) for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing molecules")
)