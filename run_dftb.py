import os
import glob
from multiprocessing import Pool
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor

# Define directories
xyz_dir = "xyz_files"
input_dir = "molecules"
slater_koster_dir = "3ob-3-1"  # Adjust if your Slater-Koster files are elsewhere

# Ensure input directory exists
os.makedirs(input_dir, exist_ok=True)

def get_element_angular_momentum(element):
    """Assign MaxAngularMomentum based on element group."""
    s_elements = {"H", "Li", "Na", "K", "Rb", "Cs", "Fr", "He", "Ne", "Ar", "Kr", "Xe", "Rn"}
    p_elements = {"Be", "Mg", "Ca", "Sr", "Ba", "Ra", "B", "Al", "Ga", "In", "Tl", "C", "Si", "Ge", "Sn", "Pb", "N", "P", "As", "Sb", "Bi", "O", "Se", "Te", "Po"}
    d_elements = {"F", "Cl", "Br", "I", "At", "S",
                  "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
                  "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
                  "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Ac", "Th", "Pa", "U", "Np", "Pu"}
    if element in s_elements:
        return "s"
    elif element in p_elements:
        return "p"
    elif element in d_elements:
        return "d"
    else:
        print(f"Warning: Unknown element {element}, assigning default MaxAngularMomentum 'p'")
        return "p"

def read_elements(xyz_file):
    """Read unique elements from an .xyz file."""
    elements = set()
    try:
        with open(xyz_file, "r") as f:
            lines = f.readlines()
            natoms = int(lines[0].strip())
            for line in lines[2:2 + natoms]:
                element = line.split()[0]
                elements.add(element)
    except (IndexError, ValueError) as e:
        print(f"Error reading {xyz_file}: {e}")
        return None
    return elements

def map_elements(xyz_files):
    """Scan all .xyz files to collect unique elements and their files using threads."""
    element_files = {}
    def process_file(xyz_file):
        elements = read_elements(xyz_file)
        if elements:
            return (xyz_file, elements)
        return (xyz_file, None)
    
    with ThreadPoolExecutor() as executor:
        results = list(executor.map(process_file, xyz_files))
    
    for xyz_file, elements in results:
        if elements:
            for element in elements:
                if element not in element_files:
                    element_files[element] = []
                element_files[element].append(xyz_file)
    
    angular_momentum = {element: get_element_angular_momentum(element) for element in sorted(element_files.keys())}
    
    print("Elements found in .xyz files and their MaxAngularMomentum:")
    for element, momentum in angular_momentum.items():
        print(f"  {element}: {momentum} (found in {len(element_files[element])} files, e.g., {element_files[element][0]})")
    
    return angular_momentum, element_files

def create_hsd_file(xyz_file, elements, base_name, job_dir, angular_momentum):
    """Create a dftb_in.hsd file in the job-specific directory."""
    hsd_content = f"""Geometry = xyzFormat {{
   <<< "../../{xyz_dir}/{base_name}.xyz"
}}
 
Driver = GeometryOptimization {{
   Optimizer = Rational {{}}
   MaxSteps = 10000
   OutputPrefix = "{base_name}"
   Convergence {{ GradElem = 1E-4 }}
}} 
        
Hamiltonian = DFTB {{
   SCC = Yes
   SCCTolerance = 1e-8
   MaxSCCIterations = 1000

   MaxAngularMomentum = {{
"""
    for element in sorted(elements):
        if element in angular_momentum:
            hsd_content += f'      {element} = "{angular_momentum[element]}"\n'
        else:
            print(f"Warning: Element {element} in {xyz_file} not in angular_momentum dictionary")
    hsd_content += f"""   }}
   SlaterKosterFiles = Type2FileNames {{
      Prefix = "../../{slater_koster_dir}/"
      Separator = "-"
      Suffix = ".skf"
      LowerCaseTypeName = No
   }}
}}
 
ParserOptions {{
   ParserVersion = 14
}}
"""
    hsd_file = os.path.join(job_dir, "dftb_in.hsd")
    with open(hsd_file, "w") as f:
        f.write(hsd_content)
    return hsd_file

def run_dftb(args):
    """Run DFTB+ for a single .xyz file in a unique subdirectory."""
    xyz_file, angular_momentum = args
    try:
        base_name = Path(xyz_file).stem
        elements = read_elements(xyz_file)
        if not elements:
            print(f"Skipping {xyz_file} due to invalid format")
            return

        job_dir = os.path.join(input_dir, base_name)
        os.makedirs(job_dir, exist_ok=True)

        hsd_file = create_hsd_file(xyz_file, elements, base_name, job_dir, angular_momentum)
        
        dftb_command = ["dftb+", hsd_file]
        result = subprocess.run(
            dftb_command,
            cwd=job_dir,
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            print(f"Successfully completed DFTB+ for {xyz_file}")
        else:
            print(f"Error running DFTB+ for {xyz_file}:\n{result.stderr}")
            
        # Optional: Clean up job directory
        # shutil.rmtree(job_dir)
            
    except Exception as e:
        print(f"Exception processing {xyz_file}: {e}")

def main():
    print(f"Number of CPU cores available: {os.cpu_count()}")

    xyz_files = glob.glob(os.path.join(xyz_dir, "*.xyz"))
    if not xyz_files:
        print(f"No .xyz files found in {xyz_dir}")
        return

    print(f"Found {len(xyz_files)} .xyz files to process")
    
    angular_momentum, _ = map_elements(xyz_files)
    
    # Run DFTB+ in parallel with explicit process count
    with Pool(processes=os.cpu_count()) as pool:
        pool.map(run_dftb, [(xyz_file, angular_momentum) for xyz_file in xyz_files])
    print("All DFTB+ jobs completed")

if __name__ == "__main__":
    main()