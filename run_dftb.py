import os
import glob
import subprocess
from multiprocessing import Pool
from pathlib import Path
from datetime import datetime
import time

# Define directories
XYZ_DIR = "xyz_files"
MOLECULES_DIR = "molecules"
SK_DIR = "3ob-3-1"  # Adjust if Slater-Koster files are elsewhere
LOG_FILE = "dftb_pipeline.log"

# Ensure molecules directory exists
os.makedirs(MOLECULES_DIR, exist_ok=True)

# Initialize log
with open(LOG_FILE, "w") as log:
    log.write(f"Starting DFTB+ pipeline at {datetime.now()}\n")

# Find dftb+ executable
DFTBPLUS_PATH = os.environ.get("DFTBPLUS_PATH", None)
print(DFTBPLUS_PATH)
if DFTBPLUS_PATH and os.path.isfile(DFTBPLUS_PATH):
    dftbplus_cmd = [DFTBPLUS_PATH]
else:
    # Try common locations and psi4conda path
    possible_paths = [
        "/home/jgduarte/psi4conda/bin/dftb+",  # Your specific path
        "/usr/local/bin/dftb+",
        "/usr/bin/dftb+",
        os.path.expanduser("~/bin/dftb+"),
        os.path.expanduser("~/dftbplus/bin/dftb+")
    ]
    dftbplus_cmd = ["dftb+"]
    for path in possible_paths:
        if os.path.isfile(path):
            dftbplus_cmd = [path]
            break
    # Verify dftb+ is accessible
    try:
        subprocess.run(dftbplus_cmd + ["--version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        with open(LOG_FILE, "a") as log:
            log.write("Error: dftb+ not found in PATH or common locations. Set DFTBPLUS_PATH environment variable or ensure /home/jgduarte/psi4conda/bin/dftb+ is accessible.\n")
        raise SystemExit(1)

def get_angular_momentum(element):
    """Get MaxAngularMomentum for an element."""
    s_elements = {"H", "Li", "Na", "K", "Rb", "Cs", "Fr", "He", "Ne", "Ar", "Kr", "Xe", "Rn"}
    p_elements = {"Be", "Mg", "Ca", "Sr", "Ba", "Ra", "B", "Al", "Ga", "In", "Tl", "C", "Si", "Ge", "Sn", "Pb", "N", "P", "As", "Sb", "Bi", "O", "Se", "Te", "Po"}
    d_elements = {"F", "Cl", "Br", "I", "At", "S", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Ac", "Th", "Pa", "U", "Np", "Pu", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"}
    
    if element in s_elements:
        return "s"
    elif element in p_elements:
        return "p"
    elif element in d_elements:
        return "d"
    else:
        with open(LOG_FILE, "a") as log:
            log.write(f"Warning: Unknown element {element}, using default 'p' at {datetime.now()}\n")
        return "p"

def process_xyz(xyz_file):
    """Process a single .xyz file."""
    try:
        base_name = Path(xyz_file).stem
        job_dir = os.path.join(MOLECULES_DIR, base_name)
        hsd_file = os.path.join(job_dir, "dftb_in.hsd")
        job_log = os.path.join(job_dir, "process.log")

        # Validate input file
        if not os.access(xyz_file, os.R_OK):
            with open(LOG_FILE, "a") as log:
                log.write(f"Error: Cannot read {xyz_file} at {datetime.now()}\n")
            return

        # Create job directory with retry
        for attempt in range(1, 4):
            try:
                os.makedirs(job_dir, exist_ok=True)
                break
            except OSError:
                with open(LOG_FILE, "a") as log:
                    log.write(f"Warning: Failed to create {job_dir} (attempt {attempt}) at {datetime.now()}\n")
                time.sleep(1)
        if not os.path.isdir(job_dir):
            with open(LOG_FILE, "a") as log:
                log.write(f"Error: Cannot create {job_dir} after 3 attempts at {datetime.now()}\n")
            return

        # Create process.log explicitly
        try:
            with open(job_log, "w"):
                pass
        except OSError:
            with open(LOG_FILE, "a") as log:
                log.write(f"Error: Cannot create {job_log} at {datetime.now()}\n")
            return

        # Extract unique elements from .xyz file
        elements = set()
        try:
            with open(xyz_file, "r") as f:
                lines = f.readlines()
                natoms = int(lines[0].strip())
                for line in lines[2:2 + natoms]:
                    if line.strip() and len(line.split()) > 0:
                        elements.add(line.split()[0])
        except (ValueError, IndexError, OSError):
            with open(LOG_FILE, "a") as log:
                log.write(f"Error: No elements found in {xyz_file} at {datetime.now()}\n")
            with open(job_log, "a") as log:
                log.write(f"Error: No elements found in {xyz_file} at {datetime.now()}\n")
            return

        if not elements:
            with open(LOG_FILE, "a") as log:
                log.write(f"Error: No elements found in {xyz_file} at {datetime.now()}\n")
            with open(job_log, "a") as log:
                log.write(f"Error: No elements found in {xyz_file} at {datetime.now()}\n")
            return

        # Generate dftb_in.hsd
        hsd_content = f"""Geometry = xyzFormat {{
   <<< '../../{XYZ_DIR}/{base_name}.xyz'
}}

Driver = GeometryOptimization {{
   Optimizer = Rational {{}}
   MaxSteps = 10000
   OutputPrefix = '{base_name}'
   Convergence {{ GradElem = 1E-4 }}
}}

Hamiltonian = DFTB {{
   SCC = Yes
   SCCTolerance = 1e-9
   MaxSCCIterations = 1000

   MaxAngularMomentum = {{
"""
        for element in sorted(elements):
            momentum = get_angular_momentum(element)
            hsd_content += f'      {element} = "{momentum}"\n'
        hsd_content += f"""   }}
   SlaterKosterFiles = Type2FileNames {{
      Prefix = '../../{SK_DIR}/'
      Separator = '-'
      Suffix = '.skf'
      LowerCaseTypeName = No
   }}
}}

Analysis = {{
  Polarisability = {{
    Static = Yes
    }}
}}

ParserOptions {{
   ParserVersion = 14
}}
"""
        try:
            with open(hsd_file, "w") as f:
                f.write(hsd_content)
        except OSError:
            with open(LOG_FILE, "a") as log:
                log.write(f"Error: Cannot create {hsd_file} at {datetime.now()}\n")
            with open(job_log, "a") as log:
                log.write(f"Error: Cannot create {hsd_file} at {datetime.now()}\n")
            return

        # Run DFTB+ in job directory
        try:
            os.chdir(job_dir)
            result = subprocess.run(
                dftbplus_cmd + [hsd_file],
                capture_output=True,
                text=True,
                timeout=600  # 10 minutes
            )
            if result.returncode == 0:
                with open(LOG_FILE, "a") as log:
                    log.write(f"Successfully completed DFTB+ for {xyz_file} at {datetime.now()}\n")
                with open(job_log, "a") as log:
                    log.write(f"Successfully completed DFTB+ for {xyz_file} at {datetime.now()}\n")
            else:
                with open(LOG_FILE, "a") as log:
                    log.write(f"Error running DFTB+ for {xyz_file} at {datetime.now()}\n")
                    log.write(f"DFTB+ output for {xyz_file}:\n{result.stderr}\n")
                with open(job_log, "a") as log:
                    log.write(f"Error running DFTB+ for {xyz_file} at {datetime.now()}\n")
                    log.write(f"DFTB+ output for {xyz_file}:\n{result.stderr}\n")
        except (subprocess.TimeoutExpired, subprocess.SubprocessError, OSError) as e:
            with open(LOG_FILE, "a") as log:
                log.write(f"Error: Failed to run DFTB+ for {xyz_file} at {datetime.now()}: {str(e)}\n")
            with open(job_log, "a") as log:
                log.write(f"Error: Failed to run DFTB+ for {xyz_file} at {datetime.now()}: {str(e)}\n")
        finally:
            os.chdir(os.path.dirname(os.path.abspath(__file__)) or ".")
    except Exception as e:
        with open(LOG_FILE, "a") as log:
            log.write(f"Exception processing {xyz_file} at {datetime.now()}: {str(e)}\n")

def main():
    print(f"Number of CPU cores available: {os.cpu_count()}")
    
    # Find all .xyz files
    xyz_files = glob.glob(os.path.join(XYZ_DIR, "*.xyz"))
    if not xyz_files:
        with open(LOG_FILE, "a") as log:
            log.write(f"Error: No .xyz files found in {XYZ_DIR} at {datetime.now()}\n")
        print(f"Error: No .xyz files found in {XYZ_DIR} at {datetime.now()}")
        return

    print(f"Found {len(xyz_files)} .xyz files to process at {datetime.now()}")
    with open(LOG_FILE, "a") as log:
        log.write(f"Found {len(xyz_files)} .xyz files to process at {datetime.now()}\n")

    # Run DFTB+ in parallel with 20 processes
    with Pool(processes=20) as pool:
        pool.map(process_xyz, xyz_files)

    with open(LOG_FILE, "a") as log:
        log.write(f"All DFTB+ jobs completed at {datetime.now()}\n")
    print(f"All DFTB+ jobs completed at {datetime.now()}")

if __name__ == "__main__":
    main()