import os
import numpy as np
from joblib import Parallel, delayed
import multiprocessing

def read_polarizability(file_path):
    """Read electric polarizability from a detailed.out file and return as a 3x3 matrix."""
    polarizability = []
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if "Electric polarisability (a.u.)" in line:
                    for j in range(i + 1, i + 4):
                        values = [float(x) for x in lines[j].split() if x]
                        polarizability.append(values)
                    break
        if len(polarizability) == 3:
            return np.array(polarizability)
        else:
            print(f"Warning: Could not parse polarizability in {file_path}")
            return None
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def compute_trace(matrix):
    """Compute the trace of a 3x3 matrix."""
    if matrix is not None and matrix.shape == (3, 3):
        return np.trace(matrix)*0.1481847/3
    return None

def write_trace_file(subfolder_path, trace):
    """Write the trace result to a trace.txt file in the subfolder."""
    output_file = os.path.join(subfolder_path, "trace.txt")
    try:
        with open(output_file, 'w') as f:
            f.write("File Path, Polarizability Trace (a.u.)\n")
            rel_path = os.path.basename(subfolder_path) + "/detailed.out"
            f.write(f"{rel_path}, {trace:.6f}\n")
        print(f"Results written to {output_file}")
    except Exception as e:
        print(f"Error writing to {output_file}: {e}")

def process_subfolder(subfolder_path, root_folder):
    """Process detailed.out file in a single subfolder and write trace.txt."""
    results = []
    for file in os.listdir(subfolder_path):
        if file == "detailed.out":
            file_path = os.path.join(subfolder_path, file)
            polarizability_matrix = read_polarizability(file_path)
            if polarizability_matrix is not None:
                trace = compute_trace(polarizability_matrix)
                if trace is not None:
                    write_trace_file(subfolder_path, trace)
                    rel_path = os.path.relpath(file_path, root_folder)
                    results.append((rel_path, trace))
    return results

def process_molecules_folder(root_folder):
    """Process all subfolders in parallel using joblib."""
    subfolders = [os.path.join(root_folder, d) for d in os.listdir(root_folder) 
                  if os.path.isdir(os.path.join(root_folder, d))]
    
    num_cores = multiprocessing.cpu_count()
    results = Parallel(n_jobs=num_cores)(
        delayed(process_subfolder)(subfolder, root_folder) for subfolder in subfolders
    )
    
    return [item for sublist in results for item in sublist]

def main():
    molecules_folder = "molecules"
    if not os.path.exists(molecules_folder):
        print(f"Error: The folder '{molecules_folder}' does not exist.")
        return
    
    results = process_molecules_folder(molecules_folder)
    if not results:
        print("No valid polarizability data found in the molecules folder.")

if __name__ == "__main__":
    main()