#!/bin/bash

# Define directories
XYZ_DIR="xyz_files"
MOLECULES_DIR="molecules"
SK_DIR="3ob-3-1"  # Adjust if Slater-Koster files are elsewhere
LOG_FILE="dftb_pipeline.log"

# Ensure molecules directory exists
mkdir -p "$MOLECULES_DIR" || { echo "Error: Cannot create $MOLECULES_DIR" >&2; exit 1; }

# Initialize log
echo "Starting DFTB+ pipeline at $(date)" > "$LOG_FILE"

# Check for dftb+ executable
if ! command -v dftb+ > /dev/null; then
    echo "Error: dftb+ not found in PATH" | tee -a "$LOG_FILE"
    exit 1
fi

# Function to get MaxAngularMomentum for an element
get_angular_momentum() {
    local element=$1
    case $element in
        H|Li|Na|K|Rb|Cs|Fr|He|Ne|Ar|Kr|Xe|Rn)
            echo "s" ;;
        Be|Mg|Ca|Sr|Ba|Ra|B|Al|Ga|In|Tl|C|Si|Ge|Sn|Pb|N|P|As|Sb|Bi|O|Se|Te|Po)
            echo "p" ;;
        F|Cl|Br|I|At|S|Sc|Ti|V|Cr|Mn|Fe|Co|Ni|Cu|Zn|Y|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Hf|Ta|W|Re|Os|Ir|Pt|Au|Hg|Ac|Th|Pa|U|Np|Pu|Cm|Bk|Cf|Es|Fm|Md|No|Lr)
            echo "d" ;;
        *)
            echo "Warning: Unknown element $element, using default 'p'" >> "$LOG_FILE"
            echo "p" ;;
    esac
}

# Function to process a single .xyz file
process_xyz() {
    local xyz_file="$1"
    local base_name=$(basename "$xyz_file" .xyz)
    local job_dir="$MOLECULES_DIR/$base_name"
    local hsd_file="$job_dir/dftb_in.hsd"
    local job_log="$job_dir/process.log"

    # Validate input file
    if [ ! -r "$xyz_file" ]; then
        echo "Error: Cannot read $xyz_file at $(date)" >> "$LOG_FILE"
        return 1
    fi

    # Create job directory with retry
    for attempt in {1..3}; do
        mkdir -p "$job_dir" 2>/dev/null && break
        echo "Warning: Failed to create $job_dir (attempt $attempt) at $(date)" >> "$LOG_FILE"
        sleep 1
    done
    if [ ! -d "$job_dir" ]; then
        echo "Error: Cannot create $job_dir after 3 attempts at $(date)" >> "$LOG_FILE"
        return 1
    fi

    # Create process.log explicitly
    touch "$job_log" 2>/dev/null || { echo "Error: Cannot create $job_log at $(date)" >> "$LOG_FILE"; return 1; }

    # Extract unique elements from .xyz file
    elements=$(awk 'NR>2 && NF>0 {print $1}' "$xyz_file" 2>/dev/null | sort -u)
    if [ -z "$elements" ]; then
        echo "Error: No elements found in $xyz_file at $(date)" | tee -a "$job_log" "$LOG_FILE"
        return 1
    fi

    # Generate dftb_in.hsd
    {
        echo "Geometry = xyzFormat {"
        echo "   <<< '../../$XYZ_DIR/$base_name.xyz'"
        echo "}"
        echo ""
        echo "Driver = GeometryOptimization {"
        echo "   Optimizer = Rational {}"
        echo "   MaxSteps = 10000"
        echo "   OutputPrefix = '$base_name'"
        echo "   Convergence { GradElem = 1E-4 }"
        echo "}"
        echo ""
        echo "Hamiltonian = DFTB {"
        echo "   SCC = Yes"
        echo "   SCCTolerance = 1e-8"
        echo "   MaxSCCIterations = 1000"
        echo ""
        echo "   MaxAngularMomentum = {"
        for element in $elements; do
            momentum=$(get_angular_momentum "$element")
            echo "      $element = \"$momentum\""
        done
        echo "   }"
        echo "   SlaterKosterFiles = Type2FileNames {"
        echo "      Prefix = '../../$SK_DIR/'"
        echo "      Separator = '-'"
        echo "      Suffix = '.skf'"
        echo "      LowerCaseTypeName = No"
        echo "   }"
        echo "}"
        echo ""
        echo "ParserOptions {"
        echo "   ParserVersion = 14"
        echo "}"
    } > "$hsd_file" 2>/dev/null || { echo "Error: Cannot create $hsd_file at $(date)" | tee -a "$job_log" "$LOG_FILE"; return 1; }

    # Run DFTB+ in job directory
    if cd "$job_dir" 2>/dev/null; then
        timeout 10m dftb+ "$hsd_file" > "$base_name.out" 2>&1
        if [ $? -eq 0 ]; then
            echo "Successfully completed DFTB+ for $xyz_file at $(date)" | tee -a "$job_log" "$LOG_FILE"
        else
            echo "Error running DFTB+ for $xyz_file at $(date)" | tee -a "$job_log" "$LOG_FILE"
            echo "DFTB+ output for $xyz_file:" >> "$LOG_FILE"
            cat "$base_name.out" | tee -a "$job_log" "$LOG_FILE"
        fi
        cd - > /dev/null
    else
        echo "Error: Cannot cd to $job_dir at $(date)" | tee -a "$job_log" "$LOG_FILE"
        return 1
    fi
}

# Export functions and variables for parallel
export -f process_xyz
export -f get_angular_momentum
export XYZ_DIR MOLECULES_DIR SK_DIR LOG_FILE

# Find all .xyz files
mapfile -t xyz_files < <(find "$XYZ_DIR" -maxdepth 1 -type f -name "*.xyz")
if [ ${#xyz_files[@]} -eq 0 ]; then
    echo "Error: No .xyz files found in $XYZ_DIR at $(date)" | tee -a "$LOG_FILE"
    exit 1
fi
echo "Found ${#xyz_files[@]} .xyz files to process at $(date)" | tee -a "$LOG_FILE"

# Check for GNU parallel
if command -v parallel > /dev/null; then
    echo "Using GNU parallel with 20 jobs at $(date)" | tee -a "$LOG_FILE"
    parallel --env process_xyz,get_angular_momentum,XYZ_DIR,MOLECULES_DIR,SK_DIR,LOG_FILE -j 20 --eta process_xyz ::: "${xyz_files[@]}"
else
    echo "GNU parallel not found, using background processes with 20 jobs at $(date)" | tee -a "$LOG_FILE"
    max_jobs=20
    job_count=0
    for xyz_file in "${xyz_files[@]}"; do
        process_xyz "$xyz_file" &
        ((job_count++))
        if [ $job_count -ge $max_jobs ]; then
            wait -n
            ((job_count--))
        fi
    done
    wait
fi

echo "All DFTB+ jobs completed at $(date)" | tee -a "$LOG_FILE"