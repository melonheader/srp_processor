#!/usr/bin/env bash
# --------------------------------------------------------------------------------
# SRP Processor: Customizable pipeline for processing SRP sequencing runs
# --------------------------------------------------------------------------------

set -euo pipefail

# Default values
OUTDIR=""
PROCESSOR=""
CLIPTYP=""
AUXDIR="auxiliary"
GPREFIX=""
NCORES=4

usage() {
    cat <<EOF
Usage: $0 -i <input> -p <processor> -g <genome_prefix> [options]
Options:
  -i <input>            Input can be:
                         - Path to SraRunTable.txt or its directory
                         - File containing SRA accessions (one per line)
                         - Space-separated SRA accessions (e.g., "SRR123 SRR456")
                         - Path/glob pattern for FASTQ files (e.g., "data/*.fastq.gz")
  -p <processor>        Processor name (used to find auxiliary scripts; mandatory)
  -g <genome_prefix>    Genome prefix for mapping (mandatory)
  -c <clip type>        CLIP type: i (iCLIP) or e (eCLIP) [default: empty]
  -o <output_directory> Output directory (default: same as input directory)
  -a <auxiliary_dir>    Auxiliary directory (default: "auxiliary")
  -n <ncores>           Number of cores (default: 4)
EOF
    exit 1
}

# Parse command-line options
while getopts "i:o:p:g:c:a:n:" opt; do
    case "$opt" in
        i) MINPUT="$OPTARG" ;; # Input (various formats supported)
        o) OUTDIR="$OPTARG" ;; # Output directory
        p) PROCESSOR="$OPTARG" ;; # Processor name
        g) GPREFIX="$OPTARG" ;; # Genome prefix
        c) CLIPTYP="$OPTARG" ;; # CLIP type (i or e)
        a) AUXDIR="$OPTARG" ;; # Auxiliary directory
        n) NCORES="$OPTARG" ;; # Number of cores
        *) usage ;;
    esac
done

# Check mandatory arguments
if [ -z "${MINPUT:-}" ] || [ -z "${PROCESSOR:-}" ] || [ -z "${GPREFIX:-}" ]; then
    usage
fi

# Function to parse different input types
parse_input() {
    local input="$1"
    local outdir="$2"
    
    # Initialize empty array for run IDs
    RUN_IDS=()
    unset FASTQ_FILES
    
    # Case 1: Input is SraRunTable.txt
    if [[ "$input" == *"SraRunTable.txt" ]] && [[ -f "$input" ]]; then
        echo "Detected SraRunTable.txt input"
        SRT_PATH="$input"
        MINPUT_DIR=$(dirname "$input")
        # Use existing parser
        mapfile -t RUN_IDS < <(python3 "$AUXDIR/parse_srt.py" -s "$SRT_PATH" -o "$outdir")
        return 0
    
    # Case 2: Input is a directory containing SraRunTable.txt
    elif [[ -d "$input" ]] && [[ -f "$input/SraRunTable.txt" ]]; then
        echo "Detected directory with SraRunTable.txt"
        SRT_PATH="$input/SraRunTable.txt"
        MINPUT_DIR="$input"
        # Use existing parser
        mapfile -t RUN_IDS < <(python3 "$AUXDIR/parse_srt.py" -s "$SRT_PATH" -o "$outdir")
        return 0
    
    # Case 3: Input is a file with SRA accessions (one per line)
    elif [[ -f "$input" ]] && grep -q -E "^SRR|^ERR|^DRR" "$input"; then
        echo "Detected file with SRA accessions"
        MINPUT_DIR=$(dirname "$input")
        # Read accessions from file
        mapfile -t RUN_IDS < <(grep -E "^SRR|^ERR|^DRR" "$input" | tr -d '\r')
        return 0
    
    # Case 4: Input is a space-separated list of SRA accessions
    elif [[ "$input" =~ ^(SRR|ERR|DRR)[0-9]+(\ +(SRR|ERR|DRR)[0-9]+)*$ ]]; then
        echo "Detected space-separated SRA accessions"
        MINPUT_DIR="$PWD"
        # Split the input string into array
        read -ra RUN_IDS <<< "$input"
        return 0
    
    # Case 5: Input is a glob pattern for FASTQ files
    elif compgen -G "$input" > /dev/null; then
        echo "Detected FASTQ file pattern"
        # Get first file's directory for default MINPUT_DIR
        MINPUT_DIR=$(dirname "$(echo "$input" | cut -d' ' -f1)")
        # Get the file paths using globbing
        FASTQ_FILES=($input)
        echo "Found ${#FASTQ_FILES[@]} FASTQ files"
        
        # Extract run IDs from filenames if possible
        for file in "${FASTQ_FILES[@]}"; do
            base=$(basename "$file")
            # Try to extract SRA ID from filename
            if [[ "$base" =~ ^(SRR|ERR|DRR)[0-9]+ ]]; then
                id=${BASH_REMATCH[0]}
                # Only add if not already in array
                if [[ ! " ${RUN_IDS[*]} " =~ " $id " ]]; then
                    RUN_IDS+=("$id")
                fi
            fi
        done
        return 0
    
    else
        echo "Error: Could not parse input format: $input" >&2
        return 1
    fi
}

# Function: Download missing runs
download_missing_runs() {
    # Skip download if we're using direct FASTQ files
    if [[ -n "${FASTQ_FILES:-}" ]]; then
        echo "Using provided FASTQ files directly, skipping download"
        return 0
    fi
    
    local missing=()
    for RUN in "${RUN_IDS[@]}"; do
        # Check if any file starting with the run ID exists in MINPUT_DIR
        if ! compgen -G "$MINPUT_DIR/${RUN}*" > /dev/null; then
            missing+=("$RUN")
        fi
    done

    if [ "${#missing[@]}" -gt 0 ]; then
        echo "The following runs are missing and will be downloaded:"
        for RUN in "${missing[@]}"; do
            echo "Downloading $RUN..."
            # Download run
            fasterq-dump "$RUN" --split-3 --progress --outdir "$MINPUT_DIR"
            # Compress FASTQ files
            pigz -p "$NCORES" "$MINPUT_DIR/${RUN}"*.fastq
        done
    fi
}

# Function: Detect mate pairs
detect_mate_pairs() {
    MATE_1=()
    MATE_2=()
    
    # If we have direct FASTQ files
    if [[ -n "${FASTQ_FILES:-}" ]]; then
        for file in "${FASTQ_FILES[@]}"; do
            fname=$(basename "$file")
            if [[ "$fname" =~ (_1\.|R1\.|_1$|R1$) ]]; then
                MATE_1+=("$fname")
            elif [[ "$fname" =~ (_2\.|R2\.|_2$|R2$) ]]; then
                MATE_2+=("$fname")
            else
                # Assume single-end if no mate pattern found
                MATE_1+=("$fname")
            fi
        done
    else
        # Original functionality for SRA accessions
        local run mfile
        for run in "${RUN_IDS[@]}"; do
            # Use globbing to detect mate files for the run
            for mfile in "$MINPUT_DIR"/${run}*; do
                if [[ -f "$mfile" ]]; then
                    fname=$(basename "$mfile")
                    if [[ "$fname" =~ (_1\.|R1\.) ]]; then
                        MATE_1+=("$fname")
                    elif [[ "$fname" =~ (_2\.|R2\.) ]]; then
                        MATE_2+=("$fname")
                    fi
                fi
            done
        done
    fi

    # Provide feedback about mate detection
    if [ "${#MATE_1[@]}" -gt 0 ] && [ "${#MATE_2[@]}" -eq 0 ]; then
        echo "No second mate files detected; treating data as single-end."
    elif [ "${#MATE_1[@]}" -gt 0 ] && [ "${#MATE_2[@]}" -gt 0 ] && [ "${#MATE_1[@]}" -ne "${#MATE_2[@]}" ]; then
        echo "Warning: Mismatch in the number of first and second mate files. Check input data."
    fi
}

# --- MAIN EXECUTION ---
# Parse input and determine RUN_IDS
parse_input "$MINPUT" "${OUTDIR:-}" || usage

# Set OUTDIR if not specified
OUTDIR=${OUTDIR:-"$MINPUT_DIR"}
mkdir -p "$OUTDIR"

download_missing_runs
detect_mate_pairs

# Pre-process: check if the preprocess script exists and run it
PREPROCESS_SCRIPT="$AUXDIR/${PROCESSOR}_preprocess.sh"
if [ ! -x "$PREPROCESS_SCRIPT" ]; then
    echo "Error: Preprocess script '$PREPROCESS_SCRIPT' not found or not executable." >&2
    exit 1
fi

bash "$PREPROCESS_SCRIPT" \
    --mate1 "${MATE_1[@]}" --mate2 "${MATE_2[@]}" \
    --outdir "$OUTDIR" --indir "$MINPUT_DIR" \
    --auxdir "$AUXDIR" --ncores "$NCORES"

# Map & postprocess: check if the mapping script exists and run it
MAP_SCRIPT="$AUXDIR/${PROCESSOR}_map.sh"
if [ ! -x "$MAP_SCRIPT" ]; then
    echo "Error: Mapping script '$MAP_SCRIPT' not found or not executable." >&2
    exit 1
fi

bash "$MAP_SCRIPT" \
    --mate1 "${MATE_1[@]}" --mate2 "${MATE_2[@]}" \
    --outdir "$OUTDIR" --indir "$MINPUT_DIR" \
    --genome "$GPREFIX" --cliptyp "$CLIPTYP" --auxdir "$AUXDIR" --ncores "$NCORES"

echo "SRP processing completed successfully."