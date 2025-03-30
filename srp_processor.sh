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
Usage: $0 -i <input_path> -p <processor> -g <genome_prefix> [options]
Options:
  -i <input_path>       Path to SraRunTable.txt or its directory (mandatory)
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
        i) MINPUT="$OPTARG" ;; # Path to SraRunTable.txt or its directory
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

# Resolve input and output directories, and locate SraRunTable.txt
if [ "$(basename "$MINPUT")" = "SraRunTable.txt" ]; then
    SRT_PATH="$MINPUT"
    MINPUT_DIR=$(dirname "$MINPUT")
    OUTDIR=${OUTDIR:-"$MINPUT_DIR"}
elif [ -f "$MINPUT/SraRunTable.txt" ]; then 
    SRT_PATH="$MINPUT/SraRunTable.txt"
    MINPUT_DIR="$MINPUT"
    OUTDIR=${OUTDIR:-"$MINPUT_DIR"}
else
    echo "Error: SraRunTable.txt not found at specified location" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

# Function: Download missing runs
download_missing_runs() {
    # Get list of RUN IDs using the parse_srt.py script
    mapfile -t RUN_IDS < <(python3 "$AUXDIR/parse_srt.py" -s "$SRT_PATH" -o "$OUTDIR")
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
            # Correct the fasterq-dump command; adjust -M flag as needed.
            fasterq-dump "$RUN" --split-3 --progress --outdir "$MINPUT_DIR"
            # Compress FASTQ files using a more specific glob pattern (e.g., *.fastq)
            pigz -p "$NCORES" "$MINPUT_DIR/${RUN}"*.fastq
        done
    fi
}

# Function: Detect mate pairs
detect_mate_pairs() {
    local run mfile
    MATE_1=()
    MATE_2=()
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

    # Provide feedback about mate detection
    if [ "${#MATE_1[@]}" -gt 0 ] && [ "${#MATE_2[@]}" -eq 0 ]; then
        echo "No second mate files detected; treating data as single-end."
    elif [ "${#MATE_1[@]}" -gt 0 ] && [ "${#MATE_2[@]}" -gt 0 ] && [ "${#MATE_1[@]}" -ne "${#MATE_2[@]}" ]; then
        echo "Warning: Mismatch in the number of first and second mate files. Check input data."
    fi
}

# --- MAIN EXECUTION ---
download_missing_runs
# Re-read RUN_IDS in case new files were downloaded
mapfile -t RUN_IDS < <(python3 "$AUXDIR/parse_srt.py" -s "$SRT_PATH" -o "$OUTDIR")
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