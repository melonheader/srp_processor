#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# CLIP-seq Preprocessing Script
# 
# Purpose: Prepare CLIP-seq fastq files for mapping by:
#   1. Trimming adapters using cutadapt
#   2. Extracting UMI barcodes and adding to read IDs
#
# Based on guidelines from: https://pureclip.readthedocs.io/en/latest/GettingStarted/preprocessing.html
# ----------------------------------------------------------------------------

# Fail on any error
set -e

# Initialize variables
FMATES=()
SMATES=()
DIROUT=""
DIRIN=""
DIRAUX=""
NCORES=4

# Help message
show_help() {
    echo "Usage: $0 --mate1 FILE1 FILE2 ... --mate2 FILE1 FILE2 ... --outdir OUTPUT_DIR --indir INPUT_DIR --auxdir AUX_DIR [--ncores CORES]"
    echo ""
    echo "Options:"
    echo "  --mate1       First mate (R1) fastq files"
    echo "  --mate2       Second mate (R2) fastq files (optional)"
    echo "  --outdir      Output directory"
    echo "  --indir       Input directory containing raw fastq files"
    echo "  --auxdir      Directory containing auxiliary scripts"
    echo "  --ncores      Number of CPU cores to use (default: 4)"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mate1)
            shift
            while [[ "$1" != --* ]] && [[ $# -gt 0 ]]; do
                FMATES+=("$1")
                shift
            done
            ;;
        --mate2)
            shift
            while [[ "$1" != --* ]] && [[ $# -gt 0 ]]; do
                SMATES+=("$1")
                shift
            done
            ;;
        --outdir)
            DIROUT="$2"
            shift 2
            ;;
        --indir)
            DIRIN="$2"
            shift 2
            ;;
        --auxdir)
            DIRAUX="$2"
            shift 2
            ;;
        --ncores)
            NCORES="$2"
            shift 2
            ;;
        --help|-h)
            show_help
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            ;;
    esac
done

# Validate required arguments
if [[ ${#FMATES[@]} -eq 0 ]] || [[ -z "$DIROUT" ]] || [[ -z "$DIRIN" ]] || [[ -z "$DIRAUX" ]]; then
    echo "Error: Missing required arguments."
    show_help
fi

# Create output directory
if [[ ! -d "$DIROUT/processed" ]]; then
    mkdir -p "$DIROUT/processed"
fi
DIROUT="$DIROUT/processed"

# Log execution
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

# ----------------------------------------------------------------------------
# Adapter trimming function
# Uses cutadapt to trim adapters from reads
# ----------------------------------------------------------------------------
trim_adapters() {
    local FMATE="" SMATE="" FEXT="" DIRIN="" DIROUT=""
    
    # Parse function arguments
    while [[ "$#" -gt 0 ]]; do
        case "$1" in
            --fm)
                FMATE="$2"
                shift 2
                ;;
            --sm)
                SMATE="$2"
                shift 2
                ;;
            --di)
                DIRIN="$2"
                shift 2
                ;;
            --do)
                DIROUT="$2"
                shift 2
                ;;
            *)
                echo "Error: Unknown option for trim_adapters: $1" >&2
                return 1
                ;;
        esac
    done
    
    # Validate required arguments
    if [[ -z "$FMATE" ]]; then
        echo "Error: FMATE (--fm) is required for trim_adapters" >&2
        return 1
    fi
    
    # Extract root sample name and file extension
    ROOT=$(echo "$FMATE" | sed 's/_\(R1\|1\).*//')
    FEXT=$(find "$DIRIN" -maxdepth 1 -type f -name "$FMATE*" | head -n 1)
    FEXT=".fa${FEXT##*.fa}"
    
    # Process paired-end reads
    if [[ -n "$SMATE" ]]; then
        # Skip if output already exists
        if [[ -f "$DIROUT/$FMATE"_trimmed.fastq.gz ]] && [[ -f "$DIROUT/$SMATE"_trimmed.fastq.gz ]]; then
            log "Skipping trimming $ROOT mates; Already preprocessed"
            return 0
        fi
        
        # First cutadapt pass
        log "Running cutadapt on $FMATE $SMATE (pass 1)"
        cutadapt \
            --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 \
            --cores "$NCORES" --report full --json "$DIROUT/${ROOT}_trimReport_1pass.txt" \
            --discard-untrimmed \
            -a "NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" \
            -g "CTTCCGATCTACAAGTT" -g "CTTCCGATCTTGGTCCT" \
            -A "AACTTGTAGATCGGA" -A "AGGACCAAGATCGGA" -A "ACTTGTAGATCGGAA" -A "GGACCAAGATCGGAA" \
            -A "CTTGTAGATCGGAAG" -A "GACCAAGATCGGAAG" -A "TTGTAGATCGGAAGA" -A "ACCAAGATCGGAAGA" \
            -A "TGTAGATCGGAAGAG" -A "CCAAGATCGGAAGAG" -A "GTAGATCGGAAGAGC" -A "CAAGATCGGAAGAGC" \
            -A "TAGATCGGAAGAGCG" -A "AAGATCGGAAGAGCG" -A "AGATCGGAAGAGCGT" -A "GATCGGAAGAGCGTC" \
            -A "ATCGGAAGAGCGTCG" -A "TCGGAAGAGCGTCGT" -A "CGGAAGAGCGTCGTG" -A "GGAAGAGCGTCGTGT" \
            -o "$DIROUT/${FMATE}_trimmed1.fastq" -p "$DIROUT/${SMATE}_trimmed1.fastq" \
            "$DIRIN/${FMATE}${FEXT}" "$DIRIN/${SMATE}${FEXT}" > /dev/null
        
        # Second cutadapt pass
        log "Running cutadapt on $FMATE $SMATE (pass 2)"
        cutadapt \
            --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 \
            --cores "$NCORES" --report full --json "$DIROUT/${ROOT}_trimReport_2pass.txt" \
            --discard-untrimmed \
            -A "AACTTGTAGATCGGA" -A "AGGACCAAGATCGGA" -A "ACTTGTAGATCGGAA" -A "GGACCAAGATCGGAA" \
            -A "CTTGTAGATCGGAAG" -A "GACCAAGATCGGAAG" -A "TTGTAGATCGGAAGA" -A "ACCAAGATCGGAAGA" \
            -A "TGTAGATCGGAAGAG" -A "CCAAGATCGGAAGAG" -A "GTAGATCGGAAGAGC" -A "CAAGATCGGAAGAGC" \
            -A "TAGATCGGAAGAGCG" -A "AAGATCGGAAGAGCG" -A "AGATCGGAAGAGCGT" -A "GATCGGAAGAGCGTC" \
            -A "ATCGGAAGAGCGTCG" -A "TCGGAAGAGCGTCGT" -A "CGGAAGAGCGTCGTG" -A "GGAAGAGCGTCGTGT" \
            -o "$DIROUT/${FMATE}_trimmed.fastq" -p "$DIROUT/${SMATE}_trimmed.fastq" \
            "$DIROUT/${FMATE}_trimmed1.fastq" "$DIROUT/${SMATE}_trimmed1.fastq" > /dev/null
        
        # Clean up intermediate files and compress output
        rm "$DIROUT/${FMATE}_trimmed1.fastq" "$DIROUT/${SMATE}_trimmed1.fastq"
        log "Compressing output files"
        pigz -p "$NCORES" "$DIROUT/${FMATE}_trimmed.fastq" "$DIROUT/${SMATE}_trimmed.fastq"
        
    # Process single-end reads
    else
        # Skip if output already exists
        if [[ -f "$DIROUT/${FMATE}_trimmed.fastq.gz" ]]; then
            log "Skipping trimming $ROOT; Already preprocessed"
            return 0
        fi
        
        log "Running cutadapt on $FMATE"
        cutadapt \
            --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 \
            --cores "$NCORES" --report full --json "$DIROUT/${ROOT}_trimReport.txt" \
            --discard-untrimmed \
            -a "NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" \
            -g "CTTCCGATCTACAAGTT" -g "CTTCCGATCTTGGTCCT" \
            -o "$DIROUT/${FMATE}_trimmed.fastq" \
            "$DIRIN/${FMATE}${FEXT}" > /dev/null
        
        # Compress output
        log "Compressing output file"
        pigz -p "$NCORES" "$DIROUT/${FMATE}_trimmed.fastq"
    fi
}

# ----------------------------------------------------------------------------
# Extract UMI barcodes and add to read IDs
# ----------------------------------------------------------------------------
extract_barcodes() {
    local FMATE="$1"
    local SMATE="$2"
    
    # Extract root sample name
    ROOT=$(echo "$FMATE" | sed 's/_\(R1\|1\).*//')
    
    # Check if files already processed
    if [[ -z "$SMATE" ]]; then
        if [[ -f "$DIROUT/${FMATE}_trimmed_bc.fastq.gz" ]]; then
            log "Skipping UMI extraction for $ROOT; UMI barcodes are already in read IDs"
            return 0
        fi
    else
        if [[ -f "$DIROUT/${FMATE}_trimmed_bc.fastq.gz" ]] && [[ -f "$DIROUT/${SMATE}_trimmed_bc.fastq.gz" ]]; then
            log "Skipping UMI extraction for $ROOT mates; UMI barcodes are already in read IDs"
            return 0
        fi
    fi
    
    log "Moving UMI barcodes to read IDs for $ROOT"
    
    # Process mates
    if [[ -n "$SMATE" ]]; then
        for MATE in "$FMATE" "$SMATE"; do
            log "Processing $MATE"
            zcat "$DIROUT/${MATE}_trimmed.fastq.gz" | awk -v l=10 \
                'BEGIN{OFS=FS=" "} 
                 substr($1, 1, 1) == "@" {print "@" substr($1, (l+3), 500) "_" substr($1, 2, l) " " $2}; 
                 substr($1, 1, 1) != "@" {print}' | \
                pigz -p "$NCORES" > "$DIROUT/${MATE}_trimmed_bc.fastq.gz"
        done
    else
        zcat "$DIROUT/${FMATE}_trimmed.fastq.gz" | awk -v l=10 \
            'BEGIN{OFS=FS=" "} 
             substr($1, 1, 1) == "@" {print "@" substr($1, (l+3), 500) "_" substr($1, 2, l) " " $2}; 
             substr($1, 1, 1) != "@" {print}' | \
            pigz -p "$NCORES" > "$DIROUT/${FMATE}_trimmed_bc.fastq.gz"
    fi
}

# ----------------------------------------------------------------------------
# Main execution
# ----------------------------------------------------------------------------
# Activate conda environment
log "Activating conda environment: clip_preprocess"
eval "$(conda shell.bash hook)"
conda activate clip_preprocess

log "Starting CLIP-seq preprocessing"
log "Found ${#FMATES[@]} samples to process"

# Process each sample
for IDX in "${!FMATES[@]}"; do
    log "Processing sample ${IDX+1}/${#FMATES[@]}: ${FMATES[$IDX]}"
    
    # Trim adapters
    trim_adapters --fm "${FMATES[$IDX]}" --sm "${SMATES[$IDX]}" --di "$DIRIN" --do "$DIROUT"
    
    # Extract UMI barcodes
    extract_barcodes "${FMATES[$IDX]}" "${SMATES[$IDX]}"
    
    log "Completed processing ${FMATES[$IDX]}"
done

log "Preprocessing completed successfully"