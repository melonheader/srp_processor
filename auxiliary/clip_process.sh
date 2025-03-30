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
GPREFIX=""
CLIPTYP=""
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
        --genome)
            GPREFIX="$2"
            shift 2
            ;;
        --cliptyp)
            CLIPTYP="$2"
            shift 2
            ;;
        --auxdir)
            DIRAUX="$2"
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

# --------------------------------------------------------------------------------
# {MAIN}
if [ ! -d "$DIROUT/mapped" ]; then
    mkdir -p "$DIROUT"/mapped
fi
DIROUT="$DIROUT"/mapped
if [ -d "$DIRIN/processed" ]; then
    DIRIN="$DIRIN/processed"
fi
GIDXLOC="$DIRAUX"/Annotations/Gencode/Indices/star/"$GPREFIX"

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
# Map reads with STAR
# ----------------------------------------------------------------------------
STARMAP() {
    local FMATE="" SMATE="" GIDX="" DIROUT="" NCORES="4" 
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
            --gi)
                GIDX="$2"
                shift 2
                ;;
            --o)
                DIROUT="$2"
                shift 2
                ;;
            --nc)
                NCORES="$2"
                shift 2
                ;;
            *)
                echo "Unknown option: $1" >&2
                return 1
                ;;
        esac
    done
    # shellcheck disable=SC2001
    ROOT=$(echo "$FMATE" | sed 's/_\(R1\|1\).*//')
    # Create a temporary directory for STAR output
    STAR_TMP="$DIROUT"/"$ROOT"_star_tmp
    
    if [ -n "$SMATE" ]; then
        if [ ! -f "$DIROUT/$ROOT.bam" ]; then
            echo mapping "$FMATE" "$SMATE"
            
            # Create temporary directory if it doesn't exist
            if [ ! -d "$STAR_TMP" ]; then
                mkdir -p "$STAR_TMP"
            fi
            
            # Run STAR for paired-end reads
            STAR --runThreadN "$NCORES" \
                --genomeDir "$GIDX" \
                --readFilesIn "$DIRIN"/"$FMATE"_trimmed_bc.fastq.gz "$DIRIN"/"$SMATE"_trimmed_bc.fastq.gz \
                --readFilesCommand zcat \
                --outFileNamePrefix "$STAR_TMP"/ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattributes Standard \
                --outFilterMultimapNmax 1 \
                --outFilterMismatchNmax 2 \
                --alignIntronMin 20 \
                --alignIntronMax 1000000 \
                --outFilterScoreMinOverLread 0.2 \
                --outFilterMatchNminOverLread 0.2 \
                --limitBAMsortRAM 10000000000 \
                --quantMode GeneCounts \
                --twopassMode Basic
            
            wait
            # Move and rename the output files
            mv "$STAR_TMP"/Aligned.sortedByCoord.out.bam "$DIROUT"/"$ROOT".bam
            
            # Create a directory for additional STAR outputs if needed
            if [ ! -d "$DIROUT"/"$ROOT" ]; then
                mkdir "$DIROUT"/"$ROOT"
            fi
            
            # Move other STAR output files to the directory
            mv "$STAR_TMP"/* "$DIROUT"/"$ROOT"/
            
            # Remove temporary directory
            rmdir "$STAR_TMP"
            
        elif [ ! -f "$DIROUT"/"$ROOT"_filtered.bam.bai ]; then
            echo filtering and indexing "$ROOT".bam
            SAMFILT --run "$ROOT" --out "$DIROUT" --paired
        else
            echo Output bam file detected. Skipping "$FMATE" "$SMATE"
            return 0
        fi
    else
        if [ ! -f "$DIROUT/$ROOT.bam" ]; then
            echo mapping "$FMATE"
            
            # Create temporary directory if it doesn't exist
            if [ ! -d "$STAR_TMP" ]; then
                mkdir -p "$STAR_TMP"
            fi
            
            # Run STAR for single-end reads
            STAR --runThreadN "$NCORES" \
                --genomeDir "$GIDX" \
                --readFilesIn "$DIRIN"/"$FMATE"_trimmed_bc.fastq.gz \
                --readFilesCommand zcat \
                --outFileNamePrefix "$STAR_TMP"/ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattributes Standard \
                --outFilterMultimapNmax 1 \
                --outFilterMismatchNmax 2 \
                --alignIntronMin 20 \
                --alignIntronMax 1000000 \
                --outFilterScoreMinOverLread 0.2 \
                --outFilterMatchNminOverLread 0.2 \
                --limitBAMsortRAM 10000000000 \
                --quantMode GeneCounts \
                --twopassMode Basic
            
            wait
            # Move and rename the output files
            mv "$STAR_TMP"/Aligned.sortedByCoord.out.bam "$DIROUT"/"$ROOT".bam
            
            # Create a directory for additional STAR outputs if needed
            if [ ! -d "$DIROUT"/"$ROOT" ]; then
                mkdir "$DIROUT"/"$ROOT"
            fi
            
            # Move other STAR output files to the directory
            mv "$STAR_TMP"/* "$DIROUT"/"$ROOT"/
            
            # Remove temporary directory
            rmdir "$STAR_TMP"
            
        elif [ ! -f "$DIROUT"/"$ROOT"_filtered.bam.bai ]; then
            echo filtering and indexing "$ROOT".bam
            SAMFILT --run "$ROOT" --out "$DIROUT"
        else
            echo Output bam file detected. Skipping "$FMATE"
            return 0
        fi
    fi
}

# ----------------------------------------------------------------------------
# Filter to keep the reads aligneing only to the main chromosomes, disregard contigs
# ----------------------------------------------------------------------------
## Filter to keep the reads aligneing only to the main chromosomes, disregard contigs
## deduplicate. keep only the R2 for eCLIP and R1 for iCLIP
SAMFILT() {
    local ROOT="" DIROUT="" PAIRED=""
    while [[ "$#" -gt 0 ]]; do
        case "$1" in
            --run)
                ROOT="$2"
                shift 2
                ;;
            --out)
                DIROUT="$2"
                shift 2
                ;;
            --paired)
                PAIRED=1
                shift 1
                ;;
            *)
                echo "Unknown option: $1" >&2
                return 1
                ;;
        esac
    done
    if [ ! -f "$DIROUT"/"$ROOT".bam.bai ]; then
        samtools index "$DIROUT"/"$ROOT".bam
    fi
    if [ ! -f "$DIROUT"/"$ROOT"_filtered.bam ]; then
        if [ -n "$PAIRED" ]; then
            samtools view -hb -f 2 -o "$DIROUT"/"$ROOT"_filtered.bam "$DIROUT"/"$ROOT".bam \
                chr1:1 chr2:1 chr3:1 chr4:1 chr5:1 chr6:1 chr7:1 chr8:1 chr9:1 chr10:1 chr11:1 chr12:1 \
                chr13:1 chr14:1 chr15:1 chr16:1 chr17:1 chr18:1 chr19:1 chr20:1 chr21:1 chr22:1 chrX:1 chrY:1
        else
            samtools view -hb -o "$DIROUT"/"$ROOT"_filtered.bam "$DIROUT"/"$ROOT".bam \
                chr1:1 chr2:1 chr3:1 chr4:1 chr5:1 chr6:1 chr7:1 chr8:1 chr9:1 chr10:1 chr11:1 chr12:1 \
                chr13:1 chr14:1 chr15:1 chr16:1 chr17:1 chr18:1 chr19:1 chr20:1 chr21:1 chr22:1 chrX:1 chrY:1
        fi 
    fi
    if [ ! -f "$DIROUT"/"$ROOT"_filtered.bam.bai ]; then
        samtools index "$DIROUT"/"$ROOT"_filtered.bam
    fi
}

# ----------------------------------------------------------------------------
# deduplicate. keep only the R2 for eCLIP and R1 for iCLIP
# ----------------------------------------------------------------------------
DEDUP() {
    local REPARR=() GROUP="" PAIRED="" DIROUT="" REPS=()
    while [[ "$#" -gt 0 ]]; do
        case "$1" in
            --input_arr)
                shift
                while [[ "$1" != --* ]] && [[ $# -gt 0 ]]; do
                    REPARR+=("$1")
                    shift
                done
                ;;
            --bam_group)
                GROUP="$2"
                shift 2
                ;;
            --p)
                PAIRED=1
                shift 1
                ;;
            --out)
                DIROUT="$2"
                shift 2
                ;;
            *)
                echo "Unknown option: $1" >&2
                return 1
                ;;
        esac
    done
    if [ -n "$PAIRED" ]; then
        if [ ! -f "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam ]; then
            for IDX in "${!REPARR[@]}"; do
                echo deduplicating "${REPARR[$IDX]}"
                umi_tools dedup -I "$DIROUT"/"${REPARR[$IDX]}"_filtered.bam --paired \
                    -S "$DIROUT"/"${REPARR[$IDX]}"_"$GROUP"_rep$(("$IDX" + 1))_filtered_dedup.bam
            done
            REPS=()
            for IDX in "${!REPARR[@]}"; do
                REPS+=("$DIROUT"/"${REPARR[$IDX]}"_"$GROUP"_rep$(("$IDX" + 1))_filtered_dedup.bam)
            done
            echo merging "${REPARR[@]}"
            samtools merge -f "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam "${REPS[@]}"
        elif [ ! -f "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam.bai ]; then
            echo indexing merged "${REPARR[@]}"
            samtools index "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam
        else
            echo merged and deduplicated bam detected, skipping step for "${REPARR[@]}"
        fi
        if [ "$CLIPTYP" == "e" ]; then
            if [ ! -f "$DIROUT"/"$GROUP"_merged.R2_filtered_dedup.bam ]; then
                echo "Extracting R2 reads for eCLIP data"
                samtools view -hb -f 130 "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam -o "$DIROUT"/"$GROUP"_merged.R2_filtered_dedup.bam
            fi
            if [ ! -f "$DIROUT"/"$GROUP"_merged.R2_filtered_dedup_pos.bw ] && [ ! -f "$DIROUT"/"$GROUP"_merged.R2_filtered_dedup_neg.bw ]; then
                echo "Generating bigwigs for eCLIP data"
                bamCoverage \
                -b "$DIROUT"/"$GROUP"_merged.R2_filtered_dedup.bam \
                --filterRNAstrand forward --normalizeUsing None \
                -o "$DIROUT"/"$GROUP"_merged.R2_filtered_dedup_pos.bw; wait
                bamCoverage \
                -b "$DIROUT"/"$GROUP"_merged.R2_filtered_dedup.bam \
                --filterRNAstrand reverse --normalizeUsing None \
                -o "$DIROUT"/"$GROUP"_merged.R2_filtered_dedup_neg.bw
            fi
        elif [ "$CLIPTYP" == "i" ]; then
            if [ ! -f "$DIROUT"/"$GROUP"_merged.R1_filtered_dedup.bam ]; then
                echo "Extracting R1 reads for iCLIP data"
                samtools view -hb -f 66 "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam -o "$DIROUT"/"$GROUP"_merged.R1_filtered_dedup.bam
            fi
            if [ ! -f "$DIROUT"/"$GROUP"_merged.R1_filtered_dedup_pos.bw ] && [ ! -f "$DIROUT"/"$GROUP"_merged.R1_filtered_dedup_neg.bw ]; then
                echo "Generating bigwigs for iCLIP data"
                bamCoverage \
                -b "$DIROUT"/"$GROUP"_merged.R1_filtered_dedup.bam \
                --filterRNAstrand forward --normalizeUsing None \
                -o "$DIROUT"/"$GROUP"_merged.R1_filtered_dedup_pos.bw; wait
                bamCoverage \
                -b "$DIROUT"/"$GROUP"_merged.R1_filtered_dedup.bam \
                --filterRNAstrand reverse --normalizeUsing None \
                -o "$DIROUT"/"$GROUP"_merged.R1_filtered_dedup_neg.bw
            fi
        fi
    else 
        if [ ! -f "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam ]; then
            for IDX in "${!REPARR[@]}"; do
                echo deduplicating "${REPARR[$IDX]}"
                umi_tools dedup -I "$DIROUT"/"${REPARR[$IDX]}"_filtered.bam \
                    -S "$DIROUT"/"${REPARR[$IDX]}"_"$GROUP"_rep$(("$IDX" + 1))_filtered_dedup.bam \
                    -L "$DIROUT"/"${REPARR[$IDX]}"_"$GROUP"_rep$(("$IDX" + 1))_filtered.dedupLog \
                    -E "$DIROUT"/"${REPARR[$IDX]}"_"$GROUP"_rep$(("$IDX" + 1))_filtered.dedupErr
            done
            REPS=()
            for IDX in "${!REPARR[@]}"; do
                REPS+=("$DIROUT"/"${REPARR[$IDX]}"_"$GROUP"_rep$(("$IDX" + 1))_filtered_dedup.bam)
            done
            echo merging "${REPARR[@]}"
            samtools merge -f "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam "${REPS[@]}"  
        elif [ ! -f "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam.bai ]; then
            echo indexing merged "${REPARR[@]}"
            samtools index "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam
        else
            echo merged and deduplicated bam detected, skipping step for "${REPARR[@]}"
        fi
        if [ ! -f "$DIROUT"/"$GROUP"_merged_filtered_dedup_pos.bw ] && [ ! -f "$DIROUT"/"$GROUP"_merged_filtered_dedup_neg.bw ]; then
            echo generating bigwigs
            bamCoverage \
            -b "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam \
            --filterRNAstrand forward --normalizeUsing None \
            -o "$DIROUT"/"$GROUP"_merged_filtered_dedup_pos.bw; wait
            bamCoverage \
            -b "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam \
            --filterRNAstrand reverse --normalizeUsing None \
            -o "$DIROUT"/"$GROUP"_merged_filtered_dedup_neg.bw
        fi
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