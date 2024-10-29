#! /usr/env/bash
# --------------------------------------------------------------------------------
# {INFO}
# following the guidlines reported at https://pureclip.readthedocs.io/en/latest/GettingStarted/preprocessing.html
# --------------------------------------------------------------------------------
# {SETUP}
set -e
FMATES=()
DIROUT=""
DIRIN=""
DIRAUX=""
GPREFIX=""
CLIPTYP=""
NCORES=4
HELP() {
    echo "Usage: $0 --mate1 FILE1 FILE2 ... --mate2  FILE1 FILE2 ... --outdir OUTPUT_DIR --indir INPUT_DIR --genome GPREFIX --cliptyp CLIPTYP --auxdir AUX_DIR"
    exit 1
}
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
        --ncores)
            NCORES="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            HELP
            ;;
    esac
done
if [[ ${#FMATES[@]} -eq 0 ]] || [[ -z "$DIROUT" ]] || [[ -z "$DIRIN" ]] || [[ -z "$DIRAUX" ]] || [[ -z "$GPREFIX" ]]; then
    HELP
fi
# --------------------------------------------------------------------------------
# {MAIN}
if [ ! -d "$DIROUT/mapped" ]; then
    mkdir -p "$DIROUT"/mapped
fi
DIROUT="$DIROUT"/mapped
if [ -d "$DIRIN/processed" ]; then
    DIRIN="$DIRIN/processed"
fi
GIDXLOC="$DIRAUX"/Annotations/Gencode/Indices/bowtie2/"$GPREFIX"
# -------------------------------------
# Function definitions
TH2MAP() {
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
    if [ -n "$SMATE" ]; then
        if [ ! -f "$DIROUT/$ROOT.bam" ]; then
            echo mapping "$FMATE" "$SMATE"
            tophat2 -p "$NCORES" \
                --library-type fr-secondstrand \
                --min-intron-length 20 \
                --max-intron-length 1000000 \
                --read-mismatches 2 \
                --read-edit-dist 2 \
                --min-anchor-length 8 \
                --max-multihits 1 \
                --segment-mismatches 3 \
                --segment-length 25 \
                --no-coverage-search \
                --output-dir "$DIROUT"/  \
                "$GIDX" \
                "$DIRIN"/"$FMATE"_trimmed_bc.fastq.gz "$DIRIN/$SMATE"_trimmed_bc.fastq.gz
            wait
            mv "$DIROUT"/accepted_hits.bam "$DIROUT"/"$ROOT".bam
            mv "$DIROUT"/unmapped.bam "$DIROUT"/"$ROOT".bam
            for REP in deletions.bed insertions.bed junctions.bed align_summary.txt prep_reads.info logs; do
                if [ ! -d "$DIROUT"/"$ROOT" ]; then
                    mkdir "$DIROUT"/"$ROOT"
                fi
                mv "$DIROUT"/"$REP" "$DIROUT"/"$ROOT"/"$REP"
            done
        elif [ ! -f "$DIROUT"/"$ROOT"_filtered.bam.bai ]; then
            echo filtering and indexing "$ROOT".bam
            SAMFILT --run "$ROOT" --out "$DIROUT"
        else
            echo Output bam file detected. Skipping "$FMATE" "$SMATE"
            return 0
        fi
    else
        if [ ! -f "$DIROUT/$ROOT.bam" ]; then
            echo mapping "$FMATE"
            tophat2 -p "$NCORES" \
                --library-type fr-secondstrand \
                --min-intron-length 20 \
                --max-intron-length 1000000 \
                --read-mismatches 2 \
                --read-edit-dist 2 \
                --min-anchor-length 8 \
                --max-multihits 1 \
                --segment-mismatches 3 \
                --segment-length 20 \
                --no-coverage-search \
                --output-dir "$DIROUT"/  \
                "$GIDX" \
                "$DIRIN"/"$FMATE"_trimmed_bc.fastq.gz
            wait
            mv "$DIROUT"/accepted_hits.bam "$DIROUT"/"$ROOT".bam
            mv "$DIROUT"/unmapped.bam "$DIROUT"/"$ROOT"_unmapped.bam
            for REP in deletions.bed insertions.bed junctions.bed align_summary.txt prep_reads.info logs; do
                if [ ! -d "$DIROUT"/"$ROOT" ]; then
                    mkdir "$DIROUT"/"$ROOT"
                fi
                mv "$DIROUT"/"$REP" "$DIROUT"/"$ROOT"/"$REP"
            done
        elif [ ! -f "$DIROUT"/"$ROOT"_filtered.bam.bai ]; then
            echo filtering and indexing "$ROOT".bam
            SAMFILT --run "$ROOT" --out "$DIROUT"
        else
            echo Output bam file detected. Skipping "$FMATE"
            return 0
        fi
    fi
}
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
DEDUP() {
    local REPARR=() GROUP="" PAIRED="" DIROUT="" 
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
            for IDX in "${!REPARR[@]}"; do
                REPS+=("$DIROUT"/"${REPARR[$IDX]}"_"$GROUP"_rep$(("$IDX" + 1))_filtered_dedup.bam)
            done
            echo merging "${REPARR[@]}"
            samtools merge -f "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam "${REPS[@]}"
        elif [ ! -f "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam.bai ]; then
            echo indexing merged "${REPARR[@]}"
            samtools index "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam "${REPS[@]}"
        else
            echo merged and deduplicated bam detected, skipping step for "${REPARR[@]}"
        fi
        if [ "$CLIPTYP" == "e" ]; then
            samtools view -hb -f 130 "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam -o "$DIROUT"/"$GROUP"_merged.R2_filtered_dedup.bam
            if [ ! -f "$DIROUT"/"$GROUP"_merged_filtered.R2_dedup_pos.bw ] && [ ! -f "$DIROUT"/"$GROUP"_merged.R2_filtered_dedup_neg.bw ]; then
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
            samtools view -hb -f 66 "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam -o "$DIROUT"/"$GROUP"_merged.R1_filtered_dedup.bam
            if [ ! -f "$DIROUT"/"$GROUP"_merged_filtered.R1_dedup_pos.bw ] && [ ! -f "$DIROUT"/"$GROUP"_merged.R1_filtered_dedup_neg.bw ]; then
                echo generating bigwigs
                bamCoverage \
                -b "$DIROUT"/"$GROUP"_merged.R1_filtered_dedup.bam \
                --filterRNAstrand forward --normalizeUsing None \
                -o "$DIROUT"/"$GROUP"_merged.R1_filtered_dedup_pos.bw; wait
                bamCoverage \
                -b "$DIROUT"/"$GROUP"_merged.R12_filtered_dedup.bam \
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
            for IDX in "${!REPARR[@]}"; do
                REPS+=("$DIROUT"/"${REPARR[$IDX]}"_"$GROUP"_rep$(("$IDX" + 1))_filtered_dedup.bam)
            done
            echo merging "${REPARR[@]}"
            samtools merge -f "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam "${REPS[@]}"  
        elif [ ! -f "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam.bai ]; then
            echo indexing merged "${REPARR[@]}"
            samtools index "$DIROUT"/"$GROUP"_merged_filtered_dedup.bam "${REPS[@]}"
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
# --------------------------------------------------------------------------------
# {EXECUTE}
eval "$(conda shell.bash hook)"
conda activate clip_map
echo "-----------------------------------------------------------------------"
echo mapping with TopHat2
for IDX in "${!FMATES[@]}"; do
    echo "------------------------------------------------"
    TH2MAP --fm "${FMATES[$IDX]}" --sm "${SMATES[$IDX]}" --gi "$GIDXLOC" --o "$DIROUT" --nc "$NCORES"; wait
done
conda deactivate
conda activate clip_dedup
echo "-----------------------------------------------------------------------"
EXPCSV=$(dirname "$DIROUT")/INFO_EXP.csv
eval "$(python3 "$DIRAUX"/split_runs.py "$EXPCSV")"
echo postprocessing INPUT "${INPUT[@]}"
DEDUP --input_arr "${INPUT[@]}" --bam_group "INPUT" --out "$DIROUT"
echo postprocessing CLIP "${CLIP[@]}"
DEDUP --input_arr "${CLIP[@]}" --bam_group "CLIP" --out "$DIROUT"