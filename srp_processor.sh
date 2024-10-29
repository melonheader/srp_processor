#! /usr/env/bash
# --------------------------------------------------------------------------------
# {INFO}
# Customisable processor-skeleton for SRP-type sequencing runs

# --------------------------------------------------------------------------------
# {SETUP}
OUTDIR=""
PROCESSOR=""
CLIPTYP=""
AUXDIR="auxiliary"
NCORES=4
HELP() {
    echo "Usage: $0 -i <input_path> [-o <output_directory>] [-p <processor>] [-a <auxiliary_directory>] [-g <genome_prefix>]"
    exit 1
}
while getopts "i:o:p:g:c:a:n:" opt; do
    case "$opt" in
        i) MINPUT="$OPTARG" ;;         # Path to SraRunTable.txt or its directory
        o) OUTDIR="$OPTARG" ;;         # Output directory
        p) PROCESSOR="$OPTARG" ;;      # Processor
        g) GPREFIX="$OPTARG" ;;        # Genome prefix
        c) CLIPTYP="$OPTARG" ;;        # i or e (iCLIP or eCLIP)
        a) AUXDIR="$OPTARG" ;;         # Auxiliary directory (default: .)
        n) NCORES="$OPTARG" ;;         # N-cores (default: 4)
        *) HELP ;;                     # Show usage if an unknown option is provided
    esac
done
if [ $OPTIND -eq 1 ]; then
    HELP
fi
##
if [ "$(basename "$MINPUT")" = "SraRunTable.txt" ]; then
    SRT_PATH="$MINPUT"
    MINPUT=$(dirname "$MINPUT")
    if [ -z "$OUTDIR" ]; then
        OUTDIR="$MINPUT"
    fi
elif [ -f "$MINPUT/SraRunTable.txt" ]; then 
    SRT_PATH="$MINPUT/SraRunTable.txt"
    if [ -z "$OUTDIR" ]; then
        OUTDIR=$MINPUT
    fi
else
    echo "Error: SraRunTable.txt is not found at the specified location" >&2
    exit 1
fi
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi
# --------------------------------------------------------------------------------
# {MAIN}
mapfile -t RUN_IDS < <(python3 "$AUXDIR"/parse_srt.py -s "$SRT_PATH" -o "$OUTDIR")
MIS_RUNS=()
for RUN in "${RUN_IDS[@]}"; do
    if ! ls "$MINPUT/$RUN"* >/dev/null 2>&1; then
        MIS_RUNS+=("$RUN")
    fi
done
if [ ${#MIS_RUNS[@]} -gt 0 ]; then
    echo "Runs without corresponding files deteced:"
    for RUN in "${MIS_RUNS[@]}"; do
        echo "$RUN; downloading..."
        fasterq-dump "$RUN" --split-3 --progress --outdir -M "$MINPUT"
        ## For complementarity with cutadapt make sure to flag --split-3 instead of --split-files https://github.com/marcelm/cutadapt/issues/197
        pigz -p "$NCORES" "$RUN"*
    done
fi
MATE_1=()
MATE_2=()
for RUN in "${RUN_IDS[@]}"; do
    MATES=()
    for MATE in "$MINPUT"/*; do
        if [[ -f "$MATE" && "$(basename "$MATE")" == "$RUN"* ]]; then
            MATES+=("$(basename "$MATE")")
        fi
    done
    for MATE in "${MATES[@]}"; do
        if [[ "$MATE" == *_1.* || "$MATE" == *R1.* ]]; then
            MATE_1+=("${MATE%%.fa*}")
        elif [[ "$MATE" == *_2.* || "$MATE" == *R2.* ]]; then
            MATE_2+=("${MATE%%.fa*}")
        fi
    done
done
if [ ${#MATE_1[@]} -gt 1 ] && [ ${#MATE_2[@]} -lt 1 ]; then
    echo "No files for second mates detected. Considering the data is single-end."
elif [ ${#MATE_1[@]} -gt 1 ] && [ ${#MATE_2[@]} -gt 1 ] && [ ${#MATE_1[@]} !=  ${#MATE_2[@]} ]; then
    echo "Not all run have both mates provided. Check the input data."
fi
# --------------------------------------------------------------------------------
# {EXECUTE}
# pre-process
bash "$AUXDIR"/"$PROCESSOR"_preprocess.sh \
    --mate1 "${MATE_1[@]}" --mate2 "${MATE_2[@]}" \
    --outdir "$OUTDIR" --indir "$MINPUT" \
    --auxdir "$AUXDIR" --ncores "$NCORES"; wait
# map & postprocess
bash "$AUXDIR"/"$PROCESSOR"_map.sh \
    --mate1 "${MATE_1[@]}" --mate2 "${MATE_2[@]}" \
    --outdir "$OUTDIR" --indir "$MINPUT" \
    --genome "$GPREFIX" --cliptyp "$CLIPTYP" --auxdir "$AUXDIR" --ncores "$NCORES"; wait