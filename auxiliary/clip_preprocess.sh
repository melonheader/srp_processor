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
NCORES=4
HELP() {
    echo "Usage: $0 --mate1 FILE1 FILE2 ... --mate2  FILE1 FILE2 ... --outdir OUTPUT_DIR --indir INPUT_DIR --auxdir AUX_DIR"
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
if [[ ${#FMATES[@]} -eq 0 ]] || [[ -z "$DIROUT" ]] || [[ -z "$DIRIN" ]] || [[ -z "$DIRAUX" ]]; then
    HELP
fi
# --------------------------------------------------------------------------------
# {MAIN}
if [ ! -d "$DIROUT/processed" ]; then
    mkdir -p "$DIROUT"/processed
fi
DIROUT="$DIROUT"/processed
# -------------------------------------
# Function definitions
TRIMAD() {
    local FMATE="" SMATE="" FEXT="" DIRIN="" DIROUT=""
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
                echo "Unknown option: $1" >&2
                return 1
                ;;
        esac
    done
    # Check if FMATE was provided
    if [[ -z "$FMATE" ]]; then
        echo "Error: FMATE (--fm) at least on read file is required" >&2
        return 1
    fi
    # shellcheck disable=SC2001
    ROOT=$(echo "$FMATE" | sed 's/_\(R1\|1\).*//')
    FEXT=$(find "$DIRIN" -maxdepth 1 -type f -name "$FMATE*" | head -n 1)
    FEXT=".fa${FEXT##*.fa}"
    # -a --> 3' end of the first mate; -g --> 5' end of the first mate; -A --> 3' end of the second mate
    ## -e (--error-rate) maximum error rate; -O minlength overlap between read and adapter; --times Remove up to COUNT adapters from each read -m
    ## prepare a command for cutadapt
    if [[ -n "$SMATE" ]]; then
        if [ -f "$DIROUT/$FMATE"_trimmed.fastq.gz ] && [ -f "$DIROUT"/"$SMATE"_trimmed.fastq.gz ]; then
            echo "Skipping trimming $ROOT mates; Already preprocessed"
            return 0
        fi
        echo "running cutadapt on $FMATE $SMATE; first pass"
        cutadapt \
            --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 --cores "$NCORES" --report full --json "$DIROUT"/"$ROOT"_trimReport_1pass.txt \
            --discard-untrimmed \
            -a "NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" \
            -g "CTTCCGATCTACAAGTT" -g "CTTCCGATCTTGGTCCT" \
            -A "AACTTGTAGATCGGA" -A "AGGACCAAGATCGGA" -A "ACTTGTAGATCGGAA" -A "GGACCAAGATCGGAA" \
            -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA \
            -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC \
            -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC \
            -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT \
            -o "$DIROUT/$FMATE"_trimmed1.fastq -p "$DIROUT/$SMATE"_trimmed1.fastq \
            "$DIRIN"/"$FMATE$FEXT" "$DIRIN"/"$SMATE$FEXT" > /dev/null
        wait
        echo "running cutadapt on $FMATE $SMATE; second pass"
        cutadapt \
            --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 --cores "$NCORES" --report full --json "$DIROUT"/"$ROOT"_trimReport_2pass.txt \
            --discard-untrimmed \
            -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA \
            -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA \
            -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC \
            -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC \
            -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT \
            -o "$DIROUT/$FMATE"_trimmed.fastq -p "$DIROUT/$SMATE"_trimmed.fastq \
            "$DIROUT/$FMATE"_trimmed1.fastq "$DIROUT/$SMATE"_trimmed1.fastq > /dev/null
        wait
        rm "$DIROUT/$FMATE"_trimmed1.fastq "$DIROUT/$SMATE"_trimmed1.fastq
        pigz -p "$NCORES" "$DIROUT/$FMATE"_trimmed.fastq "$DIROUT/$SMATE"_trimmed.fastq
        #python3 "$DIRAUX"/parse_trimreport.py -j "$DIROUT"/"$ROOT"_trimReport_1pass.txt -c "$DIROUT"/"$ROOT"_trimReport_1pass.csv
        #python3 "$DIRAUX"/parse_trimreport.py -j "$DIROUT"/"$ROOT"_trimReport_2pass.txt -c "$DIROUT"/"$ROOT"_trimReport_2pass.csv
    else
        if [ -f "$DIROUT/$FMATE"_trimmed.fastq.gz ]; then
            echo "Skipping trimming $ROOT; Already preprocessed"
            return 0
        fi
        echo "running cutadapt on $FMATE"
        cutadapt \
            --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 --cores "$NCORES" --report full --json "$DIROUT"/"$ROOT"_trimReport.txt \
            --discard-untrimmed \
            -a "NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" \
            -g "CTTCCGATCTACAAGTT" -g "CTTCCGATCTTGGTCCT" \
            -o "$DIROUT/$FMATE"_trimmed.fastq \
            "$DIRIN"/"$FMATE$FEXT" > /dev/null
        wait
        pigz -p "$NCORES" "$DIROUT/$FMATE"_trimmed.fastq
        #python3 "$DIRAUX"/parse_trimreport.py -j "$DIROUT"/"$ROOT"_trimReport.txt -c "$DIROUT"/"$ROOT"_trimReport.csv
    fi
}
GETBARCODES() {
    local FMATE=$1
    local SMATE=$2
    # shellcheck disable=SC2001
    ROOT=$(echo "$FMATE" | sed 's/_\(R1\|1\).*//')
    if [[ -z "$SMATE" ]]; then
        if [ -f "$DIROUT"/"$FMATE"_trimmed_bc.fastq.gz ]; then
            echo "Skipping UMI barcoding step for $ROOT; UMI barcodes are already in read IDs"
            return 0
        fi
    else
        if [ -f "$DIROUT"/"$FMATE"_trimmed_bc.fastq.gz ] && [ -f "$DIROUT"/"$SMATE"_trimmed_bc.fastq.gz ]; then
            echo "Skipping UMI barcoding step for $ROOT mates; UMI barcodes are already in read IDs"
            return 0
        fi
    fi
    echo "moving UMI barcodes to the read ID..."
    if [[ -n "$SMATE" ]]; then
        for MATE in "$FMATE" "$SMATE"; do
            zcat "$DIROUT"/"$MATE"_trimmed.fastq.gz | awk -v l=10 \
                'BEGIN{OFS=FS=" "} substr($1, 1, 1) == "@" {print "@" substr($1, (l+3), 500) "_" substr($1, 2, l) " " $2 }; substr($1, 1, 1) != "@" {print}; ' |\
                pigz -p "$NCORES" > "$DIROUT"/"$MATE"_trimmed_bc.fastq.gz
        done
    else
            zcat "$DIROUT"/"$FMATE"_trimmed.fastq.gz | awk -v l=10 \
                'BEGIN{OFS=FS=" "} substr($1, 1, 1) == "@" {print "@" substr($1, (l+3), 500) "_" substr($1, 2, l) " " $2 }; substr($1, 1, 1) != "@" {print}; ' |\
                pigz -p "$NCORES" > "$DIROUT"/"$FMATE"_trimmed_bc.fastq.gz
    fi
}
# --------------------------------------------------------------------------------
# {EXECUTE}
eval "$(conda shell.bash hook)"
conda activate clip_preprocess
echo "-----------------------------------------------------------------------"
echo preprocessing
for IDX in "${!FMATES[@]}"; do
    echo "------------------------------------------------"
    TRIMAD --fm "${FMATES[$IDX]}" --sm "${SMATES[$IDX]}" --di "$DIRIN" --do "$DIROUT"; wait
    GETBARCODES "${FMATES[$IDX]}" "${SMATES[$IDX]}" ; wait
done