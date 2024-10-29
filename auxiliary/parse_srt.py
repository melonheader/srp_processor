import csv
import os
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Summarize Cutadapt JSON report into a CSV file.')
    parser.add_argument('-s', '--sra-path', type=str, required=True, help='Path to the SraRunTable')
    parser.add_argument('-o', '--out-dir', type=str, required=True, help='Directory to write parsed entries from SraRunTable')
    return parser.parse_args()

def get_col_indices(header, col_names):
    indices = []
    header_lc = [col.lower() for col in header]
    for name in col_names:
        name_lc = name.lower()
        for i, col in enumerate(header_lc):
            if name_lc in col:
                indices.append(i)
                break
    return indices
def extract_cols(input_csv, output_csv, col_indices):
    if not col_indices:
        print(f"No valid columns found for {output_csv}, skipping extraction.")
        return []
    
    run_ids = []
    with open(input_csv, newline='') as csvfile, open(output_csv, 'w', newline='') as outcsv:
        reader = csv.reader(csvfile)
        writer = csv.writer(outcsv)
        for i, row in enumerate(reader):
            selected_row = [row[idx] for idx in col_indices]
            if i == 0 or (col_indices and row[col_indices[0]] != 'Run'):
                writer.writerow(selected_row)
            if i > 0 and col_indices[0] < len(row):
                run_ids.append(row[col_indices[0]])
    return run_ids
def parse_srt(input_csv, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    info_id_cols = ["Run", "BioProject", "BioSample", "SampleName", "GEOaccession", "SRAstudy", "Experiment"]
    info_exp_cols = ["Run", "AssayType", "Antibody", "LibrarySource", "LibraryLayout", "LibrarySelection", "Instrument", "Bases", "AvgSpotLen"]
    info_bio_cols = ["Run", "Organism", "strain", "genotype", "source_name", "cell_type", "cell_line"]
    with open(input_csv, newline='') as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader)
        info_id_indices = get_col_indices(header, info_id_cols)
        info_exp_indices = get_col_indices(header, info_exp_cols)
        info_bio_indices = get_col_indices(header, info_bio_cols)
    run_ids = extract_cols(input_csv, os.path.join(output_dir, "INFO_ID.csv"), info_id_indices)
    extract_cols(input_csv, os.path.join(output_dir, "INFO_EXP.csv"), info_exp_indices)
    extract_cols(input_csv, os.path.join(output_dir, "INFO_BIO.csv"), info_bio_indices)
    for run_id in run_ids:
        print(run_id)
if __name__ == "__main__":
    args = parse_args()
    parse_srt(args.sra_path, args.out_dir)