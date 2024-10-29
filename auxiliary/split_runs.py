import csv
import sys

input_runs = []
clip_runs = []

with open(sys.argv[1], newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if row['antibody'] == 'no antibody':
            input_runs.append(row['Run'])
        else:
            clip_runs.append(row['Run'])
print(f"INPUT=({' '.join(input_runs)})")
print(f"CLIP=({' '.join(clip_runs)})")