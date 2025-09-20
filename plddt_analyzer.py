#!/usr/bin/env python3

#how to run: python3 plddt_analyzer.py -i /work/10391/kaileeshlipak/ls6/AF3/models_0918/export_cif/ -o /work/10391/kaileeshlipak/ls6/AF3/analysis_results/

import csv
import statistics
from collections import defaultdict

def main():
    # Read the detailed CSV file
    data = defaultdict(lambda: defaultdict(dict))
    
    with open('/work/10391/kaileeshlipak/ls6/AF3/analysis_results/plddt_analysis_detailed.csv', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            filename = row['filename']
            median_plddt = float(row['median_plddt'])
            
            # Parse the filename to extract run code, seq, and seed
            if 'run_' in filename:
                parts = filename.split('_')
                run_code = parts[1]  # e.g., '2625777'
                seq_num = parts[2]   # e.g., 'seq1'
                
                # Extract seed info
                if 'seed-' in filename:
                    seed_part = filename.split('seed-')[1].split('_')[0]
                    sample_part = filename.split('sample-')[1].split('_')[0] if 'sample-' in filename else '0'
                    seed_key = f's{sample_part}'
                elif 'main' in filename:
                    seed_key = 'main'
                else:
                    seed_key = 'model'
                
                data[run_code][seq_num][seed_key] = median_plddt
    
    # Print the matrix
    print(f'{"Run":<8} {"Seq":<5} | Seed Models (median pLDDT)')
    print('-' * 80)
    
    for run_code in sorted(data.keys()):
        for seq_num in sorted(data[run_code].keys()):
            seeds = data[run_code][seq_num]
            
            # Create seed columns
            seed_strs = []
            all_values = []
            
            for seed_key in sorted(seeds.keys()):
                plddt = seeds[seed_key]
                seed_strs.append(f'{seed_key}:{plddt:.1f}')
                all_values.append(plddt)
            
            # Calculate consistency
            if len(all_values) > 1:
                std = statistics.stdev(all_values)
                avg = statistics.mean(all_values)
                consistency = f'avg:{avg:.1f},std:{std:.1f}'
            else:
                consistency = f'single:{all_values[0]:.1f}'
            
            print(f'{run_code:<8} {seq_num:<5} | {" | ".join(seed_strs)} ({consistency})')

if __name__ == "__main__":
    main()
