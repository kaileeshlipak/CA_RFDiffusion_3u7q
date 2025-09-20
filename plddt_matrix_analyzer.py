#!/usr/bin/env python3
import csv
import statistics
import argparse
import sys
from collections import defaultdict
from pathlib import Path

def analyze_plddt_matrix(csv_file_path):
    """
    Analyze pLDDT data from CSV and create a matrix view.
    
    Args:
        csv_file_path: Path to the CSV file containing pLDDT analysis results
    """
    # Read the detailed CSV file
    data = defaultdict(lambda: defaultdict(dict))
    
    if not Path(csv_file_path).exists():
        print(f"Error: CSV file not found: {csv_file_path}")
        return
    
    try:
        with open(csv_file_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                filename = row['filename']
                
                # Check if median_plddt column exists, fallback to mean_plddt
                if 'median_plddt' in row and row['median_plddt'] != 'NA':
                    plddt_value = float(row['median_plddt'])
                    plddt_type = "median"
                elif 'mean_plddt' in row and row['mean_plddt'] != 'NA':
                    plddt_value = float(row['mean_plddt'])
                    plddt_type = "mean"
                else:
                    print(f"Warning: No valid pLDDT data for {filename}")
                    continue
                
                # Parse the filename to extract run code, seq, and seed
                # Handle both formats: run_XXXX_seqX and ligandmpnn_XXXX_seqX
                if 'run_' in filename:
                    parts = filename.split('_')
                    run_code = parts[1]  # e.g., '2625777'
                    seq_num = parts[2]   # e.g., 'seq1'
                elif 'ligandmpnn_' in filename:
                    parts = filename.split('_')
                    run_code = parts[1]  # e.g., '2633259'
                    seq_num = parts[2]   # e.g., 'seq2'
                else:
                    print(f"Warning: Could not parse filename format: {filename}")
                    continue
                
                # Extract seed info
                if 'seed-' in filename:
                    seed_part = filename.split('seed-')[1].split('_')[0]
                    sample_part = filename.split('sample-')[1].split('_')[0] if 'sample-' in filename else '0'
                    seed_key = f's{sample_part}'
                elif 'main' in filename:
                    seed_key = 'main'
                else:
                    seed_key = 'model'
                
                data[run_code][seq_num][seed_key] = plddt_value
    
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return
    
    if not data:
        print("No data found in CSV file")
        return
    
    # Print the matrix
    print(f'{"Run":<8} {"Seq":<5} | Seed Models ({plddt_type} pLDDT)')
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

def main():
    parser = argparse.ArgumentParser(
        description="Analyze pLDDT results and create matrix view",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python plddt_matrix_analyzer.py /path/to/plddt_analysis.csv
    python plddt_matrix_analyzer.py -f /work/10391/kaileeshlipak/ls6/AF3/analysis_results/0920morn/plddt_analysis.csv
        """
    )
    
    parser.add_argument('csv_file', nargs='?', 
                       help='Path to CSV file containing pLDDT analysis results')
    parser.add_argument('-f', '--file', 
                       help='Path to CSV file (alternative to positional argument)')
    
    args = parser.parse_args()
    
    # Determine which CSV file to use
    csv_file = args.csv_file or args.file
    
    if not csv_file:
        print("Error: Please specify a CSV file path")
        print("Usage: python plddt_matrix_analyzer.py <csv_file>")
        print("   or: python plddt_matrix_analyzer.py -f <csv_file>")
        sys.exit(1)
    
    analyze_plddt_matrix(csv_file)

if __name__ == "__main__":
    main()
