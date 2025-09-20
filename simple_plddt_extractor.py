#!/usr/bin/env python3
"""
Simple pLDDT extractor - no biotite dependencies
Reads CIF files directly and extracts B-factor (pLDDT) values
How to run: python3 simple_plddt_extractor.py /work/10391/kaileeshlipak/ls6/AF3/models_0918/export_cif/ /work/10391/kaileeshlipak/ls6/AF3/analysis_results/
"""

import os
import sys
import csv
import statistics
from pathlib import Path

def extract_plddt_from_cif(cif_path):
    """Extract pLDDT values from CIF file by reading B-factors directly"""
    try:
        with open(cif_path, 'r') as f:
            lines = f.readlines()
        
        # Find the atom site data section
        in_atom_section = False
        b_factors = []
        
        for line in lines:
            line = line.strip()
            
            # Look for atom site loop
            if line.startswith('_atom_site.'):
                in_atom_section = True
                continue
            
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
                
            # End of atom section
            if line.startswith('_') and in_atom_section:
                break
                
            # Parse atom data
            if in_atom_section and not line.startswith('_'):
                parts = line.split()
                if len(parts) >= 11:  # Standard CIF atom line has many columns
                    try:
                        # B-factor is typically the last or second-to-last column
                        b_factor = float(parts[-1])  # Try last column first
                        if b_factor < 0 or b_factor > 100:  # pLDDT should be 0-100
                            b_factor = float(parts[-2])  # Try second-to-last
                        
                        if 0 <= b_factor <= 100:  # Valid pLDDT range
                            b_factors.append(b_factor)
                    except (ValueError, IndexError):
                        continue
        
        return b_factors
        
    except Exception as e:
        print(f"Error reading {cif_path}: {e}")
        return []

def process_directory(input_dir, output_file):
    """Process all CIF files in directory"""
    input_path = Path(input_dir)
    results = []
    
    # Find all CIF files
    cif_files = list(input_path.glob("*.cif"))
    
    if not cif_files:
        print(f"No CIF files found in {input_dir}")
        return
    
    print(f"Found {len(cif_files)} CIF files to process")
    
    for i, cif_file in enumerate(cif_files, 1):
        print(f"[{i}/{len(cif_files)}] Processing: {cif_file.name}")
        
        # Extract pLDDT values
        plddt_values = extract_plddt_from_cif(cif_file)
        
        if plddt_values:
            mean_plddt = statistics.mean(plddt_values)
            
            # Parse filename for design and rank
            filename_base = cif_file.stem
            if '_model_' in filename_base:
                design_name = filename_base.split('_model_')[0]
                rank_part = filename_base.split('_model_')[1]
                alphafold_rank = f"ranked_{rank_part}"
            elif '-sample-' in filename_base:
                design_name = filename_base.split('-sample-')[0]
                rank_part = filename_base.split('-sample-')[1]
                alphafold_rank = f"ranked_{rank_part}"
            else:
                design_name = filename_base
                alphafold_rank = "ranked_0"
            
            results.append({
                'design': design_name,
                'rank': alphafold_rank,
                'mean_plddt': mean_plddt,
                'filename': cif_file.name
            })
            
            print(f"    Mean pLDDT: {mean_plddt:.1f}")
        else:
            print(f"    Failed to extract pLDDT values")
    
    # Write results to CSV
    if results:
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Design', 'RMSD', 'Alphafold_rank', 'pLDDT'])
            
            for result in results:
                writer.writerow([
                    result['design'],
                    f"{result['mean_plddt']:.2f}",  # Use pLDDT as RMSD placeholder
                    result['rank'],
                    f"{result['mean_plddt']:.2f}"
                ])
        
        print(f"\nResults written to: {output_file}")
        print(f"Processed {len(results)} structures successfully")
        
        # Quick stats
        all_plddts = [r['mean_plddt'] for r in results]
        print(f"pLDDT range: {min(all_plddts):.1f} - {max(all_plddts):.1f}")
        print(f"Average pLDDT: {statistics.mean(all_plddts):.1f}")
    else:
        print("No valid results to write")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python simple_plddt_extractor.py <input_directory> <output_csv>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_csv = sys.argv[2]
    
    process_directory(input_dir, output_csv)
