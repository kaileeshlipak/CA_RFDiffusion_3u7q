#!/usr/bin/env python3
"""
Extract specific AF3 runs and sequence numbers, collecting all seed variants
Usage: python extract_specific_runs.py <af3_output_base_dir> <output_dir>

Edit the TARGETS list below to specify which runs and sequences you want
"""

import os
import sys
import shutil
from pathlib import Path

# EDIT THIS LIST: Add your (run_code, sequence_number) pairs
TARGETS = [
    # Parsed from your list:
    ("2624679", 1),
    ("2624684", 1),
    ("2625215", 3),
    ("2625215", 4),
    ("2625216", 3),
    ("2625216", 4),
    ("2625218", 3),
    ("2627566", 1),
    ("2627566", 3),
    ("2627569", 2),
    ("2627570", 3),
    ("2628842", 2),
    ("2628848", 1),
    ("2628848", 2),
    ("2628916", 1),
    ("2628916", 2),
    ("2630297", 2),
    ("2619816", 1),
    ("2619708", 1),
    ("2620033", 2),
    ("2620989", 2),
    ("2624509", 2),
    ("2624626", 2),
    ("2625777", 1),
    ("2625778", 1),
    ("2625778", 2),
]

def find_structure_path(af3_base, run_code, seq_num):
    """
    Search for structure path in various possible locations
    """
    structure_name = f"ligandmpnn_{run_code}_seq{seq_num}_0"
    
    # Method 1: Look for the directory structure
    for path in af3_base.rglob(structure_name):
        if path.is_dir():
            return path
    
    # Method 2: Look for individual CIF files and infer the structure directory
    cif_name = f"ligandmpnn_{run_code}_seq{seq_num}_0_model.cif"
    for cif_path in af3_base.rglob(cif_name):
        if cif_path.is_file():
            # Check if there's a proper structure directory, or if it's just a standalone CIF
            parent_dir = cif_path.parent
            if parent_dir.name == structure_name:
                return parent_dir
            else:
                # Create a temporary structure info for standalone CIFs
                return {"standalone_cif": cif_path}
    
    return None

def extract_specific_runs(af3_base_dir, output_dir):
    """
    Extract all seed variants for specific run codes and sequence numbers
    """
    af3_base = Path(af3_base_dir)
    output_base = Path(output_dir)
    output_base.mkdir(parents=True, exist_ok=True)
    
    found_count = 0
    missing_count = 0
    total_seeds = 0
    
    print(f"Searching for {len(TARGETS)} specific run/sequence combinations...")
    print(f"Base directory: {af3_base_dir}")
    print(f"Output directory: {output_dir}")
    
    for run_code, seq_num in TARGETS:
        print(f"\n=== Processing Run {run_code}, Sequence {seq_num} ===")
        
        # Search for this structure
        result = find_structure_path(af3_base, run_code, seq_num)
        
        if not result:
            print(f"  Structure not found for run {run_code}, seq {seq_num}")
            missing_count += 1
            continue
        
        # Handle standalone CIF files differently from full directory structures
        if isinstance(result, dict) and "standalone_cif" in result:
            print(f"  Found standalone CIF at: {result['standalone_cif']}")
            found_count += 1
            
            # Create output subdirectory for this run/sequence
            run_output_dir = output_base / f"run_{run_code}_seq{seq_num}"
            run_output_dir.mkdir(exist_ok=True)
            
            # Copy the standalone CIF file
            cif_file = result["standalone_cif"]
            shutil.copy2(cif_file, run_output_dir / f"run_{run_code}_seq{seq_num}_main_model.cif")
            print(f"  Copied standalone model (no seed variants available)")
            
        else:
            # Full directory structure with seed variants
            structure_path = result
            print(f"  Found structure at: {structure_path}")
            found_count += 1
            
            # Create output subdirectory for this run/sequence
            run_output_dir = output_base / f"run_{run_code}_seq{seq_num}"
            run_output_dir.mkdir(exist_ok=True)
            
            # Copy the main model (highest ranked)
            main_model = structure_path / f"ligandmpnn_{run_code}_seq{seq_num}_0_model.cif"
            if main_model.exists():
                shutil.copy2(main_model, run_output_dir / f"run_{run_code}_seq{seq_num}_main_model.cif")
                print(f"  Copied main model")
            
            # Copy all confidence files
            confidence_files = [
                f"ligandmpnn_{run_code}_seq{seq_num}_0_confidences.json",
                f"ligandmpnn_{run_code}_seq{seq_num}_0_summary_confidences.json",
                f"ligandmpnn_{run_code}_seq{seq_num}_0_data.json"
            ]
            
            for conf_file in confidence_files:
                source_file = structure_path / conf_file
                if source_file.exists():
                    shutil.copy2(source_file, run_output_dir / f"run_{run_code}_seq{seq_num}_{conf_file.split('_')[-1]}")
            
            # Find and copy all seed variant models
            seed_dirs = [d for d in structure_path.iterdir() if d.is_dir() and d.name.startswith("seed-")]
            
            if not seed_dirs:
                print(f"  No seed directories found")
                continue
            
            print(f"  Found {len(seed_dirs)} seed variants:")
            
            for seed_dir in sorted(seed_dirs):
                seed_model = seed_dir / "model.cif"
                seed_conf = seed_dir / "confidences.json"
                
                if seed_model.exists():
                    # Extract seed and sample numbers from directory name
                    seed_name = seed_dir.name  # e.g., "seed-837_sample-2"
                    
                    # Copy the seed model
                    output_name = f"run_{run_code}_seq{seq_num}_{seed_name}_model.cif"
                    shutil.copy2(seed_model, run_output_dir / output_name)
                    
                    # Copy seed confidence if available
                    if seed_conf.exists():
                        conf_output_name = f"run_{run_code}_seq{seq_num}_{seed_name}_confidences.json"
                        shutil.copy2(seed_conf, run_output_dir / conf_output_name)
                    
                    print(f"    • {seed_name} -> {output_name}")
                    total_seeds += 1
    
    # Summary
    print(f"\n{'='*60}")
    print(f"EXTRACTION SUMMARY")
    print(f"{'='*60}")
    print(f"Target combinations: {len(TARGETS)}")
    print(f"Found: {found_count}")
    print(f"Missing: {missing_count}")
    print(f"Total seed models extracted: {total_seeds}")
    print(f"Output directory: {output_dir}")
    
    if missing_count > 0:
        print(f"\n{missing_count} combinations were not found.")
        print("This could be because:")
        print("- The AF3 job failed or didn't complete")
        print("- The sequence number doesn't exist for that run")
        print("- The files are named differently than expected")

def main():
    if len(sys.argv) != 3:
        print("Usage: python extract_specific_runs.py <af3_output_base_dir> <output_dir>")
        print("\nExample:")
        print("  python extract_specific_runs.py /work/user/AF3/AF3_output /work/user/selected_structures")
        print("\nNote: Edit the TARGETS list in the script to specify which runs you want")
        print(f"Current targets: {len(TARGETS)} combinations")
        for run_code, seq_num in TARGETS[:5]:  # Show first 5
            print(f"  - Run {run_code}, Sequence {seq_num}")
        if len(TARGETS) > 5:
            print(f"  ... and {len(TARGETS) - 5} more")
        sys.exit(1)
    
    af3_base_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    if not os.path.exists(af3_base_dir):
        print(f"Error: Base directory '{af3_base_dir}' does not exist")
        sys.exit(1)
    
    extract_specific_runs(af3_base_dir, output_dir)

if __name__ == "__main__":
    main()
