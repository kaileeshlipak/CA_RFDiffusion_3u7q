#!/usr/bin/env python3
"""
Simple script to find all model.cif files and rename them based on their path
Usage: python curate_af3_run.py <output_dir> <work_dir>
"""

import os
import sys
import shutil
from pathlib import Path

def curate_af3_run(outdir, workdir):
    """
    Find all model.cif files and rename them with path info
    """
    # Create output directory if it doesn't exist
    os.makedirs(outdir, exist_ok=True)
    
    workdir_path = Path(workdir)
    counter = 0
    
    # Find ALL model.cif files recursively
    for file_path in workdir_path.rglob("model.cif"):
        counter += 1
        
        # Create a descriptive name from the path
        # Get the parent directory names to create unique identifier
        parent_dirs = file_path.parts[-3:]  # Get last 3 directory levels
        path_name = "_".join(parent_dirs[:-1])  # Exclude the filename itself
        
        # Clean up the name (remove problematic characters)
        clean_name = path_name.replace("/", "_").replace(" ", "_").replace("*", "star")
        
        # Create new filename
        new_filename = f"{clean_name}_{counter:03d}.cif"
        output_path = os.path.join(outdir, new_filename)
        
        print(f"Found: {file_path}")
        print(f"  -> Copying to: {new_filename}")
        
        # Copy file to output directory
        try:
            shutil.copy2(file_path, output_path)
        except Exception as e:
            print(f"  Error copying: {e}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python curate_af3_run.py <output_dir> <work_dir>")
        print("Example: python curate_af3_run.py /path/to/output /path/to/af3_results")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    work_dir = sys.argv[2]
    
    if not os.path.exists(work_dir):
        print(f"Error: Work directory '{work_dir}' does not exist")
        sys.exit(1)
    
    print(f"Finding all model.cif files in: {work_dir}")
    print(f"Output directory: {output_dir}")
    
    curate_af3_run(output_dir, work_dir)
    
    print("Done!")

if __name__ == "__main__":
    main()
