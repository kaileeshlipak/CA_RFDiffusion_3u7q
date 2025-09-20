#!/usr/bin/env python3
"""
Simplified script for analyzing pLDDT scores from Alphafold3 results
Only calculates pLDDT statistics without RMSD calculations

Usage:
    python plddt_analyzer.py -t target.cif
    python plddt_analyzer.py -i /path/to/structures/ -o ./results/
    
Or use as a module:
    from plddt_analyzer import analyze_plddt
    result = analyze_plddt(structure_path)
"""

import os
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Optional

from biotite.structure.io.pdb import PDBFile
import biotite.structure.io.pdbx as pdbx
import numpy as np


def find_femoco_center(structure):
    """
    Find the center of FeMoco (ICS residue, CX atom) in the structure.
    
    Args:
        structure: AtomArray to search
        
    Returns:
        Coordinates of FeMoco center or None if not found
    """
    # Look for FeMoco: resn ICS and name CX
    femoco_mask = (structure.res_name == "ICS") & (structure.atom_name == "CX")
    
    if np.any(femoco_mask):
        femoco_atom = structure[femoco_mask][0]
        return femoco_atom.coord
    else:
        return None


def get_residues_within_distance(structure, center_coord, distance_cutoff=10.0):
    """
    Get all residues within specified distance of a center point.
    
    Args:
        structure: AtomArray to search
        center_coord: Center coordinates
        distance_cutoff: Distance cutoff in Angstroms
        
    Returns:
        AtomArray containing all atoms of residues within distance
    """
    from biotite.structure import distance
    
    # Calculate distances from center to all atoms
    distances = np.array([distance(atom, center_coord) for atom in structure])
    
    # Find atoms within distance
    within_distance_mask = distances <= distance_cutoff
    
    # Get residue IDs of atoms within distance
    residues_within = np.unique(structure.res_id[within_distance_mask])
    
    # Select all atoms of those residues (byres selection)
    byres_mask = np.isin(structure.res_id, residues_within)
    
    return structure[byres_mask]


def analyze_plddt(structure_path: str, distance_threshold: float = 10.0) -> Dict:
    """
    Analyze pLDDT scores for a single structure.
    
    Args:
        structure_path: Path to structure file (PDB/CIF)
        distance_threshold: Distance threshold for FeMoco site analysis
        
    Returns:
        Dictionary with pLDDT analysis results
    """
    results = {
        'structure_path': structure_path,
        'mean_plddt': None,
        'median_plddt': None,
        'std_plddt': None,
        'min_plddt': None,
        'max_plddt': None,
        'femoco_site_mean_plddt': None,
        'femoco_site_median_plddt': None,
        'femoco_site_std_plddt': None,
        'num_femoco_site_residues': None,
        'error': None
    }
    
    try:
        # Load structure and b-factors (pLDDT values)
        if not os.path.exists(structure_path):
            raise FileNotFoundError(f"Structure not found: {structure_path}")
            
        if structure_path.endswith('.pdb'):
            target_file = PDBFile.read(structure_path)
            target_pdb = target_file.get_structure()[0]
            b_factors = target_file.get_b_factor()[0]
        elif structure_path.endswith('.cif'):
            target_pdb = pdbx.get_structure(
                pdbx.CIFFile.read(structure_path),
                extra_fields=['b_factor']
            )[0]
            b_factors = [float(atom.b_factor) for atom in target_pdb]
        else:
            raise ValueError("Structure must be .pdb or .cif format")
        
        # Overall pLDDT statistics
        results['mean_plddt'] = np.mean(b_factors)
        results['median_plddt'] = np.median(b_factors)
        results['std_plddt'] = np.std(b_factors)
        results['min_plddt'] = np.min(b_factors)
        results['max_plddt'] = np.max(b_factors)
        
        # FeMoco site pLDDT analysis (if FeMoco is present)
        femoco_center = find_femoco_center(target_pdb)
        if femoco_center is not None:
            femoco_site_residues = get_residues_within_distance(target_pdb, femoco_center, distance_threshold)
            femoco_site_indices = np.unique(np.where(np.isin(target_pdb.res_id, femoco_site_residues.res_id))[0])
            femoco_site_plddt = [b_factors[i] for i in femoco_site_indices if i < len(b_factors)]
            
            if femoco_site_plddt:
                results['femoco_site_mean_plddt'] = np.mean(femoco_site_plddt)
                results['femoco_site_median_plddt'] = np.median(femoco_site_plddt)
                results['femoco_site_std_plddt'] = np.std(femoco_site_plddt)
                results['num_femoco_site_residues'] = len(np.unique(femoco_site_residues.res_id))
        
        print(f"Analysis completed for {Path(structure_path).name}")
        print(f"  Overall mean pLDDT: {results['mean_plddt']:.2f}")
        print(f"  Overall median pLDDT: {results['median_plddt']:.2f}")
        print(f"  pLDDT range: {results['min_plddt']:.1f} - {results['max_plddt']:.1f}")
        
        if results['femoco_site_mean_plddt'] is not None:
            print(f"  FeMoco site mean pLDDT: {results['femoco_site_mean_plddt']:.2f}")
            print(f"  FeMoco site residues: {results['num_femoco_site_residues']}")
        else:
            print(f"  No FeMoco found in structure")
        
    except Exception as e:
        results['error'] = str(e)
        print(f"Error analyzing {structure_path}: {e}")
    
    return results


def batch_analyze_directory(input_directory: str, 
                           distance_threshold: float = 10.0,
                           output_directory: Optional[str] = None) -> List[Dict]:
    """
    Analyze pLDDT scores for all structure files in a directory.
    
    Args:
        input_directory: Directory containing structure files
        distance_threshold: Distance threshold for FeMoco site analysis
        output_directory: Directory to save results
        
    Returns:
        List of analysis results for each structure
    """
    input_dir = Path(input_directory)
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_directory}")
    
    # Find all structure files
    structure_files = []
    for ext in ['*.cif', '*.pdb']:
        structure_files.extend(input_dir.glob(ext))
    
    if not structure_files:
        print(f"No structure files (.cif or .pdb) found in {input_directory}")
        return []
    
    print(f"Found {len(structure_files)} structure files to process")
    
    # Setup output directory and CSV file
    if output_directory:
        output_dir = Path(output_directory)
        output_dir.mkdir(exist_ok=True, parents=True)
        
        csv_file = output_dir / "plddt_analysis.csv"
        header = "filename,mean_plddt,median_plddt,std_plddt,min_plddt,max_plddt,femoco_site_mean_plddt,femoco_site_median_plddt,femoco_site_std_plddt,num_femoco_site_residues"
        
        with open(csv_file, 'w') as f:
            f.write(header + "\n")
    
    # Process each structure
    all_results = []
    
    for i, structure_file in enumerate(structure_files, 1):
        print(f"\n[{i}/{len(structure_files)}] Processing: {structure_file.name}")
        
        try:
            # Analyze pLDDT
            result = analyze_plddt(str(structure_file), distance_threshold)
            result['filename'] = structure_file.name
            all_results.append(result)
            
            # Write to CSV
            if output_directory:
                row_data = [
                    structure_file.stem,
                    f"{result.get('mean_plddt'):.2f}" if result.get('mean_plddt') is not None else 'NA',
                    f"{result.get('median_plddt'):.2f}" if result.get('median_plddt') is not None else 'NA',
                    f"{result.get('std_plddt'):.2f}" if result.get('std_plddt') is not None else 'NA',
                    f"{result.get('min_plddt'):.1f}" if result.get('min_plddt') is not None else 'NA',
                    f"{result.get('max_plddt'):.1f}" if result.get('max_plddt') is not None else 'NA',
                    f"{result.get('femoco_site_mean_plddt'):.2f}" if result.get('femoco_site_mean_plddt') is not None else 'NA',
                    f"{result.get('femoco_site_median_plddt'):.2f}" if result.get('femoco_site_median_plddt') is not None else 'NA',
                    f"{result.get('femoco_site_std_plddt'):.2f}" if result.get('femoco_site_std_plddt') is not None else 'NA',
                    str(result.get('num_femoco_site_residues')) if result.get('num_femoco_site_residues') is not None else 'NA'
                ]
                
                with open(csv_file, 'a') as f:
                    f.write(",".join(row_data) + "\n")
            
        except Exception as e:
            print(f"    ERROR: {e}")
            error_result = {
                'filename': structure_file.name,
                'error': str(e)
            }
            all_results.append(error_result)
            
            # Write error to CSV
            if output_directory:
                error_row = f"{structure_file.stem},ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR"
                with open(csv_file, 'a') as f:
                    f.write(error_row + "\n")
    
    # Print final summary
    print(f"\n{'='*60}")
    print(f"PLDDT ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print(f"Total structures processed: {len(structure_files)}")
    
    successful = [r for r in all_results if r.get('mean_plddt') is not None]
    failed = [r for r in all_results if r.get('mean_plddt') is None]
    
    print(f"Successful analyses: {len(successful)}")
    print(f"Failed analyses: {len(failed)}")
    
    if successful:
        mean_plddts = [r['mean_plddt'] for r in successful]
        print(f"Mean pLDDT range: {min(mean_plddts):.2f} - {max(mean_plddts):.2f}")
        print(f"Average mean pLDDT: {np.mean(mean_plddts):.2f}")
        
        # Summary for FeMoco sites
        femoco_results = [r for r in successful if r.get('femoco_site_mean_plddt') is not None]
        if femoco_results:
            femoco_plddts = [r['femoco_site_mean_plddt'] for r in femoco_results]
            print(f"FeMoco site pLDDT range: {min(femoco_plddts):.2f} - {max(femoco_plddts):.2f}")
            print(f"Average FeMoco site pLDDT: {np.mean(femoco_plddts):.2f}")
    
    if output_directory:
        print(f"Results saved to: {csv_file}")
    
    return all_results


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description="Analyze pLDDT scores from AlphaFold structures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Single structure analysis
    python plddt_analyzer.py -t alphafold_result.cif
    
    # Batch analysis of directory
    python plddt_analyzer.py -i /path/to/alphafold_results/ -o ./analysis_output/
    
    # Custom distance threshold for FeMoco site
    python plddt_analyzer.py -i ./models/ --threshold 8.0 -o ./results/
        """
    )
    
    # Make target and input mutually exclusive
    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument("-t", "--target", type=str,
                             help="Single target structure file (PDB or CIF format)")
    target_group.add_argument("-i", "--input", type=str,
                             help="Input directory containing structure files to analyze")
    
    parser.add_argument("--threshold", type=float, default=10.0,
                       help="Distance threshold for FeMoco site analysis (default: 10.0 Å)")
    parser.add_argument("-o", "--output", type=str,
                       help="Output directory for CSV results")
    
    args = parser.parse_args()
    
    # Batch processing mode
    if args.input:
        if not args.output:
            print("Warning: No output directory specified. Results will only be printed to console.")
            print("Use -o/--output to save results to CSV file.")
        
        results = batch_analyze_directory(
            args.input,
            args.threshold,
            args.output
        )
        return results
    
    # Single structure mode
    else:
        result = analyze_plddt(args.target, args.threshold)
        
        # Handle output for single structure
        if args.output:
            output_dir = Path(args.output)
            output_dir.mkdir(exist_ok=True, parents=True)
            
            # Write results to CSV
            csv_file = output_dir / "single_plddt_analysis.csv"
            header = "filename,mean_plddt,median_plddt,std_plddt,min_plddt,max_plddt,femoco_site_mean_plddt,femoco_site_median_plddt,femoco_site_std_plddt,num_femoco_site_residues"
            
            target_name = Path(args.target).stem
            row_data = [
                target_name,
                f"{result.get('mean_plddt'):.2f}" if result.get('mean_plddt') is not None else 'NA',
                f"{result.get('median_plddt'):.2f}" if result.get('median_plddt') is not None else 'NA',
                f"{result.get('std_plddt'):.2f}" if result.get('std_plddt') is not None else 'NA',
                f"{result.get('min_plddt'):.1f}" if result.get('min_plddt') is not None else 'NA',
                f"{result.get('max_plddt'):.1f}" if result.get('max_plddt') is not None else 'NA',
                f"{result.get('femoco_site_mean_plddt'):.2f}" if result.get('femoco_site_mean_plddt') is not None else 'NA',
                f"{result.get('femoco_site_median_plddt'):.2f}" if result.get('femoco_site_median_plddt') is not None else 'NA',
                f"{result.get('femoco_site_std_plddt'):.2f}" if result.get('femoco_site_std_plddt') is not None else 'NA',
                str(result.get('num_femoco_site_residues')) if result.get('num_femoco_site_residues') is not None else 'NA'
            ]
            
            with open(csv_file, 'w') as f:
                f.write(header + "\n")
                f.write(",".join(row_data) + "\n")
            
            print(f"Results written to {csv_file}")
        
        return result


if __name__ == "__main__":
    main()
