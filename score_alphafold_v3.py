#!/usr/bin/env python3
"""
Revised script for plotting Alphafold3 results
Fixed and improved version with proper active site analysis

Usage:
    python score_alphafold_v3.py -r reference.pdb -t target.cif -a "198,202,203"
    
Or use as a module:
    from score_alphafold_v2 import analyze_structure
    result = analyze_structure(ref_path, target_path, active_site_residues)
"""

import os
import sys
import argparse
from pathlib import Path
from typing import Tuple, Dict, List, Optional, Union

import biotite
from biotite.structure.io.pdb import PDBFile
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io as bsio
from biotite.structure import AtomArray, distance
import numpy as np


def find_nearest_ca(ref_pdb: AtomArray, target_pdb: AtomArray, 
                   residue_number: int, threshold: float = 5.0) -> Tuple[Union[int, str], Union[str, str]]:
    """
    Find the nearest C-alpha atom in the target structure
    to the specified residue in the reference structure.
    
    Args:
        ref_pdb: Reference structure AtomArray
        target_pdb: Target structure AtomArray  
        residue_number: Residue number in reference structure
        threshold: Maximum distance threshold in Angstroms
        
    Returns:
        Tuple of (nearest_residue_id, nearest_residue_name) or ('skip', 'skip')
    """
    # Get specific residue C-alpha atom in reference PDB
    ref_residue_mask = (ref_pdb.res_id == residue_number) & (ref_pdb.atom_name == "CA")
    if not np.any(ref_residue_mask):
        print(f'Cannot find residue {residue_number} in reference structure')
        return 'skip', 'skip'
        
    ref_residue_ca = ref_pdb[ref_residue_mask][0]
    specified_resname = ref_residue_ca.res_name

    # Get all target C-alpha atoms of the same residue type
    target_ca_mask = (target_pdb.atom_name == "CA") & (target_pdb.res_name == specified_resname)
    target_ca = target_pdb[target_ca_mask]
    
    if len(target_ca) == 0:
        print(f'Cannot find residue type {specified_resname} in target structure')
        return 'skip', 'skip'

    # Calculate distances
    distances = [distance(ca_atom, ref_residue_ca) for ca_atom in target_ca]
    
    # Create mapping and find nearest
    residue_distance_map = {
        (atom.res_id, atom.res_name): dist 
        for atom, dist in zip(target_ca, distances)
    }
    
    # Sort by distance
    sorted_residues = sorted(residue_distance_map.items(), key=lambda x: x[1])
    
    if not sorted_residues:
        print(f'No matching residues found for {residue_number}')
        return 'skip', 'skip'
        
    nearest_res_id, nearest_res_name = sorted_residues[0][0]
    nearest_distance = sorted_residues[0][1]
    
    if nearest_distance > threshold:
        print(f'Residue {residue_number} is {nearest_distance:.2f} Å away from closest match. Omitting.')
        return 'skip', 'skip'
        
    return nearest_res_id, nearest_res_name


def map_active_site_residues(ref_pdb: AtomArray, target_pdb: AtomArray, 
                           residue_list: List[int], threshold: float = 5.0) -> Dict[str, Dict[str, Union[int, str]]]:
    """
    Map active site residues from reference to target structure based on C-alpha proximity.
    
    Args:
        ref_pdb: Reference structure
        target_pdb: Target structure  
        residue_list: List of residue numbers from reference structure
        threshold: Distance threshold for mapping
        
    Returns:
        Dictionary mapping reference residues to target residues
    """
    residue_map = {}
    
    for residue_num in residue_list:
        # Get reference residue info
        ref_mask = ref_pdb.res_id == residue_num
        if not np.any(ref_mask):
            print(f'Warning: Residue {residue_num} not found in reference structure')
            continue
            
        ref_residue_name = ref_pdb[ref_mask][0].res_name
        
        # Find nearest match in target
        nearest_id, nearest_name = find_nearest_ca(ref_pdb, target_pdb, residue_num, threshold)
        
        if nearest_id != 'skip':
            key = f'{ref_residue_name}_{residue_num}'
            residue_map[key] = {
                'target_id': nearest_id, 
                'target_name': nearest_name,
                'ref_id': residue_num,
                'ref_name': ref_residue_name
            }
            
    return residue_map


def superimpose_structures(ref_pdb: AtomArray, target_pdb: AtomArray) -> Tuple[AtomArray, float]:
    """
    Superimpose target structure onto reference structure.
    
    Args:
        ref_pdb: Reference structure
        target_pdb: Target structure
        
    Returns:
        Tuple of (superimposed_target_structure, rmsd)
    """
    try:
        # Use biotite's homolog superimposition
        oriented_target, transform = biotite.structure.superimpose_homologs(
            ref_pdb, target_pdb, return_rmsd=False
        )[:2]
        
        # Calculate RMSD on CA atoms
        ref_ca = ref_pdb[ref_pdb.atom_name == "CA"]
        target_ca = oriented_target[oriented_target.atom_name == "CA"]
        
        # For RMSD calculation, we need to align by sequence or use matched residues
        # Here we'll use a simple approach - this might need refinement
        min_len = min(len(ref_ca), len(target_ca))
        if min_len > 0:
            rmsd_val = biotite.structure.rmsd(ref_ca[:min_len], target_ca[:min_len])
        else:
            rmsd_val = float('inf')
            
        return oriented_target, rmsd_val
        
    except Exception as e:
        print(f"Error in superimposition: {e}")
        return target_pdb, float('inf')


def calculate_active_site_rmsd(ref_pdb: AtomArray, target_pdb: AtomArray, 
                              residue_map: Dict[str, Dict]) -> float:
    """
    Calculate RMSD for active site residues only.
    
    Args:
        ref_pdb: Reference structure
        target_pdb: Superimposed target structure
        residue_map: Mapping between reference and target residues
        
    Returns:
        Active site RMSD
    """
    if not residue_map:
        return float('nan')
        
    ref_coords = []
    target_coords = []
    
    for key, mapping in residue_map.items():
        ref_id = mapping['ref_id']
        target_id = mapping['target_id']
        
        # Get CA atoms for both structures
        ref_ca_mask = (ref_pdb.res_id == ref_id) & (ref_pdb.atom_name == "CA")
        target_ca_mask = (target_pdb.res_id == target_id) & (target_pdb.atom_name == "CA")
        
        if np.any(ref_ca_mask) and np.any(target_ca_mask):
            ref_coords.append(ref_pdb[ref_ca_mask][0].coord)
            target_coords.append(target_pdb[target_ca_mask][0].coord)
    
    if len(ref_coords) < 3:  # Need at least 3 points for meaningful RMSD
        return float('nan')
        
    ref_coords = np.array(ref_coords)
    target_coords = np.array(target_coords)
    
    # Calculate RMSD
    diff = ref_coords - target_coords
    rmsd_val = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    return rmsd_val


def analyze_structure(ref_structure_path: str, target_structure_path: str, 
                     active_site_residues: Optional[str] = None,
                     distance_threshold: float = 5.0) -> Dict:
    """
    Independent function to analyze a single structure comparison.
    
    Args:
        ref_structure_path: Path to reference structure (PDB)
        target_structure_path: Path to target structure (PDB/CIF)
        active_site_residues: Comma-separated residue numbers (e.g., "198,202,203")
        distance_threshold: Distance threshold for residue mapping
        
    Returns:
        Dictionary with analysis results
    """
    results = {
        'reference_path': ref_structure_path,
        'target_path': target_structure_path,
        'global_rmsd': None,
        'active_site_rmsd': None,
        'mean_plddt': None,
        'active_site_mean_plddt': None,
        'active_site_median_plddt': None,
        'active_site_std_plddt': None,
        'residue_mapping': None,
        'error': None
    }
    
    try:
        # Load reference structure
        if not os.path.exists(ref_structure_path):
            raise FileNotFoundError(f"Reference structure not found: {ref_structure_path}")
            
        ref_pdb = PDBFile.read(ref_structure_path).get_structure()[0]
        
        # Load target structure
        if not os.path.exists(target_structure_path):
            raise FileNotFoundError(f"Target structure not found: {target_structure_path}")
            
        if target_structure_path.endswith('.pdb'):
            target_file = PDBFile.read(target_structure_path)
            target_pdb = target_file.get_structure()[0]
            b_factors = target_file.get_b_factor()[0]
        elif target_structure_path.endswith('.cif'):
            target_pdb = pdbx.get_structure(
                pdbx.CIFFile.read(target_structure_path),
                extra_fields=['b_factor']
            )[0]
            b_factors = [float(atom.b_factor) for atom in target_pdb]
        else:
            raise ValueError("Target structure must be .pdb or .cif format")
        
        # Superimpose structures
        oriented_target, global_rmsd = superimpose_structures(ref_pdb, target_pdb)
        results['global_rmsd'] = global_rmsd
        
        # Calculate mean pLDDT
        results['mean_plddt'] = np.mean(b_factors)
        
        # Active site analysis
        if active_site_residues:
            residue_list = [int(x.strip()) for x in active_site_residues.split(',')]
            
            # Map active site residues
            residue_map = map_active_site_residues(
                ref_pdb, oriented_target, residue_list, distance_threshold
            )
            results['residue_mapping'] = residue_map
            
            # Calculate active site RMSD
            active_rmsd = calculate_active_site_rmsd(ref_pdb, oriented_target, residue_map)
            results['active_site_rmsd'] = active_rmsd
            
            # Calculate active site pLDDT statistics
            if residue_map:
                active_site_plddt = []
                for key, mapping in residue_map.items():
                    target_id = mapping['target_id']
                    target_mask = oriented_target.res_id == target_id
                    if np.any(target_mask):
                        # Get pLDDT values for this residue
                        residue_indices = np.where(target_mask)[0]
                        residue_plddt = [b_factors[i] for i in residue_indices]
                        active_site_plddt.extend(residue_plddt)
                
                if active_site_plddt:
                    results['active_site_mean_plddt'] = np.mean(active_site_plddt)
                    results['active_site_median_plddt'] = np.median(active_site_plddt)
                    results['active_site_std_plddt'] = np.std(active_site_plddt)
        
        print(f"Analysis completed successfully!")
        print(f"Global RMSD: {global_rmsd:.3f} Å")
        print(f"Mean pLDDT: {results['mean_plddt']:.2f}")
        
        if active_site_residues and results['active_site_rmsd'] is not None:
            print(f"Active site RMSD: {results['active_site_rmsd']:.3f} Å")
            print(f"Active site mean pLDDT: {results['active_site_mean_plddt']:.2f}")
        
    except Exception as e:
        results['error'] = str(e)
        print(f"Error in analysis: {e}")
    
    return results


def write_structure(structure: AtomArray, filename: str):
    """Write structure to PDB file."""
    try:
        pdb_file = PDBFile()
        pdb_file.set_structure(structure)
        pdb_file.write(filename)
        print(f"Structure written to {filename}")
    except Exception as e:
        print(f"Error writing structure: {e}")


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description="Analyze AlphaFold structures against reference",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python score_alphafold_v2.py -r reference.pdb -t alphafold_result.cif
    python score_alphafold_v2.py -r ref.pdb -t target.pdb -a "198,202,203,206"
    python score_alphafold_v2.py -r ref.pdb -t target.cif -a "198,202" --threshold 8.0
        """
    )
    
    parser.add_argument("-r", "--reference", type=str, required=True,
                       help="Reference structure file (PDB format)")
    parser.add_argument("-t", "--target", type=str, required=True,
                       help="Target structure file (PDB or CIF format)")
    parser.add_argument("-a", "--active_site", type=str,
                       help="Comma-separated active site residue numbers from reference structure")
    parser.add_argument("--threshold", type=float, default=5.0,
                       help="Distance threshold for residue mapping (default: 5.0 Å)")
    parser.add_argument("-o", "--output", type=str,
                       help="Output directory for aligned structures and results")
    
    args = parser.parse_args()
    
    # Run analysis
    results = analyze_structure(
        args.reference, 
        args.target, 
        args.active_site,
        args.threshold
    )
    
    # Handle output
    if args.output:
        output_dir = Path(args.output)
        output_dir.mkdir(exist_ok=True, parents=True)
        
        # Write results to CSV
        csv_file = output_dir / "analysis_results.csv"
        header = "target,global_rmsd,mean_plddt,active_site_rmsd,active_site_mean_plddt,active_site_median_plddt,active_site_std_plddt"
        
        target_name = Path(args.target).stem
        row = f"{target_name},{results.get('global_rmsd', 'NA')},{results.get('mean_plddt', 'NA')},"
        row += f"{results.get('active_site_rmsd', 'NA')},{results.get('active_site_mean_plddt', 'NA')},"
        row += f"{results.get('active_site_median_plddt', 'NA')},{results.get('active_site_std_plddt', 'NA')}"
        
        with open(csv_file, 'w') as f:
            f.write(header + "\n")
            f.write(row + "\n")
        
        print(f"Results written to {csv_file}")
    
    return results


if __name__ == "__main__":
    main()
