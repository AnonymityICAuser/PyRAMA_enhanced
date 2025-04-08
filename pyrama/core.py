# init.py
from __future__ import division, print_function

import math
import os
import sys
import argparse

import matplotlib.pyplot as plt
import numpy as np
from Bio import PDB
from matplotlib import colors

from config import RAMA_PREFERENCES, DEFAULT_OUTPUT_DIR

RAMA_PREF_VALUES = None


def _cache_RAMA_PREF_VALUES():
    f_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    RAMA_PREF_VALUES = {}
    for key, val in RAMA_PREFERENCES.items():
        RAMA_PREF_VALUES[key] = np.full((360, 360), 0, dtype=np.float64)
        with open(os.path.join(f_path, val["file"])) as fn:
            for line in fn:
                if line.startswith("#"):
                    continue
                else:
                    x = int(float(line.split()[1]))
                    y = int(float(line.split()[0]))
                    RAMA_PREF_VALUES[key][x + 180][y + 180] \
                        = RAMA_PREF_VALUES[key][x + 179][y + 179] \
                        = RAMA_PREF_VALUES[key][x + 179][y + 180] \
                        = RAMA_PREF_VALUES[key][x + 180][y + 179] \
                        = float(line.split()[2]) 
    return RAMA_PREF_VALUES


def calc_ramachandran(file_name_list):
    """
    Main calculation and plotting definition
    :param file_name_list: List of PDB files to plot
    :return: Dictionary of normal and outlier angles by chain, and a report of outliers
    """
    global RAMA_PREF_VALUES

    if RAMA_PREF_VALUES is None:
        RAMA_PREF_VALUES = _cache_RAMA_PREF_VALUES()

    # Results organized by chain
    results = {}
    # Outlier report
    outlier_report = []
    
    # Calculate the torsion angle of the inputs
    for inp in file_name_list:
        if not os.path.isfile(inp):
            continue
            
        structure = PDB.PDBParser().get_structure('input_structure', inp)
        for model in structure:
            for chain in model:
                chain_id = chain.id
                
                # Initialize data structures for this chain
                if chain_id not in results:
                    results[chain_id] = {
                        "normals": {},
                        "outliers": {}
                    }
                    for key in RAMA_PREFERENCES.keys():
                        results[chain_id]["normals"][key] = {"x": [], "y": [], "residues": []}
                        results[chain_id]["outliers"][key] = {"x": [], "y": [], "residues": []}
                
                polypeptides = PDB.PPBuilder().build_peptides(chain)
                for poly_index, poly in enumerate(polypeptides):
                    phi_psi = poly.get_phi_psi_list()
                    for res_index, residue in enumerate(poly):
                        res_name = "{}".format(residue.resname)
                        res_num = residue.id[1]
                        phi, psi = phi_psi[res_index]
                        if phi and psi:
                            if res_index + 1 < len(poly) and str(poly[res_index + 1].resname) == "PRO":
                                aa_type = "PRE-PRO"
                            elif res_name == "PRO":
                                aa_type = "PRO"
                            elif res_name == "GLY":
                                aa_type = "GLY"
                            else:
                                aa_type = "General"
                                
                            phi_deg = math.degrees(phi)
                            psi_deg = math.degrees(psi)
                            
                            # Check if this is an outlier
                            is_outlier = RAMA_PREF_VALUES[aa_type][int(psi_deg) + 180][int(phi_deg) + 180] < \
                                    RAMA_PREFERENCES[aa_type]["bounds"][1]
                                    
                            residue_info = {
                                "name": res_name,
                                "number": res_num,
                                "phi": phi_deg,
                                "psi": psi_deg
                            }
                            
                            if is_outlier:
                                results[chain_id]["outliers"][aa_type]["x"].append(phi_deg)
                                results[chain_id]["outliers"][aa_type]["y"].append(psi_deg)
                                results[chain_id]["outliers"][aa_type]["residues"].append(residue_info)
                                
                                # Add to outlier report
                                outlier_report.append({
                                    "file": os.path.basename(inp),
                                    "chain": chain_id,
                                    "residue_name": res_name,
                                    "residue_number": res_num,
                                    "phi": phi_deg,
                                    "psi": psi_deg,
                                    "type": aa_type
                                })
                            else:
                                results[chain_id]["normals"][aa_type]["x"].append(phi_deg)
                                results[chain_id]["normals"][aa_type]["y"].append(psi_deg)
                                results[chain_id]["normals"][aa_type]["residues"].append(residue_info)
                                
    return results, outlier_report


def plot_ramachandran(results, output_dir=None, show_plots=True):
    """
    Plot Ramachandran plots for each chain separately
    :param results: Dictionary with chain-specific data
    :param output_dir: Directory to save plot files (if None, won't save)
    :param show_plots: Whether to display plots (True) or just save them (False)
    :return: List of paths to saved plot files
    """
    global RAMA_PREF_VALUES
    if RAMA_PREF_VALUES is None:
        RAMA_PREF_VALUES = _cache_RAMA_PREF_VALUES()
    
    saved_files = []
    
    # Create a separate figure for each chain
    for chain_id, chain_data in results.items():
        plt.figure(figsize=(12, 10))
        plt.suptitle(f"Ramachandran Plot - Chain {chain_id}")
        
        normals = chain_data["normals"]
        outliers = chain_data["outliers"]
        
        for idx, (key, val) in enumerate(sorted(RAMA_PREFERENCES.items(), key=lambda x: x[0].lower())):
            plt.subplot(2, 2, idx + 1)
            plt.title(key)
            plt.imshow(RAMA_PREF_VALUES[key], cmap=RAMA_PREFERENCES[key]["cmap"],
                    norm=colors.BoundaryNorm(RAMA_PREFERENCES[key]["bounds"], RAMA_PREFERENCES[key]["cmap"].N),
                    extent=(-180, 180, 180, -180))
            
            # Plot normal points
            if normals[key]["x"]:
                plt.scatter(normals[key]["x"], normals[key]["y"], label="Normal")
            
            # Plot outliers
            if outliers[key]["x"]:
                plt.scatter(outliers[key]["x"], outliers[key]["y"], color="red", label="Outlier")
                
                # Add annotations for outliers
                for i, (x, y, res) in enumerate(zip(outliers[key]["x"], outliers[key]["y"], outliers[key]["residues"])):
                    plt.annotate(f"{res['name']}{res['number']}", 
                               (x, y),
                               textcoords="offset points",
                               xytext=(0, 5),
                               ha='center')
            
            plt.xlim([-180, 180])
            plt.ylim([-180, 180])
            plt.plot([-180, 180], [0, 0], color="black")
            plt.plot([0, 0], [-180, 180], color="black")
            plt.locator_params(axis='x', nbins=7)
            plt.xlabel(r'$\phi$')
            plt.ylabel(r'$\psi$')
            plt.grid()
            
            if normals[key]["x"] or outliers[key]["x"]:
                plt.legend()

        plt.tight_layout()
        
        # Save plot if output directory is provided
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            plot_file = os.path.join(output_dir, f"ramachandran_chain_{chain_id}.png")
            plt.savefig(plot_file, dpi=300)
            saved_files.append(plot_file)
        
        if show_plots:
            plt.show()
        else:
            plt.close()
    
    return saved_files


def generate_outlier_report(outlier_report):
    """
    Generate a formatted report of all outlier residues
    :param outlier_report: List of outlier dictionaries
    :return: Formatted text report
    """
    if not outlier_report:
        return "No outliers found."
        
    report = "Ramachandran Plot Outlier Report\n"
    report += "==============================\n\n"
    
    # Sort by file, then chain, then residue number
    outlier_report.sort(key=lambda x: (x["file"], x["chain"], x["residue_number"]))
    
    current_file = None
    current_chain = None
    
    for outlier in outlier_report:
        # Print headers when file or chain changes
        if outlier["file"] != current_file:
            current_file = outlier["file"]
            report += f"\nFile: {current_file}\n"
            report += "-" * (len(current_file) + 6) + "\n"
            current_chain = None
            
        if outlier["chain"] != current_chain:
            current_chain = outlier["chain"]
            report += f"\n  Chain {current_chain}:\n"
        
        # Print outlier details
        report += f"    {outlier['residue_name']} {outlier['residue_number']}: "
        report += f"φ = {outlier['phi']:.2f}°, ψ = {outlier['psi']:.2f}° "
        report += f"(Category: {outlier['type']})\n"
    
    return report


def analyze_protein_structure(file_name_list, output_dir=None, show_plots=True):
    """
    Analyze protein structures and generate reports and plots
    :param file_name_list: List of PDB files to analyze
    :param output_dir: Directory to save output files (if None, won't save)
    :param show_plots: Whether to display plots (True) or just save them (False)
    :return: Outlier report as string, paths to saved files
    """
    # Calculate Ramachandran angles and identify outliers
    results, outlier_report = calc_ramachandran(file_name_list)
    
    # Generate outlier report
    report = generate_outlier_report(outlier_report)
    
    saved_files = []
    
    # Save report to file if output directory is provided
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        report_file = os.path.join(output_dir, "outlier_report.txt")
        with open(report_file, "w") as f:
            f.write(report)
        saved_files.append(report_file)
    
    # Generate plots
    plot_files = plot_ramachandran(results, output_dir, show_plots)
    saved_files.extend(plot_files)
    
    return report, saved_files


def main():
    """Command line interface for the Ramachandran plot tool"""
    parser = argparse.ArgumentParser(description="Analyze protein structures and generate Ramachandran plots")
    
    parser.add_argument("pdb_files", nargs="+", help="PDB file(s) to analyze")
    parser.add_argument("-o", "--output-dir", default=DEFAULT_OUTPUT_DIR, 
                        help=f"Directory to save output files (default: {DEFAULT_OUTPUT_DIR})")
    parser.add_argument("--no-show", action="store_true", 
                        help="Don't display plots, just save them")
    parser.add_argument("--report-only", action="store_true",
                        help="Generate only the outlier report, no plots")
    
    args = parser.parse_args()
    
    # Validate input files
    valid_files = []
    for file_path in args.pdb_files:
        if not os.path.isfile(file_path):
            print(f"Warning: File not found: {file_path}", file=sys.stderr)
        else:
            valid_files.append(file_path)
    
    if not valid_files:
        print("Error: No valid PDB files provided.", file=sys.stderr)
        sys.exit(1)
    
    # Calculate Ramachandran angles and identify outliers
    results, outlier_report = calc_ramachandran(valid_files)
    
    # Generate and print outlier report
    report = generate_outlier_report(outlier_report)
    print(report)
    
    # Save report to file
    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, "outlier_report.txt")
    with open(report_path, "w") as f:
        f.write(report)
    print(f"Outlier report saved to: {report_path}")
    
    # Generate plots if requested
    if not args.report_only:
        plot_files = plot_ramachandran(results, args.output_dir, not args.no_show)
        if plot_files:
            print(f"Plot files saved to: {args.output_dir}")
            for plot_file in plot_files:
                print(f"  - {os.path.basename(plot_file)}")


# Allow direct execution as a script
if __name__ == "__main__":
    main()