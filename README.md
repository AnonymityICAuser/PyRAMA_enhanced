## **Enhanced PyRAMA: Ramachandran Plot Generator**
Python3 implementation of the Ramachandran plot. 

This is a new version of old PyRAMA which could process complexs with multiple chains, for my own ICA work. You could check each chain, with outliear report now. So now we can easily generate a clear ramachandran plot for chains.

Example: one of the chains of a complex:

![image](https://github.com/user-attachments/assets/f8e54cde-a9ca-4bfb-a2d6-17306d022f0c)


The PyRAMA tool is a Python-based utility designed for generating Ramachandran plots, which are widely used in structural biology to analyze backbone dihedral angles (ϕ, ψ) in protein structures. This enhanced version, **PyRAMA_multi**, introduces advanced functionality to handle protein complexes with multiple chains, generate clear visualizations, and report outliers effectively.

---

### **Technical Features**
1. **Multi-chain Analysis**:
   - The tool now supports processing protein complexes with multiple chains, enabling chain-specific analyses and visualizations.
   - Users can generate Ramachandran plots for individual chains, which are especially useful for large structures like ribosomes or multi-domain proteins.

2. **Outlier Detection**:
   - Automatically identifies and reports residues with outlier ϕ and ψ angles based on predefined preference files for different amino acid categories (General, Glycine, Proline, Pre-Proline).

3. **Customizable Visualization**:
   - The tool generates high-quality plots using **matplotlib**, with color-coded regions and labels for normal and outlier residues.
   - Annotates outliers directly on the plot for easy identification.

4. **Configurable Data Sources**:
   - Uses preference data files (e.g., `pref_general.data`, `pref_glycine.data`) to define valid ϕ/ψ regions for different residue types.
   - The configuration is customizable via the `config.py` file.

5. **Output Flexibility**:
   - Generates:
     - **Plots**: Saved as high-resolution PNG files.
     - **Outlier Reports**: Detailed text reports summarizing outlier residues by chain.
   - Output directory is configurable using `DEFAULT_OUTPUT_DIR` in `config.py`.

6. **Command-line Interface (CLI)**:
   - Provides a user-friendly CLI for analyzing one or more PDB files with options to customize output behavior.

---

### **Installation Instructions**
1. **Clone the Repository**:
   ```bash
   git clone https://github.com/gerdos/PyRAMA_enhanced.git
   cd PyRAMA_enhanced
   ```

2. **Install Dependencies**:
   Ensure Python 3.x is installed on your system, then install the required libraries:
   ```bash
   pip install -r requirements.txt
   ```
   Alternatively, install dependencies manually:
   ```bash
   pip install biopython numpy matplotlib
   ```

3. **Set Up the Tool**: (May not work in current version)
   To make the tool executable, install it as a Python package:
   ```bash
   python setup.py install
   ```

4. **Update Configuration**:
   - Modify the `DEFAULT_OUTPUT_DIR` in `pyrama/config.py` to set your preferred output directory for saving results.
   - Ensure data files (e.g., `pref_general.data`) are located in the `pyrama/data` directory.

---

### **Usage Strategies**
#### **1. Basic Command-line Usage**
To analyze a single PDB file:
```bash
python pyrama/core.py target.pdb
```

To analyze multiple PDB files:
```bash
python pyrama/core.py file1.pdb file2.pdb
```

Options:
- **`-o` / `--output-dir`**: Specify a custom directory to save results.
- **`--no-show`**: Suppress plot display; useful for batch processing.
- **`--report-only`**: Generate only the outlier report without generating plots.

Example:
```bash
python pyrama/core.py my_structure.pdb -o ./results --no-show
```

#### **2. Output Files**
- **Plots**:
  - Saved as PNG files in the specified output directory.
  - Files are named as `ramachandran_chain_<chain_id>.png`.
- **Outlier Report**:
  - A text file (`outlier_report.txt`) summarizing all outliers across chains and residues, including their dihedral angles and classifications.

#### **3. Customizing Visualization**
Edit the `RAMA_PREFERENCES` dictionary in `config.py` to:
- Adjust color schemes (via `cmap`).
- Update bounds for preference regions.

Example:
```python
"General": {
    "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
    "bounds": [0, 0.001, 0.02, 1],
},
```

#### **4. Outlier Inspection**
Use the outlier report to investigate specific residues:
- Outliers are categorized as **General**, **GLY**, **PRO**, or **PRE-PRO**.
- Each entry includes:
  - Residue name and number.
  - Dihedral angles (ϕ, ψ).
  - Classification (e.g., General, Glycine, etc.).

Example:
```
File: my_structure.pdb
-----------------------

  Chain A:
    GLY 45: φ = -60.12°, ψ = 180.45° (Category: GLY)
    PRO 78: φ = -70.34°, ψ = 150.67° (Category: PRO)
```

#### **5. Integrating into Workflows**
- **Batch Processing**:
  Incorporate the tool into larger pipelines by invoking it programmatically. For example:
  ```python
  from pyrama.core import analyze_protein_structure

  pdb_files = ["file1.pdb", "file2.pdb"]
  output_dir = "./batch_results"
  report, saved_files = analyze_protein_structure(pdb_files, output_dir, show_plots=False)
  print(report)
  ```

- **Automated Outlier Detection**:
  Use the `generate_outlier_report` function to extract only outlier data for further analysis.


### **Future Enhancements**
This enhanced PyRAMA tool is well-suited for current bioinformatics tasks, but future work could include:
- **Support for non-standard residues**: Extend the tool to handle modified residues or ligands.
- **Interactive Plots**: Use tools like Plotly for interactive visualizations.
- **Integration with Structural Databases**: Automatically fetch and analyze structures from PDB.

For questions or suggestions, feel free to contact the developer. This version would request an update to pypi.
