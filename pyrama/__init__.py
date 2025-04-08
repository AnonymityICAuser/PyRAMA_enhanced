# ramachandran_tool/__init__.py
from .core import (
    calc_ramachandran,
    plot_ramachandran,
    generate_outlier_report,
    analyze_protein_structure
)

__all__ = [
    'calc_ramachandran',
    'plot_ramachandran',
    'generate_outlier_report',
    'analyze_protein_structure'
]