# Defect_calc_tool

Python scripts to conveniently manage structure files for defect calculations in the scheme of density functional theory calculations with VASP. These scripts use atom labeling consistent with [VESTA](https://jp-minerals.org/vesta/en/), where each atom in POSCAR file is labeled with its element and number such as Si1, O1, and H1.

## Requirements
While some scripts are coded only with standard python packages (like *numpy, pandas, and sys*), the other scripts rely on *Structure* module of [pymatgen](https://pymatgen.org/).

## List of scripts

* coordinate.py
* compare_POSCAR.py
* neighbors.py
* perturb.py
* selective_dynamics.py
* defect_creation.py