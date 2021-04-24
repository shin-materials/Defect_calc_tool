# Defect_calc_tool

Python scripts to conveniently manage structure files for defect calculations in the scheme of density functional theory calculations with VASP. These scripts use atom labeling consistent with [VESTA](https://jp-minerals.org/vesta/en/), where each atom in POSCAR file is labeled with its element and number such as Si1, O1, and H1.

## Requirements
While some scripts are coded only with standard python packages (like *numpy, pandas, and sys*), the other scripts rely on *Structure* module of [pymatgen](https://pymatgen.org/).

## List of scripts
The description of each script assumes POSCAR file of SiO<sub>2</sub> (alpha-quartz). In total 9 atoms (3 silicon atoms and 6 oxygen atoms) are included, which are labeld as Si1, Si2, Si3, O1, O2, O3, O4, O5, and O6.

* coordinate.py
	* get the line number in POSCAR and coordinates of atoms.

```
$ python coordinate.py POSCAR Si1 O1
   Atom label | Line # |     x       y       z
   ───────────┼────────┼─────────────────────────
          Si1 |     8  |  0.4782  0.0000  0.6667
           O1 |    11  |  0.4178  0.2526  0.7987
```

* compare_POSCAR.py


* neighbors.py
* perturb.py
* selective_dynamics.py
* defect_creation.py