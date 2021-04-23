"""
## Command:
python selective_dynamics.py [POSCAR_filename] [optional:atom labels to get neighbors] [optional:radius of search]

[POSCAR_filename]:
    POSCAR file that the user want to get coordinates
[optional:atom labels to get neighbors]:
    follows VESTA-like labels. Separated by spaces
    ex) Si1 Si2 Si47 O28
[optional:radius of search]:
    ex) r=2.5 
    Default is 2.5 ang.
"""

import sys
import numpy as np
from pymatgen.core import Structure
import pandas as pd

sys.argv=['test','CONTCAR','Si1','r=2']
#print(sys.argv[0])

radius=None
for i in range(len(sys.argv)):
    if 'r=' in sys.argv[i]:
        radius=float(sys.argv[i].split('=')[1])
        # Once parsed, remove from the sys.argv
        sys.argv.remove(sys.argv[i])
if radius == None:
    radius=2.5

def create_df(pmg_struct):
	"""
	input:
	pmg_struct: pymatgen structure. Structure containing organic molecule to rotate

	output:
	dataframe: Pandas DataFrame 
		with columns=['site_index','atom_label','pmg_site','element']
	"""
	n_atom_count_dict=dict()
	dataframe=pd.DataFrame(columns=['site_index','atom_label','pmg_site','element'])
	for i in range(0,pmg_struct.num_sites):
	    # Update label for each element
	    if pmg_struct[i].specie in list(n_atom_count_dict.keys()):
	        n_atom_count_dict.update({pmg_struct[i].specie:n_atom_count_dict[pmg_struct[i].specie]+1})
	    else:
	        n_atom_count_dict.update({pmg_struct[i].specie:1})
	        
	    label='{0}{1}'.format(pmg_struct.species[i], n_atom_count_dict[pmg_struct[i].specie])
	    # Append a site to the data frame
	    # If this part is costly, maybe using pd.concat would be faster (not sure yet)
	    dataframe= dataframe.append({'site_index':i, \
	                'atom_label':'{0}'.format(label), \
	                'pmg_site':pmg_struct.sites[i],\
	                'element':str((pmg_struct.sites[i]).specie)},ignore_index=True)

	return dataframe

##############################################################################
########## Read POSCAR file and indexing each line with atom label
##############################################################################
if len(sys.argv) == 1:
    print('POSCAR is not provided')
    sys.exit()

##############################################################################
########## Ask user for atoms to check
##############################################################################

# If atoms are not states in the command lines
if len(sys.argv) == 2:
    print('Of which atoms do you want to get neighbors?')
    print("Separate with spaces: ex) 'Si1 Si2 O1 O3'")
    atom_input=input()
    entry_list=atom_input.split()

elif len(sys.argv) >= 3:
    entry_list=sys.argv[2:]

##############################################################################
########## Additional treatment in case user put in coordinates
##############################################################################
atom_list=[]
flag_add_to_last=0
for i,entry in enumerate(entry_list):
    if flag_add_to_last ==0:
        atom_list.append(entry)
    else: 
        atom_list[-1]=atom_list[-1]+entry
    
    if ',' in entry[-1]:
        flag_add_to_last=1
    else:
        flag_add_to_last=0

if len(atom_list)==0:
    print('No atom labels are provided')
    sys.exit()

##############################################################################
########## Find neighbors & print
##############################################################################

struct = Structure.from_file(sys.argv[1])

SG = struct.get_space_group_info()[0]
df=create_df(struct)

lattice_vector = np.array([struct.lattice.a,struct.lattice.b,struct.lattice.c])

print("   Atom label |     x       y       z    | Distance ")
print("   ───────────┼──────────────────────────┼──────────")
# list of (3,) numpy arrays
coordinate_list=[]
for atom_label in atom_list:
    # If label is atom label like 'Si1'
    if atom_label[0].isalpha():
        A1=df[df['atom_label']==atom_label]['pmg_site'].iloc[0]
        temp_array=A1.frac_coords
        print("   {0:>5}      | {1: 5.4f} {2: 5.4f} {3: 5.4f}  |  ------ ".format(atom_label,temp_array[0],temp_array[1],temp_array[2]))
        neighbors=struct.get_neighbors(A1,r=radius)
        
        for A2 in neighbors:
            A2.to_unit_cell(in_place=True)
            temp_array=A2.frac_coords
            # temp_list = lines[index_dict[atom_label]].split()
            # coordinate_list.append(np.array([float(i) for i in temp_list]))
            # temp_array=np.array([float(i) for i in temp_list])
            print("       └{0:>5} | {1: 5.4f} {2: 5.4f} {3: 5.4f}  |  {4:6.4f}"
                  .format(df[df['pmg_site']==A2]['atom_label'].iloc[0],
                          temp_array[0],temp_array[1],temp_array[2], A1.distance(A2) ))
    # If entry is coordinate
    else:
        array=np.array(eval(atom_label))
        print("     coords   | {0: 5.4f} {1: 5.4f} {2: 5.4f}  |  ------ ".format(array[0],array[1],array[2]))
        neighbors=struct.get_sites_in_sphere(array,r=radius)
        for A2_tuple in neighbors:
            A2=A2_tuple[0]
            A2.to_unit_cell(in_place=True)
            temp_array=A2.frac_coords
            # temp_list = lines[index_dict[atom_label]].split()
            # coordinate_list.append(np.array([float(i) for i in temp_list]))
            # temp_array=np.array([float(i) for i in temp_list])
            print("       └{0:>5} | {1: 5.4f} {2: 5.4f} {3: 5.4f}  |  {4:6.4f}"
                  .format(df[df['pmg_site']==A2]['atom_label'].iloc[0],
                          temp_array[0],temp_array[1],temp_array[2], A2_tuple[1] ))



A1=df[df['atom_label']=='Si1']['pmg_site'].iloc[0]
neighbors=struct.get_neighbors(A1,r=radius)
for A2 in neighbors:
    A2.to_unit_cell(in_place=True)
    A2_label=df[df['pmg_site']==A2]['atom_label'].iloc[0]
    random_vector=np.random.ranf([1,3])-0.5
    random_vector=0.15/A1.distance(A2)*random_vector[0]/np.linalg.norm(random_vector)
    A2.frac_coords = A2.frac_coords + random_vector/lattice_vector
    temp_array=A2.frac_coords
    print("       └{0:>5} | {1: 5.4f} {2: 5.4f} {3: 5.4f}  |  {4:6.4f}"
                  .format(A2_label,
                          temp_array[0],temp_array[1],temp_array[2], A1.distance(A2) ))
    
    
    