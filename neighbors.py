### Command:
# python selective_dynamics.py [POSCAR_filename] [optional:atom labels to get neighbors] [optional:radius of search]
# input_filename: POSCAR file that the user want to get coordinates
# atoms: follows VESTA-like labels. Separated by spaces
# ex) Si1 Si2 Si47 O28
# radius of search: provide liek 'r=2.5'. Default is 2.5 ang.
# 

import sys
import numpy as np
from pymatgen.core import Structure
import pandas as pd

#sys.argv=['test','CONTCAR','Si1','Si2','O1','r=2']
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
    atom_list=atom_input.split()

elif len(sys.argv) >= 3:
    atom_list=sys.argv[2:]

if len(atom_list)==0:
    print('No atom labels are provided')
    sys.exit()


##############################################################################
########## Find neighbors & print
##############################################################################

struct = Structure.from_file(sys.argv[1])
df=create_df(struct)

print("   Atom label |     x       y       z    | Bond length ")
print("   ───────────┼──────────────────────────|─────────────")
# list of (3,) numpy arrays
coordinate_list=[]
for atom_label in atom_list:
    A1=df[df['atom_label']==atom_label]['pmg_site'].iloc[0]
    temp_array=A1.frac_coords
    print("   {0:>5}      |  {1:6.4f}  {2:6.4f}  {3:6.4f}  |  ------ ".format(atom_label,temp_array[0],temp_array[1],temp_array[2]))
    neighbors=struct.get_neighbors(A1,r=radius)
    for A2 in neighbors:
        A2.to_unit_cell(in_place=True)
        temp_array=A2.frac_coords
        # temp_list = lines[index_dict[atom_label]].split()
        # coordinate_list.append(np.array([float(i) for i in temp_list]))
        # temp_array=np.array([float(i) for i in temp_list])
        print("       └{0:>5} |  {1:6.4f}  {2:6.4f}  {3:6.4f}  |  {4:6.4f}"
              .format(df[df['pmg_site']==A2]['atom_label'].iloc[0],
                      temp_array[0],temp_array[1],temp_array[2], A1.distance(A2) ))


