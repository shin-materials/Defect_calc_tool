# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 19:52:58 2021

@author: yongjin

## Command:
python selective_dynamics.py [POSCAR_filename1]  [POSCAR_filename2] [optional:threshold]

[POSCAR_filename1]: POSCAR file
[POSCAR_filename2]: POSCAR file
    ** Note that these two files should have equivalent lattice vectors
[optional:threshold]: threshold of displacement to print
    ex) 't=0.01'
    ** Default value is 0.01 angstrom
"""

import sys
import numpy as np
from pymatgen.core import Structure
import pandas as pd

#sys.argv=['test','72_small_void_neutral_2nd_relaxed.vasp','72_small_void_positive_2nd_relaxed.vasp','Si2','O1','r=0.01']
#print(sys.argv[0])

##############################################################################
########## Function definition
##############################################################################
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
########## Read POSCAR file 
##############################################################################
if len(sys.argv) <= 2:
    print('POSCAR is not provided')
    sys.exit()

# Load structure
struct1 = Structure.from_file(sys.argv[1])
struct2 = Structure.from_file(sys.argv[2])
### Threshold can be either given from user
### Or it would be 0.01 Ang
threshold=None
for i in range(len(sys.argv)):
    if '=' in sys.argv[i]:
        try:
            threshold=float(sys.argv[i].split('=')[1])
        except:
            print("Maybe you should remove equal sign (=) from structure filenames")
            sys.exit()
        # Once parsed, remove from the sys.argv
        sys.argv.remove(sys.argv[i])
if threshold == None:
    threshold=0.01

# Lattice of two structures should be equivalent
if struct1.lattice != struct2.lattice:
    print("Lattice vectors of two POSCAR files are not equivalent")
    sys.exit()

##############################################################################
########## Data analysis
##############################################################################

# Create DataFrame
df=create_df(struct1)
df1=create_df(struct1)
df2=create_df(struct2)

# Make displacement data
full_atom_list=list(df1['atom_label'])
displacement_list=[]
for atom_label in full_atom_list:
    A1_st1=df1[df1['atom_label']==atom_label]['pmg_site'].iloc[0]
    A1_st2=df2[df2['atom_label']==atom_label]['pmg_site'].iloc[0]
    displacement_list.append(A1_st1.distance(A1_st2))
    
# Insert Displ column
df.insert(2, "Displ", displacement_list, True)
# sort
df=df.sort_values(by=['Displ'],ascending=False)
# cut by threshold
df=df[df['Displ']>threshold]

##############################################################################
########## Printing
##############################################################################

### Filename treatment
filename1=sys.argv[1]
filename2=sys.argv[2]
# Align with respect to longer filename
filename1=str(filename1).ljust(max(len(filename1),len(filename2)))
filename2=str(filename2).ljust(max(len(filename1),len(filename2)))

# space for coordinate column
namespace=24
# To print, I need quotient and remainder
Q=len(filename1)//namespace
remain=len(filename1)%namespace

### Printing header
for i in range(Q):
    print('              |              | {0:<24s} | {1:<24s} |'
          .format(filename1[i*namespace:(i+1)*namespace],filename2[i*namespace:(i+1)*namespace]))
print('              | Thres.= {0:.2f} | {1:<24s} | {2:<24s} |'
          .format(threshold,filename1[Q*namespace:],filename2[Q*namespace:]))
print("   Atom label | Displ. (Ang) |     x       y       z    |     x       y       z    |")
print("   ───────────┼──────────────┼──────────────────────────┼──────────────────────────┤")


### Printing items
for atom_label in list(df['atom_label']):
    # Calculate displacement. This is always positive
    disp=df[df['atom_label']==atom_label]['Displ'].iloc[0]
    A1_st1=df1[df1['atom_label']==atom_label]['pmg_site'].iloc[0]
    A1_st2=df2[df2['atom_label']==atom_label]['pmg_site'].iloc[0]
    coord1=A1_st1.frac_coords
    coord2=A1_st2.frac_coords
    print("        {0:>5} |  {1:8.4f}    ".format(atom_label,disp) +
          "| {0: 5.4f} {1: 5.4f} {2: 5.4f}  ".format(coord1[0],coord1[1],coord1[2]) + 
          "| {0: 5.4f} {1: 5.4f} {2: 5.4f}  |".format(coord2[0],coord2[1],coord2[2]))
