"""
## Command:
python selective_dynamics.py [POSCAR_filename] [optional:atom labels to get coordinates]

[POSCAR_filename]:
    POSCAR file that the user want to get coordinates
[optional:atom labels to get coordinates]:
    follows VESTA-like labels. Separated by spaces
    ex) Si1 Si2 Si47 O28
"""
import sys
import numpy as np


##############################################################################
########## Read POSCAR file and indexing each line with atom label
##############################################################################
if len(sys.argv) == 1:
    print('POSCAR is not provided')
    sys.exit()

struct_file=open(sys.argv[1],'r')
lines=struct_file.readlines()

#### Calculate index, including header 
num_atom_dict=dict()
element_list=lines[5].split()
for i in range(len(element_list)):
    num_atom_dict[element_list[i]]=int((lines[6].split())[i])

if lines[7][0]=='s' or lines[7][0]=='S':
    # the vasp file is already tagged with 'selective dynamics'
    flag_selective=1
else:
    # the vasp file is pristine. Should add 'Selective dynamics'
    flag_selective=0

#### INDEXING INDIVIDUAL ATOMS #####
### Si1 --> 9, Si2 -->10, etc
index_dict=dict()
# index_counter to read/write in sequence from original vasp
index_counter = 8 + flag_selective
for i, element in enumerate(element_list):
    for j in range(num_atom_dict[element]):
        label_name=element+str(j+1)      
        index_dict[label_name] = index_counter
        index_counter += 1

##############################################################################
########## Ask user for tag and atoms to apply
##############################################################################

# If atoms are not states in the command lines
if len(sys.argv) == 2:
    print('Which atom do you want to get coordinate?')
    print("Separate with spaces: ex) 'Si1 Si2 O1 O3'")
    atom_input=input()
    atom_list=atom_input.split()
elif len(sys.argv) >= 3:
    atom_list=sys.argv[2:]
    

if len(atom_list)==0:
    print('No atom labels are provided')
    sys.exit()


##############################################################################
########## Printing
##############################################################################

print("   Atom label |     x       y       z   ")
print("   ───────────┼─────────────────────────")
# list of (3,) numpy arrays
coordinate_list=[]
for atom_label in atom_list:
    temp_list = lines[index_dict[atom_label]].split()
    coordinate_list.append(np.array([float(i) for i in temp_list]))
    temp_array=np.array([float(i) for i in temp_list])
    print("   {0:>10} |  {1:6.4f}  {2:6.4f}  {3:6.4f}".format(atom_label,temp_array[0],temp_array[1],temp_array[2]))


