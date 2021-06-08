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
from pymatgen.core import Structure, Element
import pandas as pd
#from pymatgen.util.coord import pbc_shortest_vectors

#sys.argv=['test','PERTURBED_initial_H-.vasp','0,0,0','r=3']
#sys.argv=['test','CONTCAR','Si1','r=3.5']
#print(sys.argv[0])

radius=None
for i in range(len(sys.argv)):
    if 'r=' in sys.argv[i]:
        radius=float(sys.argv[i].split('=')[1])
        # Once parsed, remove from the sys.argv
        sys.argv.remove(sys.argv[i])
if radius == None:
    radius=3.5

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
    print('Of which single atom do you want to perturb the position of neighbors?')
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

original_SG = struct.get_space_group_info(symprec=1e-2,angle_tolerance=0.1)[0]
df=create_df(struct)

lattice_vector = np.array([struct.lattice.a,struct.lattice.b,struct.lattice.c])

print("   Atom label |     x       y       z    | Distance ")
print("   ───────────┼──────────────────────────┼──────────")
print("    Before    | Space group: {0:<8}    |".format(original_SG))
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
            print("      └ {0:<6}| {1: 5.4f} {2: 5.4f} {3: 5.4f}  |  {4:6.4f}"
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
            print("      └ {0:<6}| {1: 5.4f} {2: 5.4f} {3: 5.4f}  |  {4:6.4f}"
                  .format(df[df['pmg_site']==A2]['atom_label'].iloc[0],
                          temp_array[0],temp_array[1],temp_array[2], A2_tuple[1] ))
            

print_list=[]
# If label is atom label like 'Si1'
if atom_list[0][0].isalpha():
    A1_label=atom_list[0]
    A1=df[df['atom_label']==A1_label]['pmg_site'].iloc[0]
    neighbors=struct.get_neighbors(A1,r=radius)
    for A2 in neighbors:
        A2.to_unit_cell(in_place=True)
        A2_label=df[df['pmg_site']==A2]['atom_label'].iloc[0]
        # random numbers
        random_vector=np.random.ranf([3,])-0.5
        # coordination dimention is (,3)
        scaling_factor = min(0.1,0.1/A1.distance(A2) )
        random_vector = scaling_factor*random_vector/np.linalg.norm(random_vector)
        # perturbation: diaplace by 0.1 \AA for atoms 1.0 \AA apart from the position
        A2.frac_coords = A2.frac_coords + random_vector/lattice_vector
        temp_array=A2.frac_coords
        struct.sites[df[df['atom_label']==A2_label]['site_index'].iloc[0]].frac_coords=temp_array
        print_list.append("      └ {0:<6}| {1: 5.4f} {2: 5.4f} {3: 5.4f}  |  {4:6.4f}"
                      .format(A2_label,
                              temp_array[0],temp_array[1],temp_array[2], A1.distance(A2) ))
        # print("      └ {0:<6}| {1: 5.4f} {2: 5.4f} {3: 5.4f}  |  {4:6.4f}"
        #               .format(A2_label,
        #                       temp_array[0],temp_array[1],temp_array[2], A1.distance(A2) ))
# If entry is coordinate
else:        
    position_array=np.array(eval(atom_label))
    neighbors=struct.get_sites_in_sphere(position_array,r=radius)
    for i, A2_tuple in enumerate(neighbors):
        A2=A2_tuple[0]
        A2.to_unit_cell(in_place=True)
        A2_label=df[df['pmg_site']==A2]['atom_label'].iloc[0]
        random_vector=np.random.ranf([3,])-0.5
        # coordination dimention is (,3)
        # scaling factor forces it to be less than 0.1 ang.
        # + 0.001 is to prevent division by zero
        scaling_factor=min(0.1, 0.1/(0.001 + A2.distance_from_point(position_array*lattice_vector)))
        random_vector = scaling_factor * random_vector/np.linalg.norm(random_vector)
        A2.frac_coords = A2.frac_coords + random_vector/lattice_vector
        temp_array=A2.frac_coords
        struct.sites[df[df['atom_label']==A2_label]['site_index'].iloc[0]].frac_coords=temp_array
        new_distance = (struct.get_sites_in_sphere(position_array,r=radius))[i][1]
        print_list.append("      └ {0:<6}| {1: 5.4f} {2: 5.4f} {3: 5.4f}  |  {4:6.4f}"
              .format(A2_label, temp_array[0], temp_array[1], temp_array[2],
                      new_distance))
        # print("      └ {0:<6}| {1: 5.4f} {2: 5.4f} {3: 5.4f}  |  {4:6.4f}"
        #       .format(A2_label, temp_array[0], temp_array[1], temp_array[2],
        #               new_distance))
        

new_SG = struct.get_space_group_info(symprec=1e-2,angle_tolerance=0.1)[0]
print("   ───────────┼──────────────────────────┼──────────")
print("    After     | Space group: {0:<8}    |".format(new_SG))

if atom_list[0][0].isalpha():
    A1=df[df['atom_label']==atom_label]['pmg_site'].iloc[0]
    temp_array=A1.frac_coords
    print("   {0:>5}      | {1: 5.4f} {2: 5.4f} {3: 5.4f}  |  ------ ".format(atom_label,temp_array[0],temp_array[1],temp_array[2]))
else:
    print("     coords   | {0: 5.4f} {1: 5.4f} {2: 5.4f}  |  ------ ".format(position_array[0],position_array[1],position_array[2]))
for i in print_list:
    print(i)
    
##############################################################################
########## Write POSCAR
##############################################################################  
# update DataFrame
label_df=create_df(struct)

# Idenfity lattice vectors
lattice = struct.lattice.matrix
# Get element of each atoms
species_list=struct.species
# Remove repeated entry in the element list
reduced_species=[str(i) for n, i in enumerate(species_list) if i not in species_list[:n]]

if sys.argv[1][-5::]=='.vasp':
    filename = sys.argv[1][:-5]+'_PERTURBED'+sys.argv[1][-5::]
else:
    filename = sys.argv[1]+'_PERTURBED'
    
out_file=open(filename,'w')
out_file.write("Perturbed POSCAR file from {0}\n".format(sys.argv[1])) #first comment line
out_file.write("1.0\n") # scale 
# Print lattice part
for i in range(np.shape(lattice)[0]):
    out_file.write("{0:20.10f} {1:20.10f} {2:20.10f}\n".format(lattice[i,0],lattice[i,1],lattice[i,2]))
# Print elements
out_file.write("  "+" ".join('%4s' % entry for entry in reduced_species))
out_file.write("\n")
# Print the number of atoms for each element
num_each_element=[species_list.count(Element(i)) for i in reduced_species]
out_file.write("  "+" ".join('%4d' % entry for entry in num_each_element))
out_file.write("\n")

out_file.write("Direct\n")
for element in reduced_species:
    site_list=label_df[label_df['element']==element]['pmg_site']
    for site in site_list:
        out_file.write("  "+"        ".join('%.10f' % entry for entry in site.frac_coords)+'\n')
out_file.close()

print("")
print('  Perturbed structure was written as {0}'.format(filename))
