"""
## Command:
python selective_dynamics.py [input_filename] [output_filename:optional]

"""
import sys
import numpy as np

### temporary_line
sys.argv=[None,'SelDy_TTT_All.vasp']

##############################################################################
########## Read POSCAR file and indexing each line with atom label
##############################################################################
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
print('Which tag do you want to apply for selective dynamics?')
print('ex) TTT, FFF, FFT')
SD_tag=input()
defect_input=input()

print('For which atoms do you want to apply this tag?')
print("Separate with spaces: ex) 'Si1 Si2 O1 O3'")
print("For all atoms, write 'All' ")
atom_input=input()

if atom_input == 'All':
    apply_atom_list = [atom_input]
    apply_index_list = []
else:
    apply_atom_list=atom_input.split()
    apply_index_list=[index_dict[i] for i in apply_atom_list]








##############################################################################
########## Write output POSCAR file
##############################################################################
if len(sys.argv) >= 3:
    out_filename=sys.argv[2]
else:
    out_filename='SelDy'+'_'+SD_tag+'_'+"_".join('%s' % entry for entry in apply_atom_list)+'.vasp'
out_file=open(out_filename,'w')

# Write header + lattice part
# Header is up to 7th line if 'selective dynamics' is not stated
# up to 8th if 'selective dynamics' is already stated
for i in range(0,8+flag_selective):
    out_file.write(lines[i])

# state 'selective dynamics' if not already stated
if flag_selective == 0:
    out_file.write("Selective dynamics \n")
    
for i in range(8+flag_selective,max(index_dict.values())):
    for j in (lines[i].split())[0:3]:
            out_file.write('  '+j)
            
    #out_file.write('  ')
    if i in apply_index_list or atom_input=='All' or atom_input=='all':
        # for j in (lines[i].split)[0:3]:
        #     out_file.write('  '+j)
        out_file.write('   {0}   {1}   {2}'.format(SD_tag[0],SD_tag[1],SD_tag[2]))
    elif flag_selective == 1:
        tag_list=lines[i].split()[3:6]
        out_file.write('   {0}   {1}   {2}'.format(tag_list[0],tag_list[1],tag_list[2]))
    else:
        # Default is to be TTT
        out_file.write('   T   T   T')
    out_file.write('\n') # end of line
# end of file
out_file.close()


print("Following POSCAR is written: {0}".format(out_filename))

