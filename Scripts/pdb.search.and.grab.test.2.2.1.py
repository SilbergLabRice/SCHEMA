"""Update 2.2.1
-Eliminated excess print commands
-Updated Comments"
"""

import re
import os
import sys
"""Note: Alignment file must have no line breaks in the middle of each sequence."""

#Returns raw coordinates of all atoms and strand breaks
def get_coordinates(pdb_file):
    
    open_pdb_file = open(pdb_file, 'rU')
    read_open_pdb_file = open_pdb_file.read()
    
    raw_coordinates = re.findall(r'ATOM\s+(\d+)\s*(\w+).*(\w\w\w)\s\w\s+(\d+)\s+(-*\w+\.\w+)\s+(-*\w+\.\w+)\s+(-*\w+\.\w+)', read_open_pdb_file)
    
    global strand_breaks
    strand_breaks = re.findall(r'TER\s+(\d+)\s+\w+\s\w+', read_open_pdb_file) 

    global coordinates
    coordinates= []
    
    #Selects only the first strand of a polymer for analysis
    for i in range(0, int(strand_breaks[0])-1):
        coordinates.append(raw_coordinates[i])
    
    return coordinates 

# Takes alignment file and list of residue identities of chimera parents and given structure sequence    
def get_alignment_matrix(mafft_file):
    
    open_mafft_file = open(mafft_file, 'rU')
    read_open_mafft_file = open_mafft_file.read()
    
    raw_align = re.findall(r'>(\w+)\s+([\w-]+)', read_open_mafft_file)
    
    global align
    align = []
    for s in raw_align:
       align.append([s[0]] + list(s[1]))


# Gets the true length of the sequence structure- this serves mainly to exclude holes at the end
def get_structure_true_length(matrix_align):
    struc_seq = matrix_align[0] #[1:] because must exclude name
    n = 1
    while struc_seq[-n] == "-":
        n +=1
    
    global true_length
    true_length = len(matrix_align[0][:-n])
    """print "True!"
    print true_length"""
# The converse of the previous function- this helps me set up indices later that will exclude holes in the beginning    
def get_true_start(matrix_align):
    struc_seq = matrix_align[0]
    n = 0
    while struc_seq[n+1] == "-":
        n +=1
        
    global true_start
    true_start = n - 1

# These dictionaries are needed for removing the nonbreakable interactions
# The first takes in atom number (position in chain) and outputs amino acid residue number (position in chain)
# The second does the opposite
# I subtract one a bunch here because the indices of a list begin with zero 
# and so I need to move the value of all atom and residue values down by one

def atom_into_res_dict_maker(coords):
    global atom_into_res_dict
    atom_into_res_dict = {} #Setting up dictionary that will take in atom numbers and output amino acid number


    for item in coords:
        #print item[0]
        atom_into_res_dict[int(item[0])-1] = int(item[3])-1 #Subtracting one to put in terms of list indices which start with 0
        


def res_into_atoms_dict_maker(coords):
    global res_into_atoms_dict
    res_into_atoms_dict = {} #Setting up dictionary that will take in amino acid number and output list of atoms
    pdb_count = int(coords[0][3])-1
    q= pdb_count
    seq_count = true_start # Refers to dictionary entry defined by position of residues in my alignment file
    #print "True Start:" + str(seq_count)
    for item in coords:
        if int(item[3]) != int(seq_count):
            """print "No Match"
            print "Item[3] AKA residue number:" + item[3]
            print "Seq count:" + str(seq_count)"""
            pdb_count = int(item[3])
            seq_count += abs(int(item[3])-int(seq_count)) #adds the difference of the alignment file residue index and the stated pdb file residue number
        else:
            print "Match"
        try:
            res_into_atoms_dict[seq_count-1].append(int(item[0])-1)
        except KeyError:
            res_into_atoms_dict[seq_count-1] = [int(item[0])-1]

# Calculates distances between all atoms
def distance(coords):
    global d
    #Creates a square matrix with sidelength as long as coordinates with elements that are the distances between the atoms that represent the indices of the matrix
    d=[ [ ((float(x[4]) - float(y[4]))**2 + (float(x[5]) - float(y[5]))**2 + (float(x[6]) - float(y[6]))**2)**0.5 for x in coords] for y in coords]

# Takes d distance matrix and makes a matrix of identical size with all entries equal to 0
# In the positions with distances less than 4.5 angstroms, it replaces the 0 with 1
def contact_matrix(distances):
    global contacts
    print "MEEE! -1 coordinates:" + str(coordinates[-1][0])
    contacts = [[0 for count in range(int(coordinates[-1][0]))] for count in range (int(coordinates[-1][0]))] 
    #Above creates an n by n matrix of zeroes where n is number of atoms in the chain. Retrieved from the first item of the last list in coordinates. Used the int() to convert to integer.
    a =0
    ind_x=0
    ind_y=0
    bad_errors=0
    for x in distances:
        for y in x:
            if float(y) < 4.5:
                if contacts[ind_x][ind_y] != 0:
                    bad_errors +=1
                contacts[ind_x][ind_y] = 1
                a += 1
            else:
                if contacts[ind_x][ind_y] != 0:
                    bad_errors +=1
                contacts[ind_x][ind_y] = 0
                a += 1
            ind_y+=1
         
        ind_x+=1
        ind_y=0
  
    #print contacts

    
    print str(a) + "--result"
    all = 0
    global zero
    global one
    zero = 0
    one = 0

    for x in contacts:
        for y in x:
            all += 1    
            if y == 0:
                zero +=1
            elif y == 1:
                one +=1
                
    print str(zero) + "--zeroes"
    print str(one) + "--ones"
    print str(bad_errors) + "--errors"
    print "Len of contacts:" + str(len(contacts))
 
# Converts the atom contact matrix to a residue contact matrix using the established dictionaries   
def atoms_into_res_contact(atoms_contact, matrix_align):
    
    #Creates an n by n matrix of zeroes where n is number of residues in the structures' chain. 
    #Retrieved from the first item of the last list in coordinates. Used the int() to convert to integer.
    global res_contacts
    res_contacts = [[0 for count in range(int(true_length))] for count in range (int(true_length))] 
    
    all_res=[]
    for res in matrix_align[0][true_start:]:
        all_res.append(res)
    print "RESS"
    print all_res
    
    #Will loop over all residues
    ind_atom=0
    print true_start
    print atom_into_res_dict
    print "Res_into_atoms_dict"
    print res_into_atoms_dict
    print "len"
    print len(res_contacts)
    for res in range(true_start, true_length): #Iterates over an index the size of the residue in the chain
        try:
            for atom in res_into_atoms_dict[res]: #Takes the index and returns a list of atoms from the res_into_atoms_dic
                print "Res:" + str(res)
                print "Atom:" + str(atom)
                for interactions in contacts[int(atom)]:
                    if interactions == 1:
                        #print "Yes res:" + str(res)
                        #print "Yes second res:" + str(atom_into_res_dict[ind_atom])
                        res_contacts[res][int(atom_into_res_dict[ind_atom])] = 1 
                        # Replaces the 0 with 1 at cross section of two interacting residues.
                        # Takes in res number and the output of the partner atom's number (given by ind_atom) into atom_into_res_dict
                    ind_atom += 1           
                ind_atom = 0   
        
        except KeyError:
            continue
    
    
    global res_zero
    global res_one
    res_zero = 0
    res_one = 0

    for x in res_contacts:
        for y in x:    
            if y == 0:
                res_zero +=1
            elif y == 1:
                res_one +=1
    
    #print res_contacts
    print str(res_zero) + "--res_zeroes"
    print str(res_one) + "--res_ones"
    return res_contacts

# This function will remove contacts for residues that don't match in the given structure and its accepted sequence    
def remove_struc_seq_mismatches(matrix_align, matrix_contacts):
    global cleaned_matrix
    cleaned_matrix = []
    cleaned_matrix = res_contacts 

    struct  = matrix_align[0][1:true_length+1]
    seq_one = matrix_align[1][1:true_length+1]
    ind=0
    for res in seq_one:
        print #"Checking %s @ %s" % (str(res), ind)
        #If sequences don't match or if the structure is blank at the residue in question, start on path to give 0 for all those residue's interactions 
        try:
        
            if res != struct[ind]:
                print "*****Mismatch: %s @ %s\n" % (str(res), ind) #all those res atoms get zero in contact matrix
                
                #This SHOULD remove both the column and row of the mismatched residue's interactions
                
                #Removes column interactions
                int_ind=0
                for ele in cleaned_matrix[ind]:

                    cleaned_matrix[ind][int_ind] = 0
                    int_ind +=1                
                
                #Removes row interactions
                for resi in cleaned_matrix:
                    resi[ind] = 0
                
            else:
                print "No Mismatch\n"
                
        except KeyError:
            print "No residue structure\n"
        
        ind+=1            
        print ind
 
    #print cleaned_matrix
    
    global zero_cleaned
    global one_cleaned
    zero_cleaned = 0
    one_cleaned = 0

    for x in cleaned_matrix:
        for y in x:   
            if y == 0:
                zero_cleaned +=1
            elif y == 1:
                one_cleaned +=1
                
    print str(zero_cleaned) + "--zeroes_cleaned"
    print str(one_cleaned) + "--ones_cleaned"

# This removes non-breakable interactions in the contact matrix
# Non breakable interactions are those in which one or more of the atoms appears in both parents
def remove_non_breakables(matrix_align, matrix_contacts):
    
    seq_one = matrix_align[1][1:true_length]
    seq_two = matrix_align[2][1:true_length]
    
    ind=0
    
    #setting up distinct breakable residues matrix so I will have both to work with
    global breakable_matrix
    breakable_matrix = matrix_contacts

    for res in seq_one:
        #If sequences match, start on path to give 0 for all those residue's atoms
        if res == seq_two[ind]:
            #print "*****Match: %s @ %s\n" % (str(res), ind) #all those res atoms get zero in contact matrix
            
            #This SHOULD remove both the column and row of the matched residue's interactions
                
            #Removes column interactions
            int_ind=0
            for ele in breakable_matrix[ind]:
                breakable_matrix[ind][int_ind] = 0
                int_ind +=1                
                
            #Removes row interactions
            for resi in breakable_matrix:
                resi[ind] = 0
        else:
            print "No Match\n"
        ind+=1
    
      
    all = 0
    zero_after = 0
    one_after = 0
    for x in breakable_matrix:
        for y in x:
            all += 1    
            if y == 0:
                zero_after +=1
            elif y == 1:
                one_after +=1
    print str(zero) + "--zeroes"
    print str(one) + "--ones"
    print str(res_zero) + "--res_zeroes"
    print str(res_one) + "--res_ones"            
    print str(zero_cleaned) + "--zeroes_after_mismatch"
    print str(one_cleaned) + "--ones_after_mismatch"    
    print str(zero_after) + "--zeroes_after_match"
    print str(one_after) + "--ones_after_match"   

# This eliminates half of the matrix so there is no double counting
def top_half_grabber(scrubbed_matrix):

    global final_matrix
    final_matrix= scrubbed_matrix
    
    ind_x=0 # x is the x coordinate of the imaginary matrix
    ind_y=0 # y is the y coordinate of the imaginary matrix, I believe

    
    for x in final_matrix:
        for y in x:
            if ind_y >= ind_x:
                """print "Ind_x:" + str(ind_x)
                print "Ind_y:" + str(ind_y)
                print "Value:" + str(final_matrix[ind_x][ind_y])"""
                
                final_matrix[ind_x][ind_y] = 0

            ind_y+=1         
        ind_x+=1
        ind_y=0 # The y has to be reset
    
    #print final_matrix[72]
    
    zero_final = 0
    one_final = 0
    for x in final_matrix:
        for y in x:   
            if y == 0:
                zero_final +=1
            elif y == 1:
                one_final +=1
    print str(zero_final) + "--zeroes_final"
    print str(one_final) + "--ones_final"

def disruption_counter(crossover_string):

    open_crossover_file = open(crossover_string, 'rU')
    read_open_crossover_file = open_crossover_file.read()
    
    #Grabs the whole file which should only be A's and B's without spacing of any kind
    crossover = re.findall(r'\w+', read_open_crossover_file)
    crossover_list = list(''.join(crossover)) #converts 1 element list into string and then makes a list of each character
    
    count_letter = "A"
    start = -1 # start at negative one so start will represent the index of my first "B"
    
    for letter in range(0, len(crossover_list)):
        start +=1
        count_letter = crossover_list[start]
        if count_letter == "B":
            break
    
    print "Start:" + str(start)
        
    end = -1 # start at negative one so start will represent the index of my first "B"
    
    for letter in range(0, len(crossover_list)):
        end +=1
        count_letter = crossover_list[end]
        if count_letter == "C":
            break
    
    print "End:" + str(end)
    print crossover_list[start:end]
    
    #Now, finally for the final count
    disruption_count = 0
    it_count= 0
    
    # Column Count
    for x in final_matrix[start:end]:
        for y in x[:start]:# I have to control for overlap by limiting how far down the column count will go
            disruption_count += y
            it_count += 1
            #print "Value:" + str(y)
            
    #print "Disruptions:" + str(disruption_count)
    print end-start    
    
    # Row Count
    for x in final_matrix[end:]:#I have to avoid counting the interactions that are carried over by the crossover
        for y in x[start:end]:
            disruption_count +=y
            it_count += 1

    #print "Iterations:" + str(it_count)
    print "Disruptions:" + str(disruption_count)

def main():
    get_coordinates(sys.argv[1])
    get_alignment_matrix(sys.argv[2])
    get_structure_true_length(align)
    get_true_start(align)
    res_into_atoms_dict_maker(coordinates)
    atom_into_res_dict_maker(coordinates)
    
    
    distance(coordinates)
    contact_matrix(d)
    
    atoms_into_res_contact(contacts, align)
   
    remove_struc_seq_mismatches(align, contacts)
    remove_non_breakables(align, cleaned_matrix)
    top_half_grabber(breakable_matrix)
    disruption_counter(sys.argv[3])
if __name__ == '__main__':
    main()