import re
import os
import sys
"""Note: Alignment file must have no line breaks in the middle of each sequence."""

trans_dict = { "ALA": "A", "ARG": "R", "ASP": "D", "ASN": "N", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}

#Returns raw coordinates of all atoms and strand breaks
def get_coordinates(pdb_file):
    
    open_pdb_file = open(pdb_file, 'rU')
    read_open_pdb_file = open_pdb_file.read()
    
    raw_coordinates = re.findall(r'ATOM\s+(\d+)\s*(\w+).*(\w\w\w)\s\w\s+(\d+)\s+(-*\w+\.\w+)\s+(-*\w+\.\w+)\s+(-*\w+\.\w+)', read_open_pdb_file)
    
    global strand_breaks
    strand_breaks = re.findall(r'TER\s+(\d+)\s+\w+\s\w+', read_open_pdb_file) 

    global coordinates
    coordinates= []
    
    #Selects only the first strand of a multimer for analysis
    for i in range(0, int(strand_breaks[0])-1):
        coordinates.append(raw_coordinates[i])
    
    return coordinates 

def get_protein_name(pdb_file):
    
    open_pdb_file = open(pdb_file, 'rU')
    read_open_pdb_file = open_pdb_file.read()
    name = re.findall(r'\w\w\w\-\d\d\s*([\w\d][\w\d][\w\d][\w\d])', read_open_pdb_file)
    open_pdb_file.close()
    return name
 
   
def num_into_res_dict_maker(coords):
    global num_into_res_dict
    num_into_res_dict = {} #Setting up dictionary that will take in atom numbers and output amino acid type at that position
    for item in coords:
        num_into_res_dict[int(item[3])] = item[2]



def unique_res_number_maker(coords):
    all_res_numbers = []
    
    for x in coordinates:
        all_res_numbers.append(x[3])
    
    #print all_res_numbers
    global unique_res_numbers
    unique_res_numbers = []
    
    pop_ind = 0
    for num in all_res_numbers:
        if num not in unique_res_numbers:
            #unique_res_numbers.append(all_res_numbers.pop
            #print all_res_numbers[pop_ind]
            unique_res_numbers.append(all_res_numbers.pop(pop_ind))
        pop_ind += 1
    
        
def text_writer(nomine):
    seq_list = []
    for num in unique_res_numbers:
        seq_list.append(trans_dict[num_into_res_dict[int(num)]]) #Takes the res number, returns res type, and then converts to one letter code
    #print seq_list
    seq = "".join(seq_list)
    #print seq
    text = ">" + nomine + "\n" + seq
    #print text
    return text
        
def main():
    get_coordinates(sys.argv[1])
    num_into_res_dict_maker(coordinates)
    unique_res_number_maker(coordinates)
    nombre = get_protein_name(sys.argv[1])[0]
    fasta = text_writer(nombre)
    print "NOMBRE:" + nombre
    
    try:
        if sys.argv[2]:
            if sys.argv[2] == "file":
                outf = open(nombre + '_seq.txt', 'w')
                outf.write(fasta)
                outf.close()
            else: 
                print fasta
                print "Enter 'file' in third system argument to write"
    except IndexError:
        print fasta
        print "Enter 'file' in third system argument to create a file" + "\n"

if __name__ == '__main__':
    main()