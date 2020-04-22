"""
Start is defined as how many residues are "past" (e.g. say 2 if you want 2 A's in the beginning) 
and end represents the last crossover residue (the last B)
system inputs: (python) (function_this) (length) (crossover_start) (crossover_end) (name)
"""

import sys

def cross_over_file_maker(length, start, end):

    list = []

    for a in range(0, int(start)):
        list.append("A")    

    for b in range(int(start), int(end)):
        list.append("B")
        
    for c in range(int(end), int(length)):
        list.append("C") 
    
    sequence = "".join(list)
    return sequence

def main():

    crossover = cross_over_file_maker(sys.argv[1], sys.argv[2], sys.argv[3])
    print crossover
    print sys.argv[2]
    
    outf = open(str(sys.argv[4]) + ".txt", 'w')
    outf.write(crossover)
    outf.close()
    
if __name__ == '__main__':
  main()