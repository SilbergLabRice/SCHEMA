# SCHEMA
SCHEMA algorithm for protein structure-guided recombination 

The pdb.search.and.grab.test.2.2.1.py takes in three data sets: 
pdb structure file, an alignment file, and a crossover file.

All three should be in .txt format.

The pdb structure can be downloaded in this format from the 
pdb website.

The alignment file must be in fasta format and not have white
spaces between residues. It must have an alignment of the structure
sequence and the two parent chimera sequences.  The included
sequence_grabber.py script can obtain the sequence from the pdb
file. I often use MAFFT, selecting FASTA format, to make the 
alignments. 
Here is an example:

\>AAAA 
TWYRAPSQIL
\>seq_one 
TWYRAPSQID 
\>seq_three 
ARDVOARDVO

The crossover file can be made via cross_maker.py. Instructions
for use are included in comments within the cross_maker.py file.
You will have the option to give the length of residues for your
protein, the start and ending indices of the crossover, and name
for the file.  
General inputs: (python) (function_this) (length) (crossover_start) (crossover_end) (name)

For example:
"python cross_maker.py 10 2 8 T1" in the command line will yield

AABBBBBBCC

As a file named T1.

General notes:
This script will generate disruption in the form of E.  It can
tolerate holes in sequence structure, as long as the pdb file 
follows the general convention of making its atom numbering
sequentially continuous even when there are gaps in structure
(I.e. the missing "atoms" are not counted in the file.  If the
final atom # is 1000, then there will be 1000 atom entries.)

It is generally best to take the sequence from the pdb file 
using the sequence_grabber.py script, rather than downloading
the sequence file from the pdb website.  Older pdb structures
will have discrepancies between listed structure sequences and
the sequences that are actually present in the structure file.

Older pdb files may have unforeseen complications and require
more scrutiny.

How it works:

1. The pdb file is read.  For each atom, the atom's number (index in chain),
the atom type, the residue in which the atom resides (residue index),
the parent residue's type, and the atom's coordinates are extracted.
This information is contained in the list of lists known as coordinatess.

2. The align file is read.  The residues from the the structure sequence
and two parent sequences are sorted into a list of lists.

3.  True beginning and true end of chain are extracted from align file.
This is a small step that accounts for the holes at the beginning and end
of structure sequences and helps set up indices for later steps.

4. Atom to res and res to atom dictionaries are made.  Using the coordinates
list of lists, two dictionaries are made.  One that takes in atom number
(position in chain) and gives out residue number (residue position in chain)
and one that takes in residue number and gives out a list of the atom numbers
it contains.

5. The distances between all atoms are calculated, using coordinates, and
sorted into matrix "d".

6. Matrix d is used to make another matrix, contacts.  Contacts is n x n
where n is the number of atoms in the chain.  All of contacts elements are
0 initially, and replaced with 1 if the two atoms which represent the indices
of that element have a distance less than 4.5 A.  Contacts therefore represents
a binary accounting of all atom contacts less than 4.5 angstroms.

7. Matrix contacts is converted into res_contacts.  What this means is that another
n x n matrix is created.  This time n represents the number of residues in 
structure chain.  This matrix is also populated only with 0's.  Again the 
res_contacts matrix is iterated over and for each element whose indices,
which represent residues, have contacting atoms the 0 is replaced with a 1.
This is done using dictionaries created in 3 and the contacts matrix

8. Residues that don't match between given sequence and sequences structures 
are accounted for. All of these residues' columns and rows are given a 0 in
the res_contacts matrix.

9. Non-breakable interactions are accounted for. Non-breakable interactions
are defined as those which are identical in both parents.  This step is done
using the align file list of lists. All of these residues' columns and rows are 
given a 0 in the res_contacts matrix. 

10. The bottom half of the matrix is removed.  I don't want to double count
interactions, so in this step I remove the bottom left half of the matrix
(below the diagonal) by giving all elements in that region 0.

11. Disruptions are counted.  Using the crossover file, the regions in which 
the crossover occurs are counted in the residue contacts matrix and added up.
