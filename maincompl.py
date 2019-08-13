
from modules_arg import *
from functions import *

if verbose: 
    sys.stderr.write("If no input (-i) provided, current directory will be used\n")
    sys.stderr.write("If no output (-o) provided, resulting pdb files will be created in current directory\n")
    sys.stderr.write("If no distance (-d) provided, interactions between chains will be considered if chains are 6Ä or less apart\n\n")

##Apply the function that filters the input and create a list with the correct input files (only PDB)
pdbfiles=list(filterinput(folder,verbose,outfolder))

##Parse the pdb files and filter those that have 2 interacting chains 
struc, seq_ids=parser_filter(pdbfiles,verbose,folder, outfolder,in_dist)

##Superimposition and reconstruction of the complex 

# This function returns a list of tuples
# with all the possible combinations of chains to superpose
refs = tuple_superposition(seq_ids, verbose)
if verbose: 
	sys.stderr.write("%s\n\n" % refs)

# And this function is a generator of all the possible models
result= Model(struc, refs, verbose)

# This loop starts a counter and for each
# element of the result generator
# it initializes the PDBIO object and starts a structure defined 
# from the function, and saves it into a pdb named Model+the counter value
c=1
if verbose: 
	sys.stderr.write("The generator ended, now creating the .pdb files: \n\n")
for element in list(result):
	save= PDBIO()
	save.set_structure(element.structure)
	save.save(outfolder+"/"+"Model"+str(c)+".pdb")
	if verbose: 
		sys.stderr.write("File %s was created!!!\n\n" % ("Model"+str(c)+".pdb"))
	c +=1