from Bio.PDB import * 
from glob import glob
import sys  
import os
import Bio.PDB as pdb
import random
import string
ppb = pdb.Polypeptide.PPBuilder() 
import shutil 
import argparse
from modules_arg import *
import string
import random
import re 
from Bio.pairwise2 import align



def filterinput(inputdir,verbose,output):
	"""Function that filters the input and creates a list with the correct input files (only PDB)."""
	c=0
	if inputdir:
		path =os.path.dirname(os.path.realpath(__file__))
		pdbfolder= inputdir
		for filename in os.listdir(pdbfolder): 
			if filename.endswith(".pdb"):
				c=c+1
				yield filename
			else:
				if verbose:
					sys.stderr.write("The file %s is not in a PDB format\n" % filename)
	elif os.path.isfile(inputdir):
		yield inputdir

	else: ##if a directory is not provided it takes the current directory
		path = "./"
		for files in os.listdir(path):
			if files.endswith(".pdb"):
				c=c+1
				yield files
	if verbose:
		sys.stderr.write("Number of PDB files found in the input folder: %d  \n" % c)

	# #Knows if the output folder has been already created or not, if not the program creates the folder. 
	if output:
		if not os.path.isdir(path+'/'+output):
			os.mkdir(path+'/'+output)
	 	#elif:
	 	#	sys.stderr.write("The output folder %s already exist \n" % output)
	else:
		raise IOError("No output folder provided")


def chains_interaction(structure,in_dist):
		"""Function that filters the files that contain chains that are interacting """
		chains = list(structure.get_chains())

		## We are looking for files with 2 chains
		if len(chains) == 2:	
			atoms=[]

			# obtain info of each sequence.
			for chain in chains:
				atoms.append(chain.get_atoms())			
			#look if the two chains are close to interact.	
			CA1= list(filter(lambda x: x.get_id()=="CA", atoms[0]))
			CA2= list(filter(lambda x: x.get_id()=="CA", atoms[1]))
			ns= NeighborSearch(CA1)
			interactions =[]
			for atom in CA2:
				interact = ns.search(atom.get_coord(), in_dist) 
				interactions.extend([tuple(sorted([str(atom.get_parent().resname), str(x.get_parent().resname)])) for x in interact])
			if interactions:
				return interactions

def parser_filter(pdbfiles,verbose,folder,outfolder,in_dist):
	"""Parse the pdb files, filter using the chains_interaction function and dna_prot_filter. So it's going to filter 
	the pdb files that only have  two interacting chains. Also is going to remove the atomic information of nucleotides (DNA) 
	from protien and create pdb files without the dna information, in order to create only a protein model."""
	f=[]
	s={}	
	seq_ids={} 
	used_letters=[]
	path =os.path.dirname(os.path.realpath(__file__))
	
	if verbose:
			sys.stderr.write("Parsing the pdb files, returning the structure objects of each of them\n\n")
			sys.stderr.write("Selecting files with two interacting chains\n")
	for file in pdbfiles:
		ID=file.strip(".pdb")
		fulldir_path=(path+'/'+folder+'/'+file)
		#PDBParser is a class and get.structure is an instance method of the PDBParser class.Return the structure so we can know the chains.
		structure=PDBParser(QUIET=True).get_structure(ID, fulldir_path)
		#Only the files that have two chains interacting and that chains have interaction will be used:
		if verbose:
			sys.stderr.write("Looking if %s contains two interacting chains\n" % file)
		interactions=chains_interaction(structure,in_dist)
		if interactions:
			#From the files that have intercation (the ones that we are interested in) we are going to look if they have nucleic acids in their sequences in order to creat 'clean' pdb files.
			dna_prot_filter(file,verbose,folder,outfolder)
			if verbose:
				sys.stderr.write("File %s is selected\n\n" % file)
			
			
	

	#We have all the clean files without DNA information or extra information in  filtered_pdb_files inside the output file 
	for file in os.listdir(path+'/'+outfolder+'/'+'filtered_pdb_files'):
		if file.endswith(".pdb"):
			ID=file.strip(".pdb")
			fulldir_path=(path+'/'+outfolder+'/'+'filtered_pdb_files'+'/'+file)
		structure=PDBParser(QUIET=True).get_structure(ID, fulldir_path)
		s[ID]=structure
		c=0
		chains = list(structure.get_chains())
		for chain in chains:
				c+=1
				seq_ids[ID+str(c)] = chain
	if verbose:
		if seq_ids:
			sys.stderr.write("Files processed correctly, look for filtered files in filtered_pdb_files directory\n\n")
			sys.stderr.write("Starting to reconstruct the macrocomplex...\n")
	return s, seq_ids
	
	

def dna_prot_filter(file,verbose,folder,outfolder):
	""" This function creates a folder (filtered_pdb_files) inside the output folder, inside this folder are the filtered files (without 
	information for nucleic acids), these files will be used later to generate the macrocomplex"""
	path =os.path.dirname(os.path.realpath(__file__))
	if not os.path.exists(path+'/'+outfolder+'/'+'filtered_pdb_files'):
		os.mkdir (path+'/'+outfolder+'/'+'filtered_pdb_files')

	#Create output files
	
	name_out=outfolder+'/'+'filtered_pdb_files'+'/'+'filtered_'+file
	new_pdb= open(name_out,'w')
			
	file=path+'/'+folder+'/'+file
	file = open(file, 'r')

	found = False
	for line in file:
		if line.startswith('ATOM'):
			line = line.strip()
			search = re.search(r'[A-Z]*\s*\d{1,5}\s*[A-Z]{1,3}[0-9]?\'?\s*([A-Z]{1,3}).*', line)
			residue = search.group(1)
			if residue in dna_in_pdb:
				while (found==False):
					if verbose:
						sys.stderr.write("The file %s contains nucleic acids \n" % filename)
						found = True
			else:
				new_pdb.write("%s\n" % (line))
	if verbose:
		sys.stderr.write("""   ...Filtering DNA, RNA or other non-protein elements... \n""")
dna_in_pdb = ["DC", "DG", "DT", "DA"]	

def Align(chain1, chain2):
	"""Function that performes the pairwise alignment between the sequences of two chains"""
	for residue in ppb.build_peptides(chain1):
		seq1= residue.get_sequence()
	for residue in ppb.build_peptides(chain2):
		seq2= residue.get_sequence()

	alignment=pairwise2.align.globalxx(seq1, seq2)
	return alignment	


def format_alignment(align1, align2, score, begin, end):
	"""Format the alignment prettily into a string.Since Biopython 1.71 identical matches are shown with a pipe character, mismatches as a dot, and gaps as a space.
	Note that spaces are also used at the start/end of a local
	alignment.
	Prior releases just used the pipe character to indicate the
	aligned region (matches, mismatches and gaps).
	"""
	s = []
	s.append("%s\n" % align1)
	s.append(" " * begin)
	for a, b in zip(align1[begin:end], align2[begin:end]):
		if a == b:
			s.append("|")  # match
		elif a == "-" or b == "-":
			s.append(" ")  # gap
		else:
			s.append(".")  # mismatch
	s.append("\n")
	s.append("%s\n" % align2)
	s.append("  Score=%g\n" % int(score))
	## We assume that we can use either chain lengths 
	## because a difference in chain lengths would also lead
	## to mismatches, so our filter will also be achieved
	score2= int(score)/len(align1)
	s.append("\nIdentity=%f\n" % score2)
	
	return s, score2


def new_chain_id(current_chain):
	alphabet= string.ascii_letters

	#obtain next ID.
	new_id= alphabet[current_chain]
	try:
		return new_id
	except Exception:
		new_id= alphabet[current_chain+1]
		return new_id


dna_in_pdb = ["DC", "DG", "DT", "DA"]

def tuple_superposition(seq_ids,verbose):
	"""
	It is a list that is going to contain a tuple with the combination 
	of the chains that we want to superpose (chains that have an identity > 0,95),
	without repeating the combinations.
	"""
	si_chains=[] 
	chain_counter=0

	for filename, chain in seq_ids.items():
		for filename2, chain2 in seq_ids.items():
			if filename[:-1] != filename2[:-1]:
				aligment=Align(chain, chain2)
				for AA in aligment:
					al, identity=(format_alignment(*AA))
					if identity > 0.95:
						if (chain, chain2) and (chain2, chain) not in si_chains:
							si_chains.append((chain, chain2))

	
		chain.id = new_chain_id(chain_counter)
		chain_counter += 1
	if verbose:
		sys.stderr.write(""" All possible combination of chains are detected:\n""")
	return si_chains
def find_common_atoms(chain_fix, chain_mov):
	'''
	Function that takes 2 chain objects
	previously filtered and considered the same in 
	different .pdb files. Thenreturns two lists
	of atoms, the ones that will be fixed and that 
	will be moved when superimposed
	'''
	for residue in ppb.build_peptides(chain_fix):
		seq_fix= residue.get_sequence()
	for residue in ppb.build_peptides(chain_mov):
		seq_mov= residue.get_sequence()

	atoms_fix = list()
	atoms_mov = list()
	alignment = align.globalxx(seq_fix, seq_mov) #align sequences
	count_fix = 0
	count_move = 0
	for i in range(len(alignment[0][0])):
		if alignment[0][0][i] == '-':
			count_fix -= 1
		elif alignment[0][1][i] == '-':
			count_move -= 1
		elif alignment[0][0][i] != '-' and alignment[0][1][i] != '-':
			atoms_fix.extend(chain_fix.get_atoms())
			atoms_mov.extend(chain_mov.get_atoms())
		
	return atoms_fix, atoms_mov

def get_chain_to_move(chain_fix, chain_mov):
	'''
	Function that will move the chain of the mov pdb file
	'''
	sup = Superimposer()
	# print ("inputs: ", chain_fix, chain_mov)
	# atoms_fix, atoms_mov = find_common_atoms(chain_fix, chain_mov)
	fix2, mov2 = find_common_atoms(chain_fix, chain_mov)
	if len(mov2) < len(fix2):
		fix2= fix2[:len(mov2)]
	elif len(mov2) > len(fix2):
		mov2= mov2[:len(fix2)]
	# sup.set_atoms(atoms_fix, atoms_mov)
	sup.set_atoms(fix2,mov2)
	chain_mov.parent.transform(rot = sup.rotran[0], tran = sup.rotran[1])
	# print ("chain_mov parent: ", list(chain_mov.parent.parent.get_chains()))

	structure = chain_mov.parent.parent
	if len(list(structure.get_chains())) == 2:
		structure2 = structure.copy()
		# print (id(structure2), id(structure))
		# exit()
	
	for chain in structure2.get_chains():
		if chain.id != chain_mov.id:
			# print ("hi", list(structure2.get_chains()))
			moving_chain= chain
			# print ("mov_chain: ", moving_chain)

	return moving_chain

def clashes(structure, new_chain):
	"""Function that filters the files that contain chains that are interacting """
	chains = list(structure.get_chains())

	new_atoms= []
	atoms=[]
	#obtain info of each sequence.
	for chain in chains:
		if chain == new_chain:
			new_atoms.append(chain.get_atoms())
		else:
			atoms.append(chain.get_atoms())			
	#look if the two chains are close to interact.	
	
	CA1= list(filter(lambda x: x.get_id()=="CA", list(new_atoms[0])))
	
	ns= NeighborSearch(CA1)
	
	for list_a in atoms:
		interactions =[]
		CA2 = list(filter(lambda x: x.get_id()=="CA", list(list_a) ))
		
		for atom in CA2:
			interact = ns.search(atom.get_coord(), 0.4) 
			interactions.extend([tuple(sorted([atom.get_parent(), x.get_parent()])) for x in interact])
			
		if interactions:
			inter_chains= (interactions[0][0].get_parent(), interactions[0][1].get_parent())
			return True
			

def Model(files, tup_list,verbose):
	"""
	This function will add one new chain to the model
	for each file found, as each file will contain 2 chains, 
	and the second one will be use as guide to superpose it 
	in the model and know the correct coordinates
	"""
	
	# Some variables to make less iterations over the structures
	used_files=files.copy()
	used_tuples=tup_list.copy()
	
	# We need 2 indented loops to compare 2 structures from different files
	for file1, s1 in files.items():
		del used_files[file1] #We delete files from the copy, 
							  #so that the second iteration won't compare files with itself. 
		
		for file2, s2 in used_files.items():
			
			#The tuple list contains all the combinations previously
			# selected to superpose because of their identity
			for tup in tup_list:
				
				# So, if both elements of the tuple are in the files structures:
				if tup[0] in s1.get_chains() and tup[1] in s2.get_chains():
					
					# Define the chain chosen to stay fixed
					# And the one from the file we will move
					fix= tup[0]
					mov= tup[1]

					# And pass them as arguments for the superimpose 
					# function, that will return the chain we have to add
					# to the new model						
					moving_chain = get_chain_to_move(fix,mov)
					
					# Create a new model (PDBIO object)
					new_model = pdb.PDBIO()
					new_model.set_structure(s1)
					
					#If the chain is not already in the model 
					# (because PDBIO will complain if we try to add
					# chains with same id)
					if moving_chain not in list(new_model.structure.get_chains()):
						# Add the chain returned by the superimposer
						new_model.structure[0].add(moving_chain)
						
						# And calculate if there are clashes between the added 
						# chain and the previous model
						ns = clashes(new_model.structure, moving_chain)
						
						# The previous function returns True if it finds clashes
						# so: 
						if not ns: 
							# If the new chain of the model doesn't add clash:
							yield new_model
														
						else:
							# if the new chain added clash, delete it and return to 
							# the previous model
							new_model.structure[0].__delitem__(moving_chain.get_id())					