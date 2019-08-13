##MODLUES:
from Bio.PDB import * 
from glob import glob
import sys  
import os
import Bio.PDB as pdb 
import random
import string
ppb = pdb.Polypeptide.PPBuilder() 
import shutil 
from Bio import pairwise2
import argparse
import string
import random
from functions import *


parser = argparse.ArgumentParser(description="Given a set of interacting pairs (prot-prot), reconstruct the complete macrocomplex")

parser.add_argument('-i', '--input',
                        dest="inputdir",
                        action="store",
                        default="./",
                        help="Name of the directory that contains the files (in PDB format) with the chain pairs that will be used to obtain the protein complex")

parser.add_argument('-o', '--output',
                        dest="outputdir",
                        action="store",
                        default="./",
                        help="Directory that contains the set of outputs that are obtained (pdb files with the resulting macrocomplexes and other outputs)")

parser.add_argument('-v', '--verbose',
                        dest="verbose",
                        action="store_true",
                        default=None,
                        help="Using this option the user can follow step by step what the program is doing (Print log in stderr)")

parser.add_argument('-d', '--distance',
                        dest="inter_dist",
                        type=int,
                        action="store",
                        default=6,
                        help="Interaction distance of the chains of the input pdb files. Used to filter the input. The default value is 6")


options = parser.parse_args()


folder=options.inputdir
outfolder=options.outputdir
verbose=options.verbose
in_dist=options.inter_dist
path=os.path.dirname(os.path.realpath(__file__))
