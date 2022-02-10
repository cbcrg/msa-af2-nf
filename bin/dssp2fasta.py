#!/usr/bin/env python2

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBList
import os
import glob

# a list of PDB IDs to download
pdb_list = glob.glob("*[0-9].pdb")

# parse structure
p = PDBParser()
for i in pdb_list:
	prefix = os.path.splitext(os.path.basename(i))[0]
	struct_name = os.path.splitext(os.path.basename(i))[0].replace("_","/")
	structure = p.get_structure(prefix, i)
	# calculate DSSP
	dssp = DSSP(structure[0],i)
	# extract sequence and secondary structure from the DSSP tuple
	sequence = ''
	sec_structure = ''
	for z in range(len(dssp)):
		a_key = dssp.keys()[z]
		sequence += dssp[a_key][1]
		sec_structure += dssp[a_key][2]

	#
    # The DSSP codes for secondary structure used here are:
    # =====     ====
    # Code      Structure
    # =====     ====
    # H         Alpha helix (4-12)
    # B         Isolated beta-bridge residue
    # E         Strand
    # G         3-10 helix
    # I         Pi helix
    # T         Turn
    # S         Bend
    # -         None
    # =====     ====
    #
    # if desired, convert DSSP's 8-state assignments into 3-state [C - coil, E - extended (beta-strand), H - helix]
	sec_structure = sec_structure.replace('-', 'C')
	sec_structure = sec_structure.replace('I', 'C')
	sec_structure = sec_structure.replace('T', 'C')
	sec_structure = sec_structure.replace('S', 'C')
	sec_structure = sec_structure.replace('G', 'H')
	sec_structure = sec_structure.replace('B', 'E')

	f1=open(prefix+'.dssp', 'w+')
	f1.write('>'+struct_name+"\n"+sec_structure+"\n")
	f1.close()
