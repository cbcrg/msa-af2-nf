#!/usr/bin/env python
import sys
import os

msa={}
msa_file=sys.argv[1]
with open(msa_file) as f:
	for i in f:
		i=i.rstrip()
		if i and i[0]=='>':
			ID=i[1:]
			continue
		msa[ID]=msa.get(ID,'')+i

#for k in sorted(msa.iterkeys()):
for k in sorted(iter(msa.keys())):
	print(">{}\n{}".format(k,msa[k]))

