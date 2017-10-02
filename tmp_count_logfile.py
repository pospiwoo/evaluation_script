import os, time, operator, logging, copy , sys
from collections import defaultdict

infile = open(sys.argv[1],'r')
d_hash = {}
for line in infile:
	data = line.strip().split('\t')
	d_hash[data[1]] = 1

for i in d_hash:
	print i
print len(d_hash)

