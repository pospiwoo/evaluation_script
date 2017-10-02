import os, time, operator, logging, copy , sys
from datetime import timedelta
from collections import defaultdict
#from scipy import stats
#import numpy as np

cov_threshold = 3

cnt_all = 0.0
cnt_unmasked = 0.0
sum_SNP = 0.0
infile = sys.argv[1]
infile = open(infile,'r')
for line in infile:
	line = line.strip()
	if line.startswith('>'):
		continue
	elif line.startswith('u') or line.startswith('-') or line.startswith('+'):
		indicator = line.strip()
		continue
	elif line.startswith('N') or line.startswith('T') or line.startswith('C') or line.startswith('G') or line.startswith('A'):
		continue
	else:
		cov = []
		for j in xrange(0,len(line)):
			cnt_all += 1.0
			if int(line[j]) >= cov_threshold:
				cnt_unmasked += 1.0
				if indicator[j] == '*':
					sum_SNP += 1.0
		continue


#>KRASex2;17.944_2;3Spotters_cnt:1;gene_uniq_sixmer:1;mutation_search_cnt:0;2Spotters_cnt:4;3Spotters:GCGTAG;2Spotters:CANNGT,GCTCNN,NNTGGT,GANNGC;Mutation3Spotters:
#NNNNNNNNNNNNNNNNNNNNNNNNGCGTAGNNNNNNNNNNNNNNNNNNNN
#uuuuuuuuuuuuuuuuuuuuuuuu------uuuuuuuuuuuuuuuuuuuu
#00000000000000000000000011111100000000000000000000

print 'cnt_all', cnt_all
print 'sequenced percentage', float(cnt_unmasked) / float(cnt_all) * 100
print 'bp error rate', float(sum_SNP) / float(cnt_unmasked) * 100
