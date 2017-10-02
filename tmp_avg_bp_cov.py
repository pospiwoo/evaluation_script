import os, time, operator, logging, copy , sys
from datetime import timedelta
from collections import defaultdict
#from scipy import stats
#import numpy as np

infile_dir = sys.argv[1]
file_list = os.listdir(infile_dir)
global_bp_cov_sum = 0.0
global_bp_cov_div = 0.
global_masked_cnt = 0.0
global_masked_div = 0.0
for i in file_list:
	if (i.startswith('3_avg_') or i.startswith('out_5_')) and i.endswith('.txt'):
		inFile = open(os.path.join(infile_dir,i),'r')
		for line in inFile:
			if line.find('\tref\t') > -1:
				data = line.strip().split('\tref\t')
				global_bp_cov_sum += float(data[1])
				global_bp_cov_div += 1.0
	elif (i.startswith('1_assembled_reads') or i.startswith('1_')) and i.endswith('.fa'):
		inFile1 = open(os.path.join(infile_dir,i),'r')
		for line in inFile1:
			if line.startswith('u') or line.startswith('-') or line.startswith('*'):
				for i in xrange(0,len(line)):
					if line[i] == 'u':
						global_masked_cnt += 1.0
						global_masked_div += 1.0
					else:
						global_masked_div += 1.0

print 'global bp coverage:', global_bp_cov_sum / global_bp_cov_div
print 'masked percentage:', global_masked_cnt / global_masked_div * 100











