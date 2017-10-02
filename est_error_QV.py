# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 22:04:20 2016

@author: pospiwoo
"""
import os, math
from os.path import dirname, abspath

### Path class ###
class Path:
    def __init__(self):
        self.script_dir = dirname(__file__)
        self.parent_dir = dirname(dirname(abspath(__file__)))
        self.input_dir = os.path.join(self.parent_dir, 'Output')
        self.output_dir = self.script_dir

        
### Error Estimation class for artificial SNP ###
class ErrorEstimation:
    def __init__(self, fasta):
        self.fasta = fasta
        self.cnt_all_mu = 0
        self.total_correct_molecule = 0
        self.cnt_total_molecule = 0
        self.total_bp = 0
        self.total_snp = 0
        self.cnt_masked = 0
        
        self.bp_level_error_sum = {}
        self.bp_level_error_cnt = {}
        self.bp_level_error = {}
        self.bp_masked_cnt = {}

        self.curr_target = ''
        self.script_dir = dirname(__file__)
        self.parent_dir = dirname(dirname(abspath(__file__)))
        self.input_dir = os.path.join(self.parent_dir, 'Input')
        self.output_dir = os.path.join(self.parent_dir, 'Output')
        #self.output_dir = os.path.join(self.output_dir, 'guided_SNPs_test')
        

    def countCorrectMolecues(self, read_str):
        check_found = False
        if read_str.startswith('A') or read_str.startswith('T') or \
            read_str.startswith('C') or read_str.startswith('G'):
            self.cnt_total_molecule += 1
            for i in xrange(0, len(self.error_dic[self.curr_target])):
                if read_str == self.error_dic[self.curr_target][i][0]:
                    self.error_dic[self.curr_target][i][2] += 1
                    check_found = True
                    self.total_correct_molecule += 1
            if check_found == False:
                self.error_dic[self.curr_target][-1][2] += 1
        
        
    def countCorrectSNPs(self, read_str):
        if read_str.startswith('A') or read_str.startswith('T') or \
            read_str.startswith('C') or read_str.startswith('G') or \
            read_str.startswith('N'):
            pos_dic_list = self.mu_pos_dic[self.curr_target]
            for i in xrange(0, len(pos_dic_list)):
                if pos_dic_list[i][0] != 'SNP':
                    continue
                snp_pos = pos_dic_list[i][1]
                ref_base = pos_dic_list[i][2]
                mu_base = pos_dic_list[i][3]
#                mu_name = pos_dic_list[i][4]
                if read_str[snp_pos] == 'N':
                    pos_dic_list[i][5] += 1
                    #print snp_pos, read_str[snp_pos], read_str
                elif read_str[snp_pos] == ref_base:
                    pos_dic_list[i][6] += 1
                elif read_str[snp_pos] == mu_base:
                    pos_dic_list[i][7] += 1
                else: # wrong snp called
                    pos_dic_list[i][8] += 1
        

    def addQVErrors(self, line):        
		for i in xrange(0,len(self.bp_level_error_sum[self.curr_target])):
			if len(line) != len(self.bp_level_error_sum[self.curr_target]):
#				print 'wrong:', self.curr_target
#				print ' ', self.bp_level_error_sum[self.curr_target]
#				print ' ', line
				continue
			if line[i] == '-' :
				self.bp_level_error_cnt[self.curr_target][i] += 1.0
				continue
			elif line[i] == 'u' :
				self.bp_masked_cnt[self.curr_target][i] += 1.0
				continue
			elif line[i] == '*' or \
				line[i] == '=' or \
				line[i] == '+':
				self.bp_level_error_sum[self.curr_target][i] += 1.0
				self.bp_level_error_cnt[self.curr_target][i] += 1.0

    def Process(self, out_file_name):
        with open(self.fasta,'r') as assm_file:
            for line in assm_file:
                line = line.strip()
                if line.startswith('>'):
                    self.curr_target = line.split(';')[0].replace('>','')
                    continue
                elif line.startswith('-') or line.startswith('u') or line.startswith('*'):
                    if self.curr_target in self.bp_level_error_sum:
						self.addQVErrors(line)
                    else:
						self.bp_level_error_sum[self.curr_target] = [0.0]*len(line)
						self.bp_level_error_cnt[self.curr_target] = [0.0]*len(line)
						self.bp_level_error[self.curr_target] = [0.0]*len(line)
						self.bp_masked_cnt[self.curr_target] = [0.0]*len(line)
						self.addQVErrors(line)
                    self.cnt_masked += line.count('u')
                    self.cnt_all_mu += line.count('*')
                    self.cnt_all_mu += line.count('=')
                    self.cnt_all_mu += line.count('+')
                    self.total_bp += len(line)
                    continue                    

        oFile = open(out_file_name,'w')
        for g in self.bp_level_error_cnt:
            oFile.write(g+'_phred\t')
#            oFile.write(str(self.bp_level_error_sum[g])+'\t')
#            oFile.write(str(self.bp_level_error_cnt[g])+'\t')
#            print '1', self.bp_level_error_sum[g]
#            print '2', self.bp_level_error_cnt[g]
            for i in xrange(3,len(self.bp_level_error_cnt[g])-3):
                self.bp_level_error[g][i] = self.bp_level_error_sum[g][i] / self.bp_level_error_cnt[g][i]
                phred = 999999
                if self.bp_level_error_cnt[g][i] == 0.0:
                    phred = 0
                elif self.bp_level_error[g][i] == 0.0 and self.bp_level_error_cnt[g][i] > 0.0:
                    phred = 60
                    #print self.bp_level_error_sum[g][i] , self.bp_level_error_cnt[g][i]
                else:
                    phred = -10.0 * math.log(float(self.bp_level_error[g][i]), 10)
                oFile.write(str(phred)+'\t')
            oFile.write('\n')
            oFile.write(g+'_error_rate\t')
            for i in xrange(3,len(self.bp_level_error[g])-3):
                    oFile.write(str(self.bp_level_error[g][i])+'\t')
            oFile.write('\n')
            oFile.write(g+'_wrong_cnt\t')
            for i in xrange(3,len(self.bp_level_error_sum[g])-3):
                    oFile.write(str(self.bp_level_error_sum[g][i])+'\t')
            oFile.write('\n')
            oFile.write(g+'_all_called_bp\t')
            for i in xrange(3,len(self.bp_level_error_cnt[g])-3):
                    oFile.write(str(self.bp_level_error_cnt[g][i])+'\t')
            oFile.write('\n')
            oFile.write(g+'_masked_bp\t')
            for i in xrange(3,len(self.bp_masked_cnt[g])-3):
                    oFile.write(str(self.bp_masked_cnt[g][i])+'\t')
            oFile.write('\n')
            oFile.write('\n')

        bp_e_rate_g = float(self.cnt_all_mu) / (float(self.total_bp) - float(self.cnt_masked))
        if bp_e_rate_g == 0.0:
            phred = 60
        else:
            phred = -10.0 * math.log(float(bp_e_rate_g), 10)
        oFile.write('Global BP error rate'+'\t'+str(bp_e_rate_g)+'\n')
        oFile.write('phred'+'\t'+str(phred)+'\n')
        oFile.write('base called'+'\t'+str(float(self.total_bp-self.cnt_masked)/self.total_bp*100)+'\n')
        print 'Global BP error rate', bp_e_rate_g
        print 'phred', phred
        print 'base called:', float(self.total_bp-self.cnt_masked)/self.total_bp*100, '%'
        oFile.close()


if __name__=='__main__':
    path = Path()
    fa_file_name = os.path.join(path.input_dir, 'out_1_assm_reads.fa')
    Error_inst = ErrorEstimation(fa_file_name)
    Error_inst.Process()
    