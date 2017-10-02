# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 22:04:20 2016

@author: pospiwoo
"""
import os
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
    def __init__(self, fasta, seqlog_txt):
        self.fasta = fasta
        self.seqlog_txt = seqlog_txt
        self.total_correct_molecule = 0
        self.cnt_total_molecule = 0
        self.total_bp = 0
        self.total_snp = 0
        self.cnt_total_bp_g = 0
        self.cnt_masked_g = 0
        self.cnt_unmasked_g = 0
        self.cnt_correct_g = 0
        self.cnt_wrong_g = 0
        self.cnt_wrong_MTM = 0
        self.cnt_all_MTM = 0
        self.parse_from_file = True
        self.curr_grid_pos = -1
        self.seqlog_dic = {}
        self.curr_seq = ''
        self.curr_read_str = ''	
        self.curr_call_type = ''
        self.match_gene_dic = {}
        self.wrong_grid_pos = 3
		
        self.turn_on_indels = False # False True
        self.cov_thres = 3
        self.cov_thres = 2
        self.cov_thres = -1

        self.curr_target = ''
        self.script_dir = dirname(__file__)
        self.parent_dir = dirname(dirname(abspath(__file__)))
        self.input_dir = os.path.join(self.parent_dir, 'Input')
        self.output_dir = os.path.join(self.parent_dir, 'Output')
        #self.output_dir = os.path.join(self.output_dir, 'guided_SNPs_test')
        
        self.error_dic = {}
        self.mu_pos_dic = {}
        if self.parse_from_file:
			self.parseMutations()
        else:
			self.hardcodeMutations()
       
        
    def parseMutations(self):
        mutation_dir = os.path.join(self.parent_dir, 'Input', 'mutation')
        mu_file_name = os.path.join(mutation_dir, 'predefined_mutations.txt')
        mu_file_name = os.path.join(mutation_dir, 'predefined_mutations_random_sim_2.txt')
        mu_file_name = os.path.join(mutation_dir, 'predefined_mutations_random_sim_1.txt')
        mu_file_name = os.path.join(mutation_dir, 'predefined_mutations_random_sim.txt')
        mu_file_name = os.path.join(mutation_dir, 'predefined_mutations_random_sim_4.txt')
        with open(mu_file_name, 'r') as inFile:
            for line in inFile:
                if line.startswith('#') or line.strip() == '':
                    continue
                #PICK3CA:A:G:chr?:26:SNP
                data = line.strip().split(':')
                gene = data[0]
                ref = data[1]
                mu = data[2]
                chr_str = data[3]
                pos = data[4]
                mu_type = data[5]
                tmp_list = [mu_type, int(pos), ref, mu, line.strip(), 0, 0, 0, 0]
                if gene in self.mu_pos_dic:
                    self.mu_pos_dic[gene].append(tmp_list)
                else:
                    self.mu_pos_dic[gene] = []
                    self.mu_pos_dic[gene].append(tmp_list)


    def countCorrect_IN_DEL(self, read_str, read_str_0):
        if read_str.startswith('-') or read_str.startswith('u') or \
            read_str.startswith('*') or read_str.startswith('=') or \
            read_str.startswith('+') or read_str.startswith('^'):
            if self.curr_target not in self.mu_pos_dic:
				return
            number_list = []
            for k in xrange(0, len(read_str_0)):
                number_list.append(int(read_str_0[k]))
            pos_dic_list = self.mu_pos_dic[self.curr_target]

            for i in xrange(0, len(pos_dic_list)):
                if pos_dic_list[i][0] == 'SNP':
                    continue
                elif pos_dic_list[i][0] == 'DEL':
                    snp_pos = pos_dic_list[i][1]
#                    mu_name = pos_dic_list[i][4]
                    if number_list[snp_pos] <= self.cov_thres or read_str[snp_pos] == 'u':
                        pos_dic_list[i][5] += 1
                        #print snp_pos, read_str[snp_pos], read_str
                    elif read_str[snp_pos] == '-':
                        pos_dic_list[i][8] += 1
                    elif read_str[snp_pos] == '=':
                        pos_dic_list[i][6] += 1
#                    else: # wrong snp called
#                        pos_dic_list[i][8] += 1
                elif pos_dic_list[i][0] == 'INS':
                    snp_pos = pos_dic_list[i][1]
#                    mu_name = pos_dic_list[i][4]
                    if number_list[snp_pos] <= self.cov_thres or read_str[snp_pos] == 'u':
                        pos_dic_list[i][5] += 1
                        #print snp_pos, read_str[snp_pos], read_str
                    elif read_str[snp_pos] == '-':
                        pos_dic_list[i][8] += 1
                    elif read_str[snp_pos] == '+':
                        pos_dic_list[i][6] += 1
#                    else: # wrong snp called
#                        pos_dic_list[i][8] += 1


    def calcPercentage(self, n_masked, ref, correct_mu, wrong_mu):
        correct = float(ref)
        divide = float(ref)+float(correct_mu)+float(wrong_mu)
        divide += 0.000001
        float_percentage = correct / divide * 100
        return float_percentage
                
        
    def calcPercentage_old(self, n_masked, ref, correct_mu, wrong_mu):
        correct = float(correct_mu)
        divide = float(ref)+float(correct_mu)+float(wrong_mu)
        divide += 0.000001
        float_percentage = correct / divide * 100
        return float_percentage
        
        
    def countCorrectSNPs(self, read_str):
        if self.curr_grid_pos in self.seqlog_dic:
            self.curr_seq = self.seqlog_dic[self.curr_grid_pos]
        else:
            self.curr_seq = ''
            return
        if read_str.startswith('A') or read_str.startswith('T') or \
            read_str.startswith('C') or read_str.startswith('G') or \
            read_str.startswith('N'):
            if self.curr_target not in self.mu_pos_dic:
				return
            pos_dic_list = self.mu_pos_dic[self.curr_target]
            for i in xrange(0, len(pos_dic_list)):
                if pos_dic_list[i][0] == 'DEL':
                    continue
                elif pos_dic_list[i][0] == 'SNP':
                    snp_pos = pos_dic_list[i][1]
                    ref_base = pos_dic_list[i][2]
                    mu_base = pos_dic_list[i][3]
#                    mu_name = pos_dic_list[i][4]
                    if read_str[snp_pos] == 'N':
                        pos_dic_list[i][5] += 1
                        #print snp_pos, read_str[snp_pos], read_str
                    elif read_str[snp_pos] == ref_base:
                        pos_dic_list[i][6] += 1
                    elif read_str[snp_pos] == mu_base:
                        pos_dic_list[i][7] += 1
                    else: # wrong snp called
                        pos_dic_list[i][8] += 1
        
        
    def countCorrectSNPsSIM(self, read_str):
        if self.curr_grid_pos in self.seqlog_dic:
            self.curr_seq = self.seqlog_dic[self.curr_grid_pos]
        else:
            self.curr_seq = ''
            return
        if read_str.startswith('A') or read_str.startswith('T') or \
            read_str.startswith('C') or read_str.startswith('G') or \
            read_str.startswith('N'):
            if self.curr_target not in self.mu_pos_dic:
				return
            
            if len(read_str) != len(self.curr_seq):
                print self.curr_grid_pos
                print read_str
                print self.curr_seq
                print ''
                return
            
            pos_dic_list = self.mu_pos_dic[self.curr_target]
            for i in xrange(0, len(pos_dic_list)):
                if pos_dic_list[i][0] == 'SNP':
                    snp_pos = pos_dic_list[i][1]
                    ref_base = pos_dic_list[i][2]
                    mu_base = pos_dic_list[i][3]
                    if read_str[snp_pos] == 'N':
                        pos_dic_list[i][5] += 1
                        #print snp_pos, read_str[snp_pos], read_str
                    elif read_str[snp_pos] == self.curr_seq[snp_pos]:
                        pos_dic_list[i][6] += 1
                    else: # wrong
                        pos_dic_list[i][8] += 1
        
        
    def countCorrectSNPsSIM_cov(self, read_str):
        if self.curr_grid_pos in self.seqlog_dic:
            self.curr_seq = self.seqlog_dic[self.curr_grid_pos]
        else:
            self.curr_seq = ''
            return
        if read_str.startswith('A') or read_str.startswith('T') or \
            read_str.startswith('C') or read_str.startswith('G') or \
            read_str.startswith('N'):
            self.curr_read_str = read_str
        else:
            if self.curr_target not in self.mu_pos_dic:
				return
            if len(read_str) != len(self.curr_seq):
#                print self.curr_grid_pos
#                print read_str
#                print self.curr_seq
#                print ''
                return
            number_list = []
            for k in xrange(0, len(read_str)):
                number_list.append(int(read_str[k]))
            pos_dic_list = self.mu_pos_dic[self.curr_target]
            for i in xrange(0, len(pos_dic_list)):
                if pos_dic_list[i][0] == 'SNP':
                    snp_pos = pos_dic_list[i][1]
                    ref_base = pos_dic_list[i][2]
                    mu_base = pos_dic_list[i][3]
                    if number_list[snp_pos] <= self.cov_thres or self.curr_read_str[snp_pos] == 'N':
                        pos_dic_list[i][5] += 1
                        #print number_list[i], snp_pos, self.curr_read_str[snp_pos], self.curr_read_str
                        #print snp_pos, self.curr_read_str[snp_pos], self.curr_read_str
                    elif self.curr_read_str[snp_pos] == self.curr_seq[snp_pos]:
                        pos_dic_list[i][6] += 1
                    else: # wrong
                        pos_dic_list[i][8] += 1
        
        
    def error_est_QV(self, read_str):
        if read_str.startswith('A') or read_str.startswith('T') or \
            read_str.startswith('C') or read_str.startswith('G') or \
            read_str.startswith('N'):
            if self.curr_target not in self.mu_pos_dic:
				return            
#            self.cnt_all_MTM += 1
        
            if len(read_str) != len(self.curr_seq):
#                self.cnt_wrong_MTM += 1
#                print self.curr_grid_pos
#                print read_str
#                print self.curr_seq
#                print ''
                return
            
            for pos in xrange(0, len(read_str)):
                    self.cnt_total_bp_g += 1
                    if read_str[pos] == 'N':
                        self.cnt_masked_g += 1
                    elif read_str[pos] == self.curr_seq[pos]:
                        self.cnt_correct_g += 1
                        self.cnt_unmasked_g += 1
                    else: # wrong
                        self.cnt_wrong_g += 1
                        self.cnt_unmasked_g += 1
        
        
    def error_est_QV_cov(self, read_str):
        if read_str.startswith('A') or read_str.startswith('T') or \
            read_str.startswith('C') or read_str.startswith('G') or \
            read_str.startswith('N'):
            self.curr_read_str = read_str
        else:
            if self.curr_target not in self.mu_pos_dic:
				return

            if len(read_str) != len(self.curr_seq):
                print self.curr_grid_pos
                print read_str
                print self.curr_seq
                print ''
                return

            number_list = []
            for k in xrange(0, len(read_str)):
                number_list.append(int(read_str[k]))
            for pos in xrange(0, len(read_str)):
                    self.cnt_total_bp_g += 1
                    if number_list[pos] <= self.cov_thres or read_str[pos] == 'N':
                        self.cnt_masked_g += 1
                    elif read_str[pos] == self.curr_seq[pos]:
                        self.cnt_correct_g += 1
                        self.cnt_unmasked_g += 1
                    else: # wrong
                        self.cnt_wrong_g += 1
                        self.cnt_unmasked_g += 1


    def Process(self):
        with open(self.fasta,'r') as assm_file, open(self.seqlog_txt,'r') as seqlog_file:
            for line in seqlog_file:
                line = line.strip().split()
                pos = line[0].split('_')[-1]
                seq = line[-1].strip()
                #@Seq_GridPos_15,469.17  ATGACAAAGAAAGCTATATAAGATATTATTTTATTTTACAGAGTAACAGACTAGCTAGAGACAATGAATTAAGGGAAAATGACAAAGAACAGCTCAAAGCAATTTCTACACGAGATCCTCTCTCTGAAATCACTGAGCAGGAGAAAGATT

                self.seqlog_dic[pos] = seq
            
            for line in assm_file:
                line = line.strip()
                if line.startswith('>'):
                    #print line.split(';')
                    self.curr_target = line.split(';')[0].replace('>','')                    
                    if self.wrong_grid_pos == -1:                    
                        self.curr_grid_pos = ''
                        self.wrong_grid_pos = line.split(';')[1]
                    else:
                        self.curr_grid_pos = self.wrong_grid_pos
                        self.wrong_grid_pos = line.split(';')[1]

                    continue
                elif line.startswith('-') or line.startswith('u') or line.startswith('*'):
                    #self.cnt_masked += line.count('u')
                    if self.turn_on_indels:
                        next_line = next(assm_file).strip()
                        self.countCorrect_IN_DEL(line, next_line)
                    continue
                else:
                    if self.curr_grid_pos == '':
                        continue
                    if not self.turn_on_indels:
                        #self.countCorrectSNPsSIM(line)
                        #self.error_est_QV(line)
                        self.countCorrectSNPsSIM_cov(line)
                        self.error_est_QV_cov(line)
                    

        
        sum_error_list = [0.0, 0.0, 0.0, 0.0]
        for i in self.mu_pos_dic:
            for j in xrange(0,len(self.mu_pos_dic[i])):
                if self.mu_pos_dic[i][j][5] + self.mu_pos_dic[i][j][6] +self.mu_pos_dic[i][j][7] +self.mu_pos_dic[i][j][8] == 0.0:
					continue
                sum_error_list[0] += self.mu_pos_dic[i][j][5]
                sum_error_list[1] += self.mu_pos_dic[i][j][6]
                sum_error_list[2] += self.mu_pos_dic[i][j][7]
                sum_error_list[3] += self.mu_pos_dic[i][j][8]
                corr_percent = self.calcPercentage(self.mu_pos_dic[i][j][5], self.mu_pos_dic[i][j][6], self.mu_pos_dic[i][j][7], self.mu_pos_dic[i][j][8])
                #print i, j, self.mu_pos_dic[i][j], "%.2f"%(corr_percent), '%'
                print     j, self.mu_pos_dic[i][j], "%.2f"%(corr_percent), '%'

        if not self.turn_on_indels:
            corr_percent = self.calcPercentage(sum_error_list[0], sum_error_list[1], sum_error_list[2], sum_error_list[3])
            print "%.2f"%(corr_percent), '%'

            print ""
            print 'masked_percentage', float(self.cnt_masked_g) / float(self.cnt_total_bp_g) * 100
            print 'base_pair_error_percentage', float(self.cnt_wrong_g) / float(self.cnt_unmasked_g) * 100

#            print ""
#            print 'wrong_MTM_percentage', float(self.cnt_wrong_MTM) / float(self.cnt_all_MTM) * 100
			





if __name__=='__main__':
    path = Path()
    fa_file_name = os.path.join(path.input_dir, 'out_1_assm_reads.fa')
    Error_inst = ErrorEstimation(fa_file_name)
    Error_inst.Process()
    