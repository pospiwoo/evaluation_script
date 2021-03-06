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
    def __init__(self, fasta):
        self.fasta = fasta
        self.total_correct_molecule = 0
        self.cnt_total_molecule = 0
        self.total_bp = 0
        self.total_snp = 0
        self.cnt_masked = 0
        self.parse_from_file = True

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
        mu_file_name = os.path.join(mutation_dir, 'predefined_mutations_run083_artificial_100_test.txt')
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
        
#    def countCorrectMolecues(self, read_str):
#        check_found = False
#        if read_str.startswith('A') or read_str.startswith('T') or \
#            read_str.startswith('C') or read_str.startswith('G'):
#            self.cnt_total_molecule += 1
#            self.total_bp += len(read_str)
#            for i in xrange(0, len(self.error_dic[self.curr_target])):
#                if read_str == self.error_dic[self.curr_target][i][0]:
#                    self.error_dic[self.curr_target][i][2] += 1
#                    check_found = True
#                    self.total_correct_molecule += 1
#            if check_found == False:
#                self.error_dic[self.curr_target][-1][2] += 1
#        if read_str.startswith('-') or read_str.startswith('u') or \
#            read_str.startswith('*') or read_str.startswith('+') or \
#            read_str.startswith('=') or read_str.startswith('^'):
        
        
    def countCorrectSNPs(self, read_str):
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
        
        
    def countCorrect_IN_DEL(self, read_str):
        if read_str.startswith('-') or read_str.startswith('u') or \
            read_str.startswith('*') or read_str.startswith('=') or \
            read_str.startswith('+') or read_str.startswith('^'):
            if self.curr_target not in self.mu_pos_dic:
				return
            pos_dic_list = self.mu_pos_dic[self.curr_target]
            for i in xrange(0, len(pos_dic_list)):
                if pos_dic_list[i][0] == 'SNP':
                    continue
                elif pos_dic_list[i][0] == 'DEL':
                    snp_pos = pos_dic_list[i][1]
                    ref_base = pos_dic_list[i][2]
                    mu_base = pos_dic_list[i][3]
#                    mu_name = pos_dic_list[i][4]
                    if read_str[snp_pos] == 'u':
                        pos_dic_list[i][5] += 1
                        #print snp_pos, read_str[snp_pos], read_str
                    elif read_str[snp_pos] == '-':
                        pos_dic_list[i][6] += 1
                    elif read_str[snp_pos] == '=':
                        pos_dic_list[i][7] += 1
#                    else: # wrong snp called
#                        pos_dic_list[i][8] += 1
                elif pos_dic_list[i][0] == 'INS':
                    snp_pos = pos_dic_list[i][1]
                    ref_base = pos_dic_list[i][2]
                    mu_base = pos_dic_list[i][3]
#                    mu_name = pos_dic_list[i][4]
                    if read_str[snp_pos] == 'u':
                        pos_dic_list[i][5] += 1
                        #print snp_pos, read_str[snp_pos], read_str
                    elif read_str[snp_pos] == '-':
                        pos_dic_list[i][6] += 1
                    elif read_str[snp_pos] == '+':
                        pos_dic_list[i][7] += 1
#                    else: # wrong snp called
#                        pos_dic_list[i][8] += 1
                
        
    def calcPercentage(self, n_masked, ref, correct_mu, wrong_mu):
        correct = float(correct_mu)
        divide = float(ref)+float(correct_mu)+float(wrong_mu)
        divide += 0.000001
        float_percentage = correct / divide * 100
        return float_percentage
        
    def Process(self):
        with open(self.fasta,'r') as assm_file:
            for line in assm_file:
                line = line.strip()
                if line.startswith('>'):
                    self.curr_target = line.split(';')[0].replace('>','')
                    continue
                elif line.startswith('-') or line.startswith('u'):
                    self.cnt_masked += line.count('u')
                    self.countCorrect_IN_DEL(line)
                    continue
                else:
                    #self.countCorrectMolecues(line)
                    self.countCorrectSNPs(line)
                    
        #adjusted = self.cnt_bp_error - self.total_correct_molecule
        #self.total_bp = self.cnt_total_molecule * self.len_mole
#        print 'total molecules:', self.cnt_total_molecule
#        print 'total_correct_molecule'
#        print 'molecule error rate:', float(self.cnt_total_molecule - self.total_correct_molecule) / self.cnt_total_molecule * 100 , '%'
#        print 'total bp:', self.total_bp
#        print 'base called:', float(self.total_bp-self.cnt_masked)/self.total_bp*100, '%'
        
#        for i in self.error_dic:
#            for j in xrange(0,len(self.error_dic[i])):
#                print i, j, self.error_dic[i][j]
#        print ''

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

        corr_percent = self.calcPercentage(sum_error_list[0], sum_error_list[1], sum_error_list[2], sum_error_list[3])
        print "%.2f"%(corr_percent), '%'


    def hardcodeMutations(self):
			self.error_dic['BRAFex15'] = []
			self.error_dic['BRAFex15'].append(['AATAGGTGATTTTGGTCTAGCTACAGTGAAATCTCGATGGAGTGGGTCCC', 'ref', 0] )
			self.error_dic['BRAFex15'].append(['AATAGGTGATTTTGGTCTAGCTACAGAGAAATCTCGATGGAGTGGGTCCC', 'COSM476', 0] )
			self.error_dic['BRAFex15'].append(['AATAGGTGATTTTGGTCTAGCTACAAAGAAATCTCGATGGAGTGGGTCCC', 'COSM473', 0] )
			self.error_dic['BRAFex15'].append(['N', 'N', 0] )
			
			self.error_dic['EGFRex18'] = []
			self.error_dic['EGFRex18'].append(['CTGAATTCAAAAAGATCAAAGTGCTGGGCTCCGGTGCGTTCGGCACGGTG', 'ref', 0] )
			self.error_dic['EGFRex18'].append(['CTGAATTCAAAAAGATCAAAGTGCTGAGCTCCGGTGCGTTCGGCACGGTG', 'COSM6252', 0] )
			self.error_dic['EGFRex18'].append(['N', 'N', 0] )
			
			self.error_dic['EGFRex20'] = []
			self.error_dic['EGFRex20'].append(['CCTCACCTCCACCGTGCAGCTCATCACGCAGCTCATGCCCTTCGGCTGCC', 'ref', 0] )
			self.error_dic['EGFRex20'].append(['CCTCACCTCCACCGTGCAGCTCATCATGCAGCTCATGCCCTTCGGCTGCC', 'COSM6240', 0] )
			self.error_dic['EGFRex20'].append(['N', 'N', 0] )
			
			self.error_dic['FLT3'] = []
			self.error_dic['FLT3'].append(['TGTGACTTTGGATTGGCTCGAGATATCATGAGTGATTCCAACTATGTTGT', 'ref', 0] )
			self.error_dic['FLT3'].append(['TGTGACTTTGGATTGGCTCGAGATATGAGTGATTCCAACTATGTTGT', 'COSM797', 0] )
			self.error_dic['FLT3'].append(['N', 'N', 0] )
			
			self.error_dic['IDH1'] = []
			self.error_dic['IDH1'].append(['ATGGGTAAAACCTATCATCATAGGTCGTCATGCTTATGGGGATCAAGTAA', 'ref', 0] )
			self.error_dic['IDH1'].append(['ATGGGTAAAACCTATCATCATAGGTTGTCATGCTTATGGGGATCAAGTAA', 'COSM28747', 0] )
			self.error_dic['IDH1'].append(['N', 'N', 0] )
			
			self.error_dic['JAK2'] = []
			self.error_dic['JAK2'].append(['ATTTGGTTTTAAATTATGGAGTATGTGTCTGTGGAGACGAGAGTAAGTAA', 'ref', 0] )
			self.error_dic['JAK2'].append(['ATTTGGTTTTAAATTATGGAGTATGTTTCTGTGGAGACGAGAGTAAGTAA', 'COSM12600', 0] )
			self.error_dic['JAK2'].append(['N', 'N', 0] )
			
			self.error_dic['KRASex2'] = []
			self.error_dic['KRASex2'].append(['AACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATA', 'ref', 0] )
			self.error_dic['KRASex2'].append(['AACTTGTGGTAGTTGGAGCTCGTGGCGTAGGCAAGAGTGCCTTGACGATA', 'COSM518', 0] )
			self.error_dic['KRASex2'].append(['AACTTGTGGTAGTTGGAGCTGCTGGCGTAGGCAAGAGTGCCTTGACGATA', 'COSM522', 0] )
			self.error_dic['KRASex2'].append(['AACTTGTGGTAGTTGGAGCTGGTGACGTAGGCAAGAGTGCCTTGACGATA', 'COSM532', 0] )
			self.error_dic['KRASex2'].append(['N', 'N', 0] )
			
			self.error_dic['MEK1'] = []
			self.error_dic['MEK1'].append(['GCAGGTTCTGCATGAGTGCAACTCTCCGTACATCGTGGGCTTCTATGGTG', 'ref', 0] )
			self.error_dic['MEK1'].append(['GCAGGTTCTGCATGAGTGCAACTCTCTGTACATCGTGGGCTTCTATGGTG', 'COSM1315861', 0] )
			self.error_dic['MEK1'].append(['N', 'N', 0] )
			
			self.error_dic['NOTCH1'] = []
			self.error_dic['NOTCH1'].append(['TTCCTGCGGGAGCTCAGCCGCGTGCTGCACACCAACGTGGTCTTCAAGCG', 'ref', 0] )
			self.error_dic['NOTCH1'].append(['TTCCTGCGGGAGCTCAGCCGCGTGCCGCACACCAACGTGGTCTTCAAGCG', 'COSM1155771', 0] )
			self.error_dic['NOTCH1'].append(['N', 'N', 0] )
			
			self.error_dic['NRAS'] = []
			self.error_dic['NRAS'].append(['GTTGGACATACTGGATACAGCTGGACAAGAAGAGTACAGTGCCATGAGAG', 'ref', 0] )
			self.error_dic['NRAS'].append(['GTTGGACATACTGGATACAGCTGGAAAAGAAGAGTACAGTGCCATGAGAG', 'COSM580', 0] )
			self.error_dic['NRAS'].append(['N', 'N', 0] )
			
			self.error_dic['PICK3CA'] = []
			self.error_dic['PICK3CA'].append(['TTTCATGAAACAAATGAATGATGCACATCATGGTGGCTGGACAACAAAAA', 'ref', 0] )
			self.error_dic['PICK3CA'].append(['TTTCATGAAACAAATGAATGATGCACGTCATGGTGGCTGGACAACAAAAA', 'COSM775', 0] )
			self.error_dic['PICK3CA'].append(['N', 'N', 0] )
			
			
			self.mu_pos_dic['BRAFex15'] = []
			self.mu_pos_dic['BRAFex15'].append(['SNP', 26, 'T', 'A', 'COSM476', 0, 0, 0, 0] ) #BRAFex15:T:A:chr?:26:SNP
			self.mu_pos_dic['BRAFex15'].append(['SUBS', 25, 'GT', 'AA', 'COSM473', 0, 0, 0, 0] ) #BRAFex15:GT:AA:chr?:25:SNP

			self.mu_pos_dic['EGFRex18'] = []
			self.mu_pos_dic['EGFRex18'].append(['SNP', 26, 'G', 'A', 'COSM6252', 0, 0, 0, 0] ) #EGFRex18:G:A:chr?:26:SNP
			
			self.mu_pos_dic['EGFRex20'] = []
			self.mu_pos_dic['EGFRex20'].append(['SNP', 26, 'C', 'T', 'COSM6240', 0, 0, 0, 0] ) #EGFRex20:C:T:chr?:26:SNP

			self.mu_pos_dic['FLT3'] = []
			self.mu_pos_dic['FLT3'].append(['DEL', 26, 'CAT', '', 'COSM797', 0, 0, 0, 0] ) #FLT3:CAT::chr?:26:DEL
			
			self.mu_pos_dic['IDH1'] = []
			self.mu_pos_dic['IDH1'].append(['SNP', 25, 'C', 'T', 'COSM28747', 0, 0, 0, 0] ) #IDH1:C:T:chr?:25:SNP

			self.mu_pos_dic['JAK2'] = []
			self.mu_pos_dic['JAK2'].append(['SNP', 26, 'G', 'T', 'COSM12600', 0, 0, 0, 0] ) #JAK2:G:T:chr?:26:SNP

			self.mu_pos_dic['KRASex2'] = []
			self.mu_pos_dic['KRASex2'].append(['SNP', 20, 'G', 'C', 'COSM518', 0, 0, 0, 0] ) #KRASex2:G:C:chr?:20:SNP
			self.mu_pos_dic['KRASex2'].append(['SNP', 21, 'G', 'C', 'COSM522', 0, 0, 0, 0] ) #KRASex2:G:C:chr?:21:SNP
			self.mu_pos_dic['KRASex2'].append(['SNP', 24, 'G', 'A', 'COSM532', 0, 0, 0, 0] ) #KRASex2:G:A:chr?:24:SNP

			self.mu_pos_dic['MEK1'] = []
			self.mu_pos_dic['MEK1'].append(['SNP', 26, 'C', 'T', 'COSM1315861', 0, 0, 0, 0] ) #MEK1:C:T:chr?:26:SNP

			self.mu_pos_dic['NOTCH1'] = []
			self.mu_pos_dic['NOTCH1'].append(['SNP', 25, 'T', 'C', 'COSM1155771', 0, 0, 0, 0] ) #NOTCH1:T:C:chr?:25:SNP
			
			self.mu_pos_dic['NRAS'] = []
			self.mu_pos_dic['NRAS'].append(['SNP', 25, 'C', 'A', 'COSM580', 0, 0, 0, 0] ) #NRAS:C:A:chr?:25:SNP

			self.mu_pos_dic['PICK3CA'] = []
			self.mu_pos_dic['PICK3CA'].append(['SNP', 26, 'A', 'G', 'COSM775', 0, 0, 0, 0] ) #PICK3CA:A:G:chr?:26:SNP


if __name__=='__main__':
    path = Path()
    fa_file_name = os.path.join(path.input_dir, 'out_1_assm_reads.fa')
    Error_inst = ErrorEstimation(fa_file_name)
    Error_inst.Process()
    