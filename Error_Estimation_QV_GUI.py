# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 09:54:59 2016

@author: Sunghee Woo
"""
import os, Tkinter, tkFileDialog #, time,sys
#import est_error_modules as Error_est
#import est_error_from_predefined_mutations as Error_est_barcode
import est_error_QV as Error_est_barcode


class Error_Estimation_GUI(Tkinter.Frame):
    def __init__(self, root):
        self.path = Error_est.Path()
        self.script_dir = self.path.script_dir
        self.parent_dir = self.path.parent_dir
        self.input_dir = os.path.join(self.parent_dir, 'Output')
        self.target_fa_input_dir = os.path.join(self.parent_dir, 'Input')
        self.output_dir = self.path.script_dir
        self.fasta = '' # Result FASTA file
        self.target_fa = '' # target sequence FASTA file: 'target_sequeces.fa'
        self.output_file_name = ''
        Tkinter.Frame.__init__(self, root)

        # define buttons
        button_row_ind = 0
        Tkinter.Label(self, text='Select result FASTA file').grid(row=button_row_ind,column=0)
        Tkinter.Button(self, text='Select result FASTA file', command=self.ask_fasta_filename).grid(row=button_row_ind,column=1)
        
        button_row_ind += 1
        Tkinter.Label(self, text='Select target FASTA file').grid(row=button_row_ind,column=0)
        Tkinter.Button(self, text='Select target FASTA file', command=self.ask_target_seq_file_name).grid(row=button_row_ind,column=1)
        
        button_row_ind += 1
        Tkinter.Label(self, text='Output file name').grid(row=button_row_ind,column=0)
        self.output_file_name = Tkinter.StringVar()
        Tkinter.Entry(self, textvariable=self.output_file_name).grid(row=button_row_ind,column=1)
        
        button_row_ind += 1
        Tkinter.Button(self, text ='Run', command = self.submit).grid(row=button_row_ind,column=0)

    def ask_fasta_filename(self):
        self.fasta = tkFileDialog.askopenfilename(title="Select result FASTA file", initialdir=self.input_dir)
        Tkinter.Label(self, text=os.path.basename(self.fasta)).grid(row=0,column=2)

    def ask_target_seq_file_name(self):
        self.target_fa = tkFileDialog.askopenfilename(title="Select target FASTA file", initialdir=self.target_fa_input_dir)
        Tkinter.Label(self, text=os.path.basename(self.target_fa)).grid(row=1,column=2)

    def submit(self):
        self.output_file_name = os.path.join(self.output_dir, str(self.output_file_name.get())+'.txt')
        error_inst.Process(self.output_file_name)

        
        
        


class Error_Estimation_GUI_barcode(Tkinter.Frame):
    def __init__(self, root):
        self.path = Error_est_barcode.Path()
        self.script_dir = self.path.script_dir
        self.parent_dir = self.path.parent_dir
        self.input_dir = os.path.join(self.parent_dir, 'Output')
        self.target_fa_input_dir = os.path.join(self.parent_dir, 'Input')
        self.output_dir = self.path.script_dir
        self.fasta = '' # Result FASTA file
        self.target_fa = '' # target sequence FASTA file: 'target_sequeces.fa'
        self.output_file_name = ''
        Tkinter.Frame.__init__(self, root)

        # define buttons
        button_row_ind = 0
        Tkinter.Label(self, text='Select result FASTA file').grid(row=button_row_ind,column=0)
        Tkinter.Button(self, text='Select result FASTA file', command=self.ask_fasta_filename).grid(row=button_row_ind,column=1)
        
        button_row_ind += 1
        Tkinter.Label(self, text='Select target FASTA file').grid(row=button_row_ind,column=0)
        Tkinter.Button(self, text='Select target FASTA file', command=self.ask_target_seq_file_name).grid(row=button_row_ind,column=1)
        
        button_row_ind += 1
        Tkinter.Label(self, text='Output file name').grid(row=button_row_ind,column=0)
        self.output_file_name = Tkinter.StringVar()
        Tkinter.Entry(self, textvariable=self.output_file_name).grid(row=button_row_ind,column=1)
        
        button_row_ind += 1
        Tkinter.Button(self, text ='Run', command = self.submit).grid(row=button_row_ind,column=0)

    def ask_fasta_filename(self):
        self.fasta = tkFileDialog.askopenfilename(title="Select result FASTA file", initialdir=self.input_dir)
        Tkinter.Label(self, text=os.path.basename(self.fasta)).grid(row=0,column=2)

    def ask_target_seq_file_name(self):
        self.target_fa = tkFileDialog.askopenfilename(title="Select target FASTA file", initialdir=self.target_fa_input_dir)
        Tkinter.Label(self, text=os.path.basename(self.target_fa)).grid(row=1,column=2)

    def submit(self):
        self.output_file_name = os.path.join(self.output_dir, str(self.output_file_name.get())+'.txt')
        #error_inst = Error_est_barcode.ErrorEstimation_old(self.fasta)
        error_inst = Error_est_barcode.ErrorEstimation(self.fasta)
        error_inst.Process(self.output_file_name)

        
        
        
        

if __name__=='__main__':
    root = Tkinter.Tk()
    #Error_Estimation_GUI(root).pack()
    Error_Estimation_GUI_barcode(root).pack()
    root.wm_title("Hyb&Seq Software Evaluator")
    root.mainloop()
    



########################## Deprecated ##########################
        #Tkinter.Label(self, text='Running').grid(row=3,column=0)
        #Tkinter.Label(self, text='Finished').grid(row=3,column=0)