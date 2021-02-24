# ********************************************************* #
#  This is a Test Inteface for JLAB-CSSA Beam-Beam code.    #
#  ---                                                      #
#  River Huang and Vasiliy Morozov                          #
#  2018-2019                                                #
# ********************************************************* #


import os
import tkinter as tk
from tkinter import *
from tkinter.scrolledtext import ScrolledText
import sys, time, math, cmath
import os.path
import numpy as np
from scipy import special
from string import *
from casa_elegant_functions import *


class Btn_def():
    """ Button Function """
    def load(self, filename, contents):
        """ Open File """
        try:
            with open(filename) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass

    def skip_edit(self):
        os._exit(0)

    def save_and_continue(self, filename, contents):
        try:
            with open(filename, 'w') as file:
                file.write(contents.get('1.0', END))
        except FileNotFoundError:
            pass
        os._exit(0)


# main windows (parameter window)
class ele_Editor(tk.Tk):
    
    def __init__(self):
        super().__init__()
        # self.pack() # if using tk.Frame, keep this line !!!
        self.title('Check and Edit Beam-2:  Elegant-continue-command')
        self.geometry('660x600')        
        
        # Self interface
        self.setupUI()
        
    def setupUI(self):

        contents = ScrolledText()
        contents.pack(side=BOTTOM, expand=True, fill=BOTH)
        
        current_path = os.getcwd()
        os.chdir("Elegant_Process")
        parameter = input_parameters()
        os.chdir(current_path)
        
        filename = 'Elegant_Process/ele_command/' + parameter[1][4] + '.ele'
        generate_beam_continue_ele(filename, parameter[1])
        
        try:
            with open(filename) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass
        
        Btn = Btn_def()
        btn1_name = 'Re-load ele-file:    ' + parameter[1][4] + '.ele'
        # btn1 = Button(self, text='Re-load ELE-command', bg="orange", command=lambda: Btn.load(filename, contents)).pack(side=LEFT)
        btn1 = Button(self, text=btn1_name, bg="orange", command=lambda: Btn.load(filename, contents)).pack(side=LEFT)
        btn2 = Button(self, text='Save and continue', bg="yellow", command=lambda: Btn.save_and_continue(filename, contents)).pack(side=RIGHT)
        btn3 = Button(self, text='Skip and continue', bg="lightgreen", command=lambda: Btn.skip_edit()).pack(side=RIGHT)
        
        return
        
        
if __name__ == '__main__':
    app = ele_Editor()
    app.mainloop()

