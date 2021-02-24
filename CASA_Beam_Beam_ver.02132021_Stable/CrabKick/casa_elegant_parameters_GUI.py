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
class ParameterEditor(tk.Tk):
    
    def __init__(self):
        super().__init__()
        # self.pack() # if using tk.Frame, keep this line !!!
        self.title("Check and edit the parameters for Elegant")
        self.geometry('760x460')        
        
        # Self interface
        self.setupUI()
        
    def setupUI(self):

        contents = ScrolledText()
        contents.pack(side=BOTTOM, expand=True, fill=BOTH)

        filename = 'Elegant_Process/beam_parameters/Beam_parameters.casa'
        
        try:
            with open(filename) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass
        
        Btn = Btn_def()
        btn1 = Button(self, text='Re-load Parameters', bg="orange", command=lambda: Btn.load(filename, contents)).pack(side=LEFT)
        btn2 = Button(self, text='Save and continue', bg="yellow", command=lambda: Btn.save_and_continue(filename, contents)).pack(side=RIGHT)
        btn3 = Button(self, text='Skip and continue', bg="lightgreen", command=lambda: Btn.skip_edit()).pack(side=RIGHT)
        
        return
        
        
if __name__ == '__main__':
    app = ParameterEditor()
    app.mainloop()

