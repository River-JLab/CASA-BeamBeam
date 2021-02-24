# ********************************************************* #
#  This is the Test Inteface for JLAB-CSSA Beam-Beam code.  #
#  ---                                                      #
#  River Huang and Vasiliy Morozov                          #
#  2020                                                     #
# ********************************************************* #


import numpy as np
import os
import tkinter as tk
from tkinter import *
from tkinter.scrolledtext import ScrolledText

class Btn_def():

    def load_luminosity(self, Readme, contents):
        try:
            with open(Readme) as file:
                contents.delete('1.0', END)
                # contents.insert(INSERT, '(4.3321047482379751e+31, 5.6539877697226666e+31, 5.6539877697226666e+31, 'Not Availavle')'+'\n')
                contents.insert(INSERT, '\n')
                # contents.insert(INSERT, '---------------------------------------------------------------------------------------------------'+'\n')
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass

    def exit_luminosity(self):
        sys.exit(0)

           
# main windows (parameter window)
class luminosity_run(tk.Tk):
    
    def __init__(self):
        super().__init__()
        # self.pack() # if using tk.Frame, keep this line !!!
        self.title("Luminosity via CASA Beam-Beam GUI")
        self.geometry('960x820')        
        
        # Self interface
        self.setupUI()
        
    def setupUI(self):
    
        contents = ScrolledText()
        contents.pack(side=BOTTOM, expand=True, fill=BOTH)

        Readme = 'Readme_Luminosity_Data_Instruction.casa'
        Readme_Luminosity = 'File_of_Output/Luminosity.casa'

        try:
            with open(Readme) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass
        
        Btn = Btn_def()
        btn1 = Button(self, text='Show Luminosity', bg="lightblue", command=lambda: Btn.load_luminosity(Readme_Luminosity, contents)).pack(side=LEFT)
        btn2 = Button(self, text='Exit', bg="lightgreen", command=lambda: Btn.exit_luminosity()).pack(side=RIGHT)
        
        return
        
        
if __name__ == '__main__':
    app = luminosity_run()
    app.mainloop()
