# ********************************************************* #
#  This is a Test Inteface for JLAB-CSSA Beam-Beam code.    #
#  ---                                                      #
#  River Huang and Vasiliy Morozov                          #
#  2018-2019                                                #
# ********************************************************* #


import numpy as np
import os
import tkinter as tk
from tkinter import *
from tkinter.scrolledtext import ScrolledText

class Btn_def():

    def load_sigma_beam1(self, Readme, contents):
        try:
            with open(Readme) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, '  Beam1_Sigma_x(m),        Beam1_Sigma_y(m),        Beam1_Sigma_z(m)'+'\n')
                contents.insert(INSERT, '---------------------------------------------------------------------------'+'\n')
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass

    def load_sigma_beam2(self, Readme, contents):
        try:
            with open(Readme) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, '  Beam2_Sigma_x(m),        Beam2_Sigma_y(m),        Beam2_Sigma_z(m)'+'\n')
                contents.insert(INSERT, '---------------------------------------------------------------------------'+'\n')
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass

    def load_luminosity(self, Readme, contents):
        try:
            with open(Readme) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, '\n')
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass

    def exit_Crab(self):
        sys.exit(0)

           
# main windows (parameter window)
class Crab_run(tk.Tk):
    
    def __init__(self):
        super().__init__()
        # self.pack() # if using tk.Frame, keep this line !!!
        self.title("Crab Kick via CASA Beam-Beam GUI")
        self.geometry('880x820')        
        
        # Self interface
        self.setupUI()
        
    def setupUI(self):
    
        contents = ScrolledText()
        contents.pack(side=BOTTOM, expand=True, fill=BOTH)

        Readme = 'Readme_CrabKick_Data_Instruction.casa'
        Readme_sigma_beam1 = 'Files_of_Output/Sigma_beam1.casa'
        Readme_sigma_beam2 = 'Files_of_Output/Sigma_beam2.casa'
        Readme_luminosity  = 'Files_of_Output/Luminosity.casa'

        try:
            with open(Readme) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass
        
        Btn = Btn_def()
        btn1 = Button(self, text='Show RMS-zise of Beam-1', bg="lightblue", command=lambda: Btn.load_sigma_beam1(Readme_sigma_beam1, contents)).pack(side=LEFT)
        btn2 = Button(self, text='Show RMS-zise of Beam-2', bg="lightblue", command=lambda: Btn.load_sigma_beam2(Readme_sigma_beam2, contents)).pack(side=LEFT)
        btn2 = Button(self, text='Show Luminosity', bg="lightblue", command=lambda: Btn.load_luminosity(Readme_luminosity, contents)).pack(side=LEFT)
        btn4 = Button(self, text='Exit', bg="lightgreen", command=lambda: Btn.exit_Crab()).pack(side=RIGHT)
        
        return
        
        
if __name__ == '__main__':
    app = Crab_run()
    app.mainloop()
