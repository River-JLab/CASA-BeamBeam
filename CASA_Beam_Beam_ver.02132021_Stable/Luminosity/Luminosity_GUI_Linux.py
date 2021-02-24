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


class Btn_def():
    """ Button Function """
    # def save(self, filename, contents):
        # """ Save File """
        # try:
            # with open(filename, 'w') as file:
                # file.write(contents.get('1.0', END))
        # except FileNotFoundError:
            # pass

    def load(self, filename, contents):
        """ Open File """
        try:
            with open(filename) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass

    def Luminosity_CAL(self, parameter, filename, contents):
        try:
            with open(parameter, 'w') as file01:
                file01.write(contents.get('1.0', END))
        except FileNotFoundError:
            pass

        os.system("python3 Luminosity_Asking_sigma_GUI.py")
        os.system("python3 Luminosity_Calculation_Linux.py")
        os.system("python3 Luminosity_Result_GUI.py")
        os._exit(0)
           
# main windows (parameter window)
class LuminosityParameterEditor(tk.Tk):
    
    def __init__(self):
        super().__init__()
        # self.pack() # if using tk.Frame, keep this line !!!
        self.title("Please edit and check the parameters, then Click 'Calculate Luminosity'")
        self.geometry('760x450')        
        
        # Self interface
        self.setupUI()
        
    def setupUI(self):

        contents = ScrolledText()
        contents.pack(side=BOTTOM, expand=True, fill=BOTH)

        input1 = 'File_of_Input/Parameters_For_Luminosity.casa'
        output = 'File_of_Output/Luminosity.casa'
        
        try:
            with open(input1) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass
        
        Btn = Btn_def()
        btn1 = Button(self, text='Re-Load Parameters for Calculating Luminosity', bg="lightblue", command=lambda: Btn.load(input1, contents)).pack(side=LEFT)
        
        btn3 = Button(self, text="Calculate Luminosity", bg="lightgreen", command=lambda: Btn.Luminosity_CAL(input1, output, contents)).pack(side=RIGHT)
        
        return
        
        
if __name__ == '__main__':
    app = LuminosityParameterEditor()
    app.mainloop()