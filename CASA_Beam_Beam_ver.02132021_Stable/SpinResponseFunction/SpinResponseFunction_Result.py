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

    def exit_SRF(self):
        sys.exit(0)

           
# main windows (parameter window)
class SRF_run(tk.Tk):
    
    def __init__(self):
        super().__init__()
        # self.pack() # if using tk.Frame, keep this line !!!
        self.title("SRF via CASA Beam-Beam GUI")
        self.geometry('360x820')        
        
        # Self interface
        self.setupUI()
        
    def setupUI(self):
    
        contents = ScrolledText()
        contents.pack(side=BOTTOM, expand=True, fill=BOTH)

        Readme = 'Files_of_Output/SpinFunctionOutput.casa'

        try:
            with open(Readme) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, '   z(m),              |F|'+'\n')
                contents.insert(INSERT, '----------------------------------'+'\n')
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass
        
        Btn = Btn_def()
        btn1 = Button(self, text='Exit', bg="lightblue", command=lambda: Btn.exit_SRF()).pack(side=TOP)
        
        return
        
        
if __name__ == '__main__':
    app = SRF_run()
    app.mainloop()
