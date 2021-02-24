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

    def NO_exit(self):
        os.system('rm -f File_of_Input/s*')
        sys.exit(0)

    def YES_exit(self):
        sys.exit(0)

           
# main windows (parameter window)
class Ask_sigma(tk.Tk):
    
    def __init__(self):
        super().__init__()
        # self.pack() # if using tk.Frame, keep this line !!!
        self.title("Checking if the RMS-files exist")
        self.geometry('660x460')        
        
        # Self interface
        self.setupUI()
        
    def setupUI(self):
    
        contents = ScrolledText()
        contents.pack(side=BOTTOM, expand=True, fill=BOTH)

        contents.delete('1.0', END)
        contents.insert(INSERT, '\n')
        contents.insert(INSERT, '\n')
        contents.insert(INSERT, '\n')
        contents.insert(INSERT, '    ************************************************************************'+'\n')
        contents.insert(INSERT, '    *                                                                      *'+'\n')
        contents.insert(INSERT, '    *   Do you have RMS-files and put it in the folder "File_of_Input" ?   *'+'\n')
        contents.insert(INSERT, '    *                                                                      *'+'\n')
        contents.insert(INSERT, '    *   If YES, CASA BeamBeam will use your RMS-data.                      *'+'\n')
        contents.insert(INSERT, '    *   IF NO,  CASA BeamBeam will generate RMS-data.                      *'+'\n')
        contents.insert(INSERT, '    *                                                                      *'+'\n')
        contents.insert(INSERT, '    ************************************************************************'+'\n')

        Btn = Btn_def()
        btn1 = Button(self, text='YES', bg="lightblue", command=lambda: Btn.YES_exit()).pack(side=RIGHT)
        btn2 = Button(self, text='NO', bg="lightgreen", command=lambda: Btn.NO_exit()).pack(side=RIGHT)
        
        return
        
        
if __name__ == '__main__':
    app = Ask_sigma()
    app.mainloop()
