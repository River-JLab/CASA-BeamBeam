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

    def NO_Btn(self):
        os._exit(0)

    def YES_Btn(self):
        os._exit(1)

# main windows (parameter window)
class CrabKickParameterEditor(tk.Tk):
    
    def __init__(self):
        super().__init__()
        # self.pack() # if using tk.Frame, keep this line !!!
        self.title("Please select the initial phase coordinates")
        self.geometry('720x450')        
        
        # Self interface
        self.setupUI()
        
    def setupUI(self):

        contents = ScrolledText()
        contents.pack(side=BOTTOM, expand=True, fill=BOTH)
        
        contents.tag_config('orange', background = 'orange')
        contents.tag_config('lightgreen', background = 'lightgreen')
        contents.insert(INSERT, "\n")
        contents.insert(INSERT, "\n")
        contents.insert(INSERT, "\n")
        contents.insert(INSERT, "   CASA Beam-Beam found that the user's own initial phase-coordinates already existed.\n")
        contents.insert(INSERT, "   ___________________________________________________________________________________\n")
        contents.insert(INSERT, "\n")
        contents.insert(INSERT, "\n")
        contents.insert(INSERT, "   If mant to use CASA Beam-Beam initial phase-coordinates, click the button ")
        contents.insert(INSERT, "'Yes'", 'lightgreen')
        contents.insert(INSERT, "\n\n")
        contents.insert(INSERT, "   If mant to use the user's own initial phase-coordinates, click the button ")
        contents.insert(INSERT, "'No'", 'orange')
        contents.insert(INSERT, "\n")
        
        Btn = Btn_def()
        btn1 = Button(self, text="No.  Using user's phase-coordinates", bg="orange", command=lambda: Btn.NO_Btn()).pack(side=LEFT)
        btn2 = Button(self, text="Yes.  Using CASA Beam-Beam's phases-coordinates", bg="lightgreen", command=lambda: Btn.YES_Btn()).pack(side=RIGHT)

        return
        
        
if __name__ == '__main__':
    app = CrabKickParameterEditor()
    app.mainloop()
