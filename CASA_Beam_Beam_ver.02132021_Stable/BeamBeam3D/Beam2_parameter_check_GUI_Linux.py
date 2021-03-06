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

    def save_and_continue(self, parameter, contents):
        try:
            with open(parameter, 'w') as file:
                file.write(contents.get('1.0', END))
        except FileNotFoundError:
            pass
        os._exit(0)


# main windows (parameter window)
class BeamParameterEditor(tk.Tk):
    
    def __init__(self):
        super().__init__()
        # self.pack() # if using tk.Frame, keep this line !!!
        self.title("Beam-2:           Check and edit the parameters of Beam-2")
        self.geometry('1280x680')        
        
        # Self interface
        self.setupUI()
        
    def setupUI(self):

        contents = ScrolledText()
        contents.pack(side=BOTTOM, expand=True, fill=BOTH)

        parameter = 'BeamBeam3D/beam2.in'
        
        try:
            with open(parameter) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass
        
        Btn = Btn_def()
        btn1 = Button(self, text='Re-load Parameters of Beam-2', bg="lightgreen", command=lambda: Btn.load(parameter, contents)).pack(side=LEFT)
        btn2 = Button(self, text='Save and continue', bg="lightblue", command=lambda: Btn.save_and_continue(parameter, contents)).pack(side=RIGHT)
        btn3 = Button(self, text='Skip and continue', bg="yellow", command=lambda: Btn.skip_edit()).pack(side=RIGHT)
        
        return
        
        
if __name__ == '__main__':
    app = BeamParameterEditor()
    app.mainloop()