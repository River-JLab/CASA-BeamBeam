# ********************************************************* #
#  This GUI is used for JLAB-CSSA Beam-Beam package.        #
#  ---                                                      #
#  River Huang and Vasiliy Morozov                          #
#  2018-2020                                                #
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

    def make_clean(self):
        current_path = os.getcwd()
        os.chdir("BeamBeam3D")
        if (os.name == 'nt'):
            os.system("make -f Makefile.Windows clean")
        else:
            os.system("make -f Makefile.Linux clean")
        os.chdir(current_path)

    def make(self, parameter, contents):
        try:
            with open(parameter, 'w') as file:
                file.write(contents.get('1.0', END))
        except FileNotFoundError:
            pass

        current_path = os.getcwd()
        os.chdir("BeamBeam3D")
        if (os.name == 'nt'):
            os.system("make -f Makefile.Windows")
        else:
            os.system("make -f Makefile.Linux")
        os.chdir(current_path)
        os._exit(0)
        
    def skip_make(self):
        os._exit(0)
           
# main windows (parameter window)
class MakeEditor(tk.Tk):
    
    def __init__(self):
        super().__init__()
        # self.pack() # if using tk.Frame, keep this line !!!
        self.title("Check and edit the make file, then click 'make clean' or 'make' or 'skip' to continue")
        self.geometry('800x920')        
        
        # Self interface
        self.setupUI()
        
    def setupUI(self):

        contents = ScrolledText()
        contents.pack(side=BOTTOM, expand=True, fill=BOTH)

        if (os.name == 'nt'):
            parameter = 'BeamBeam3D/Makefile.Windows'
        else:
            parameter = 'BeamBeam3D/Makefile.Linux'
        
        try:
            with open(parameter) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass
        
        Btn = Btn_def()
        btn1 = Button(self, text='Re-load Makefile', bg="lightblue", command=lambda: Btn.load(parameter, contents)).pack(side=LEFT)
        btn3 = Button(self, text="make", bg="yellow", command=lambda: Btn.make(parameter, contents)).pack(side=RIGHT)
        btn2 = Button(self, text='make clean', bg="lightgreen", command=lambda: Btn.make_clean()).pack(side=RIGHT)
        btn4 = Button(self, text="skip", bg="lightblue", command=lambda: Btn.skip_make()).pack(side=RIGHT)
        
        return
        
        
if __name__ == '__main__':
    app = MakeEditor()
    app.mainloop()