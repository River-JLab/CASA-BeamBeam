# ********************************************************* #
#  This is a Test Inteface for JLAB-CSSA Beam-Beam code.    #
#  ---                                                      #
#  River Huang and Vasiliy Morozov                          #
#  2018-2020                                                #
# ********************************************************* #


import os
import tkinter as tk
from tkinter import *
from tkinter.scrolledtext import ScrolledText


class Btn_def():
    ''' Button Function '''
    def BB_continue(self):
        os.system('python3 casa_elegant_parameters_GUI.py')
        os.system('python3 casa_elegant_lte1_GUI.py')
        os.system('python3 casa_elegant_lte2_GUI.py')
        os.system('python3 casa_elegant_ele1_Start_GUI.py')
        os.system('python3 casa_elegant_ele1_Continue_GUI.py')
        os.system('python3 casa_elegant_ele2_Start_GUI.py')
        os.system('python3 casa_elegant_ele2_Continue_GUI.py')
        os.system('python3 CrabKick_parameter_GUI.py')
        os.system('python3 CrabKick.py')
        os.system("python3 CrabKick_Output_GUI_Linux.py")
        os._exit(0)
           
    def BB_exit(self):
        os._exit(0)
           
# main windows (parameter window)
class CrabKickParameterEditor(tk.Tk):
    
    def __init__(self):
        super().__init__()
        # self.pack() # if using tk.Frame, keep this line !!!
        self.title("JLab-CASA Beam-Beam for Crab-Kicking")
        self.geometry('580x460')        
        
        # Self interface
        self.setupUI()
        
    def setupUI(self):

        contents = ScrolledText()
        contents.pack(side=BOTTOM, expand=True, fill=BOTH)

        readme = 'Readme_CrabKick.casa'
        
        try:
            with open(readme) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass
        
        Btn = Btn_def()
        btn1 = Button(self, text="Calculate Beam-Beam interaction", bg="lightgreen", command=lambda: Btn.BB_continue()).pack(side=RIGHT)
        btn2 = Button(self, text="Exit", bg="lightblue", command=lambda: Btn.BB_exit()).pack(side=LEFT)
        return
        
        
if __name__ == '__main__':
    app = CrabKickParameterEditor()
    app.mainloop()
