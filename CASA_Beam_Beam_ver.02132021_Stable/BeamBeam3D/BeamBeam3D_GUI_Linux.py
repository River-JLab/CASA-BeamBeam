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

    def run_BB3D(self):
        current_path = os.getcwd()
        os.chdir("BeamBeam3D")
        cpu_number = os.cpu_count()
        if (os.name == 'nt'):
            if not os.path.isfile('WindowsBB3D.exe'):
                print()
                print('There is no excuted file, please make sure if it was made by ...')
                print()
                os.chdir(current_path)
                sys.exit(0)
            cpu_number = int(cpu_number / 4)
            commd = 'mpiexec -n ' + str(cpu_number) + ' WindowsBB3D'
            os.system(commd)
        else:
            if not os.path.isfile('LinuxBB3D'):
                print()
                print('There is no excuted file, please make sure if it was made by ...')
                print()
                os.chdir(current_path)
                sys.exit(0)
            cpu_number = int(cpu_number / 2)
            commd = 'mpirun -n ' + str(cpu_number) + ' ./LinuxBB3D'
            os.system(commd)

        os.chdir(current_path)
        os._exit(0)

    def exit_BB3D(self):
        os._exit(0)

           
# main windows (parameter window)
class BB3D_run(tk.Tk):
    
    def __init__(self):
        super().__init__()
        # self.pack() # if using tk.Frame, keep this line !!!
        self.title("BeamBeam3D via CASA Beam-Beam GUI")
        self.geometry('800x480')        
        
        # Self interface
        self.setupUI()
        
    def setupUI(self):

        contents = ScrolledText()
        contents.pack(side=BOTTOM, expand=True, fill=BOTH)

        Readme = 'Readme_BeamBeam3D.casa'
        try:
            with open(Readme) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass
        
        Btn = Btn_def()
        btn2 = Button(self, text='Exit', bg="lightblue", command=lambda: Btn.exit_BB3D()).pack(side=RIGHT)
        btn1 = Button(self, text='Continue to run BeamBeam3D', bg="orange", command=lambda: Btn.run_BB3D()).pack(side=RIGHT)
        
        return
        
        
if __name__ == '__main__':
    app = BB3D_run()
    app.mainloop()