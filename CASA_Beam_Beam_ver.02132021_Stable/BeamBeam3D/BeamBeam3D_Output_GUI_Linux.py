# ********************************************************* #
#  This is a Test Inteface for JLAB-CSSA Beam-Beam code.    #
#  ---                                                      #
#  River Huang and Vasiliy Morozov                          #
#  2018-2020                                                #
# ********************************************************* #


import numpy as np
import os
import tkinter as tk
from tkinter import *
from tkinter.scrolledtext import ScrolledText

import shutil
import glob

def save_CASA_data():
    if not os.path.isfile("BeamBeam3D/luminosity.data"):
        print()
        print('    Luminosity-file does not exist.')
        print()
        sys.exit(0) 
    if not os.path.isfile("BeamBeam3D/fort.24"):
        print()
        print('    RMS_x of Beam-1 does not exist.')
        print()
        sys.exit() 
    if not os.path.isfile("BeamBeam3D/fort.25"):
        print()
        print('    RMS_y of Beam-1 does not exist.')
        print()
        sys.exit() 
    if not os.path.isfile("BeamBeam3D/fort.26"):
        print()
        print('    RMS_x of Beam-1 does not exist.')
        print()
        sys.exit() 
    if not os.path.isfile("BeamBeam3D/fort.34"):
        print()
        print('    RMS_x of Beam-2 does not exist.')
        print()
        sys.exit() 
    if not os.path.isfile("BeamBeam3D/fort.35"):
        print()
        print('    RMS_y of Beam-2 does not exist.')
        print()
        sys.exit() 
    if not os.path.isfile("BeamBeam3D/fort.36"):
        print()
        print('    RMS_x of Beam-2 does not exist.')
        print()
        sys.exit() 


    def len_of_file(filename):
        with open(filename) as file_line:
            return len(list(file_line))
        file_line.close()
    
    if not os.path.isfile("BeamBeam3D/beam1.in"):
        print()
        print('   The file beam1.in does not exist.')
        print()
        sys.exit() 

    if not os.path.isfile("BeamBeam3D/beam2.in"):
        print()
        print('   The file beam2.in does not exist.')
        print()
        sys.exit() 

    Filename1="BeamBeam3D/beam1.in"
    Filename2="BeamBeam3D/beam2.in"

    File1=open(Filename1,"r")
    holder=File1.readlines()
    File1.close()

    temp = holder[1]
    temp = temp.split()
    macro_number1 = float(temp[1])
    
    temp = holder[9]
    temp = temp.split()
    number1 = float(temp[0])
    
    File2=open(Filename2,"r")
    holder=File2.readlines()
    File2.close()

    temp = holder[1]
    temp = temp.split()
    macro_number2 = float(temp[1])

    temp = holder[9]
    temp = temp.split()
    number2 = float(temp[0])

    # ratio = number1 * number2 / macro_number1 / macro_number2
    ratio = 1


    def Delete_directory(directory_name):
        for root, dirs, files in os.walk(directory_name, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        os.rmdir(directory_name)

    if not os.path.exists('BeamBeam3D_Output_Data_CASA'):
        os.mkdir('BeamBeam3D_Output_Data_CASA')
    else:
        Delete_directory('BeamBeam3D_Output_Data_CASA/')
        os.mkdir('BeamBeam3D_Output_Data_CASA')


    number_of_line_luminosity = len_of_file('BeamBeam3D/luminosity.data')
    Temp = np.loadtxt('BeamBeam3D/luminosity.data')
    if (number_of_line_luminosity == 1):
        Luminosity = Temp[4]
    else:
        Luminosity = Temp[:, 4]
    del Temp
    Luminosity = Luminosity * (ratio /10000)        # /10000 means translating meter to centimeter


    np.savetxt("BeamBeam3D_Output_Data_CASA/Luminosity_Data.casa", np.transpose(Luminosity))
    del Luminosity
    
    number_of_line_sigma_beam1 = len_of_file('BeamBeam3D/fort.24')
    Temp = np.loadtxt('BeamBeam3D/fort.24')
    if (number_of_line_sigma_beam1 == 1):
        sigma_x_of_Beam1 = Temp[2]
    else:
        sigma_x_of_Beam1 = Temp[:, 2]
    del Temp
    sigma_x_of_Beam1 = sigma_x_of_Beam1 * 100.0

    Temp = np.loadtxt('BeamBeam3D/fort.25')
    if (number_of_line_sigma_beam1 == 1):
        sigma_y_of_Beam1 = Temp[2]
    else:
        sigma_y_of_Beam1 = Temp[:, 2]
    del Temp
    sigma_y_of_Beam1 = sigma_y_of_Beam1 * 100.0

    Temp = np.loadtxt('BeamBeam3D/fort.26')
    if (number_of_line_sigma_beam1 == 1):
        sigma_z_of_Beam1 = Temp[2]
    else:
        sigma_z_of_Beam1 = Temp[:, 2]
    del Temp
    sigma_z_of_Beam1 = sigma_z_of_Beam1 * 100.0

    temp_arr = np.array([sigma_x_of_Beam1, sigma_y_of_Beam1, sigma_z_of_Beam1])
    np.savetxt("BeamBeam3D_Output_Data_CASA/Sigma_beam1.casa", temp_arr.T, fmt='%0.11f    %0.11f    %0.11f')
    del sigma_x_of_Beam1
    del sigma_y_of_Beam1
    del sigma_z_of_Beam1
    del temp_arr


    number_of_line_sigma_beam2 = len_of_file('BeamBeam3D/fort.34')
    Temp = np.loadtxt('BeamBeam3D/fort.34')
    if (number_of_line_sigma_beam2 == 1):
        sigma_x_of_Beam2 = Temp[2]
    else:
        sigma_x_of_Beam2 = Temp[:, 2]
    del Temp
    sigma_x_of_Beam2 = sigma_x_of_Beam2 * 100.0

    Temp = np.loadtxt('BeamBeam3D/fort.35')
    if (number_of_line_sigma_beam2 == 1):
        sigma_y_of_Beam2 = Temp[2]
    else:
        sigma_y_of_Beam2 = Temp[:, 2]
    del Temp
    sigma_y_of_Beam2 = sigma_y_of_Beam2 * 100.0

    Temp = np.loadtxt('BeamBeam3D/fort.36')
    if (number_of_line_sigma_beam2 == 1):
        sigma_z_of_Beam2 = Temp[2]
    else:
        sigma_z_of_Beam2 = Temp[:, 2]
    del Temp
    sigma_z_of_Beam2 = sigma_z_of_Beam2 * 100.0

    temp_arr = np.array([sigma_x_of_Beam2, sigma_y_of_Beam2, sigma_z_of_Beam2])
    np.savetxt("BeamBeam3D_Output_Data_CASA/Sigma_beam2.casa", temp_arr.T, fmt='%0.11f    %0.11f    %0.11f')
    del sigma_x_of_Beam2
    del sigma_y_of_Beam2
    del sigma_z_of_Beam2
    del temp_arr

    file_casa = glob.glob('BeamBeam3D/*.casa')
    for i in file_casa:
        try:
            shutil.move(i, 'BeamBeam3D_Output_Data_CASA/')
        except:
            print("Failed to move files: ", sys.exc_info())

    os.remove('BeamBeam3D/luminosity.data')


class Btn_def():

    def load_luminosity(self, Readme, contents):
        try:
            with open(Readme) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, 'Luminosity (1/cm^2s）'+'\n')
                contents.insert(INSERT, '------------------------'+'\n')
                
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass

    def load_sigma_beam1(self, Readme, contents):
        try:
            with open(Readme) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, 'Beam1 ' + '\u03C3' + 'x(cm)    Beam1 ' + '\u03C3' + 'y(cm)    Beam1 ' + '\u03C3' + 'z(cm)\n')
                contents.insert(INSERT, '------------------------------------------------'+'\n')
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass

    def load_sigma_beam2(self, Readme, contents):
        try:
            with open(Readme) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, 'Beam2 ' + '\u03C3' + 'x(cm)    Beam2 ' + '\u03C3' + 'y(cm)    Beam2 ' + '\u03C3' + 'z(cm)\n')
                contents.insert(INSERT, '------------------------------------------------'+'\n')
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass

    def exit_BB3D(self):
        sys.exit(0)

           
# main windows (parameter window)
class BB3D_run(tk.Tk):
    
    def __init__(self):
        super().__init__()
        # self.pack() # if using tk.Frame, keep this line !!!
        self.title("BeamBeam3D via CASA Beam-Beam GUI")
        self.geometry('720x820')        
        
        # Self interface
        self.setupUI()
        save_CASA_data()
        
    def setupUI(self):
    
        contents = ScrolledText()
        contents.pack(side=BOTTOM, expand=True, fill=BOTH)

        Readme = 'BeamBeam3D_Instructions/BeamBeam3D_Data_Instruction.casa'
        Readme_luminosity  = 'BeamBeam3D_Output_Data_CASA/Luminosity_Data.casa'
        Readme_sigma_beam1 = 'BeamBeam3D_Output_Data_CASA/Sigma_beam1.casa'
        Readme_sigma_beam2 = 'BeamBeam3D_Output_Data_CASA/Sigma_beam2.casa'

        try:
            with open(Readme) as file:
                contents.delete('1.0', END)
                contents.insert(INSERT, file.read())
        except FileNotFoundError:
            pass
        
        Btn = Btn_def()
        btn1 = Button(self, text='Show Luminosity', bg="lightgreen", command=lambda: Btn.load_luminosity(Readme_luminosity, contents)).pack(side=LEFT)
        btn2 = Button(self, text='Show RMS-size of Beam-1', bg="yellow", command=lambda: Btn.load_sigma_beam1(Readme_sigma_beam1, contents)).pack(side=LEFT)
        btn3 = Button(self, text='Show RMS-size of Beam-2', bg="lightblue", command=lambda: Btn.load_sigma_beam2(Readme_sigma_beam2, contents)).pack(side=LEFT)
        btn4 = Button(self, text='Exit', bg="orange", command=lambda: Btn.exit_BB3D()).pack(side=RIGHT)
        
        return
        
        
if __name__ == '__main__':
    app = BB3D_run()
    app.mainloop()