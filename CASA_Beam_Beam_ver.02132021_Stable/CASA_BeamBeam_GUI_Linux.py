# ********************************************************* #
#  This is a Test Inteface for JLAB-CSSA Beam-Beam code.    #
#  ---                                                      #
#  River Huang and Vasiliy Morozov                          #
#  2018-2019                                                #
# ********************************************************* #

import os
from tkinter import *
from tkinter import filedialog
import urllib.request


class App_Listed:
    def __init__(self, master):
    
        self.information = StringVar()
        self.information.set('Thanks for using JLAB-CASA Beam-Beam code package.\n\n River Huang and Vasiliy Morozov\n 2018-2020')
        
        # tkinter   index of row and column; cell located
        Luminosity_btn     = Button(master, text="Luminosity",             width=30,height=4,borderwidth=3,
                                  command=self.Luminosity)
        CrabKick_btn       = Button(master, text="Crab Kick",              width=30,height=4,borderwidth=3,
                                  command=self.CrabKick)
        BeamBeam3D_btn     = Button(master, text="BeamBeam3D",             width=30,height=4,borderwidth=3,
                                  command=self.BeamBeam3D)
        SpinResponse_btn   = Button(master, text="Spin Response Function", width=30,height=4,borderwidth=3,
                                  command=self.SpinResponse)
        OtherAPP_btn          = Button(master, text="OtherAPP",                  width=30,height=4,borderwidth=3,
                                  command=self.OtherAPP)

        Luminosity_btn.grid  (row=0,column=0,ipadx=5,ipady=3,padx=10,pady=3)
        CrabKick_btn.grid    (row=1,column=0,ipadx=5,ipady=3,padx=10,pady=3)
        BeamBeam3D_btn.grid  (row=2,column=0,ipadx=5,ipady=3,padx=10,pady=3)
        SpinResponse_btn.grid(row=3,column=0,ipadx=5,ipady=3,padx=10,pady=3)
        OtherAPP_btn.grid       (row=4,column=0,ipadx=5,ipady=3,padx=10,pady=3)
        
        # self.information.grid(row=0,column=1,rowspan=2,ipadx=200,ipady=15,padx=10,pady=0)    #,padx=20,pady=10,
        self.information_label = Label(master,textvariable=self.information,wraplength=250,fg='black',bg='lightblue',)
        self.information_label.grid(row=0,column=1,rowspan=5,ipadx=100,ipady=15,padx=10,pady=0)    #,padx=20,pady=10,

    def Luminosity(self):
        current_path = os.getcwd()
        os.chdir("Luminosity")
        os.system("python3 Luminosity_GUI_Linux.py")
        os.chdir(current_path)
        
        
    def CrabKick(self):
        current_path = os.getcwd()
        os.chdir("CrabKick")
        os.system("python3 CrabKick_GUI_Linux.py")
        # os.system("python3 CrabKick_Output_GUI_Linux.py")
        os.chdir(current_path)

    def BeamBeam3D(self):
        current_path = os.getcwd()
        os.chdir("BeamBeam3D")
        os.system("python3 make_file_GUI_Linux.py")
        os.system("python3 Beam1_parameter_check_GUI_Linux.py")
        os.system("python3 Beam2_parameter_check_GUI_Linux.py")
        os.system("python3 BeamBeam3D_GUI_Linux.py")
        os.system("python3 BeamBeam3D_Output_GUI_Linux.py")
        os.chdir(current_path)

    def SpinResponse(self):
        current_path = os.getcwd()
        os.chdir("SpinResponseFunction")
        os.system("python3 SpinResponseFunction.py")
        os.chdir(current_path)

    def OtherAPP(self):
        current_path = os.getcwd()
        os.chdir("OtherAPP")
        os.chdir(current_path)

if __name__=='__main__':

    os.system("python3 CASA_LOGO_Linux.py")

    Master_Window = Tk()
    Master_Window.title("JLAB-CSSA Beam-Beam code")
    # Master_Window.columnconfigure(0, weight=3)
    Master_Window.columnconfigure(0, weight=3)
    Master_Window.columnconfigure(1, weight=7)
    Master_Window.rowconfigure(0, weight=1)
    Master_Window.rowconfigure(1, weight=1)
    Master_Window.rowconfigure(2, weight=1)
    Master_Window.rowconfigure(3, weight=1)
    Master_Window.rowconfigure(4, weight=1)
    Master_Window.geometry('640x260') # Set up the window size 640x480
    # Master_Window.geometry('640x800') # Set up the window size 640x800
    app = App_Listed(Master_Window)
    Master_Window.mainloop()
