import time
import threading
import os
from tkinter import *

CASA_logo = Tk()
ws = CASA_logo.winfo_screenwidth() # width of the screen
hs = CASA_logo.winfo_screenheight() # height of the screen

logo_pos_w = 960 # width of the CASA-LOGO + 2
logo_pos_h = 504 # height of the CASA-LOGO + 3

logo_pos_x = ws/2 - logo_pos_w/2
logo_pos_y = hs/2 - logo_pos_h/2

# initialize tk

CASA_logo.overrideredirect(1)

# CASA_logo.title('Jefferson Lab CASA Beam-Beam Code Package')
# CASA_logo.geometry('960x504+200+200')
CASA_logo.geometry('%dx%d+%d+%d' % (logo_pos_w, logo_pos_h, logo_pos_x, logo_pos_y))
bm = PhotoImage(file = 'CASA_Beam_Beam.png')
Show_logo = Label(CASA_logo, image = bm)
Show_logo.bm = bm
Show_logo.pack()

def autoClose_window():
    # for i in range(5):
        # time.sleep(1)
    time.sleep(2) # show CASA-LOGO for about 2 seconds
    CASA_logo.destroy()

t = threading.Thread(target=autoClose_window)
t.start()

CASA_logo.mainloop()
os._exit(0)

