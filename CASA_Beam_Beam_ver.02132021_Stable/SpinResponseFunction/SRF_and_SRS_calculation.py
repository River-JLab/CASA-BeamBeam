#  This is the fuction code used for Spin Response Function
#  Vasiliy and River, 2019

import os
import sys, time, math
import cmath
import numpy as np
from scipy import *
from TwissParameter_Linux import *
from SRF_and_SRS_functions import *

parameters = parameter_input()
twiss_filename      = 'Files_of_Input/' + parameters[0] 
SRF_result_filename = 'Files_of_Output/' + parameters[1]
SRS_result_filename = 'Files_of_Output/' + parameters[2]

Starting_enargy   = float(parameters[3])
Ending_enargy     = float(parameters[4])
Energy_time_steps = int(parameters[5]) + 1

Anomalous_magnetic_moment =  float(parameters[6])

delta_phi   = float(parameters[7])
delta_y     = float(parameters[8])
alpha_snake = float(parameters[9])

Ending_enargy - Starting_enargy
if (Ending_enargy == Starting_enargy):
    Energy_time_steps = 1
    delta_energy = 0.0
else:
    delta_energy = (Ending_enargy - Starting_enargy) / (Energy_time_steps - 1)

energy = np.zeros(Energy_time_steps)
for i in range (Energy_time_steps):
    energy[i] = Starting_enargy + i * delta_energy

SRS_data = open(SRS_result_filename, mode="w")
SRS_data.write('energy(GeV)      Zero-Integer-SRS\n')
SRS_data.write('\n')
for i in range (Energy_time_steps):
    F = F_integral(twiss_filename, energy[i], alpha_snake, Anomalous_magnetic_moment, delta_phi, delta_y)

    length_F = len(F[1])
    SpinFunctionOutput = SRF_result_filename + '.' + str(energy[i]) + 'GeV'
    with open(SpinFunctionOutput, "w") as f:
        for j in range (length_F):
            temp_p = str(F[0][j]) + '       '+ str(F[1][j]) + '\n'
            f.write(temp_p)
    
    SRS_data.write(str(energy[i]) + '             ' + str(F[2]) + '\n')
    # print(energy[i], F[1])
SRS_data.close()

