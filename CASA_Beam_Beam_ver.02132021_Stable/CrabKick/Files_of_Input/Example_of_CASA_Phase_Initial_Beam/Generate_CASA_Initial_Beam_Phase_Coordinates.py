#  This code is used to Generate CASA_Initial_Beam_Phase_Coordinates
#  River and Vasiliy, 2019

import numpy as np
import math
import random

c_speed = 2.99792458e8

# *********** beam 1 ******************************************************
Energy1                 = 0.51e3        # unit: MeV (Kinetic Energy)
mass1                   = 0.51099906    # unit: MeV
gamma_1                 = Energy1 / mass1
p0_1                    = math.sqrt(gamma_1**2.0 - 1.0)  # p0 = beta*gamma  
Emittance_x_1           = 87.0e-9       # unit: m
Emittance_y_1           = 0.87e-9       # unit: m
beta_x_star_1           = 0.026         # unit: m
beta_y_star_1           = 0.009         # unit: m

sigma1_x_star           = 48.0e-6       # unit: m
sigma1_y_star           = 2.8e-6        # unit: m
sigma1_z                = 20.0e-3       # unit: m
sigma1_p                = 0.00045
sigma1_xp               = math.sqrt(Emittance_x_1 / beta_x_star_1)
sigma1_yp               = math.sqrt(Emittance_y_1 / beta_y_star_1)
macro_particle_number1  = 10000
# *************************************************************************

# *********** beam 2 ******************************************************
Energy2                 = 0.01e3        # unit: MeV (Kinetic Energy)
mass2                   = 0.51099906    # unit: MeV
gamma_2                 = Energy2 / mass2
p0_2                    = math.sqrt(gamma_2**2.0 - 1.0)  # p0 = beta*gamma  
Emittance_x_2           = 51.0e-9       # unit: m
Emittance_y_2           = 2.55e-9       # unit: m
beta_x_star_2           = 0.045         # unit: m
beta_y_star_2           = 0.003         # unit: m

sigma2_x_star           = 48.0e-6       # unit: m
sigma2_y_star           = 2.8e-6        # unit: m
sigma2_z                = 1.0e-3        # unit: m
sigma2_p                = 0.00045
sigma2_xp               = math.sqrt(Emittance_x_2 / beta_x_star_2)
sigma2_yp               = math.sqrt(Emittance_y_2 / beta_y_star_2)
macro_particle_number2  = 10000
# *************************************************************************

seed1_value = [] 
seed2_value = []
for i in range (6):
    seed1_value.append(random.randint(0, 10000))
    seed2_value.append(random.randint(0, 10000))

# ******** Generate Beam-1 Phase Coordinates ***********************************
np.random.seed(seed1_value[0])
x1  = np.random.normal(0, sigma1_x_star, macro_particle_number1) 
np.random.seed(seed1_value[1])
xp1 = np.random.normal(0, sigma1_xp, macro_particle_number1) 

np.random.seed(seed1_value[2])
y1  = np.random.normal(0, sigma1_y_star, macro_particle_number1) 
np.random.seed(seed1_value[3])
yp1 = np.random.normal(0, sigma1_yp, macro_particle_number1) 

''' **********  for CASA Beam-Beam and Elegant ****************************** '''
np.random.seed(seed1_value[4])
t1  = np.random.normal(0, sigma1_z, macro_particle_number1) / c_speed
np.random.seed(seed1_value[5])
p1  = abs(np.random.normal(p0_1, sigma1_p, macro_particle_number1))
''' ************************************************************************* '''
''' **********  for BB3D **************************************************** '''
# t1  = t1 * c_speed          # here: t means z
# p1  = p1 / p0_1 - 1.0       # here: p means pz
''' ************************************************************************* '''
# *******************************************************************************

# ******** Generate Beam-2 Phase Coordinates ************************************
np.random.seed(seed2_value[0])
x2  = np.random.normal(0, sigma2_x_star, macro_particle_number2) 
np.random.seed(seed2_value[1])
xp2 = np.random.normal(0, sigma2_xp, macro_particle_number2) 

np.random.seed(seed2_value[2])
y2  = np.random.normal(0, sigma2_y_star, macro_particle_number2) 
np.random.seed(seed2_value[3])
yp2 = np.random.normal(0, sigma2_yp, macro_particle_number2) 

''' **********  for CASA Beam-Beam and Elegant ****************************** '''
np.random.seed(seed2_value[4])
t2  = np.random.normal(0, sigma2_z, macro_particle_number2) / c_speed 
np.random.seed(seed1_value[5])
p2  = abs(np.random.normal(p0_2, sigma2_p, macro_particle_number2))
''' ************************************************************************** '''
''' **********  for BB3D **************************************************** '''
# t2  = t2 * c_speed          # here: t means z
# p2  = p2 / p0_2 - 1.0       # here: p means pz
''' ************************************************************************* '''
# ********************************************************************************

# ********* Save Beam-1 Phase-Coordinates ***********************************************
beam1_phase = np.array([x1, xp1, y1, yp1, t1, p1])
np.savetxt('CASA_Phase_Initial_Beam_1.casa', beam1_phase.T, fmt='% .12e', delimiter='  ')
# ***************************************************************************************

# ********* Save Beam-1 Phase-Coordinates ***********************************************
beam2_phase = np.array([x2, xp2, y2, yp2, t2, p2])
np.savetxt('CASA_Phase_Initial_Beam_2.casa', beam2_phase.T, fmt='% .12e', delimiter='  ')
# ***************************************************************************************



