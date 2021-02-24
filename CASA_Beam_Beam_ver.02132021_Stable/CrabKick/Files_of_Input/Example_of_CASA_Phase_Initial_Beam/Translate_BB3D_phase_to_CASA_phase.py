#  This code is used to Translate_BB3D_phase_to_CASA_phase
#  River and Vasiliy, 2019

import numpy as np
import math

c_speed = 2.99792458e8

# *********** beam 1 ******************************************************
Energy1                 = 0.51e3        # unit: MeV (Kinetic Energy)
mass1                   = 0.51099906    # unit: MeV
gamma_1                 = Energy1 / mass1
p0_1                    = math.sqrt(gamma_1**2.0 - 1.0)  # p0 = beta*gamma  
# *************************************************************************

# *********** beam 2 ******************************************************
Energy2                 = 0.01e3        # unit: MeV (Kinetic Energy)
mass2                   = 0.51099906    # unit: MeV
gamma_2                 = Energy2 / mass2
p0_2                    = math.sqrt(gamma_2**2.0 - 1.0)  # p0 = beta*gamma  
# *************************************************************************


bb3d_phase_beam1_filename = 'CASA_Phase_out_Beam_1_turn_0.casa'
bb3d_phase_beam2_filename = 'CASA_Phase_out_Beam_2_turn_0.casa'

CASA_phase_beam1_filename = 'CASA_Phase_Initial_Beam_1.casa'
CASA_phase_beam2_filename = 'CASA_Phase_Initial_Beam_2.casa'

phase1 = np.loadtxt(bb3d_phase_beam1_filename).T
phase2 = np.loadtxt(bb3d_phase_beam2_filename).T

phase1[4] = -phase1[4] / c_speed
phase2[4] = -phase2[4] / c_speed
phase1[5] = (phase1[5] + 1.0) * p0_1
phase2[5] = (phase2[5] + 1.0) * p0_2

np.savetxt(CASA_phase_beam1_filename, phase1.T, fmt='% .12e', delimiter='  ')
np.savetxt(CASA_phase_beam2_filename, phase2.T, fmt='% .12e', delimiter='  ')


