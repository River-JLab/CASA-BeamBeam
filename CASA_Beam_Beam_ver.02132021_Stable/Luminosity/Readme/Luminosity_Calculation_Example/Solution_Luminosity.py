# ********************************************************* #
#  This is a Code for JLAB-CSSA Luminosity Calculation      #
#  ---                                                      #
#  River Huang and Vasiliy Morozov                          #
#  2020                                                     #
# ********************************************************* #


import scipy.integrate

from math import sin
from math import cos
from math import tan
from math import exp
from math import pi
from math import sqrt

from scipy.special import kn

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def Nominal_Luminosity(Number_of_particle1, Number_of_particle2, Number_of_colliding_bunches, \
                         Collision_Frequency, sigma_star_of_beam1, sigma_star_of_beam2):

    sig_1xs = sigma_star_of_beam1[0]
    sig_1ys = sigma_star_of_beam1[1]

    sig_2xs = sigma_star_of_beam2[0]
    sig_2ys = sigma_star_of_beam2[1]

    sig_xs = sqrt(0.5*(sig_1xs**2 + sig_2xs**2))
    sig_ys = sqrt(0.5*(sig_1ys**2 + sig_1ys**2))
    
    L0 = 0.25e-4 * Number_of_particle1 * Number_of_particle2 * Number_of_colliding_bunches * Collision_Frequency / (pi * sig_xs * sig_ys)
    
    return L0

def Reduction_Factor_for_Numerical_Solution(sigma_star_of_beam1, sigma_star_of_beam2, \
                                            beta_star_of_beam1, beta_star_of_beam2, \
                                            Half_Crossing_Angle, Offset_defference):

    sig_1xs = sigma_star_of_beam1[0]
    sig_1ys = sigma_star_of_beam1[1]
    sig_1s  = sigma_star_of_beam1[2]

    sig_2xs = sigma_star_of_beam2[0]
    sig_2ys = sigma_star_of_beam2[1]
    sig_2s  = sigma_star_of_beam2[2]

    sig_xs = sqrt(0.5*(sig_1xs**2 + sig_2xs**2))
    sig_ys = sqrt(0.5*(sig_1ys**2 + sig_2ys**2))
    sig_s  = sqrt(0.5*(sig_1s**2  + sig_2s**2))

    beta_1xs = beta_star_of_beam1[0]
    beta_1ys = beta_star_of_beam1[1]

    beta_2xs = beta_star_of_beam2[0]
    beta_2ys = beta_star_of_beam2[1]

    phi_x = Half_Crossing_Angle[0]
    phi_y = Half_Crossing_Angle[1]

    delta_x = Offset_defference[0]
    delta_y = Offset_defference[1]

    A = 0.5 * (sig_1xs**2 / beta_1xs**2 + sig_2xs**2 / beta_2xs**2) / sig_xs**2
    B = 0.5 * (sig_1ys**2 / beta_1ys**2 + sig_2ys**2 / beta_2ys**2) / sig_ys**2

    def Numerical_lumin(s):
        a = sin(phi_x)**2 / (sig_xs**2 * (A*s**2 + 1)) + \
            sin(phi_y)**2 / (sig_ys**2 * (B*s**2 + 1)) + \
            (1 + tan(phi_x)**2 + tan(phi_y)**2) / sig_s**2

        b = delta_x * sin(phi_x) / (sig_xs**2 * (A*s**2 + 1)) + \
            delta_y * sin(phi_y) / (sig_ys**2 * (B*s**2 + 1))

        c = delta_x**2 / (4*sig_xs**2 * (A*s**2 + 1)) + \
            delta_y**2 / (4*sig_ys**2 * (B*s**2 + 1))

        lumin = exp(-a*s**2 - b* s - c) / sqrt((A*s**2 + 1) * (B*s**2 + 1))
        return lumin

    intgrl_range = 10 * (sig_1s + sig_2s)
    integral = scipy.integrate.quad(Numerical_lumin, -intgrl_range, intgrl_range)
    
    angle_part = sqrt(2 * cos(phi_x)**2 + 2 * cos(phi_y)**2 - 3 * cos(phi_x)**2 * cos(phi_y)**2) / \
                         (cos(phi_x)**2 +     cos(phi_y)**2 -     cos(phi_x)**2 * cos(phi_y)**2)

    index = sqrt(pi) * sig_s

    Reduction_Factor = angle_part / index * integral[0]

    return Reduction_Factor

def Reduction_Factor_for_Analynic_Solution_1(sigma_star_of_beam1, sigma_star_of_beam2, \
                                             beta_star_of_beam1, beta_star_of_beam2, \
                                             Half_Crossing_Angle, Offset_defference):

    sig_1xs = sigma_star_of_beam1[0]
    sig_1ys = sigma_star_of_beam1[1]
    sig_1s  = sigma_star_of_beam1[2]

    sig_2xs = sigma_star_of_beam2[0]
    sig_2ys = sigma_star_of_beam2[1]
    sig_2s  = sigma_star_of_beam2[2]

    sig_xs = sqrt(0.5*(sig_1xs**2 + sig_2xs**2))
    sig_ys = sqrt(0.5*(sig_1ys**2 + sig_2ys**2))
    sig_s  = sqrt(0.5*(sig_1s**2  + sig_2s**2))

    beta_1xs = beta_star_of_beam1[0]
    beta_1ys = beta_star_of_beam1[1]

    beta_2xs = beta_star_of_beam2[0]
    beta_2ys = beta_star_of_beam2[1]

    phi_x = Half_Crossing_Angle[0]
    phi_y = Half_Crossing_Angle[1]

    B = 0.5 * (sig_1ys**2 / beta_1ys**2 + sig_2ys**2 / beta_2ys**2) / sig_ys**2

    a = sin(phi_x)**2 / sig_xs**2 + 1 / (cos(phi_x)**2 * sig_s**2)

    angle_part = sqrt(2 * cos(phi_x)**2 + 2 * cos(phi_y)**2 - 3 * cos(phi_x)**2 * cos(phi_y)**2) / \
                         (cos(phi_x)**2 +     cos(phi_y)**2 -     cos(phi_x)**2 * cos(phi_y)**2)

    index = sqrt(pi) * sig_s

    D = 0.5 * a / B
    if (D>512.0):
        D = sqrt(1.0 / D)
        R_L = D - 0.125 * D**3 + 0.0703125 * D**5 - 0.0732421875 * D**7
    else:
        R_L = exp(D) * kn(0, D)

    Reduction_Factor = angle_part / index * R_L / sqrt(B)
    
    return Reduction_Factor
    
def Reduction_Factor_for_Analynic_Solution_2(sigma_star_of_beam1, sigma_star_of_beam2, \
                                             beta_star_of_beam1, beta_star_of_beam2, \
                                             Half_Crossing_Angle, Offset_defference):

    sig_1xs = sigma_star_of_beam1[0]
    sig_1ys = sigma_star_of_beam1[1]
    sig_1s  = sigma_star_of_beam1[2]

    sig_2xs = sigma_star_of_beam2[0]
    sig_2ys = sigma_star_of_beam2[1]
    sig_2s  = sigma_star_of_beam2[2]

    sig_xs = sqrt(0.5*(sig_1xs**2 + sig_2xs**2))
    sig_ys = sqrt(0.5*(sig_1ys**2 + sig_2ys**2))
    sig_s  = sqrt(0.5*(sig_1s**2  + sig_2s**2))

    beta_1xs = beta_star_of_beam1[0]
    beta_1ys = beta_star_of_beam1[1]

    beta_2xs = beta_star_of_beam2[0]
    beta_2ys = beta_star_of_beam2[1]

    phi_x = Half_Crossing_Angle[0]
    phi_y = Half_Crossing_Angle[1]

    delta_x = Offset_defference[0]
    delta_y = Offset_defference[1]

    A = 0.5 * (sig_1xs**2 / beta_1xs**2 + sig_2xs**2 / beta_2xs**2) / sig_xs**2
    B = 0.5 * (sig_1ys**2 / beta_1ys**2 + sig_2ys**2 / beta_2ys**2) / sig_ys**2

    a = sin(phi_x)**2 / sig_xs**2 + sin(phi_y)**2 / sig_ys**2 + (1 + tan(phi_x)**2 + tan(phi_y)**2) / sig_s**2
    b = delta_x * sin(phi_x) / sig_xs**2 + delta_y * sin(phi_y) / sig_ys**2
    c = 0.25 * (delta_x**2 / sig_xs**2 + delta_y**2 / sig_ys**2)

    angle_part = sqrt(2 * cos(phi_x)**2 + 2 * cos(phi_y)**2 - 3 * cos(phi_x)**2 * cos(phi_y)**2) / \
                         (cos(phi_x)**2 +     cos(phi_y)**2 -     cos(phi_x)**2 * cos(phi_y)**2)

    index = sqrt(a) * sig_s

    Reduction_Factor = angle_part / index * exp(0.25 * b**2 / a - c)
    
    return Reduction_Factor


def Reduction_Factor_for_Analynic_Solution_3(sigma_star_of_beam1, sigma_star_of_beam2, \
                                             beta_star_of_beam1, beta_star_of_beam2, \
                                             Half_Crossing_Angle, Offset_defference):

    sig_1xs = sigma_star_of_beam1[0]
    sig_1ys = sigma_star_of_beam1[1]
    sig_1s  = sigma_star_of_beam1[2]

    sig_2xs = sigma_star_of_beam2[0]
    sig_2ys = sigma_star_of_beam2[1]
    sig_2s  = sigma_star_of_beam2[2]

    sig_xs = sqrt(0.5*(sig_1xs**2 + sig_2xs**2))
    sig_ys = sqrt(0.5*(sig_1ys**2 + sig_2ys**2))
    sig_s  = sqrt(0.5*(sig_1s**2  + sig_2s**2))

    phi_x = Half_Crossing_Angle[0]
    phi_y = Half_Crossing_Angle[1]

    delta_x = Offset_defference[0]
    delta_y = Offset_defference[1]

    index = sqrt(1 + sig_s**2 / sig_1xs**2 * (tan(phi_x))**2 )
    exp_term = exp(-0.25 * delta_x**2 /((sig_1xs * cos(phi_x))**2 + (sig_s * sin(phi_x))**2) - 0.25 * (delta_y / sig_1ys)**2)

    Reduction_Factor = exp_term / index
    
    return Reduction_Factor


def Luminosity(L0_luminosity, reduction_factor):
    return L0_luminosity * reduction_factor


# # ***************************************
# #              Examples                 *
# # ***************************************

# Number_of_particle1 = 1.4e10
# Number_of_particle2 = 0.18e10

# Number_of_colliding_bunches = 1.0

# sig_1xs = 48.0e-6
# sig_1ys = 2.8e-6
# sig_2xs = 48.0e-6
# sig_2ys = 2.8e-6

# sig_1s  = 20.0e-3
# sig_2s  = 1.0e-3

# beta_1xs = 0.026
# beta_1ys = 0.009
# beta_2xs = 0.045
# beta_2ys = 0.003

# Collision_Frequency = 356.0e6

# phi_x = 0.0
# phi_y = 0.0

# delta_x = 0.0
# delta_y = 0.0

# # ***************************************************************************************************

# sigma_star_of_beam1 = [sig_1xs, sig_1ys, sig_1s]
# sigma_star_of_beam2 = [sig_2xs, sig_2ys, sig_2s]

# beta_star_of_beam1 = [beta_1xs, beta_1ys]
# beta_star_of_beam2 = [beta_2xs, beta_2ys]

# Half_Crossing_Angle = [phi_x, phi_y]
# Offset_defference   = [delta_x, delta_y]



# L0_luminosity = Nominal_Luminosity(Number_of_particle1, Number_of_particle2, Number_of_colliding_bunches, \
                         # Collision_Frequency, sigma_star_of_beam1, sigma_star_of_beam2)
# print('L0_luminosity        =', L0_luminosity)
# print()

# reduction_numerical = Reduction_Factor_for_Numerical_Solution(sigma_star_of_beam1, sigma_star_of_beam2, \
                                                              # beta_star_of_beam1, beta_star_of_beam2, \
                                                              # Half_Crossing_Angle, Offset_defference)
# Numerical_Luminosity = Luminosity(L0_luminosity, reduction_numerical)
# print('Numerical Reduction  =', reduction_numerical)
# print('Numerical Luminosity =', Numerical_Luminosity)
# print()

# reduction_analytic_1 = Reduction_Factor_for_Analynic_Solution_1(sigma_star_of_beam1, sigma_star_of_beam2, \
                                                              # beta_star_of_beam1, beta_star_of_beam2, \
                                                              # Half_Crossing_Angle, Offset_defference)
# Analytic_1_Luminosity = Luminosity(L0_luminosity, reduction_analytic_1)
# print('Analytic Reduction 1 =', reduction_analytic_1)
# print('Analytic1 Luminosity =', Analytic_1_Luminosity)
# print()

# reduction_analytic_2 = Reduction_Factor_for_Analynic_Solution_2(sigma_star_of_beam1, sigma_star_of_beam2, \
                                                              # beta_star_of_beam1, beta_star_of_beam2, \
                                                              # Half_Crossing_Angle, Offset_defference)
# Analytic_2_Luminosity = Luminosity(L0_luminosity, reduction_analytic_2)
# print('Analytic Reduction 2 =', reduction_analytic_2)
# print('Analytic2 Luminosity =', Analytic_2_Luminosity)
# print()

# reduction_analytic_3 = Reduction_Factor_for_Analynic_Solution_2(sigma_star_of_beam1, sigma_star_of_beam2, \
                                                              # beta_star_of_beam1, beta_star_of_beam2, \
                                                              # Half_Crossing_Angle, Offset_defference)
# Analytic_3_Luminosity = Luminosity(L0_luminosity, reduction_analytic_3)
# print('Analytic Reduction 3 =', reduction_analytic_3)
# print('Analytic3 Luminosity =', Analytic_3_Luminosity)
# print()


