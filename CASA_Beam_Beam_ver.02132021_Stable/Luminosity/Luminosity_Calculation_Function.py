# ********************************************************* #
#  This is the Code for JLAB-CSSA Luminosity Calculation    #
#  ---                                                      #
#  River Huang and Vasiliy Morozov                          #
#  2020-2021                                                     #
# ********************************************************* #

import os
import sys
import os.path
from string import *

import tkinter as tk
from tkinter import *
from tkinter.scrolledtext import ScrolledText

from math import sin
from math import cos
from math import tan
from math import exp
from math import pi
from math import sqrt

import scipy.integrate
from scipy.special import kn

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FuncFormatter
from pylab import *


def Luminosity_parameter():

    if not os.path.isfile("File_of_Input/Parameters_For_Luminosity.casa"):
        sys.exit() 
    
    # read parameters input file
    Filename_1="File_of_Input/Parameters_For_Luminosity.casa"
    File_1=open(Filename_1,"r")
    holder=File_1.readlines()
    File_1.close()
    
    temp = holder[2]
    temp = temp.split()
    particle1_mass = float(temp[0])
    particle2_mass = float(temp[1])
    
    temp = holder[3]
    temp = temp.split()
    beam1_energy = float(temp[0])
    beam2_energy = float(temp[1])
    
    temp = holder[5]
    temp = temp.split()
    beta_star_x_of_beam1 = float(temp[0])
    beta_star_x_of_beam2 = float(temp[1])
    
    temp = holder[6]
    temp = temp.split()
    beta_star_y_for_beam1 = float(temp[0])
    beta_star_y_for_beam2 = float(temp[1])
    
    temp = holder[7]
    temp = temp.split()
    Normalized_emittance_x_of_beam1 = float(temp[0])
    Normalized_emittance_x_of_beam2 = float(temp[1])
    
    temp = holder[8]
    temp = temp.split()
    Normalized_emittance_y_of_beam1 = float(temp[0])
    Normalized_emittance_y_of_beam2 = float(temp[1])
    
    temp = holder[10]
    temp = temp.split()
    Length_of_beam1 = float(temp[0])
    Length_of_beam2 = float(temp[1])
    
    temp = holder[11]
    temp = temp.split()
    Number_of_particle1 = float(temp[0])
    Number_of_particle2 = float(temp[1])
    
    temp = holder[13]
    temp = temp.split()
    Collision_Frequency = float(temp[0])
    
    temp = holder[14]
    temp = temp.split()
    Angle_x = float(temp[0])
    
    temp = holder[15]
    temp = temp.split()
    Angle_y = float(temp[0])
    
    temp = holder[16]
    temp = temp.split()
    offset_x = float(temp[0])
    
    temp = holder[17]
    temp = temp.split()
    offset_y = float(temp[0])

    temp = holder[19]
    temp = temp.split()
    Number_of_colliding_bunches = float(temp[0])

    
    return particle1_mass, particle2_mass, beam1_energy, beam2_energy, beta_star_x_of_beam1, beta_star_x_of_beam2, beta_star_y_for_beam1, beta_star_y_for_beam2, Normalized_emittance_x_of_beam1, Normalized_emittance_x_of_beam2, Normalized_emittance_y_of_beam1, Normalized_emittance_y_of_beam2, Length_of_beam1, Length_of_beam2, Number_of_particle1, Number_of_particle2, Collision_Frequency, Angle_x, Angle_y, offset_x, offset_y, Number_of_colliding_bunches

def Generate_RMS(parameters):
    particle1_mass = float(parameters[0])
    particle2_mass = float(parameters[1])
    beam1_energy = float(parameters[2])
    beam2_energy = float(parameters[3])

    beta_star_x1 = float(parameters[4])
    beta_star_x2 = float(parameters[5])
    beta_star_y1 = float(parameters[6])
    beta_star_y2 = float(parameters[7])

    Normalized_emittance_x1 = float(parameters[8])
    Normalized_emittance_x2 = float(parameters[9])
    Normalized_emittance_y1 = float(parameters[10])
    Normalized_emittance_y2 = float(parameters[11])

    Length_of_beam1 = float(parameters[12])
    Length_of_beam2 = float(parameters[13])
    
    gamma1 = beam1_energy / particle1_mass
    gamma2 = beam2_energy / particle2_mass
    
    beta_gama_beam1 = sqrt(gamma1**2.0 - 1.0)
    beta_gama_beam2 = sqrt(gamma2**2.0 - 1.0)

    sigma_x1 = sqrt(Normalized_emittance_x1 * beta_star_x1 / beta_gama_beam1)
    sigma_y1 = sqrt(Normalized_emittance_y1 * beta_star_y1 / beta_gama_beam1)
    sigma_z1 = Length_of_beam1
    sigma_x2 = sqrt(Normalized_emittance_x2 * beta_star_x2 / beta_gama_beam2)
    sigma_y2 = sqrt(Normalized_emittance_y2 * beta_star_y2 / beta_gama_beam2)
    sigma_z2 = Length_of_beam2

    SigmaData1 = str(sigma_x1) + '    ' + str(sigma_y1)  + '    ' + str(sigma_z1) + '\n'
    SigmaData2 = str(sigma_x2) + '    ' + str(sigma_y2)  + '    ' + str(sigma_z2) + '\n'

    with open('File_of_Input/sigma_beam1.casa', 'w') as output1:
        output1.write(SigmaData1)
    with open('File_of_Input/sigma_beam2.casa', 'w') as output2:
        output2.write(SigmaData2)

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

def Luminosity_Calculation(L0_luminosity, reduction_factor):
    return L0_luminosity * reduction_factor

def Luminosity_2D_Curve(Numerical_data, Analytic1_data, Analytic2_data, Analytic3_data):

    # import matplotlib.pyplot as plt
    # from mpl_toolkits.mplot3d import Axes3D
    # from matplotlib.ticker import FuncFormatter

    # from pylab import tick_params
    # from pylab import ticklabel_format

    turn_number = range(len(Numerical_data))

    font1 = {'weight': 'bold',   'size': 22, }
    font2 = {'weight': 'bold',   'size': 18, }

    plt.rcParams['figure.figsize'] = (8.8, 8.0)
    plt.rcParams.update({'font.size': 12, 'font.family': 'serif'})

    ax1 = plt.subplot(3, 1, 1)
    tick_params(top=False, bottom=True, left=True, right=False, which='both', direction='in')
    ax1.tick_params(labeltop=False, labelbottom=False, labelleft=True, labelright=False)
    ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
    # plt.ylim(-1.2, 1.2)
    plot1 = ax1.plot(turn_number, Numerical_data, '-o', markerfacecolor='none', label='Numerical Luminosity')
    plot2 = ax1.plot(turn_number, Analytic1_data, '--r', label='Analytic-1 Luminosity')
    ax1.legend(loc=1) # 1 means position of Legend at the position of upper right corner
    # plt.title('Luminosity vs Turn-Number', font1)

    ax2 = plt.subplot(3, 1, 2)
    ax2.tick_params(top=False, bottom=True, left=True, right=False, which='both', direction='in')
    ax2.tick_params(labeltop=False, labelbottom=False, labelleft=True, labelright=False)
    ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
    # plt.ylim(-1.2, 1.2)
    plot1 = ax2.plot(turn_number, Numerical_data, '-o', markerfacecolor='none', label='Numerical Luminosity')
    plot2 = ax2.plot(turn_number, Analytic2_data, '--r',label='Analytic-2 Luminosity')
    ax2.legend(loc=1) # 1 means position of Legend at the position of upper upper corner
    plt.ylabel('Luminosity $(cm^{-2}s^{-1})$', font2)

    ax3 = plt.subplot(3, 1, 3)
    ax3.tick_params(top=False, bottom=True, left=True, right=False, which='both', direction='in')
    ax3.tick_params(labeltop=False, labelbottom=True, labelleft=True, labelright=False)
    ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
    # plt.ylim(-1.2, 1.2)
    plot1 = plt.plot(turn_number, Numerical_data, '-o', markerfacecolor='none', label='Numerical Luminosity')
    plot2 = plt.plot(turn_number, Analytic3_data, '--r',label='Analytic-3 Luminosity')
    ax3.legend(loc=1) # 1 means position of Legend at the position of upper right corner
    plt.xlabel('Turn-Number', font2)

    plt.show()







