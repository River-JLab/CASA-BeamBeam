#  This is the fuction code used for testing KICK of CRIB Crossing
#  River and Vasiliy, 2018-2019

from sdds_load_save import *
from casa_elegant_functions import *

import os, sys, time, math, cmath

import numpy as np

from scipy import special
from scipy.special import kn
import scipy.integrate

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

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FuncFormatter
from pylab import *

import shutil
import glob


def Delete_directory(directory_name):
    for root, dirs, files in os.walk(directory_name, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    os.rmdir(directory_name)

def Average_of_list(lst):
    return np.mean(lst)

def Product_of_position_and_momentum(position, momentum):
    return np.multiply(position, momentum)

def product_average(position, momentum):
    return np.mean(np.multiply(position, momentum))

def square_value_in_a_list(lst):
    return np.multiply(lst, lst)

def Data_mode(input):
    Data_input = SDDS(0)
    mode = Data_input.Data_Storage_mode(input)
    return mode

def SDDS_Data_Mode(input1, input2):
    mode = Data_mode(input1)
    if (mode == 0):
        print ('Beam-1 data storage mode is', 'ASCII', '\n')
    if (mode == 1):
        print ('Beam-1 data storage mode is', 'BINARY', '\n')

    mode = Data_mode(input2)
    if (mode == 0):
        print ('Beam-2 data storage mode is', 'ASCII', '\n')
    if (mode == 1):
        print ('Beam-2 data storage mode is', 'BINARY', '\n')
    print ()

def Crab_Kick_Input_Parameters(input_file_name):
 
    file = open(input_file_name, "r")
    holder = file.readlines()
    file.close()

    temp = holder[2]
    temp = temp.split()
    Beam_name_1_in = 'Files_of_Input/' + temp[0]
    Beam_name_2_in = 'Files_of_Input/' + temp[1]

    temp = holder[3]
    temp = temp.split()
    Beam_name_1_out = 'Files_of_Output/' + temp[0]
    Beam_name_2_out = 'Files_of_Output/' + temp[1]

    temp = holder[4]
    temp = temp.split()
    Sigma_beam1_file = 'Files_of_Output/' + temp[0]
    Sigma_beam2_file = 'Files_of_Output/' + temp[1]

    temp = holder[6]
    temp = temp.split()
    Particle_1_Radius = float(temp[0])
    Particle_2_Radius = float(temp[1])

    temp = holder[7]
    temp = temp.split()
    Particle_1_Numbers = float(temp[0])
    Particle_2_Numbers = float(temp[1])

    temp = holder[8]
    temp = temp.split()
    c_speed = float(temp[0])

    temp = holder[9]
    temp = temp.split()
    case_flag = int(temp[0])

    temp = holder[11]
    temp = temp.split()
    Crossing_angle_xs = float(temp[0])

    temp = holder[12]
    temp = temp.split()
    Crossing_angle_ys = float(temp[0])

    temp = holder[13]
    temp = temp.split()
    Beam_1_tilt = float(temp[0])
    Beam_2_tilt = float(temp[1])

    temp = holder[15]
    temp = temp.split()
    Collision_Frequency = float(temp[0])

    temp = holder[16]
    temp = temp.split()
    force_sign = float(temp[0])

    temp = holder[17]
    temp = temp.split()
    number_of_turn = int(temp[0])

    temp = holder[18]
    temp = temp.split()
    Beam_1_Slices = int(temp[0])
    Beam_2_Slices = int(temp[1])

    temp = holder[20]
    temp = temp.split()
    sdds_output_frequency = int(temp[0])

    temp = holder[21]
    temp = temp.split()
    recording_phase = temp[0]

    temp = holder[22]
    temp = temp.split()
    phase_output = int(temp[0])

    # return Beam_name_1_in, \                [0]
           # Beam_name_1_out, \               [1]
           # Beam_name_2_in, \                [2]
           # Beam_name_2_out, \               [3]
           # Beam_1_Slices, \                 [4]
           # Beam_2_Slices, \                 [5]
           # Particle_1_Radius, \             [6]
           # Particle_2_Radius, \             [7]
           # Crossing_angle_xs, \             [8]
           # Crossing_angle_ys, \             [9]
           # Beam_1_tilt, \                   [10]
           # Beam_2_tilt, \                   [11]
           # c_speed, \                       [12]
           # case_flag, \                     [13]
           # Collision_Frequency, \           [14]
           # force_sign, \                    [15]
           # number_of_turn, \                [16]
           # Sigma_beam1_file, \              [17]
           # Sigma_beam2_file, \              [18]
           # sdds_output_frequency, \         [19]
           # Particle_1_Numbers, \            [20]
           # Particle_2_Numbers, \            [21]
           # recording_phase, \               [22]
           # phase_output                     [23]

    return Beam_name_1_in, \
           Beam_name_1_out, \
           Beam_name_2_in, \
           Beam_name_2_out, \
           Beam_1_Slices, \
           Beam_2_Slices, \
           Particle_1_Radius, \
           Particle_2_Radius, \
           Crossing_angle_xs, \
           Crossing_angle_ys, \
           Beam_1_tilt, \
           Beam_2_tilt, \
           c_speed, \
           case_flag, \
           Collision_Frequency, \
           force_sign, \
           number_of_turn, \
           Sigma_beam1_file, \
           Sigma_beam2_file, \
           sdds_output_frequency, \
           Particle_1_Numbers, \
           Particle_2_Numbers, \
           recording_phase, \
           phase_output

def Parameters_from_SDDS(input):

    Parameter_input = SDDS(0)
    Parameter_input.Load_Data(input)
    
    for i in range (len(Parameter_input.parameterName)):
        if (Parameter_input.parameterName[i] == 'Step'):
            p00 = Parameter_input.parameterData[i][0]
        if (Parameter_input.parameterName[i] == 'pCentral'):
            p01 = Parameter_input.parameterData[i][0]
        if (Parameter_input.parameterName[i] == 'Charge'):
            p02 = Parameter_input.parameterData[i][0]
        if (Parameter_input.parameterName[i] == 'Particles'):
            p03 = Parameter_input.parameterData[i][0]
        if (Parameter_input.parameterName[i] == 'IDSlotsPerBunch'):
            p04 = Parameter_input.parameterData[i][0]
        if (Parameter_input.parameterName[i] == 'Pass'):
            p05 = Parameter_input.parameterData[i][0]
        if (Parameter_input.parameterName[i] == 'PassLength'):
            p06 = Parameter_input.parameterData[i][0]
        if (Parameter_input.parameterName[i] == 'PassCentralTime'):
            p07 = Parameter_input.parameterData[i][0]
        if (Parameter_input.parameterName[i] == 'ElapsedTime'):
            p08 = Parameter_input.parameterData[i][0]
        if (Parameter_input.parameterName[i] == 'ElapsedCoreTime'):
            p09 = Parameter_input.parameterData[i][0]
        if (Parameter_input.parameterName[i] == 's'):
            p10 = Parameter_input.parameterData[i][0]
        if (Parameter_input.parameterName[i] == 'Description'):
            p11 = Parameter_input.parameterData[i][0]
        if (Parameter_input.parameterName[i] == 'PreviousElementName'):
            p12 = Parameter_input.parameterData[i][0]

    return p00, p01, p02, p03, p04, p05, p06, p07, p08, p09, p10, p11, p12

def Data_from_SDDS(input): 

    Data_input = SDDS(0)
    Data_input.Load_Data(input)
    
    x           = Data_input.columnData[0][0]
    xp          = Data_input.columnData[1][0]           # px/pz
    y           = Data_input.columnData[2][0]
    yp          = Data_input.columnData[3][0]           # py/pz
    t           = Data_input.columnData[4][0]
    p           = Data_input.columnData[5][0]
    dt          = Data_input.columnData[6][0]
    particle_ID = Data_input.columnData[7][0]

    return x, xp, y, yp, t, p, dt, particle_ID

def Gamma(g1, g2, angle):
    s1 = sqrt(1.0 - 1.0 / g1**2.0) * cos(angle)
    s2 = sqrt(1.0 - 1.0 / g2**2.0) * cos(angle)
    s0 = (s1 + s2) / (1 + s1 * s2)
    gamma = 1 / sqrt(1.0 - s0**2.0)
    
    # print('g1 = ', g1, 'g2 = ', g2, 'gamma = ', gamma)
    # input('>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    return gamma

def Lab_frame_to_Boost_frame(Beam_Data, Beam_parameters, Crossing_angle, Beam_Tilt, Light_speed):

    sine        = sin(Crossing_angle)
    cosine      = cos(Crossing_angle)
    tangent     = tan(Crossing_angle)
    tangent_sqr = tangent**2

    tangent_tilt = tan(Beam_Tilt)
    
    p_Central         = Beam_parameters[1]
    particle_number   = Beam_parameters[3]
    Pass_Central_Time = Beam_parameters[7]

    x_beam_bst  = []
    y_beam_bst  = []
    z_beam_bst  = []
    
    px_beam_bst = []
    py_beam_bst = []
    pz_beam_bst = []

    ID_beam_bst = [i for i in range(0, particle_number)]

    def Lab_to_Boost_append(x, xp, y, yp, t, p):
        z = (Pass_Central_Time - t) * Light_speed
        
        pz = (abs(p) - p_Central) / p_Central
        pz_plus_1 = pz + 1.0
        
        h = pz_plus_1 - sqrt(abs(pz_plus_1*pz_plus_1 - xp**2 - yp**2))     # Hamiltonian
        
        px_bst = (xp - tangent * h) / cosine
        py_bst = yp / cosine
        pz_bst = pz - tangent * xp + tangent_sqr * h
        pz__bst_plus_1 = pz_bst + 1.0

        px_beam_bst.append(px_bst)
        py_beam_bst.append(py_bst)
        pz_beam_bst.append(pz_bst)

        square_root = sqrt(abs(pz__bst_plus_1*pz__bst_plus_1 - px_bst*px_bst - py_bst*py_bst))
        
        hx_bst = px_bst / square_root
        hy_bst = py_bst / square_root
        hz_bst = 1.0 - pz__bst_plus_1 / square_root
    
        x_bst = tangent * z + (1.0 + hx_bst * sine) * x + tangent_tilt * z
        y_bst = y + hy_bst * sine * x
        z_bst = z / cosine + hz_bst * sine * x
        
        x_beam_bst.append(x_bst)
        y_beam_bst.append(y_bst)
        z_beam_bst.append(z_bst)
       
    list(map(Lab_to_Boost_append, Beam_Data[0], Beam_Data[1], Beam_Data[2], Beam_Data[3], Beam_Data[4], Beam_Data[5]))
    
    return x_beam_bst, y_beam_bst, z_beam_bst, px_beam_bst, py_beam_bst, pz_beam_bst, ID_beam_bst

def Beam_bst(xyz_pxpypz_ID_bst, number_of_Beam_Slices):
    Beam_x_bst  = [[] for row in range(number_of_Beam_Slices)]
    Beam_y_bst  = [[] for row in range(number_of_Beam_Slices)]
    Beam_z_bst  = [[] for row in range(number_of_Beam_Slices)]
    Beam_px_bst = [[] for row in range(number_of_Beam_Slices)]
    Beam_py_bst = [[] for row in range(number_of_Beam_Slices)]
    Beam_pz_bst = [[] for row in range(number_of_Beam_Slices)]
    Beam_ID_bst = [[] for row in range(number_of_Beam_Slices)]

    particle_number = len(xyz_pxpypz_ID_bst[0])

    max_temp = max(xyz_pxpypz_ID_bst[2])
    min_temp = min(xyz_pxpypz_ID_bst[2])
    delta_temp = (max_temp - min_temp) / number_of_Beam_Slices

    max_slice_ID = number_of_Beam_Slices - 1
    
    def bst_append(x, y, z, px, py, pz, ID):
        criterion = int((z - min_temp) / delta_temp)
        ID_of_Slice = max_slice_ID - min(max_slice_ID, criterion)

        Beam_x_bst[ID_of_Slice].append(x)
        Beam_y_bst[ID_of_Slice].append(y)
        Beam_z_bst[ID_of_Slice].append(z)
        Beam_px_bst[ID_of_Slice].append(px)
        Beam_py_bst[ID_of_Slice].append(py)
        Beam_pz_bst[ID_of_Slice].append(pz)
        Beam_ID_bst[ID_of_Slice].append(ID)

    list(map(bst_append, xyz_pxpypz_ID_bst[0], xyz_pxpypz_ID_bst[1], xyz_pxpypz_ID_bst[2], xyz_pxpypz_ID_bst[3], xyz_pxpypz_ID_bst[4], xyz_pxpypz_ID_bst[5], xyz_pxpypz_ID_bst[6]))

    # ******* delete the empty list *******************
    Beam_x_bst = [x for x in Beam_x_bst if x]
    Beam_y_bst = [x for x in Beam_y_bst if x]
    Beam_z_bst = [x for x in Beam_z_bst if x]
    Beam_px_bst = [x for x in Beam_px_bst if x]
    Beam_py_bst = [x for x in Beam_py_bst if x]
    Beam_pz_bst = [x for x in Beam_pz_bst if x]
    Beam_ID_bst = [x for x in Beam_ID_bst if x]
    # *************************************************

    return Beam_x_bst, Beam_y_bst, Beam_z_bst, Beam_px_bst, Beam_py_bst, Beam_pz_bst, Beam_ID_bst

def Group_Beam_Slice_in_Boost_frame(xyz_pxpypz_ID_bst, number_of_Beam_Slices):

    bst = Beam_bst(xyz_pxpypz_ID_bst, number_of_Beam_Slices)
    Beam_x_bst  = bst[0]
    Beam_y_bst  = bst[1]
    Beam_z_bst  = bst[2]
    Beam_px_bst = bst[3]
    Beam_py_bst = bst[4]
    Beam_pz_bst = bst[5]
    Beam_ID_bst = bst[6]

    cnt = 0
    while(cnt<len(Beam_x_bst)):
        if (len(Beam_x_bst[cnt])>1):
            cnt += 1
            continue
        if (cnt <= (len(Beam_x_bst)-1)/2):
            Beam_x_bst[cnt+1]  += Beam_x_bst[cnt]
            Beam_y_bst[cnt+1]  += Beam_y_bst[cnt]
            Beam_z_bst[cnt+1]  += Beam_z_bst[cnt]
            Beam_px_bst[cnt+1] += Beam_px_bst[cnt]
            Beam_py_bst[cnt+1] += Beam_py_bst[cnt]
            Beam_pz_bst[cnt+1] += Beam_pz_bst[cnt]
            Beam_ID_bst[cnt+1] += Beam_ID_bst[cnt]
            del Beam_x_bst[cnt]
            del Beam_y_bst[cnt]
            del Beam_z_bst[cnt]
            del Beam_px_bst[cnt]
            del Beam_py_bst[cnt]
            del Beam_pz_bst[cnt]
            del Beam_ID_bst[cnt]
            cnt += 1
            continue
        if (cnt > (len(Beam_x_bst)-1)/2):
            Beam_x_bst[cnt-1]  += Beam_x_bst[cnt]
            Beam_y_bst[cnt-1]  += Beam_y_bst[cnt]
            Beam_z_bst[cnt-1]  += Beam_z_bst[cnt]
            Beam_px_bst[cnt-1] += Beam_px_bst[cnt]
            Beam_py_bst[cnt-1] += Beam_py_bst[cnt]
            Beam_pz_bst[cnt-1] += Beam_pz_bst[cnt]
            Beam_ID_bst[cnt-1] += Beam_ID_bst[cnt]
            del Beam_x_bst[cnt]
            del Beam_y_bst[cnt]
            del Beam_z_bst[cnt]
            del Beam_px_bst[cnt]
            del Beam_py_bst[cnt]
            del Beam_pz_bst[cnt]
            del Beam_ID_bst[cnt]
            cnt += 1
            continue

    number_Slices = len(Beam_x_bst)

    return Beam_x_bst, Beam_y_bst, Beam_z_bst, Beam_px_bst, Beam_py_bst, Beam_pz_bst, Beam_ID_bst, number_Slices

def Slice_center(grouped_z_position_bst, number_of_Beam_Slices):
    z_center = np.zeros(number_of_Beam_Slices)
    for i in range(number_of_Beam_Slices):
        z_center[i] = Average_of_list(grouped_z_position_bst[i])  

    return z_center

def S_for_collision(grouped_z_position_bst, number_of_Beam_Slices):
    S_collision = np.zeros(number_of_Beam_Slices)
    for i in range(number_of_Beam_Slices):
        if(len(grouped_z_position_bst[i]) == 0):
            S_collision[i] = 0
        else:
            S_collision[i] = max(grouped_z_position_bst[i]) - min(grouped_z_position_bst[i])

    return S_collision

def Sigma_Calculation(x):    
    return np.std(x)

def Slice_series_when_interaction(Beam_1_Slice_number, Beam_2_Slice_number, Step_of_interaction):
    n1 = Beam_1_Slice_number
    n2 = Beam_2_Slice_number
    step = Step_of_interaction

    interaction_beam_1_slices = []
    interaction_beam_2_slices = []
    
    N1 = step - n1
    touch_number = min(step, n2 - N1, n1, n2)
    summation = step - 1
    add = max(0, N1)
    
    i = [n for n in range(0, touch_number)]
    def interaction(i):
        show = i + add
        interaction_beam_1_slices.append(summation - show)
        interaction_beam_2_slices.append(show)
    list(map(interaction, i))
        
    return interaction_beam_1_slices, interaction_beam_2_slices

def X_Y_of_boost_frame_case_1(xb, yb, pxb, pyb, S):

    S_square = S*S

    sigma_x_square = Sigma_Calculation(xb)
    sigma_x_square = sigma_x_square**2
    
    product_x_px = Product_of_position_and_momentum(xb, pxb)
    x_px_product_avg = Sigma_Calculation(product_x_px)
    
    sigma_px_square = Sigma_Calculation(pxb)
    sigma_px_square = sigma_px_square**2

    sigma_x_kick = sqrt(abs(sigma_x_square + 2.0 * x_px_product_avg * S + sigma_px_square * S_square))

    sigma_y_square = Sigma_Calculation(yb)
    sigma_y_square = sigma_y_square**2
    
    product_y_py = Product_of_position_and_momentum(yb, pyb)
    y_py_product_avg = Sigma_Calculation(product_y_py)
    
    sigma_py_square = Sigma_Calculation(pyb)
    sigma_py_square = sigma_py_square**2

    sigma_y_kick = sqrt(abs(sigma_y_square + 2.0 * y_py_product_avg * S + sigma_py_square * S_square))

    beam_X_temp = np.array(xb)
    beam_PX_temp = np.array(pxb)
    X_temp = beam_X_temp + beam_PX_temp * S
    X = X_temp.tolist()

    beam_Y_temp = np.array(yb)
    beam_PY_temp = np.array(pyb)
    Y_temp = beam_Y_temp + beam_PY_temp * S
    Y = Y_temp.tolist()

    return X, Y, sigma_x_kick, sigma_y_kick, sigma_px_square, sigma_py_square

def X_Y_of_boost_frame_case_2(xb, yb, pxb, pyb, S):

    S_square = S**2

    x_square = square_value_in_a_list(xb)
    sigma_x_square = Average_of_list(x_square)
    
    x_px_product_avg = product_average(xb, pxb)
    
    px_square = square_value_in_a_list(pxb)
    sigma_px_square = Average_of_list(px_square)

    sigma_x_kick = sqrt(abs(sigma_x_square + 2.0 * x_px_product_avg * S + sigma_px_square * S_square))

    y_square = square_value_in_a_list(yb)
    sigma_y_square = Average_of_list(y_square)
    y_py_product_avg = product_average(yb, pyb)
    py_square = square_value_in_a_list(pyb)
    sigma_py_square = Average_of_list(py_square)

    sigma_y_kick = sqrt(abs(sigma_y_square + 2.0 * y_py_product_avg * S + sigma_py_square * S_square))

    beam_X_temp = np.array(xb)
    beam_PX_temp = np.array(pxb)
    X_temp = beam_X_temp + beam_PX_temp * S
    X = X_temp.tolist()

    beam_Y_temp = np.array(yb)
    beam_PY_temp = np.array(pyb)
    Y_temp = beam_Y_temp + beam_PY_temp * S
    Y = Y_temp.tolist()

    return X, Y, sigma_x_kick, sigma_y_kick, sigma_px_square, sigma_py_square

def Fx_and_Fy(x, y, sigma_x, sigma_y):
    if (sigma_x == sigma_y):
        commom_part = x*x + y*y
        if (commom_part == 0):
            return 0, 0
        elif (sigma_x == 0):
            return 0, 0
        else:
            commom_part = (1.0 - exp(- 0.5 * commom_part / (sigma_x*sigma_x))) / commom_part
            commom_part = sqrt(2.0 * pi) * commom_part
            Fx = commom_part * x
            Fy = commom_part * y
            return commom_part * x, commom_part * y
    else:
        Fy_plus_iFX_sign = 1
        if (sigma_x < sigma_y):
            x, y = y, x
            sigma_x, sigma_y = sigma_y, sigma_x
            Fy_plus_iFX_sign = -1
    
        Fy_sign = 1.0
        if (y < 0.0):
            y = -y
            Fy_sign = -1.0
        
        sigma_x2 = sigma_x*sigma_x
        sigma_y2 = sigma_y*sigma_y
        delta_sigma2 = sigma_x2 - sigma_y2
    
        cmplx_sqrt = cmath.sqrt(2 * delta_sigma2)
        sigma_frac = sigma_y / sigma_x
    
        z1 = complex(x, y) / cmplx_sqrt
        z2 = complex(z1.real * sigma_frac, z1.imag / sigma_frac)

        w1 = special.wofz(z1)
        w2 = special.wofz(z2)
        
        exponential = exp(-0.5 * (x*x / sigma_x2 + y*y / sigma_y2))
        
        Fy_plus_iFX = sqrt(4.0 * pi) / cmplx_sqrt * (w1 - exponential * w2)

        if (Fy_plus_iFX_sign == 1):
            return Fy_plus_iFX.imag, Fy_sign * Fy_plus_iFX.real
        else:
            return Fy_sign * Fy_plus_iFX.real, Fy_plus_iFX.imag

def Boost_frame_to_Lab_frame(x_bst, y_bst, z_bst, px_bst, py_bst, pz_bst, Beam_parameters, Crossing_angle, c_speed):

    p_Central         = Beam_parameters[1]
    Pass_Central_Time = Beam_parameters[7]

    x  = []
    y  = []
    xp = []
    yp = []
    p  = []
    t  = []

    cosine  = cos(Crossing_angle)
    sine    = sin(Crossing_angle)
    tangent = tan(Crossing_angle)

    cosine_sqr = cosine**2
    tan_sqr = tangent**2
    one_over_cos = 1.0 / cosine

    def Boost_to_Lab(x_bst, y_bst, z_bst, px_bst, py_bst, pz_bst):
        pz_i_temp = pz_bst + 1.0

        h_bst = pz_i_temp - sqrt(abs(pz_i_temp**2 - px_bst**2 - py_bst**2))
        h = h_bst * cosine_sqr

        sqrt_for_h = h - pz_i_temp

        hx_bst = 2.0 * px_bst / sqrt_for_h
        hy_bst = 2.0 * py_bst / sqrt_for_h
        hz_bst = 1.0 - 2.0 * pz_i_temp / sqrt_for_h

        hx_bst_sine = hx_bst * sine
        hy_bst_sine = hy_bst * sine
        hz_bst_sine = hz_bst * sine

        hx_bst_sine_1 = hx_bst * sine + 1.0
        tan_h = tangent * h
        tan_sqr_h = tan_h * tangent

        z = (z_bst - hx_bst_sine * x_bst / hx_bst_sine_1) / ( one_over_cos - hz_bst_sine * tangent / hx_bst_sine_1)
        x.append((x_bst - tangent * z) / hx_bst_sine_1)
        y.append(y_bst - hy_bst_sine)

        px = px_bst * cosine + tan_h
        py = py_bst * cosine
        pz = pz_bst + tangent * px + tan_sqr_h

        t.append(Pass_Central_Time - z / c_speed)                 # for orginal SDDS format

        xp.append(px)
        yp.append(py)
        p.append((pz + 1) * p_Central)                      # for orginal SDDS format

    list(map(Boost_to_Lab, x_bst, y_bst, z_bst, px_bst, py_bst, pz_bst))

    return x, xp, y, yp, t, p
    
def Save_to_SDDS_after_Kick(input_filename, data_after_kick, output_filename):

    """Save the SDDS file using the SDDS class."""

    SaveData = SDDS(0)
    SaveData.Load_Data(input_filename)
    
    SaveData.columnData[0][0] = data_after_kick[0]
    SaveData.columnData[1][0] = data_after_kick[1]
    SaveData.columnData[2][0] = data_after_kick[2]
    SaveData.columnData[3][0] = data_after_kick[3]
    SaveData.columnData[4][0] = data_after_kick[4]
    SaveData.columnData[5][0] = data_after_kick[5]
    
    ''' Examples '''
    # # x_position
    # print(SaveData.columnData[0][0][index= 0,1,2,...])
    # SaveData.columnData[0][0][index= 0,1,2,...] = 1.0
    # print(SaveData.columnData[0][0][index= 0,1,2,...])
    # print()
    
    # # y_position
    # print(SaveData.columnData[2][0][index= 0,1,2,...])
    # SaveData.columnData[2][0][index= 0,1,2,...] = 2.0
    # print(SaveData.columnData[2][0][index= 0,1,2,...])
    # print()
     
    # # x_p
    # print(SaveData.columnData[1][0][index= 0,1,2,...])
    # SaveData.columnData[1][0][index= 0,1,2,...] = 4.0
    # print(SaveData.columnData[1][0][index= 0,1,2,...])
    # print()
    
    # # y_p
    # print(SaveData.columnData[3][0][index= 0,1,2,...])
    # SaveData.columnData[3][0][index= 0,1,2,...] = 5.0
    # print(SaveData.columnData[3][0][index= 0,1,2,...])
    # print()

    # # p (total momentum) z_p = sqrt(p*p/(gamma_beta*gamma_bata) - x_p*x_p - y_p*y_p)
    # print(SaveData.columnData[5][0][index= 0,1,2,...])
    # SaveData.columnData[5][0][index= 0,1,2,...] = 6.0
    # print(SaveData.columnData[5][0][index= 0,1,2,...])
    # print()

    # # dt (z_position = dt x lightspeed)
    # print(SaveData.columnData[6][0][index= 0,1,2,...])
    # SaveData.columnData[6][0][index= 0,1,2,...] = 6.0
    # print(SaveData.columnData[6][0][index= 0,1,2,...])
    # print()
          
    SaveData.Save_Data(output_filename)
     
    del SaveData

def kick_beam_1(temp_particle_1_number, temp_particle_2_number, Particle_1_Radius, one_over_gamma, ratio_particle2, X_1, X_2_avg, Y_1, Y_2_avg, sigma_x_2_kick, sigma_y_2_kick, sigma_px_1_square, sigma_py_1_square, S, force_sign, Beam_1_Slice_ID, grouped_xyz_pxpypz_bst_1):

    number_radius_over_gamma = temp_particle_2_number * Particle_1_Radius * one_over_gamma * ratio_particle2
    
    sigmax_square = sigma_x_2_kick**2
    sigmay_square = sigma_y_2_kick**2
    Gxy = - 2.0 * (sigmax_square - sigmay_square)
    ratio_sigmay_sigmax = sigma_y_2_kick / sigma_x_2_kick if (sigma_x_2_kick) else 0
    
    for j in range(temp_particle_1_number):
        xy1_temp = Fx_and_Fy(X_1[j] - X_2_avg, Y_1[j] - Y_2_avg, sigma_x_2_kick, sigma_y_2_kick)
        fx = xy1_temp[0] * number_radius_over_gamma
        fy = xy1_temp[1] * number_radius_over_gamma
        
        if (Gxy == 0):
            gx = 0
            gy = 0
        elif (sigma_x_2_kick == 0):
            X_fx_Y_fy = X_1[j] * fx + Y_1[j] * fy - 2.0
            gx = X_fx_Y_fy / Gxy
            gy = gx
        else:
            math_exp = exp(-0.5 * ( X_1[j]**2 / sigmax_square + Y_1[j]**2 / sigmay_square))
            X_fx_Y_fy = X_1[j] * fx + Y_1[j] * fy - 2.0
            gx = (X_fx_Y_fy + 2.0 * ratio_sigmay_sigmax * math_exp) / Gxy
            gy = (X_fx_Y_fy + 2.0 * math_exp / ratio_sigmay_sigmax) / Gxy
    
        g = sigma_px_1_square * gx + sigma_py_1_square * gy
        g = g * number_radius_over_gamma * S
        
        f_X = force_sign * fx
        f_Y = force_sign * fy
        px = grouped_xyz_pxpypz_bst_1[3][Beam_1_Slice_ID][j]
        py = grouped_xyz_pxpypz_bst_1[4][Beam_1_Slice_ID][j]
        grouped_xyz_pxpypz_bst_1[0][Beam_1_Slice_ID][j] += S * f_X
        grouped_xyz_pxpypz_bst_1[1][Beam_1_Slice_ID][j] += S * f_Y
        grouped_xyz_pxpypz_bst_1[3][Beam_1_Slice_ID][j] += - f_X
        grouped_xyz_pxpypz_bst_1[4][Beam_1_Slice_ID][j] += - f_Y
        grouped_xyz_pxpypz_bst_1[5][Beam_1_Slice_ID][j] += - 0.5 * f_X * (px - 0.5 * f_X) - 0.5 * f_Y * (py - 0.5 * f_Y) - g

    return grouped_xyz_pxpypz_bst_1

def kick_beam_2(temp_particle_2_number, temp_particle_1_number, Particle_2_Radius, one_over_gamma, ratio_particle1, X_2, X_1_avg, Y_2, Y_1_avg, sigma_x_1_kick, sigma_y_1_kick, sigma_px_2_square, sigma_py_2_square, S, force_sign, Beam_2_Slice_ID, grouped_xyz_pxpypz_bst_2):
    
    number_radius_over_gamma = temp_particle_1_number * Particle_2_Radius * one_over_gamma * ratio_particle1
    
    sigmax_square = sigma_x_1_kick**2
    sigmay_square = sigma_y_1_kick**2
    Gxy = - 2.0 * (sigmax_square - sigmay_square)
    ratio_sigmay_sigmax = sigma_y_1_kick / sigma_x_1_kick if (sigma_x_1_kick) else 0
    
    for j in range(temp_particle_2_number):
        xy2_temp = Fx_and_Fy(X_2[j] - X_1_avg, Y_2[j] - Y_1_avg, sigma_x_1_kick, sigma_y_1_kick)
        fx = xy2_temp[0] * number_radius_over_gamma
        fy = xy2_temp[1] * number_radius_over_gamma
        
        if (Gxy == 0):
            gx = 0
            gy = 0
        elif (sigma_x_1_kick == 0):
            X_fx_Y_fy = X_2[j] * fx + Y_2[j] * fy - 2.0
            gx = X_fx_Y_fy / Gxy
            gy = gx
        else:
            math_exp = exp(-0.5 * ( X_2[j]**2 / sigmax_square + Y_2[j]**2 / sigmay_square))
            X_fx_Y_fy = X_2[j] * fx + Y_2[j] * fy - 2.0
            gx = (X_fx_Y_fy + 2.0 * ratio_sigmay_sigmax * math_exp) / Gxy
            gy = (X_fx_Y_fy + 2.0 * math_exp / ratio_sigmay_sigmax) / Gxy
    
        g = sigma_px_2_square * gx + sigma_py_2_square * gy
        g = g * number_radius_over_gamma * S
        
        f_X = fx
        f_Y = fy
        px = grouped_xyz_pxpypz_bst_2[3][Beam_2_Slice_ID][j]
        py = grouped_xyz_pxpypz_bst_2[4][Beam_2_Slice_ID][j]
        grouped_xyz_pxpypz_bst_2[0][Beam_2_Slice_ID][j] += S * f_X
        grouped_xyz_pxpypz_bst_2[1][Beam_2_Slice_ID][j] += S * f_Y
        grouped_xyz_pxpypz_bst_2[3][Beam_2_Slice_ID][j] += - f_X
        grouped_xyz_pxpypz_bst_2[4][Beam_2_Slice_ID][j] += - f_Y
        grouped_xyz_pxpypz_bst_2[5][Beam_2_Slice_ID][j] += - 0.5 * f_X * (px - 0.5 * f_X) - 0.5 * f_Y * (py - 0.5 * f_Y) - g

    return grouped_xyz_pxpypz_bst_2    

def Crab_kick_track_coordinates(CK_parameter, i_turn):

    Beam_name_1_in          = CK_parameter[0]
    Beam_name_1_out         = CK_parameter[1]
    Beam_name_2_in          = CK_parameter[2]
    Beam_name_2_out         = CK_parameter[3]
    number_of_Beam_1_Slices = CK_parameter[4]
    number_of_Beam_2_Slices = CK_parameter[5]
    Particle_1_Radius       = CK_parameter[6]
    Particle_2_Radius       = CK_parameter[7]
    Crossing_angle          = CK_parameter[8]       # half angle for our method
    Beam_1_Tilt             = CK_parameter[10]      # beam 2 tilt angle for our method
    Beam_2_Tilt             = CK_parameter[11]      # beam 2 tilt angle for our method
    c_speed                 = CK_parameter[12]
    case_flag               = CK_parameter[13]
    force_sign              = CK_parameter[15]
    Particle_1_number       = CK_parameter[20]
    Particle_2_number       = CK_parameter[21]
    recording_phase         = CK_parameter[22]
    phase_output            = CK_parameter[23]


    # SDDS_Data_Mode(Beam_name_1_in, Beam_name_2_in)
    ''' *********************************************************************** '''

    ''' Data from SDDS ******************************************************** '''
    # ***** beam 1 *****
    Beam_1_parameters = Parameters_from_SDDS(Beam_name_1_in)
    Beam_1_Data       = np.array(Data_from_SDDS(Beam_name_1_in))
    # ***** beam 2 *****
    Beam_2_parameters = Parameters_from_SDDS(Beam_name_2_in)
    Beam_2_Data       = np.array(Data_from_SDDS(Beam_name_2_in))
    ''' *********************************************************************** '''


    # Input User's own initial Distribution: x, x', y, y', z, p *********************************************************************************
    check_initial_phase = 1
    if (i_turn == 0):
        if (os.path.isfile('Files_of_Input/CASA_Phase_Initial_Beam_1.casa') and os.path.isfile('Files_of_Input/CASA_Phase_Initial_Beam_2.casa')):
            check_initial_phase = os.system('python3 CrabKick_initial_phase_input_GUI.py')

        if (check_initial_phase == 0):
            CASA_Beam_1_data = np.loadtxt('Files_of_Input/CASA_Phase_Initial_Beam_1.casa').T
            CASA_Beam_2_data = np.loadtxt('Files_of_Input/CASA_Phase_Initial_Beam_2.casa').T

            Beam_1_Data[0] = CASA_Beam_1_data[0]
            Beam_1_Data[1] = CASA_Beam_1_data[1]
            Beam_1_Data[2] = CASA_Beam_1_data[2]
            Beam_1_Data[3] = CASA_Beam_1_data[3]
            Beam_1_Data[4] = CASA_Beam_1_data[4]
            Beam_1_Data[5] = CASA_Beam_1_data[5]

            Beam_2_Data[0] = CASA_Beam_2_data[0]
            Beam_2_Data[1] = CASA_Beam_2_data[1]
            Beam_2_Data[2] = CASA_Beam_2_data[2]
            Beam_2_Data[3] = CASA_Beam_2_data[3]
            Beam_2_Data[4] = CASA_Beam_2_data[4]
            Beam_2_Data[5] = CASA_Beam_2_data[5]
            
            del CASA_Beam_1_data
            del CASA_Beam_2_data

    # *********************************************************************************************************************************************

    g1 = sqrt(1.0 + Beam_1_parameters[1]**2.0)
    g2 = sqrt(1.0 + Beam_2_parameters[1]**2.0)

    number_of_kicking = 1
    ratio_particle1 = Particle_1_number / Beam_1_parameters[3]
    ratio_particle2 = Particle_2_number / Beam_2_parameters[3]
    # one_over_gamma   = 1.0 / Gamma(g1, g2, Crossing_angle)
    one_over_gamma_1   = 1.0 / g1
    one_over_gamma_2   = 1.0 / g2


    ''' Step 1: Lab frame to Boost frame ************************************** '''
    print('1: Lab frame to Boost frame')

    # for Beam 1 ******************************************************************
    xyz_pxpypz_ID_bst_1 = Lab_frame_to_Boost_frame(Beam_1_Data, Beam_1_parameters, Crossing_angle, Beam_1_Tilt, c_speed)
    print('   Beam 1: -- Completed.')
    grouped_xyz_pxpypz_bst_1 = Group_Beam_Slice_in_Boost_frame(xyz_pxpypz_ID_bst_1, number_of_Beam_1_Slices)
    number_of_Beam_1_Slices = int(grouped_xyz_pxpypz_bst_1[7])

    # for Beam 2 ******************************************************************
    xyz_pxpypz_ID_bst_2 = Lab_frame_to_Boost_frame(Beam_2_Data, Beam_2_parameters, Crossing_angle, Beam_2_Tilt, c_speed)
    print('   Beam 2: -- completed.\n')

    grouped_xyz_pxpypz_bst_2 = Group_Beam_Slice_in_Boost_frame(xyz_pxpypz_ID_bst_2, number_of_Beam_2_Slices)
    number_of_Beam_2_Slices = int(grouped_xyz_pxpypz_bst_2[7])

    ''' *********************************************************************** '''

    ''' Step 2: Kick in Boost frame ******************************************* '''
    print('2: Kick in Boost frame')

    # Collision_S for Beam 1 ******************************************************************
    S_of_beam_1_slices = S_for_collision(grouped_xyz_pxpypz_bst_1[2], number_of_Beam_1_Slices)
    # S_of_beam_1_slices = Slice_center(grouped_xyz_pxpypz_bst_1[2], number_of_Beam_1_Slices)

    # Collision_S for Beam 2 ******************************************************************
    S_of_beam_2_slices = S_for_collision(grouped_xyz_pxpypz_bst_2[2], number_of_Beam_2_Slices)
    # S_of_beam_2_slices = Slice_center(grouped_xyz_pxpypz_bst_2[2], number_of_Beam_2_Slices)

    dirs = 'Files_of_Output/Phase_coordinates/' + 'turn_' + str(i_turn) + '_in_Boost_Frame/'
    if not os.path.exists(dirs):
        os.makedirs(dirs)


    Total_interaction_steps = number_of_Beam_1_Slices + number_of_Beam_2_Slices - 1
    for interaction in range(Total_interaction_steps):
        Step_of_interaction = interaction + 1
        series_number = Slice_series_when_interaction(number_of_Beam_1_Slices, number_of_Beam_2_Slices, Step_of_interaction)
        for i in range(len(series_number[0])):
            Beam_1_Slice_ID = series_number[0][i]
            Beam_2_Slice_ID = series_number[1][i]

            S1 = S_of_beam_1_slices[Beam_1_Slice_ID]
            S2 = S_of_beam_2_slices[Beam_2_Slice_ID]

            # Calculate X_1, Y_1, ..... *************************************************
            xb  = grouped_xyz_pxpypz_bst_1[0][Beam_1_Slice_ID]
            yb  = grouped_xyz_pxpypz_bst_1[1][Beam_1_Slice_ID]
            pxb = grouped_xyz_pxpypz_bst_1[3][Beam_1_Slice_ID]
            pyb = grouped_xyz_pxpypz_bst_1[4][Beam_1_Slice_ID]
            if (case_flag == 1):
                Temp = X_Y_of_boost_frame_case_1(xb, yb, pxb, pyb, S2)
            if (case_flag == 2):
                Temp = X_Y_of_boost_frame_case_2(xb, yb, pxb, pyb, S2)

            X_1 = Temp[0]
            Y_1 = Temp[1]
            sigma_x_1_kick = Temp[2]
            sigma_y_1_kick = Temp[3]
            sigma_px_1_square = Temp[4]
            sigma_py_1_square = Temp[5]
            X_1_avg = Average_of_list(X_1)
            Y_1_avg = Average_of_list(Y_1)
            del Temp
            # ***************************************************************************

            # Calculate X_2, Y_2, ..... *************************************************
            xb  = grouped_xyz_pxpypz_bst_2[0][Beam_2_Slice_ID]
            yb  = grouped_xyz_pxpypz_bst_2[1][Beam_2_Slice_ID]
            pxb = grouped_xyz_pxpypz_bst_2[3][Beam_2_Slice_ID]
            pyb = grouped_xyz_pxpypz_bst_2[4][Beam_2_Slice_ID]
            if (case_flag == 1):
                Temp = X_Y_of_boost_frame_case_1(xb, yb, pxb, pyb, S1)
            if (case_flag == 2):
                Temp = X_Y_of_boost_frame_case_2(xb, yb, pxb, pyb, S1)

            X_2 = Temp[0]
            Y_2 = Temp[1]
            sigma_x_2_kick = Temp[2]
            sigma_y_2_kick = Temp[3]
            sigma_px_2_square = Temp[4]
            sigma_py_2_square = Temp[5]
            X_2_avg = Average_of_list(X_2)
            Y_2_avg = Average_of_list(Y_2)
            del Temp
            # ***************************************************************************

            temp_particle_1_number = len(X_1)
            temp_particle_2_number = len(X_2)

            # **** Starting Kicking Beam 1 ******************************************************************
            grouped_xyz_pxpypz_bst_1 = kick_beam_1(temp_particle_1_number, temp_particle_2_number, Particle_1_Radius, one_over_gamma_1, ratio_particle2, X_1, X_2_avg, Y_1, Y_2_avg, sigma_x_2_kick, sigma_y_2_kick, sigma_px_1_square, sigma_py_1_square, S2, force_sign, Beam_1_Slice_ID, grouped_xyz_pxpypz_bst_1)
            # **** Ending Kicking Beam 1 ******************************************************************** 

            # **** Starting Kicking Beam 2*******************************************************************
            grouped_xyz_pxpypz_bst_2 = kick_beam_2(temp_particle_2_number, temp_particle_1_number, Particle_2_Radius, one_over_gamma_2, ratio_particle1, X_2, X_1_avg, Y_2, Y_1_avg, sigma_x_1_kick, sigma_y_1_kick, sigma_px_2_square, sigma_py_2_square, S1, force_sign, Beam_2_Slice_ID, grouped_xyz_pxpypz_bst_2)
            # **** Ending Kicking Beam 2 ******************************************************************** 

        print('   Step of interaction : %4d' % Step_of_interaction + '  -- Completed.')


        # ********** Record all detail changes of the phase-coordinates in the collision *************************************

        if (recording_phase == 'yes'):
            tx1  = np.array(grouped_xyz_pxpypz_bst_1[0][0])
            ty1  = np.array(grouped_xyz_pxpypz_bst_1[1][0])
            tz1  = np.array(grouped_xyz_pxpypz_bst_1[2][0])
            txp1 = np.array(grouped_xyz_pxpypz_bst_1[3][0])
            typ1 = np.array(grouped_xyz_pxpypz_bst_1[4][0])
            tzp1 = np.array(grouped_xyz_pxpypz_bst_1[5][0])
            for i_iD in range(len(grouped_xyz_pxpypz_bst_1[0]) - 1):
                tx1  = np.append(tx1, grouped_xyz_pxpypz_bst_1[0][i_iD+1])
                ty1  = np.append(ty1, grouped_xyz_pxpypz_bst_1[1][i_iD+1])
                tz1  = np.append(tz1, grouped_xyz_pxpypz_bst_1[2][i_iD+1])
                txp1 = np.append(txp1, grouped_xyz_pxpypz_bst_1[3][i_iD+1])
                typ1 = np.append(typ1, grouped_xyz_pxpypz_bst_1[4][i_iD+1])
                tzp1 = np.append(tzp1, grouped_xyz_pxpypz_bst_1[5][i_iD+1])

            tx2  = np.array(grouped_xyz_pxpypz_bst_2[0][0])
            ty2  = np.array(grouped_xyz_pxpypz_bst_2[1][0])
            tz2  = np.array(grouped_xyz_pxpypz_bst_2[2][0])
            txp2 = np.array(grouped_xyz_pxpypz_bst_2[3][0])
            typ2 = np.array(grouped_xyz_pxpypz_bst_2[4][0])
            tzp2 = np.array(grouped_xyz_pxpypz_bst_2[5][0])
            for i_iD in range(len(grouped_xyz_pxpypz_bst_2[0]) - 1):
                tx2  = np.append(tx2, grouped_xyz_pxpypz_bst_2[0][i_iD+1])
                ty2  = np.append(ty2, grouped_xyz_pxpypz_bst_2[1][i_iD+1])
                tz2  = np.append(tz2, grouped_xyz_pxpypz_bst_2[2][i_iD+1])
                txp2 = np.append(txp2, grouped_xyz_pxpypz_bst_2[3][i_iD+1])
                typ2 = np.append(typ2, grouped_xyz_pxpypz_bst_2[4][i_iD+1])
                tzp2 = np.append(tzp2, grouped_xyz_pxpypz_bst_2[5][i_iD+1])


            [tx1, txp1, ty1, typ1, tz, tzp] = Boost_frame_to_Lab_frame(tx1, ty1, tz1, txp1, typ1, tzp1, Beam_1_parameters, Crossing_angle, c_speed)
            tz1  = tz  if (phase_output == 1) else tz1    # "phase_output == 1" outputs "t" otherwise outputs "z" 
            tzp1 = tzp if (phase_output == 1) else tzp1   # "phase_output == 1" outputs "p" otherwise outputs "dp/p0"

            [tx2, txp2, ty2, typ2, tz, tzp] = Boost_frame_to_Lab_frame(tx2, ty2, tz2, txp2, typ2, tzp2, Beam_2_parameters, Crossing_angle, c_speed)
            tz2  = tz  if (phase_output == 1) else tz2    # "phase_output == 1" outputs "t" otherwise outputs "z" 
            tzp2 = tzp if (phase_output == 1) else tzp2   # "phase_output == 1" outputs "p" otherwise outputs "dp/p0"

            temp_name_1 = dirs + str(interaction) + '_1.casa'
            temp_name_2 = dirs + str(interaction) + '_2.casa'

            tt_1 = np.array([tx1, txp1, ty1, typ1, tz1, tzp1])
            tt_2 = np.array([tx2, txp2, ty2, typ2, tz2, tzp2])

            np.savetxt(temp_name_1, tt_1.T, fmt='% .12e', delimiter='  ')
            np.savetxt(temp_name_2, tt_2.T, fmt='% .12e', delimiter='  ')
        # ****************************************************************************************************

        del series_number  # very import !!!!


    ''' *********************************************************************** '''


    ''' Step 3: Boost frame to Lab frame ************************************** '''
    print('')
    print('3: Boost frame to Lab frame')

    for j in range(number_of_Beam_1_Slices):
        number_inside_slice = len(grouped_xyz_pxpypz_bst_1[0][j])
        # new x,y,z-position and new px,py,pz in Boost frame (in the same sequence of orginal data)
        for k in range(number_inside_slice):
            x  = grouped_xyz_pxpypz_bst_1[0][j][k]
            y  = grouped_xyz_pxpypz_bst_1[1][j][k]
            z  = grouped_xyz_pxpypz_bst_1[2][j][k]
            px = grouped_xyz_pxpypz_bst_1[3][j][k]
            py = grouped_xyz_pxpypz_bst_1[4][j][k]
            pz = grouped_xyz_pxpypz_bst_1[5][j][k]
            ID = grouped_xyz_pxpypz_bst_1[6][j][k]

            xyz_pxpypz_ID_bst_1[0][ID] = x
            xyz_pxpypz_ID_bst_1[1][ID] = y
            xyz_pxpypz_ID_bst_1[2][ID] = z
            xyz_pxpypz_ID_bst_1[3][ID] = px
            xyz_pxpypz_ID_bst_1[4][ID] = py
            xyz_pxpypz_ID_bst_1[5][ID] = pz

    for j in range(number_of_Beam_2_Slices):
        number_inside_slice = len(grouped_xyz_pxpypz_bst_2[0][j])
        # new x,y,z-position and new px,py,pz in Boost frame (in the same sequence of orginal data)
        for k in range(number_inside_slice):
            x  = grouped_xyz_pxpypz_bst_2[0][j][k]
            y  = grouped_xyz_pxpypz_bst_2[1][j][k]
            z  = grouped_xyz_pxpypz_bst_2[2][j][k]
            px = grouped_xyz_pxpypz_bst_2[3][j][k]
            py = grouped_xyz_pxpypz_bst_2[4][j][k]
            pz = grouped_xyz_pxpypz_bst_2[5][j][k]
            ID = grouped_xyz_pxpypz_bst_2[6][j][k]

            xyz_pxpypz_ID_bst_2[0][ID] = x
            xyz_pxpypz_ID_bst_2[1][ID] = y
            xyz_pxpypz_ID_bst_2[2][ID] = z
            xyz_pxpypz_ID_bst_2[3][ID] = px
            xyz_pxpypz_ID_bst_2[4][ID] = py
            xyz_pxpypz_ID_bst_2[5][ID] = pz

    Beam_1_data_after_kick = Boost_frame_to_Lab_frame(xyz_pxpypz_ID_bst_1[0], xyz_pxpypz_ID_bst_1[1], xyz_pxpypz_ID_bst_1[2], xyz_pxpypz_ID_bst_1[3], xyz_pxpypz_ID_bst_1[4], xyz_pxpypz_ID_bst_1[5], Beam_1_parameters, Crossing_angle, c_speed)

    print('   Beam 1: -- Completed.')

    Beam_2_data_after_kick = Boost_frame_to_Lab_frame(xyz_pxpypz_ID_bst_2[0], xyz_pxpypz_ID_bst_2[1], xyz_pxpypz_ID_bst_2[2], xyz_pxpypz_ID_bst_2[3], xyz_pxpypz_ID_bst_2[4], xyz_pxpypz_ID_bst_2[5], Beam_2_parameters, Crossing_angle, c_speed)

    print('   Beam 2: -- Completed.\n')
    ''' *********************************************************************** '''


    ''' Step 4: Save Data ***************************************************** '''
    print('4: Save Data (SDDS)')

    mode = Data_mode(Beam_name_1_in)
    if (mode == 0):
        print('   Beam-1 data storage mode is ASCII')
    if (mode == 1):
        print('   Beam-1 data storage mode is BINARY')
    mode = Data_mode(Beam_name_2_in)
    if (mode == 0):
        print('   Beam-2 data storage mode is ASCII')
    if (mode == 1):
        print('   Beam-2 data storage mode is BINARY')

    print('   Data have been saved.')

    Save_to_SDDS_after_Kick(Beam_name_1_in, Beam_1_data_after_kick, Beam_name_1_out)
    Save_to_SDDS_after_Kick(Beam_name_2_in, Beam_2_data_after_kick, Beam_name_2_out)
    ''' *********************************************************************** '''


    Beam_1_final = np.array(Beam_1_data_after_kick)
    Beam_2_final = np.array(Beam_2_data_after_kick)
    Beam_1_final[4] = Beam_1_final[4] * c_speed
    Beam_2_final[4] = Beam_2_final[4] * c_speed
    np.savetxt('Files_of_Output/Phase_coordinates/' + str(i_turn) + '-beam-1-final.casa', Beam_1_final.T, fmt='% .12e', delimiter='  ')
    np.savetxt('Files_of_Output/Phase_coordinates/' + str(i_turn) + '-beam-2-final.casa', Beam_2_final.T, fmt='% .12e', delimiter='  ')



    sigma_beam1_x = Sigma_Calculation(Beam_1_data_after_kick[0])
    sigma_beam1_y = Sigma_Calculation(Beam_1_data_after_kick[2])
    sigma_beam1_z = Sigma_Calculation(Beam_1_data_after_kick[4]) * c_speed
    sigma_beam2_x = Sigma_Calculation(Beam_2_data_after_kick[0])
    sigma_beam2_y = Sigma_Calculation(Beam_2_data_after_kick[2])
    sigma_beam2_z = Sigma_Calculation(Beam_2_data_after_kick[4]) * c_speed



    Sigma_beam1 = np.array([sigma_beam1_x, sigma_beam1_y, sigma_beam1_z])
    Sigma_beam2 = np.array([sigma_beam2_x, sigma_beam2_y, sigma_beam2_z])

    x1_avg  = Average_of_list(Beam_1_data_after_kick[0])
    xp1_avg = Average_of_list(Beam_1_data_after_kick[1])
    y1_avg  = Average_of_list(Beam_1_data_after_kick[2])
    yp1_avg = Average_of_list(Beam_1_data_after_kick[3])
    beam1_phase_avg = np.array([x1_avg, xp1_avg, y1_avg, yp1_avg])

    x2_avg  = Average_of_list(Beam_2_data_after_kick[0])
    xp2_avg = Average_of_list(Beam_2_data_after_kick[1])
    y2_avg  = Average_of_list(Beam_2_data_after_kick[2])
    yp2_avg = Average_of_list(Beam_2_data_after_kick[3])
    beam2_phase_avg = np.array([x2_avg, xp2_avg, y2_avg, yp2_avg])

    gamma_beta_after_turn = np.array([Beam_1_parameters[1], Beam_2_parameters[1]])

    return Sigma_beam1, Sigma_beam2, gamma_beta_after_turn, beam1_phase_avg, beam2_phase_avg

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


