#  This is the fuction code used for dealing with Twiss-file
# this is a tentative method (just works for Vasiliy's sample-file).
#  River and Vasiliy, 2018

import cmath
import numpy as np
import math

def TwissParameter_nu_Gamma(Filename):
    count = len(open(Filename,'r').readlines()) # total number of lines of the text-file
    file = open(Filename, "r")
    holder = file.readlines()
    file.close()
    
    for i in range (count):
        temp = holder[i]
        temp = temp.split()
        if (temp):
            if (temp[1] == 'Q2'):
                nu_y = temp[3]
            elif (temp[1] == 'MASS'):
                mass = temp[3]

    return float(nu_y), float(mass)


def ParameterForSRF(Filename):
    count = len(open(Filename,'r').readlines()) # total number of lines of the text-file

    file = open(Filename, "r")
    holder = file.readlines()
    file.close()
    
    parameter_in  = []
    for i in range (count):
        temp = holder[i]
        temp = temp.split()
        if (temp):
            if ((temp[0] == '@') and (temp[1] != 'Q2') and (temp[1] != 'GAMMA')) or (temp[0] == '*') or (temp[0] == '$'):
                pass
            elif (temp[1] == 'Q2'):
                nu_y = temp[3]
            elif (temp[1] == 'GAMMA'):
                Gamma = temp[3]
            elif (temp):
                parameter_in.append(temp)
    
    zeta    = []
    S       = []
    length  = []
    angle   = []
    beta_y  = []
    alpha_y = []
    mu_y    = []
    K1_L    = []

    init_length = len(parameter_in)
    flag = 1.0

    for i in range (init_length):
        str1 = 'SNAKE1'
        str2 = 'SNAKE2'
        str3 = '"SNAKE1"'
        str4 = '"SNAKE2"'
        if (parameter_in[i][0] == str1) or (parameter_in[i][0] == str2 or parameter_in[i][0] == str3) or (parameter_in[i][0] == str4) :
            flag = - flag
        zeta.append(flag)
        S.append(parameter_in[i][2])
        length.append(parameter_in[i][3])
        angle.append(parameter_in[i][4])
        beta_y.append(float(parameter_in[i][5]))
        alpha_y.append(float(parameter_in[i][6]))
        mu_y.append(float(parameter_in[i][7]))
        K1_L.append(float(parameter_in[i][8]))

    del parameter_in

    return zeta, S, length, angle, beta_y, alpha_y, mu_y, K1_L
    

def Parameter_W_ForSRF(TwiParameters):
    length = len(TwiParameters[0])
    
    Parameter_Write = 'Files_of_Output/parameters_for_calculation.txt'
    with open(Parameter_Write, "w") as f:
        for i in range (length):
            f.write('%2.2f    %5.10f    %5.10f    %5.10f    %5.10f    %5.10f    %5.10f \n' % (float(TwiParameters[0][i]), float(TwiParameters[1][i]), float(TwiParameters[2][i]), float(TwiParameters[3][i]), float(TwiParameters[4][i]), float(TwiParameters[5][i]), float(TwiParameters[6][i])))
            # f.write('%3s  %16s  %16s  %16s  %16s  %16s  %16s \n' % (TwiParameters[0][i], TwiParameters[1][i], TwiParameters[2][i], TwiParameters[3][i], TwiParameters[4][i], TwiParameters[5][i], TwiParameters[6][i]))

