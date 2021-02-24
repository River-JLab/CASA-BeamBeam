#  This is the fuction code used for Spin Response Function
#  Vasiliy and River, 2019

import os
import sys, time, math
import cmath
import numpy as np
from scipy import *
from TwissParameter_Linux import *

def complex_array(real_array):
    number    = len(real_array)
    compx_arr = []
    for i in range(number):
        compx_arr.append(complex(real_array[i], 1.0))
    return np.array(compx_arr)


def complex_exp(theta):
    return complex(math.cos(theta), math.sin(theta))


def complex_exp_array(variablb_array):
    number =  len(variablb_array)
    complex_array = [0] * number
    for i in range(number):
        complex_array[i] = complex(math.cos(variablb_array[i]), math.sin(variablb_array[i]))

    return np.array(complex_array)

    
def fy(beta_y_array, mu_y_array):
    a = np.sqrt(beta_y_array)
    b = complex_exp_array(mu_y_array)
    c = a * b
    return np.array(c)


def fy_conjugate(fy_array):
    a = np.conjugate(fy_array)
    return np.array(a)


def fy_derivative(beta_y_array, mu_y_array, alpha_y_array):
    mu_y    = np.array(complex_exp_array(mu_y_array))
    alpha_y = np.array(complex_array(-alpha_y_array))
    beta_y  = np.sqrt(beta_y_array)
    return np.array(mu_y * alpha_y / beta_y)


def fy_conjugate_derivative(beta_y_array, mu_y_array, alpha_y_array):
    mu_y    = np.array(complex_exp_array(-mu_y_array))
    alpha_y = np.array(complex_array(alpha_y_array))
    beta_y  = np.sqrt(beta_y_array)
    return np.array(-mu_y * alpha_y / beta_y)


def psi_z(gamma, B_moment_anomaly, zeta, angle):
    number = len(zeta)
    Psi = np.zeros(number)
    
    Psi[0] = gamma * B_moment_anomaly * zeta[0] * angle[0]
    for i in range(number - 1):
        Psi[i+1] = Psi[i] + gamma * B_moment_anomaly * zeta[i+1] * angle[i+1]
    return np.array(Psi)


def complex_exp_psi(Psi_array):
	complx_exp_psi = complex_exp_array(Psi_array) 
	return np.array(complx_exp_psi)


def derivative_psi(Gamma, B_moment_anomaly, Kappa, zeta):
	d_psi = Kappa * zeta
	d_psi = array_times_number(Gamma * B_moment_anomaly, d_psi)
	return np.array(d_psi)


def Phi_z(alpha, complx_exp_psi, d_psi, zeta):
	zeta = 1.0 - zeta
	zeta = array_times_number(alpha, zeta)
	zeta = complex_exp_array(zeta)
	Phi =  zeta * complx_exp_psi * d_psi
	Phi = array_times_number(complex(0,1), Phi)
	
	return np.array(Phi)


def array_times_number(number, array):
    return np.array([i * number for i in array])


def parameter_input():
    Filename="Files_of_Input/SRF_parameters.casa"

    def len_of_file():
        with open(Filename) as file_of_para:
            return len(list(file_of_para))
        file_of_para.close()
    number_of_para = len_of_file()

    File=open(Filename,"r")
    holder=File.readlines()
    File.close()

    temp = holder[0]
    temp = temp.split()
    twiss_filename = temp[0]

    temp = holder[1]
    temp = temp.split()
    SRF_result_filename = temp[0]

    temp = holder[2]
    temp = temp.split()
    SRS_result_filename = temp[0]

    temp = holder[4]
    temp = temp.split()
    Starting_enargy = temp[0]
    Ending_enargy   = temp[1]

    temp = holder[5]
    temp = temp.split()
    Energy_time_steps = temp[0]

    temp = holder[6]
    temp = temp.split()
    Anomalous_magnetic_moment = temp[0]

    temp = holder[7]
    temp = temp.split()
    delta_phi = temp[0]

    temp = holder[8]
    temp = temp.split()
    delta_y = temp[0]

    temp = holder[9]
    temp = temp.split()
    alpha_snake = temp[0]

    return twiss_filename, SRF_result_filename, SRS_result_filename, Starting_enargy, Ending_enargy, Energy_time_steps, Anomalous_magnetic_moment, delta_phi, delta_y, alpha_snake


# Spin Response Function
def F_integral(twiss_filename, energy, alpha_snake, Anomalous_magnetic_moment, delta_phi, delta_y):
    nu_mass = TwissParameter_nu_Gamma(twiss_filename)
    ParametersForSRF = ParameterForSRF(twiss_filename)

    # Parameter_W_ForSRF(ParametersForSRF)
    
    zeta    = np.array([float(i) for i in ParametersForSRF[0]])
    S       = np.array([float(i) for i in ParametersForSRF[1]])
    length  = np.array([float(i) for i in ParametersForSRF[2]])
    angle   = np.array([float(i) for i in ParametersForSRF[3]])

    beta_y  = np.array([float(i) for i in ParametersForSRF[4]])
    alpha_y = np.array([float(i) for i in ParametersForSRF[5]])
    mu_y    = np.array([float(i) for i in ParametersForSRF[6]]) 
    mu_y_temp_num = 2.0 * math.pi
    mu_y    = array_times_number(mu_y_temp_num, mu_y)
    
    K1_L    = np.array([float(i) for i in ParametersForSRF[7]]) 
    
    nu_y = float(nu_mass[0])
    mass = float(nu_mass[1])

    gamma = float(energy) / mass
    
    alpha_snake = float(alpha_snake)
    B_moment_anomaly = float(Anomalous_magnetic_moment)
    
    delta_phi = float(delta_phi)
    delta_y   = float(delta_y)

    exp_compx_nu1 = complex_exp(2.0 * math.pi * nu_y)
    exp_compx_nu1 = exp_compx_nu1 - 1.0
    exp_compx_nu1 = 1.0 / exp_compx_nu1
    
    exp_compx_nu2 = complex_exp(-2.0 * math.pi * nu_y)
    exp_compx_nu2 = exp_compx_nu2 - 1.0
    exp_compx_nu2 = 1.0 / exp_compx_nu2
    
    fy_z = fy(beta_y, mu_y)
    fy_z_conj = fy_conjugate(fy_z)

    d_fy_z = fy_derivative(beta_y, mu_y, alpha_y)
    d_fy_z_conj = fy_conjugate(d_fy_z)
    
    len_dfy = len(d_fy_z)
    d_fy_z_avg = np.zeros(len_dfy, dtype=complex)
    d_fy_z_avg[0] = d_fy_z[0]
    for i in range(len_dfy-1): 
        d_fy_z_avg[i+1] = 0.5 * (d_fy_z[i] + d_fy_z[i+1])

    d_fy_z_conj_avg = fy_conjugate(d_fy_z_avg)

    
    psi = psi_z(gamma, B_moment_anomaly, zeta, angle)
    complx_exp_psi = complex_exp_psi(psi)

    
    G1j = np.zeros(len_dfy, dtype=complex)
    # G1j[0] = d_fy_z_avg[0] * cos(alpha_snake * ( 1 - zeta[0])) * (complx_exp_psi[0] - complx_exp_psi[0])
    for i in range(len_dfy-1): 
        G1j[i+1] = d_fy_z_avg[i+1] * cos(alpha_snake * ( 1 - zeta[i+1])) * (complx_exp_psi[i+1] - complx_exp_psi[i])

    SUM_G1j = np.zeros(len_dfy, dtype=complex)
    SUM_G1j[0] = G1j[0]
    for i in range(len_dfy-1): 
        SUM_G1j[i+1] = G1j[i+1] + SUM_G1j[i]

    G2j = np.zeros(len_dfy, dtype=complex)
    # G2j[0] = d_fy_z_conj_avg[0] * cos(alpha_snake * ( 1 - zeta[0])) * (complx_exp_psi[0] - complx_exp_psi[0])
    for i in range(len_dfy-1): 
        G2j[i+1] = d_fy_z_conj_avg[i+1] * cos(alpha_snake * ( 1 - zeta[i+1])) * (complx_exp_psi[i+1] - complx_exp_psi[i])

    SUM_G2j = np.zeros(len_dfy, dtype=complex)
    SUM_G2j[0] = G2j[0]
    for i in range(len_dfy-1): 
        SUM_G2j[i+1] = G2j[i+1] + SUM_G2j[i]


    F1j_const = exp_compx_nu1 * SUM_G1j[len_dfy-1]
    F2j_const = exp_compx_nu2 * SUM_G2j[len_dfy-1]

    F1j = np.zeros(len_dfy, dtype=complex)
    for i in range(len_dfy): 
        F1j[i] = SUM_G1j[i] + F1j_const

    F2j = np.zeros(len_dfy, dtype=complex)
    for i in range(len_dfy): 
        F2j[i] = SUM_G2j[i] + F2j_const

    F3 = np.zeros(len_dfy, dtype=complex)
    for i in range(len_dfy): 
        F3[i] = gamma * B_moment_anomaly / complex(0.0, 2.0) * ( fy_z_conj[i] * F1j[i] - fy_z[i] * F2j[i])


    F3 = np.absolute(F3)

    length_F = len(F3)

    For_dipoles = np.zeros(len_dfy)
    for i in range(len_dfy): 
        For_dipoles[i] = (angle[i] * delta_phi * 0.001)**2.0
    
    For_quadrupoles = np.zeros(len_dfy)
    for i in range(len_dfy): 
        For_quadrupoles[i] = (K1_L[i] * delta_y * 0.001)**2.0
   
    SUM_D_and_Q = np.zeros(len_dfy)
    for i in range(len_dfy): 
        SUM_D_and_Q[i] = For_dipoles[i] + For_quadrupoles[i]

    element_of_strength = np.zeros(len_dfy)
    for i in range(len_dfy): 
        element_of_strength[i] = F3[i]**2 * SUM_D_and_Q[i]

    Zero_integer_SRS = np.sum(element_of_strength)
    Zero_integer_SRS = 0.5 * math.sqrt(Zero_integer_SRS) / math.pi

    return S, F3, Zero_integer_SRS
    

