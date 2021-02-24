# ********************************************************* #
#                                                           #
#  This is a Test Inteface for JLAB-CSSA Beam-Beam code.    #
#                                                           #
#  Vasiliy Morozov and River Huang                          #
#  2018-2020                                                #
#                                                           #
# ********************************************************* #


import os, math, sys
from math import cos, sin 
from numpy import *
import random
# from numba import *

def input_parameters():
    Filename="beam_parameters/Beam_parameters.casa"

    def len_of_file():
        with open(Filename) as file_of_para:
            return len(list(file_of_para))
        file_of_para.close()
    number_of_para = len_of_file()
    
    # print('number_of_para', number_of_para)

    File=open(Filename,"r")
    holder=File.readlines()
    File.close()

    parameter1 = []
    parameter2 = []

    temp = holder[0]
    temp = temp.split()
    parameter1.append(temp[0])
    parameter2.append(temp[1])

    temp = holder[1]  # dish line ---------------------- "
    
    for i in range (8-2):
        temp = holder[i+2]
        temp = temp.split()
        parameter1.append(temp[0])
        parameter2.append(temp[1])
    
    for i in range (number_of_para-8):
        temp = holder[i+8]
        temp = temp.split()
        if temp:
            parameter1.append(float(temp[0]))
            parameter2.append(float(temp[1]))
        else:
            # print('i = ', i, temp)
            break
        
    return parameter1, parameter2

def generate_mtr(parameters):

    filename = "Elegant_Process/mtr_file/" + str(parameters[1]) + '.mtr'
    
    betax = (parameters[7])
    betay = (parameters[8])
    alphx = (parameters[9])
    alphy = (parameters[10])
    nux   = (parameters[11])
    nuy   = (parameters[12])
    nus   = (parameters[13])
    sigmz = (parameters[14])
    sigmp = (parameters[15])
    
    r11=cos(2*pi*nux)+alphx*sin(2*pi*nux)
    r12=betax*sin(2*pi*nux)
    r21=-(1+alphx**2)*sin(2*pi*nux)/betax
    r22=cos(2*pi*nux)-alphx*sin(2*pi*nux)
    
    r33=cos(2*pi*nuy)+alphy*sin(2*pi*nuy)
    r34=betay*sin(2*pi*nuy)
    r43=-(1+alphy**2)*sin(2*pi*nuy)/betay
    r44=cos(2*pi*nuy)-alphy*sin(2*pi*nuy)
    
    r55=cos(2*pi*nus)
    r56=sigmz*sin(2*pi*nus)/sigmp
    r65=-sigmp*sin(2*pi*nus)/sigmz
    r66=cos(2*pi*nus)
    
    file=open(filename,"w")
    
    # file.write("ICR: 0.0 0.0 0.0 0.0 0.0 0.0 \n")   # old version
    file.write("C: 0.0 0.0 0.0 0.0 0.0 0.0 \n")
    file.write("R1: "+str(r11)+" "+str(r12)+" 0.0 0.0 0.0 0.0 \n")
    file.write("R2: "+str(r21)+" "+str(r22)+" 0.0 0.0 0.0 0.0 \n")
    file.write("R3: 0.0 0.0 "+str(r33)+" "+str(r34)+" 0.0 0.0 \n")
    file.write("R4: 0.0 0.0 "+str(r43)+" "+str(r44)+" 0.0 0.0 \n")
    file.write("R5: 0.0 0.0 0.0 0.0 "+str(r55)+" "+str(r56)+" \n")
    file.write("R6: 0.0 0.0 0.0 0.0 "+str(r65)+" "+str(r66)+" \n")
    
    file.close()

def generate_lte(lte_filename, parameters):
    if not os.path.isfile(lte_filename):
        file=open(lte_filename,"w")
    
        sentance = 'STM: MATR, L=0.0, filename="' + "mtr_file/" + parameters[1] + '.mtr' + '"; \n'
        file.write(sentance)
        sentance = 'W1: WATCH,mode="coordinates",filename="' + "sdds_file/" + str(parameters[5]) + '.sdds' + '",interval=1; \n'
        file.write(sentance)
        sentance = 'W2: WATCH,mode="coordinates",filename="' + "sdds_file/" + str(parameters[6]) + '.sdds' + '",interval=1; \n'
        file.write(sentance)
        sentance = parameters[3] + ': LINE=(W1)\n'
        file.write(sentance)
        sentance = parameters[4] + ': LINE=(STM,W2)\n'
        file.write(sentance)
        
        file.close()

def generate_beam_start_ele(beam_start_ele_filename, parameters):
    if not os.path.isfile(beam_start_ele_filename):
        file=open(beam_start_ele_filename,"w")
    
        sentance = '&change_particle\n'
        file.write(sentance)
        sentance = '  name = ' + parameters[0] +'\n'
        file.write(sentance)
        sentance = '&end\n'
        file.write(sentance)
        sentance = '\n'
        file.write(sentance)
        sentance = '&run_setup\n'
        file.write(sentance)
        sentance = '    lattice       = ' + 'lte_file/' + parameters[2] + '.lte' +',\n'
        file.write(sentance)
        sentance = '    use_beamline  = ' + parameters[3] +',\n'
        file.write(sentance)
        sentance = '    p_central_mev = ' + str(parameters[18]) +',\n'
        file.write(sentance)
        sentance = '    random_number_seed = ' + str(random.randint(0,1000)) +',\n'
        file.write(sentance)
        sentance = '&end\n'
        file.write(sentance)
        sentance = '\n'
        file.write(sentance)
        sentance = '&run_control\n'
        file.write(sentance)
        sentance = '    n_passes = 1\n'
        file.write(sentance)
        sentance = '&end\n'
        file.write(sentance)
        sentance = '\n'
        file.write(sentance)
        sentance = '&bunched_beam\n'
        file.write(sentance)
        sentance = '    n_particles_per_bunch = ' + str(parameters[19]) +',\n'
        file.write(sentance)
        sentance = '    emit_nx = ' + str(parameters[16]) +',\n'
        file.write(sentance)
        sentance = '    emit_ny = ' + str(parameters[17]) +',\n'
        file.write(sentance)
        sentance = '    beta_x = ' + str(parameters[7]) +',\n'
        file.write(sentance)
        sentance = '    beta_y = ' + str(parameters[8]) +',\n'
        file.write(sentance)
        sentance = '    sigma_dp = ' + str(parameters[15]) +',\n'
        file.write(sentance)
        sentance = '    sigma_s = ' + str(parameters[14]) +',\n'
        file.write(sentance)
        sentance = '    distribution_type[0] = "gaussian","gaussian","gaussian"' +'\n'
        file.write(sentance)
        sentance = '    enforce_rms_values[0] = 0, 0, 1' +'\n'
        file.write(sentance)
        sentance = '    distribution_cutoff[0] = 6, 6, 6' +'\n'
        file.write(sentance)
        sentance = '&end\n'
        file.write(sentance)
        sentance = '\n'
        file.write(sentance)
        sentance = '&track &end\n'
        file.write(sentance)
        sentance = '&stop &end\n'
        file.write(sentance)

        file.close()

def generate_beam_continue_ele(beam_continue_ele_filename, parameters):
    if not os.path.isfile(beam_continue_ele_filename):
        file=open(beam_continue_ele_filename,"w")
    
        sentance = '&change_particle\n'
        file.write(sentance)
        sentance = '  name = ' + parameters[0] +'\n'
        file.write(sentance)
        sentance = '&end\n'
        file.write(sentance)
        sentance = '\n'
        file.write(sentance)
        sentance = '&run_setup\n'
        file.write(sentance)
        sentance = '    lattice       = ' + 'lte_file/' + parameters[2] + '.lte' +',\n'
        file.write(sentance)
        sentance = '    use_beamline  = ' + parameters[4] +',\n'
        file.write(sentance)
        sentance = '    p_central_mev = ' + str(parameters[18]) +',\n'
        file.write(sentance)
        sentance = '&end\n'
        file.write(sentance)
        sentance = '\n'
        file.write(sentance)
        sentance = '&run_control\n'
        file.write(sentance)
        sentance = '    n_passes = 1\n'
        file.write(sentance)
        sentance = '&end\n'
        file.write(sentance)
        sentance = '\n'
        file.write(sentance)
        sentance = '&sdds_beam\n'
        file.write(sentance)
        sentance = '    input = ' + 'sdds_file/' + parameters[5] + '.sdds' +',\n'
        file.write(sentance)
        sentance = '&end\n'
        file.write(sentance)
        sentance = '\n'
        file.write(sentance)
        sentance = '&track &end\n'
        file.write(sentance)
        sentance = '&stop &end\n'
        file.write(sentance)

        file.close()


def elegant_starting_generating_sdds(parameter):
    sentance = 'elegant ele_command/' + parameter[0][3] + '.ele'
    os.system(sentance)
    
    sentance = 'elegant ele_command/' + parameter[1][3] + '.ele'
    os.system(sentance)
    

def elegant_continue(parameter, turn_number):

    number = str(int(turn_number))

    sentance = 'elegant ele_command/' + parameter[0][4] + '.ele'
    os.system(sentance)
    
    sentance = 'sddsprocess sdds_file/' + parameter[0][6] + '.sdds' + ' -redefine,parameter,Pass,' + number +',type=long -noWarning'
    os.system(sentance)

    sentance = 'elegant ele_command/' + parameter[1][4] + '.ele'
    os.system(sentance)

    sentance = 'sddsprocess sdds_file/' + parameter[1][6] + '.sdds' + ' -redefine,parameter,Pass,' + number +',type=long -noWarning'
    os.system(sentance)





