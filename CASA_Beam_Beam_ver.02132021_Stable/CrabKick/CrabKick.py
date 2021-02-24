#  This code is used to test the kick process of CRIB crossing
#  River and Vasiliy, 2018-2020

from function_to_deal_sdds_and_crab_crossing import *

print('')

# t1 = time.clock()
t1 = time.process_time()

''' Parameters for Calculating Crab Kick Process ************************ '''
if not os.path.isfile('Files_of_input/Crab_Kick_Parameters.casa'):
    sys.exit()

CK_parameter = Crab_Kick_Input_Parameters("Files_of_input/Crab_Kick_Parameters.casa")
Beam_name_1_in        = CK_parameter[0]
Beam_name_1_out       = CK_parameter[1]
Beam_name_2_in        = CK_parameter[2]
Beam_name_2_out       = CK_parameter[3]
number_of_turn        = CK_parameter[16]
Sigma_beam1_file      = CK_parameter[17]
Sigma_beam2_file      = CK_parameter[18]
sdds_output_frequency = CK_parameter[19]

if not os.path.isfile(Sigma_beam1_file):
    pass
else:
    Sigma_beam1_file_old = Sigma_beam1_file + '_old'
    shutil.move(Sigma_beam1_file, Sigma_beam1_file_old)

if not os.path.isfile(Sigma_beam2_file):
    pass
else:
    Sigma_beam2_file_old = Sigma_beam2_file + '_old'
    shutil.move(Sigma_beam2_file, Sigma_beam2_file_old)


# Generating initial sdds-file via Elegant *****************************
current_path = os.getcwd()
os.chdir("Elegant_Process")

parameter = input_parameters()
# generate_mtr_lte_ele(parameter)
elegant_starting_generating_sdds(parameter)

# copy all sdds to other directory
file_sdds = glob.glob('sdds_file/*.sdds')
for i in file_sdds:
    try:
        shutil.copy(i, 'beam_initial/')
    except:
        print("Failed to copy file: ", sys.exc_info())

os.chdir(current_path)
# **********************************************************************


if (parameter[0][0] == 'proton'):
    particle1_mass = 938.2733       # MeV
else:
    particle1_mass = 0.51099906     # MeV

if (parameter[1][0] == 'proton'):
    particle2_mass = 938.2733
else:
    particle2_mass = 0.51099906

gamma1 = parameter[0][18] / particle1_mass
gamma2 = parameter[1][18] / particle2_mass

beta_gama_beam1 = math.sqrt(gamma1**2.0 - 1.0)
beta_gama_beam2 = math.sqrt(gamma2**2.0 - 1.0)


sigma_x1 = math.sqrt(parameter[0][16] * parameter[0][7] / beta_gama_beam1)
sigma_y1 = math.sqrt(parameter[0][17] * parameter[0][8] / beta_gama_beam1)
sigma_z1 = parameter[0][14]
sigma_x2 = math.sqrt(parameter[1][16] * parameter[1][7] / beta_gama_beam2)
sigma_y2 = math.sqrt(parameter[1][17] * parameter[1][8] / beta_gama_beam2)
sigma_z2 = parameter[1][14]


# Starting Kicking ..................................................................
Sigma_beam1 = []
Sigma_beam2 = []
Sigma_beam1.append([sigma_x1, sigma_y1, sigma_z1])
Sigma_beam2.append([sigma_x2, sigma_y2, sigma_z2])

gamma_beta = []
gamma_beta.append([beta_gama_beam1, beta_gama_beam2])

beam1_phase_avg = []
beam2_phase_avg = []


if not os.path.exists('Files_of_Output'):
    os.mkdir('Files_of_Output')
else:
    Delete_directory('Files_of_Output')
    os.mkdir('Files_of_Output')


for i in range(number_of_turn):

    # ........  Starting Elegant Process ....................
    current_path = os.getcwd()
    os.chdir("Elegant_Process")
    elegant_continue(parameter, i)
    
    Beam1_in_name  = 'sdds_file/' + parameter[0][6] + '.sdds'
    Beam1_out_name = '../' + Beam_name_1_in
    shutil.copy(Beam1_in_name, Beam1_out_name)
    
    Beam2_in_name  = 'sdds_file/' + parameter[1][6] + '.sdds'
    Beam2_out_name = '../' + Beam_name_2_in
    shutil.copy(Beam2_in_name, Beam2_out_name)

    os.chdir(current_path)
    # ........  Ending Elegant Process ......................

    # ........  CASA Crab kicking ...........................

    # Result_temp = Crab_kick(CK_parameter)                     # do not record the phase-coordinates
    Result_temp = Crab_kick_track_coordinates(CK_parameter, i)   # record all phase-coordinates of all interactions between slices for each step for one turn
    
    Sigma_beam1.append(Result_temp[0])
    Sigma_beam2.append(Result_temp[1])
    
    gamma_beta.append(Result_temp[2])
    
    beam1_phase_avg.append(Result_temp[3])
    beam2_phase_avg.append(Result_temp[4])
    
    # Crab_kick(CK_parameter)
    Beam1_in_name  = Beam_name_1_out
    Beam1_out_name = 'Elegant_Process/sdds_file/' + parameter[0][5] + '.sdds'
    shutil.copy(Beam1_in_name, Beam1_out_name)

    # sentance = 'cp ' + Beam1_in_name + ' ' + Beam1_out_name
    # os.system(sentance)

    Beam2_in_name  = Beam_name_2_out
    Beam2_out_name = 'Elegant_Process/sdds_file/' + parameter[1][5] + '.sdds'
    shutil.copy(Beam2_in_name, Beam2_out_name)

    # sentance = 'cp ' + Beam2_in_name + ' ' + Beam2_out_name
    # os.system(sentance)

    
    if (sdds_output_frequency <= 0):
        pass
    elif ((sdds_output_frequency > 0) and (i == 0)):
        str_add = '_' + str(i+1)
        rename_sdds = Beam1_in_name + str_add
        shutil.move(Beam1_in_name, rename_sdds)
        # sentance = 'mv ' + Beam1_in_name + ' ' + rename_sdds
        # os.system(sentance)
        rename_sdds = Beam2_in_name + str_add
        shutil.move(Beam2_in_name, rename_sdds)
        # sentance = 'mv ' + Beam2_in_name + ' ' + rename_sdds
        # os.system(sentance)
    elif ((sdds_output_frequency > 0) and (i == (number_of_turn-1))):
        str_add = '_' + str(i+1)
        rename_sdds = Beam1_in_name + str_add
        shutil.move(Beam1_in_name, rename_sdds)
        # sentance = 'mv ' + Beam1_in_name + ' ' + rename_sdds
        # os.system(sentance)
        rename_sdds = Beam2_in_name + str_add
        shutil.move(Beam2_in_name, rename_sdds)
        # sentance = 'mv ' + Beam2_in_name + ' ' + rename_sdds
        # os.system(sentance)
    elif (((i+1) % sdds_output_frequency == 0) and (sdds_output_frequency > 0)):
        str_add = '_' + str(i+1)
        rename_sdds = Beam1_in_name + str_add
        shutil.move(Beam1_in_name, rename_sdds)
        # sentance = 'mv ' + Beam1_in_name + ' ' + rename_sdds
        # os.system(sentance)
        rename_sdds = Beam2_in_name + str_add
        shutil.move(Beam2_in_name, rename_sdds)
        # sentance = 'mv ' + Beam2_in_name + ' ' + rename_sdds
        # os.system(sentance)


# np.savetxt('Files_of_Output/sigma_beam1.casa', Sigma_beam1)
# np.savetxt('Files_of_Output/sigma_beam2.casa', Sigma_beam2)

np.savetxt(Sigma_beam1_file, Sigma_beam1)
np.savetxt(Sigma_beam2_file, Sigma_beam2)

np.savetxt('Files_of_Output/gamma_beta.casa', gamma_beta)

np.savetxt('Files_of_Output/beam1_phase_avg.casa', beam1_phase_avg)
np.savetxt('Files_of_Output/beam2_phase_avg.casa', beam2_phase_avg)


# Calculatiing Luminosity ...........................................................
beta_star_x_for_beam1 = parameter[0][7]
beta_star_y_for_beam1 = parameter[0][8]
beta_star_x_for_beam2 = parameter[1][7]
beta_star_y_for_beam2 = parameter[1][8]

Number_of_beam1_Slices  = CK_parameter[4]
Number_of_beam2_Slices  = CK_parameter[5]
Collision_Frequency     = CK_parameter[14]
Number_of_particle1     = CK_parameter[20]
Number_of_particle2     = CK_parameter[21]

Angle_x                 = CK_parameter[8]   # half angle for our method
Angle_y                 = 0.0
offset_x                = 0.0
offset_y                = 0.0
Sigma_beam1_file        = CK_parameter[17]
Sigma_beam2_file        = CK_parameter[18]

sigma_star_of_beam1 = Sigma_beam1
sigma_star_of_beam2 = Sigma_beam2

# number_of_sigma = len(sigma_star_of_beam2)
def len_of_file():
    with open(Sigma_beam2_file) as file_of_sigma:
        return len(list(file_of_sigma))
    file_of_sigma.close()
number_of_sigma = len_of_file()

beta_star_of_beam1  = [beta_star_x_for_beam1, beta_star_y_for_beam1]
beta_star_of_beam2  = [beta_star_x_for_beam2, beta_star_y_for_beam2]

Half_Crossing_Angle = [Angle_x, Angle_y]
Offset_defference   = [offset_x, offset_y]

Number_of_colliding_bunches = 1


LuminosityData = []
if (number_of_sigma == 1):

    L0_luminosity = Nominal_Luminosity(Number_of_particle1, Number_of_particle2, Number_of_colliding_bunches, \
                         Collision_Frequency, sigma_star_of_beam1, sigma_star_of_beam2)

    reduction_numerical = Reduction_Factor_for_Numerical_Solution(sigma_star_of_beam1, sigma_star_of_beam2, \
                                                              beta_star_of_beam1, beta_star_of_beam2, \
                                                              Half_Crossing_Angle, Offset_defference)
    reduction_analytic_1 = Reduction_Factor_for_Analynic_Solution_1(sigma_star_of_beam1, sigma_star_of_beam2, \
                                                              beta_star_of_beam1, beta_star_of_beam2, \
                                                              Half_Crossing_Angle, Offset_defference)
    reduction_analytic_2 = Reduction_Factor_for_Analynic_Solution_2(sigma_star_of_beam1, sigma_star_of_beam2, \
                                                              beta_star_of_beam1, beta_star_of_beam2, \
                                                              Half_Crossing_Angle, Offset_defference)
    reduction_analytic_3 = Reduction_Factor_for_Analynic_Solution_3(sigma_star_of_beam1, sigma_star_of_beam2, \
                                                              beta_star_of_beam1, beta_star_of_beam2, \
                                                              Half_Crossing_Angle, Offset_defference)

    Numerical_Luminosity = Luminosity_Calculation(L0_luminosity, reduction_numerical)
    Analytic1_Luminosity = Luminosity_Calculation(L0_luminosity, reduction_analytic_1)
    Analytic2_Luminosity = Luminosity_Calculation(L0_luminosity, reduction_analytic_2)
    Analytic3_Luminosity = Luminosity_Calculation(L0_luminosity, reduction_analytic_3)

    LuminosityData.append([reduction_numerical, Numerical_Luminosity, Analytic1_Luminosity, Analytic2_Luminosity, Analytic3_Luminosity])

    Luminosity_output = np.array(LuminosityData)
else:
    for i in range (number_of_sigma):
        L0_luminosity = Nominal_Luminosity(Number_of_particle1, Number_of_particle2, Number_of_colliding_bunches, \
                         Collision_Frequency, sigma_star_of_beam1[i], sigma_star_of_beam2[i])

        reduction_numerical = Reduction_Factor_for_Numerical_Solution(sigma_star_of_beam1[i], sigma_star_of_beam2[i], \
                                                              beta_star_of_beam1, beta_star_of_beam2, \
                                                              Half_Crossing_Angle, Offset_defference)
        reduction_analytic_1 = Reduction_Factor_for_Analynic_Solution_1(sigma_star_of_beam1[i], sigma_star_of_beam2[i], \
                                                              beta_star_of_beam1, beta_star_of_beam2, \
                                                              Half_Crossing_Angle, Offset_defference)
        reduction_analytic_2 = Reduction_Factor_for_Analynic_Solution_2(sigma_star_of_beam1[i], sigma_star_of_beam2[i], \
                                                              beta_star_of_beam1, beta_star_of_beam2, \
                                                              Half_Crossing_Angle, Offset_defference)
        reduction_analytic_3 = Reduction_Factor_for_Analynic_Solution_3(sigma_star_of_beam1[i], sigma_star_of_beam2[i], \
                                                              beta_star_of_beam1, beta_star_of_beam2, \
                                                              Half_Crossing_Angle, Offset_defference)

        Numerical_Luminosity = Luminosity_Calculation(L0_luminosity, reduction_numerical)
        Analytic1_Luminosity = Luminosity_Calculation(L0_luminosity, reduction_analytic_1)
        Analytic2_Luminosity = Luminosity_Calculation(L0_luminosity, reduction_analytic_2)
        Analytic3_Luminosity = Luminosity_Calculation(L0_luminosity, reduction_analytic_3)

        LuminosityData.append([reduction_numerical, Numerical_Luminosity, Analytic1_Luminosity, Analytic2_Luminosity, Analytic3_Luminosity])

if (number_of_sigma > 1):
    Luminosity_output = np.array(LuminosityData)
    Luminosity_output_Curve = Luminosity_output.T
    Numerical_data = Luminosity_output_Curve[1]
    Analytic1_data = Luminosity_output_Curve[2]
    Analytic2_data = Luminosity_output_Curve[3]
    Analytic3_data = Luminosity_output_Curve[4]
    Luminosity_2D_Curve(Numerical_data, Analytic1_data, Analytic2_data, Analytic3_data)


np.savetxt('Files_of_Output/Luminosity.casa', Luminosity_output, fmt='%0.4f       %.4e   %.4e   %.4e   %.4e', header='Numerical    Numerical    Analytic-1   Analytic-2   Analytic-3\nReduction    Luminosity   Luminosity   Luminosity   Luminosity\n', comments='')


print('')
print('**************************************************************')
print('*                                                            *')
print('*       Crab Kick calculation is completed.                  *')
print('*                                                            *')
print('*       Thanks for using JLab CASA Beam-Beam.                *')
print('*                                                            *')
print('*                                                            *')
print('**************************************************************')
print('')

# t2 = time.clock()
t2 = time.process_time()

print('')
print("run time:%f s" % (t2 - t1))
print('')



