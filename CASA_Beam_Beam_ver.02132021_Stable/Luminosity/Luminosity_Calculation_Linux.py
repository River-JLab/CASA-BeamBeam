# ********************************************************* #
#  This is a Test Code for JLAB-CSSA Luminosity Calculation #
#  ---                                                      #
#  River Huang and Vasiliy Morozov                          #
#  2020                                                     #
# ********************************************************* #

from Luminosity_Calculation_Function import *

# parameters
parameters = Luminosity_parameter()

particle1_mass = float(parameters[0])
particle2_mass = float(parameters[1])
beam1_energy = float(parameters[2])
beam2_energy = float(parameters[3])

beta_star_x_for_beam1 = float(parameters[4])
beta_star_x_for_beam2 = float(parameters[5])
beta_star_y_for_beam1 = float(parameters[6])
beta_star_y_for_beam2 = float(parameters[7])

Normalized_emittance_x_of_beam1 = float(parameters[8])
Normalized_emittance_x_of_beam2 = float(parameters[9])
Normalized_emittance_y_of_beam1 = float(parameters[10])
Normalized_emittance_y_of_beam2 = float(parameters[11])

Length_of_beam1 = float(parameters[12])
Length_of_beam2 = float(parameters[13])

Number_of_particle1 = float(parameters[14])
Number_of_particle2 = float(parameters[15])

Collision_Frequency = float(parameters[16])

Angle_x = float(parameters[17])
Angle_y = float(parameters[18])

offset_x = float(parameters[19])
offset_y = float(parameters[20])

Number_of_colliding_bunches = float(parameters[21])

if not os.path.isfile("File_of_Input/sigma_beam1.casa"):
    Generate_RMS(parameters)
if not os.path.isfile("File_of_Input/sigma_beam2.casa"):
    Generate_RMS(parameters)


# number_of_sigma = len(sigma_of_beam1)
def len_of_file():
    with open('File_of_Input/sigma_beam2.casa') as file_of_sigma:
        return len(list(file_of_sigma))
    file_of_sigma.close()
number_of_sigma = len_of_file()


sigma_star_of_beam1 = np.loadtxt('File_of_Input/sigma_beam1.casa')
sigma_star_of_beam2 = np.loadtxt('File_of_Input/sigma_beam2.casa')

beta_star_of_beam1  = [beta_star_x_for_beam1, beta_star_y_for_beam1]
beta_star_of_beam2  = [beta_star_x_for_beam2, beta_star_y_for_beam2]
Half_Crossing_Angle = [Angle_x, Angle_y]
Offset_defference   = [offset_x, offset_y]

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


# filename = open('File_of_Output/Luminosity.casa', 'w')
# filename.write(' '.join(map(str, Luminosity_output)))
# filename.close()
np.savetxt('File_of_Output/Luminosity.casa', Luminosity_output, fmt='%0.4f       %.4e   %.4e   %.4e   %.4e', header='Numerical    Numerical    Analytic-1   Analytic-2   Analytic-3\nReduction    Luminosity   Luminosity   Luminosity   Luminosity\n', comments='')


print('')
print('**************************************************************')
print('*                                                            *')
print('*       Luminosity Calculation is completed. Thanks.         *')
print('*                                                            *')
print('**************************************************************')
print('')


