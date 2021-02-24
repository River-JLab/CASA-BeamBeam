#  This is the example-codes used for generating, saving and loading SDDS-format file
#  River and Vasiliy, 2018

import sys, time, math
import os.path
import cmath
import numpy as np
from sdds_load_save import *
from scipy import special

# The followings are 3 examples for generating and saving SDDS data ******************************
def Generate_Save_SDDS_Data_Example_1(output):
    '''Generate and Save a example-SDDS file using the SDDS class.'''

    data = SDDS(0)
    data.mode = 2 # ASCII mode
    
    data.description[0] = "text"
    data.description[1] = "contents"
    data.parameterName = ["Short", "Long", "Float", "Double", "String", "Character"]
    data.parameterData = [[1, 6], [2, 7], [3, 8], [4, 9], ["Five", "Ten"], ["Six", "Eleven"]]
    data.parameterDefinition = [["", "", "", "", data.SDDS_SHORT,     ""],
                                ["", "", "", "", data.SDDS_LONG,      ""],
                                ["", "", "", "", data.SDDS_FLOAT,     ""],
                                ["", "", "", "", data.SDDS_DOUBLE,    ""],
                                ["", "", "", "", data.SDDS_STRING,    ""],
                                ["", "", "", "", data.SDDS_CHARACTER, ""]]
    data.columnName = ["short_integer", "long_integer", "float", "double", "string", "character"]
    data.columnData = [[[1, 2, 3], [-1.01, -1.02, -1.03, -1.04]],
                       [[1, 2, 3], [2.01, 2.02, 2.03, 2.04]],
                       [[1, 2, 3], [-3.01, -3.02, -3.03, -3.04]],
                       [[1, 2, 3], [4.01, 4.02, 4.03, 4.04]],
                       [["row 1", "row 2", "row 3"], ["row 1", "row 2", "row 3", "row 4"]],
                       [["x", "y", "z"], ["a", "b", "c", "d"]]]
    data.columnDefinition = [["", "", "", "", data.SDDS_SHORT,    0],
                             ["", "", "", "", data.SDDS_LONG,     0],
                             ["", "", "", "", data.SDDS_FLOAT,    0],
                             ["", "", "", "", data.SDDS_DOUBLE,   0],
                             ["", "", "", "", data.SDDS_STRING,   0],
                             ["", "", "", "", data.SDDS_CHARACTER,0]]
    data.Save_Data(output)
    del data

def Generate_Save_SDDS_Data_Example_2(output):
    '''Generate and Save a example-SDDS file using the SDDS class.'''

    data = SDDS(0)
    data.mode = 2 # ASCII mode

    data.setDescription("text", "contents")
    names = ["Short", "Long", "Float", "Double", "String", "Character"]
    types = [data.SDDS_SHORT, data.SDDS_LONG, data.SDDS_FLOAT, data.SDDS_DOUBLE, data.SDDS_STRING, data.SDDS_CHARACTER]

    for i in range(6):
        data.defineSimpleParameter(names[i] + "P", types[i])
        data.defineSimpleColumn(names[i] + "C", types[i])

    parameterData = [[1, 6], [2, 7], [3, 8], [4, 9], ["Five", "Ten"], ["Six", "Eleven"]]
    for i in range(6):
        data.setParameterValueList(names[i] + "P", parameterData[i])

    columnData = [[[1, 2, 3], [-1.01, -1.02, -1.03, -1.04]],
                  [[1, 2, 3], [2.01, 2.02, 2.03, 2.04]],
                  [[1, 2, 3], [-3.01, -3.02, -3.03, -3.04]],
                  [[1, 2, 3], [4.01, 4.02, 4.03, 4.04]],
                  [["row 1", "row 2", "row 3"], ["row 1", "row 2", "row 3", "row 4"]],
                  [["x", "y", "z"], ["a", "b", "c", "d"]]]
    for i in range(6):          
        data.setColumnValueLists(names[i] + "C", columnData[i])

    data.Save_Data(output)
    del data

# Integer number was difined in SDDS class. Please see 'sdds_load_save.py'
def Generate_Save_SDDS_Data_Example_for_Elegant(output):
    '''Elegant Format: Generate and Save a SDDS file using the SDDS class.'''

    data = SDDS(0)
    data.mode = 1 # Binary mode
    data.description = ['watch-point phase space--input: bunched_beam.ele  lattice: elegant.lte',
                        'watch-point phase space']
                        
    # Parameter Details *************************************************************
    data.parameterName = [  'Step', 
                            'pCentral', 
                            'Charge', 
                            'Particles', 
                            'IDSlotsPerBunch', 
                            'Pass', 
                            'PassLength', 
                            'PassCentralTime', 
                            'ElapsedTime', 
                            'ElapsedCoreTime', 
                            's', 
                            'Description', 
                            'PreviousElementName']
    data.parameterData = [  [1],
                            [63.94634722234839], 
                            [0.0], 
                            [3], 
                            [3], 
                            [0], 
                            [169.26000000000013], 
                            [0.0], 
                            [0.0], 
                            [0.0], 
                            [0.0], 
                            [''], 
                            ['_BEG_']]
    data.parameterDefinition = [['',         '',        'Simulation step',                                 '',   3, ''      ], 
                                ['p$bcen$n', 'm$be$nc', 'Reference beta*gamma',                            '',   1, ''      ], 
                                ['',         'C',       'Beam charge',                                     '',   1, ''      ], 
                                ['',         '',        'Number of particles',                             '',   3, ''      ], 
                                ['',         '',        'Number of particle ID slots reserved to a bunch', '',   3, ''      ], 
                                ['',         '',        '',                                                '',   3, ''      ], 
                                ['',         'm',       '',                                                '',   1, ''      ], 
                                ['',         's',       '',                                                '',   1, ''      ], 
                                ['',         's',       '',                                                '',   1, ''      ], 
                                ['',         's',       '',                                                '',   1, ''      ], 
                                ['',         'm',       '',                                                '',   1, ''      ], 
                                ['',         '',        '',                                                '%s', 7, ''      ], 
                                ['',         '',        '',                                                '%s', 7, '_BEG_' ]]
    # **********************************************************************************
    
    # Data Details *********************************************************************
    data.columnName = ['x', 'xp', 'y', 'yp', 't', 'p', 'dt', 'particleID']
    x  = [1.1, 1.2, 1.3]
    xp = [2.1, 2.2, 2.3]
    y  = [3.1, 3.2, 3.3]
    yp = [4.1, 4.2, 4.3]
    t  = [5.1, 5.2, 5.3]
    p  = [6.1, 6.2, 6.3]
    dt = [7.1, 7.2, 7.3]
    particleID = [1, 2, 3]
    data.columnData = [ [x], [xp], [y], [yp], [t], [p], [dt], [particleID]]

    data.columnDefinition = [   ['', 'm',       '', '', 1, 0], 
                                ['', '',        '', '', 1, 0], 
                                ['', 'm',       '', '', 1, 0], 
                                ['', '',        '', '', 1, 0], 
                                ['', 's',       '', '', 1, 0], 
                                ['', 'm$be$nc', '', '', 1, 0], 
                                ['', 's',       '', '', 1, 0], 
                                ['', '',        '', '', 3, 0]]
    # **********************************************************************************

    # Save SDDS data *******************************************************************
    data.Save_Data(output)
    # **********************************************************************************

    return
    
''' Three Examples of Generating SDDS data'''
Generate_Save_SDDS_Data_Example_1('Example_1.sdds')
Generate_Save_SDDS_Data_Example_2('Example_2.sdds')
Generate_Save_SDDS_Data_Example_for_Elegant('Example_Elegant.sdds')
# The above are examples of generating and saving SDDS data **************************************






# Belowings show how to load SDDS data ***********************************************************
''' For proper using the parameter and variable of SDDS data, 
    we suggest you run "sddsquery xxxxxxxx.sdds"
    so that you can know the structure and format of SDDS data.
    River and Vasiliy
'''
def Load_SDDS_Data(input):

    data = SDDS(0)
    data.Load_Data(input)

    parameter_number = len(data.parameterData)
    print('parameter_number =', parameter_number)
    parameter = [[] for row in range(parameter_number)]
    for i in range(parameter_number):
        parameter[i] = data.parameterData[i][0]
        print (parameter[i])
    
    variable_number = len(data.columnData)
    print('variable_number =', variable_number)
    variable = [[] for row in range(variable_number)]
    for i in range(variable_number):
        variable[i] = data.columnData[i][0]
        print (variable[i])
    
    return parameter, variable

''' Three Examples of loading SDDS data'''
Load_SDDS_Data('Example_Elegant.sdds')
Load_SDDS_Data('Example_1.sdds')
Load_SDDS_Data('Example_2.sdds')
# The above are examples of loading SDDS data ****************************************************



























