
#  ******     This code is used for reading and writing SDDS-data
#  ******     River and Vasiliy
#  ******     2018-2020



import sddsdata, sys, time


class SDDS:
    """This class implements SDDS datasets."""

    def __init__(self, index):
        #define common SDDS definitions
        self.SDDS_VERBOSE_PrintErrors = 1
        self.SDDS_EXIT_PrintErrors = 2
        self.SDDS_CHECK_OKAY = 0
        self.SDDS_CHECK_NONEXISTENT = 1
        self.SDDS_CHECK_WRONGTYPE = 2
        self.SDDS_CHECK_WRONGUNITS = 3
        self.SDDS_LONGDOUBLE = 1
        self.SDDS_DOUBLE = 2
        self.SDDS_REAL64 = 2
        self.SDDS_FLOAT = 3
        self.SDDS_REAL32 = 3
        self.SDDS_LONG = 4
        self.SDDS_INT32 = 4
        self.SDDS_ULONG = 5
        self.SDDS_UINT32 = 5
        self.SDDS_SHORT = 6
        self.SDDS_INT16 = 6
        self.SDDS_USHORT = 7
        self.SDDS_UINT16 = 7
        self.SDDS_STRING = 8
        self.SDDS_CHARACTER = 9
        self.SDDS_NUM_TYPES = 9
        self.SDDS_BINARY = 1
        self.SDDS_ASCII = 2
        #only indexes of 0 through 19 are allowed
        if index >= 0 and index < 20:
            self.index = index
        else:
            self.index = 0
        #initialize data storage variables
        self.description = ["", ""]
        self.parameterName = []
        self.columnName = []
        self.parameterDefinition = []
        self.columnDefinition = []
        self.parameterData = []
        self.columnData = []
        self.mode = self.SDDS_ASCII

    def Data_Storage_mode(self, input):
        if sddsdata.InitializeInput(self.index, input) != 1:
            time.sleep(2)
            if sddsdata.InitializeInput(self.index, input) != 1:
                sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)
        self_data_mode = sddsdata.GetMode(self.index)

        #close SDDS file
        if sddsdata.Terminate(self.index) != 1:
            sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)

        return self_data_mode

    def Load_Data(self, input):
        """Load an SDDS file into an SDDS class."""

        try:
            #open SDDS file
            if sddsdata.InitializeInput(self.index, input) != 1:
                time.sleep(2)
                if sddsdata.InitializeInput(self.index, input) != 1:
                    sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)

            #get data storage mode (SDDS_ASCII or SDDS_BINARY)
            self.mode = sddsdata.GetMode(self.index)
            # if (self.mode == 0):
                # print ()
                # print ('Data storage mode is ', 'ASCII')
                # print ()
            # if (self.mode == 1):
                # print ()
                # print ('Data storage mode is ', 'BINARY')
                # print ()

            #get description text and contents
            self.description = sddsdata.GetDescription(self.index)
            # print ()
            # print ('Description text and contents ')
            # print (self.description)
            # print ()

            #get parameter names
            self.parameterName = sddsdata.GetParameterNames(self.index)
            # print ()
            # print ('Parameter names ')
            # print (self.parameterName)
            # print ()
            numberOfParameters = len(self.parameterName)
            # print ()
            # print ('Number of Parameters = ', numberOfParameters)
            # print ()

            #get column names
            self.columnName = sddsdata.GetColumnNames(self.index)
            # print ()
            # print ('Column names ')
            # print (self.columnName)
            # print ()
            numberOfColumns = len(self.columnName)
            # print ()
            # print ('Number of Columns = ', numberOfColumns)
            # print ()

            #get parameter definitions
            self.parameterDefinition = list(range(numberOfParameters))
            for i in range(numberOfParameters):
                self.parameterDefinition[i] = sddsdata.GetParameterDefinition(self.index, self.parameterName[i])
            # print ()
            # print ('Parameter definitions ')
            # print (self.parameterDefinition)
            # print ()

            #get column definitions
            self.columnDefinition = list(range(numberOfColumns))
            for i in range(numberOfColumns):
                self.columnDefinition[i] = sddsdata.GetColumnDefinition(self.index, self.columnName[i])
            # print ()
            # print ('Column definitions ')
            # print (self.columnDefinition)
            # print ()

            #initialize parameter and column data     !!!! this is mentioned by River and Vasiliy !!!!
            self.parameterData = list(range(numberOfParameters))
            self.columnData = list(range(numberOfColumns))
            for i in range(numberOfParameters):
                self.parameterData[i] = []
            for i in range(numberOfColumns):
                self.columnData[i] = []

            #read in SDDS data
            page = sddsdata.ReadPage(self.index)
            if page != 1:
                sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)
            while page > 0:
                for i in range(numberOfParameters):
                    self.parameterData[i].append(sddsdata.GetParameter(self.index,i))
                rows = sddsdata.RowCount(self.index);
                if rows > 0:
                    for i in range(numberOfColumns):
                           self.columnData[i].append(sddsdata.GetColumn(self.index,i))
                else:
                    for i in range(numberOfColumns):
                        self.columnData[i].append([])
                 
                page = sddsdata.ReadPage(self.index)

                # global x_position
                # x_position = sddsdata.GetColumn(0,0)
                # print (x_position[0], x_position[2], x_position[9999])
                # global y_position
                # y_position = sddsdata.GetColumn(0,2)
                # print (y_position[0], y_position[2], y_position[9999])
                # global z_position
                # z_position = sddsdata.GetColumn(0,4)
                # print (z_position[0], z_position[2], z_position[9999])
 
                # print (sddsdata.GetColumn(0,0))  # (self.index,i) = (0, 0) means all data for x position
                # print (sddsdata.GetColumn(0,2))  # (self.index,i) = (0, 2) means all data for y position
                # print (sddsdata.GetColumn(0,4))  # (self.index,i) = (0, 4) means all data for z position 
                page = sddsdata.ReadPage(self.index)

            #close SDDS file
            if sddsdata.Terminate(self.index) != 1:
                sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)

        except:
            sddsdata.PrintErrors(self.SDDS_VERBOSE_PrintErrors)
            raise
        # print ()
        # print (x_position[0], x_position[2], x_position[9999])
        # return x_position[0]
        # return x_position, y_position, z_position
        # return sddsdata.columnData[0][0], sddsdata.columnData[2][0], sddsdata.columnData[4][0]


    def Save_Data(self, output):
        """Save an SDDS class to an SDDS file."""

        try:
            #check for invalid SDDS data
            numberOfParameters = len(self.parameterName)
            numberOfColumns = len(self.columnName)
            pages = 0
            if numberOfParameters != len(self.parameterData):
                raise Exception("unmatched parameterName and parameterData")
            if numberOfColumns != len(self.columnData):
                raise Exception("unmatched columnName and columnData")
            if numberOfParameters != len(self.parameterDefinition):
                raise Exception("unmatched parameterName and parameterDefinition")
            if numberOfColumns != len(self.columnDefinition):
                raise Exception("unmatched columnName and columnDefinition")
            if numberOfParameters > 0:
                pages = len(self.parameterData[0])
            elif numberOfColumns > 0:
                pages = len(self.columnData[0])
            for i in range(numberOfParameters):
                if pages != len(self.parameterData[i]):
                    raise Exception("unequal number of pages in parameter data")
            for i in range(numberOfColumns):
                if pages != len(self.columnData[i]):
                    raise Exception("unequal number of pages in column data")
            for page in range(pages):               
                rows = 0
                if numberOfColumns > 0:
                    rows = len(self.columnData[0][page])
                for i in range(numberOfColumns):
                    if rows != len(self.columnData[i][page]):
                        raise Exception("unequal number of rows in column data")

            #open SDDS output file
            if sddsdata.InitializeOutput(self.index, self.mode, 1, self.description[0], self.description[1], output) != 1:
                sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)

            #define parameters and columns
            for i in range(numberOfParameters):
                if sddsdata.DefineParameter(self.index, self.parameterName[i],
                                         self.parameterDefinition[i][0],
                                         self.parameterDefinition[i][1],
                                         self.parameterDefinition[i][2],
                                         self.parameterDefinition[i][3],
                                         self.parameterDefinition[i][4],
                                         self.parameterDefinition[i][5]) == -1:
                      sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)
            for i in range(numberOfColumns):               
                if sddsdata.DefineColumn(self.index, self.columnName[i],
                                         self.columnDefinition[i][0],
                                         self.columnDefinition[i][1],
                                         self.columnDefinition[i][2],
                                         self.columnDefinition[i][3],
                                         self.columnDefinition[i][4],
                                         self.columnDefinition[i][5]) == -1:
                    sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)

            #write SDDS header
            if sddsdata.WriteLayout(self.index) != 1:
                sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)

            #write SDDS data
            for page in range(pages):               
                rows = 0
                if numberOfColumns > 0:
                    rows = len(self.columnData[0][page])
                if sddsdata.StartPage(self.index, rows) != 1:
                    sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)
                for i in range(numberOfParameters):
                    if sddsdata.SetParameter(self.index, i, self.parameterData[i][page]) != 1:
                        sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)
                for i in range(numberOfColumns):
                    if sddsdata.SetColumn(self.index, i, self.columnData[i][page]) != 1:
                        sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)
                if sddsdata.WritePage(self.index) != 1:
                    sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)

            #close SDDS output file
            if sddsdata.Terminate(self.index) != 1:
                sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)
        except:
            sddsdata.PrintErrors(self.SDDS_VERBOSE_PrintErrors)
            raise

    def setDescription(self, text, contents):
        self.description = [text, contents]

    def defineParameter(self, name, symbol, units, description, formatString, type, fixedValue):
        self.parameterName.append(name)
        self.parameterDefinition.append([symbol, units, description, formatString, type, fixedValue])
        self.parameterData.append([])

    def defineSimpleParameter(self, name, type):
        self.parameterName.append(name)
        self.parameterDefinition.append(["", "", "", "", type, ""])
        self.parameterData.append([])

    def defineColumn(self, name, symbol, units, description, formatString, type, fieldLength):
        self.columnName.append(name)
        self.columnDefinition.append([symbol, units, description, formatString, type, fieldLength])
        self.columnData.append([])

    def defineSimpleColumn(self, name, type):
        self.columnName.append(name)
        self.columnDefinition.append(["", "", "", "", type, 0])
        self.columnData.append([])

    def setParameterValueList(self, name, valueList):
        numberOfParameters = len(self.parameterName)
        for i in range(numberOfParameters):
            if self.parameterName[i] == name:
                self.parameterData[i] = valueList
                return
        msg = "invalid parameter name " + name
        raise Exception(msg)
    
    def setParameterValue(self, name, value, page):
        page = page - 1
        numberOfParameters = len(self.parameterName)
        for i in range(numberOfParameters):
            if self.parameterName[i] == name:
                if len(self.parameterData[i]) == page:
                    self.parameterData[i][page:] = [value]
                elif len(self.parameterData[i]) < page or page < 0:
                    msg = "invalid page " + str(page+1)
                    raise Exception(msg)
                else:
                    self.parameterData[i][page] = [value]
                return
        msg = "invalid parameter name " + name
        raise Exception(msg)

    def setColumnValueLists(self, name, valueList):
        numberOfColumns = len(self.columnName)
        for i in range(numberOfColumns):
            if self.columnName[i] == name:
                self.columnData[i] = valueList
                return
        msg = "invalid column name " + name
        raise Exception(msg)
    
    def setColumnValueList(self, name, valueList, page):
        page = page - 1
        numberOfColumns = len(self.columnName)
        for i in range(numberOfColumns):
            if self.columnName[i] == name:
                if len(self.columnData[i]) == page:
                    self.columnData[i][page:] = [valueList]
                elif len(self.columnData[i]) < page or page < 0:
                    msg = "invalid page " + str(page+1)
                    raise Exception(msg)
                else:
                    self.columnData[i][page] = [valueList]
                return
        msg = "invalid column name " + name
        raise Exception(msg)
    
    def setColumnValue(self, name, value, page, row):
        page = page - 1
        row = row - 1
        numberOfColumns = len(self.columnName)
        for i in range(numberOfColumns):
            if self.columnName[i] == name:
                if len(self.columnData[i]) == page:
                    if row == 0:
                        self.columnData[i][page:] = [[value]]
                    else:
                        msg = "invalid row " + str(row+1)
                        raise Exception(msg)
                elif len(self.columnData[i]) < page or page < 0:
                    msg = "invalid page " + str(page+1)
                    raise Exception(msg)
                else:
                    if len(self.columnData[i][page]) == row:
                        self.columnData[i][page][row:] = [value]
                    elif len(self.columnData[i][page]) < row or row < 0:
                        msg = "invalid row " + str(row+1)
                        raise Exception(msg)
                    else:
                        self.columnData[i][page][row] = [value]
                return
        msg = "invalid column name " + name
        raise Exception(msg)







