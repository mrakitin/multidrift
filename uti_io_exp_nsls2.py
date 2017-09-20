
from __future__ import print_function #Python 2.7 compatibility
from srwlib import *

import os
#from array import *

#**********************Auxiliary function to read-in and parse tabulated scan data (for HXN and maybe other beamlines)
def uti_io_read_scan_data(_file_path):

    sSep = ','
    f = open(_file_path, 'r')
    lines = f.readlines()

    nLines = len(lines)
    resData = []
    for i in range(nLines):
        curParts = lines[i].split(sSep)
        if(len(curParts) > 0): resData.append(curParts)
    
    return resData

#**********************Extracting data columnns from all scan file data
def uti_io_extract_scan_cols(_all_data, _col_names):

    allColNames = _all_data[0]
    nAllCols = len(allColNames)
    #print(allColNames)
    
    nDataPts = len(_all_data) - 1
    nCols = len(_col_names)
    resData = []
    for i in range(nCols):
        curName = _col_names[i]
        for ii in range(nAllCols):
            if(curName == allColNames[ii]):
                curDataAr = array('d', [0]*nDataPts)
                for j in range(nDataPts):
                    curDataAr[j] = float(_all_data[j + 1][ii])
                resData.append(curDataAr)
    return resData

#**********************

sFolderName = 'C:\\SR_Magnets\\InsertionDevices\\Commissioning\\NSLSII\\HXN\\MeasHXN_Nov02_2015'
sInFileNameCore = 'scan_4247'
sOutFileName = sInFileNameCore + '_spec.dat'
arColNamesToExtract = ['dcm_th', 'sclr2_ch2']

allData = uti_io_read_scan_data(os.path.join(os.getcwd(), sFolderName, sInFileNameCore + '.txt'))
#print(allData)

resDataCols = uti_io_extract_scan_cols(allData, arColNamesToExtract)
#print(resDataCols)

srwl_uti_write_data_cols(os.path.join(os.getcwd(), sFolderName, sOutFileName), resDataCols, '\t')


