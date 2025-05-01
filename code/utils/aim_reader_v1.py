"""
Script Name: aim_reader_v1.py
Description: This script:
                - Reads Scanco AIM files and converts them to numpy arrays

*Adaptation of the script by Ghislain MAQUER, University of Bern
ghislain.maquer@istb.unibe.ch
"""


import sys
import os
import time
import struct
import numpy
import vtk 
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import re
from time import strptime

def vtk2numpy(imvtk): 
    """ 
    converts a vtk array in a numpy array
    """
    dim  = imvtk.GetDimensions()
    data = imvtk.GetPointData().GetScalars()
    imnp = vtk_to_numpy(data)
    imnp = imnp.reshape(dim[2], dim[1], dim[0])
    # vtk and numpy have different array conventions
    imnp = imnp.transpose(2,1,0)
    imnp = imnp.transpose(1,0,2)
    imnp = numpy.flip(imnp,axis=0)
    return imnp

def get_AIM_ints( f ): 
    """
    reads the integer data of the AIM file to find the header length
    by Glen L. Niebur, University of Notre Dame (2010)
    """
    nheaderints = 32
    nheaderfloats = 8
    f.seek(0)
    # Read the first 128 bytes of the bunary file (each 4 byte one 32 bit integer)
    binints = f.read( nheaderints * 4 )
    # Convert bytes to integers
    header_int = struct.unpack('=32i', binints )
    return( header_int )

def AIMreader( fileINname ):
    """
    reads the AIM file
    returns infomation about image size and BMD calibration
    """
    # reads header
    # print ("     "+fileINname)
    f = open( fileINname, 'rb' )
    AIM_ints = get_AIM_ints( f )
    # checks AIM version and compression
    if int(AIM_ints[5])==16:
        # print ("     -> version 020")
        if   int(AIM_ints[10]) == 131074:   format = 'short';#          print("     -> format "+format)
        elif int(AIM_ints[10]) == 65537:    format = 'char';#           print("     -> format "+format)
        elif int(AIM_ints[10]) == 1376257:  format = 'bin compressed';# print("     -> format "+format+" not supported! Exiting!")#; exit(1)
        else:                               format = 'unknown';#        print("     -> format "+format+"! Exiting!")#; exit(1)
        header = f.read( AIM_ints[2] )
        header = header.decode("ISO-8859-1")
        header_len = len(header)+160 
        extents    = ( 0, AIM_ints[14]-1, 0, AIM_ints[15]-1, 0, AIM_ints[16]-1 )
        coords = (AIM_ints[11], AIM_ints[12], AIM_ints[13] )
    else:
        # print("     -> version 030")
        if   int(AIM_ints[17]) == 131074:   format = 'short';#          print("     -> format "+format)
        elif int(AIM_ints[17]) == 65537:    format = 'char';#           print("     -> format "+format)
        elif int(AIM_ints[17]) == 1376257:  format = 'bin compressed';# print("     -> format "+format+" not supported! Exiting!")#; exit(1)
        else:                               format = 'unknown';#        print("     -> format "+format+"! Exiting!")#; exit(1)
    
        # Reading the header
        header = f.read( AIM_ints[8] )
        header = header.decode("ISO-8859-1")
        header_len = len(header)+280
        # Dimensions of the AIM file, and coordinates of the image in space
        extents    = ( 0, AIM_ints[24]-1, 0, AIM_ints[26]-1, 0, AIM_ints[28]-1 )
        coords = (AIM_ints[18], AIM_ints[20], AIM_ints[22] )
    f.close()
    # collect data from header if existing
    header = re.sub("(?i) +"," ", header)
    header = header.split("\n")
    header.pop(0); header.pop(0); header.pop(0); header.pop(0)
    Scaling=None; Slope=None; Intercept=None
    # Extracting the required info from the AIM header and store in a dictionary
    
    for line in header: 
        if line.find("Original Creation-Date")>-1: Date  = \
            (int(line.split(" ")[2].split("-")[0]),(int(strptime(line.split(" ")[2].split("-")[1],'%b').tm_mon)),\
             int(line.split(" ")[2].split("-")[2]))
        if line.find("Orig-ISQ-Dim-p")>-1:         OrigDimp  = ( float(line.split(" ")[1]), float(line.split(" ")[2]), float(line.split(" ")[3]))
        if line.find("Orig-ISQ-Dim-um")>-1:        OrigDimum = ( float(line.split(" ")[1]))
        if line.find("Orig-GOBJ-Dim-p")>-1:        OrigDimp  = ( float(line.split(" ")[1]), float(line.split(" ")[2]), float(line.split(" ")[3]))
        if line.find("Orig-GOBJ-Dim-um")>-1:       OrigDimum = ( float(line.split(" ")[1]))
        if line.find("Patient Name")>-1:           Pat_ID    = str(line.split(" ")[-2])
        if line.find("Mu_Scaling")>-1:             Scaling   = float(line.split(" ")[-1])
        if line.find("Density: intercept")>-1:     Intercept = float(line.split(" ")[-1])
        if line.find("Density: slope")>-1:         Slope     = float(line.split(" ")[-1])
        if line.find("scale (el_size factor)")>-1: Factor    = ( float(line.split(" ")[-3]), float(line.split(" ")[-2]), float(line.split(" ")[-1]))
    
    vox = float(OrigDimum/(OrigDimp[0]*1000))
    if vox >= 0.081 and vox<= 0.083: XCT_gen = 1
    elif vox >= 0.06 and vox<= 0.061:XCT_gen = 2
    
    header_dic = {"coordinates":coords,"dimensions":(extents[1]+1,extents[3]+1,extents[5]+1),\
                  "slope":Slope,"intercept":Intercept,"XCT_gen":XCT_gen,\
                      "Mu_scaling":Scaling,"voxel":vox,"grid":OrigDimp,"Date":Date,"Patient_ID":Pat_ID}

    print(header_dic)
    # read AIM
    reader = vtk.vtkImageReader2()
    reader.SetFileName( fileINname )
    reader.SetDataByteOrderToLittleEndian()
    reader.SetFileDimensionality( 3 )
    reader.SetDataExtent( extents )
    reader.SetHeaderSize( header_len )
    if format == 'short':  reader.SetDataScalarTypeToShort()
    elif format == 'char': reader.SetDataScalarTypeToChar()
    reader.Update()
    imvtk = reader.GetOutput()
    return imvtk, header_dic

def aim2np(AIM):
    """
    Convert the AIM file to numpy array
    """
    imvtk, Header = AIMreader(AIM)
    Header
    aimnp = vtk2numpy(imvtk)
    return aimnp,Header