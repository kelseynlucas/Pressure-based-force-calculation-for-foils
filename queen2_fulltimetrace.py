# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 15:39:16 2015

@author: Kelsey





"""



#Load some built-in file-reading functions
from os import listdir
from os.path import isfile, join


def queen2_results_folders(directory):
    """
    Reads a file path (directory) - for a test type - and finds the folders 
    within the directory with queen2 results.   Returns a list of the data folders.
    
    Input
    
    -directory - file path for a test type containing folders with queen2
        data.  Note that all '\'s in the default file path must be 
        changed to '/'. No ending slash
        
    """
    
    #If have a subfolder (not a file), get name and add to list of subfolders
    parentVids = [f for f in listdir(directory) if not isfile(join(directory,f))]

    #Initialize storage for data folders with queen2 results
    queen2folders = {}

    #for each parent vid, get the queen2 results folder
    for vid in parentVids:
        
        #set filepath to the video's queen2 results folder
        path = directory + '/' + vid + '/results_3cyc'
        
        #add to the dictionary
        queen2folders[vid] = path

    #return the dictionary
    return queen2folders   
    
    
#import the dataframe package    
import pandas as pd

#note this function was copied from the queen2_phaseavg.py script, since it is
#also needed in this script
def queen2_results(folder):
    """
    Reads in a folder path where queen2 data are stored.  Retrieves the queen2
    data files and combines to a dataframe.  Applies converstion from N/m and
    Nm/m to N and Nmm.
    
    Input
    -folder - folder where queen2 data are stored.  All '\'s in the default
        file path must be changed to '/'.
    
    """
    #set file name and path for Fx, Fy, and Tz data
    PFx_file = folder + '/queen2_rodsim_dT10ms_xforce.xlsx'
    PFy_file = folder + '/queen2_rodsim_dT10ms_yforce.xlsx'
    PTz_file = folder + '/queen2_rodsim_dT10ms_torque.xlsx'
    
    #Read in Fx data
    PFx = pd.read_excel(PFx_file, sheetname = 'Sheet1', header = None)
    #Data is natively in one row.  Need it in one column, so transpose.
    PFx = PFx.T
    #Label the column something informative so it can be retrieved later
    PFx.columns = ['PFx']
    #Scale the pressure code data from N/m to m by multiplying by foil depth.
    #Sign is assigned to match the definition of +/- given by the Flapper.
    PFx['PFx'] = PFx['PFx']*-0.08
    
    #repeat above for Fy
    PFy = pd.read_excel(PFy_file, sheetname = 'Sheet1', header = None)
    PFy = PFy.T
    PFy.columns = ['PFy']
    PFy['PFy'] = PFy['PFy']*0.08
    
    #repeat above for Tz
    PTz = pd.read_excel(PTz_file, sheetname = 'Sheet1', header = None)
    PTz = PTz.T
    PTz.columns = ['PTz']
    PTz['PTz'] = PTz['PTz']*-80.
    
    #Create a time sequence
    time = range(0, len(PFx['PFx']))
    #scale the time to ms
    time = [t/100. for t in time]
    
    #Store time sequence in final dataframe
    combo = pd.DataFrame(time, columns=['time'])
    
    #Add Fx, Fy, and Tz to the final dataframe
    combo['PFx'] = PFx['PFx']
    combo['PFy'] = PFy['PFy']
    combo['PTz'] = PTz['PTz']
    
    #return the Fx, Fy, Tz, and time data for the test
    return combo
    
    

        



#make a list of all the test types
tests = ['3D_2cm','3Dmid','3Dedge','dynamic','static']

#ID the directory containing the results for ALL the tests
directory = 'G:/Nonunif_performance_PIV/Pressure code analysis/'

#for each test,
for t in tests:
    #make the full path for the test's data folder
    path = directory + t
    #ID all the folders containing queen2 results
    folders = queen2_results_folders(path)
    
    #for each folder with queen2 data,
    for f in folders.keys():
        #retrieve the queen2 data and store it in one dataframe
        data = queen2_results(folders[f])
        
        #save out the data combined into one file
        finalpath = folders[f] + '/queen2_results_fulltimetrace_axisfix.xlsx'
        data.to_excel(finalpath)
    
        #confirm done
        print 'Saved ' + str(f) + 'file as ' + finalpath.split('/')[-1]