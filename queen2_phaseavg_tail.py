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

    #Initialize storage for PIV data folders
    queen2folders = {}

    #for each parent vid, get a list of the PIV data folders
    for vid in parentVids:
        
        #set filepath to the video's queen2 results folder
        path = directory + '/' + vid + '/results_3cyc'
        
        #add to the dictionary
        queen2folders[vid] = path

    #return the dictionary
    return queen2folders   
    
    
    
import pandas as pd
    
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
    PFx_file = folder + '/queen2_rodsim_dT10ms_xforce_fixed.xlsx'
    PFy_file = folder + '/queen2_rodsim_dT10ms_yforce_fixed.xlsx'
    PTz_file = folder + '/queen2_rodsim_dT10ms_torque_fixed.xlsx'
    
    #Read in Fx data
    PFx = pd.read_excel(PFx_file, sheetname = 'Sheet1', header = None)
    #Data is natively in one row.  Need it in one column, so transpose.
    PFx = PFx.T
    #Label the column something informative so it can be retrieved later
    PFx.columns = ['PFx']
    #Scale the pressure code data from N/m to m by multiplying by foil depth.
    #Sign is assigned to match the definition of +/- given by the Flapper.
    PFx['PFx'] = PFx['PFx']*-1.
    
    #repeat above for Fy
    PFy = pd.read_excel(PFy_file, sheetname = 'Sheet1', header = None)
    PFy = PFy.T
    PFy.columns = ['PFy']
    PFy['PFy'] = PFy['PFy']*1.
    
    #repeat above for Tz
    PTz = pd.read_excel(PTz_file, sheetname = 'Sheet1', header = None)
    PTz = PTz.T
    PTz.columns = ['PTz']
    PTz['PTz'] = PTz['PTz']*-1000.
    
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
    
    
    
    

def phaseAvg_queen2(path, df):
    """
    Finds phase average and error of pressure code results.

    Input:
    -path = filepath for the results folder
    -df = dataframe with time, PFx, PFy, and PTz columns from queen2.
    
    
    Output:
    Saves out a spreadsheet with all the queen2 data phase averaged.
    """
    #load useful package
    import numpy as np
    
    #find out the type of test that this is.  In the path for the results folder, 
    #test type is the name of the 3rd subfolder (as currently organized)
    testType = path.split('/')[3]
    
    #if the test type is NOT static,
    if testType != 'static':
        
        #find out the name of the parent video - it's the name of the 4th subfolder
        vid = path.split('/')[4]
        
        #Extract the flapping frequency from the video name
        freq = float(vid.split('_')[2][0:3])
        
        #calculate the period for flapping
        period = 1/freq
        #because the pressure code samples at 100Hz, multiply by 100
        p = period*100
        #p equals the number of dataframe rows per 1 motion cycle
        
        #round p to an integer
        p = int(round(p))
    
    #in the static case:
    else:
        #The static case was divided up into 1000-frame chunks since there aren't
        #real "motion cycles"
        p = 100
        
    #Create a dictionary to store averages in
    avgs = {}
    
    #Labels to be used for the columns
    columns = ['PFx','PFy','PTz']
        
    #for each column to calculate an average for:
    for col in columns:
            
        #initialize a list where each time-point's phase-average will be stored
        averageList = []
        #initialize a list where each time-point's standard deviation will be stored
        stdevList = []
            
        #look at the ith row in a cycle
        for i in range(0,p):
            
            #initialize storage for the values at this time point in each cycle
            values = []
            
            #initialize storage for the sum of corresponding values
            total = 0
            
            #look at the ith row in the nth cycle (all in turn)
            for n in range(0, 3):
                #Because the pressure code cuts off the last couple of frames
                #during calculations, initialize a trigger to indicate we've run
                #out of rows
                out_of_rows = False
            
                #Get the value in the ith row of the nth cycle
            
                #try to get the value using the row index
                try:
                    #add the current value to the sum of corresponding values
                    total += df[col][i+n*p]
                    #add the current value to the list of corresponding values
                    values.append(df[col][i+n*p])
                    
                except:
                    #if we run out of rows in the 3rd cycle, trigger the "out_of_rows"
                    #flag, and continue
                    out_of_rows = True
                
            #if we ran out rows, only 2 values are being averaged instead of three.
            if out_of_rows:
                newAvg = total/2.
                
            #if we had 3 values, average all three.
            else:
                newAvg = total/3.

            
            #add the average to the collection
            averageList.append(newAvg)
            
            #find the standard deviation of the corresponding values
            #make a list of square errors
            sqerrors = []
                
            #calculate the standard deviation
            for v in values:
                sqerrors.append(abs(v-newAvg)**2.)
            sumsqerrors = np.sum(sqerrors)
            stdev = np.sqrt(sumsqerrors/(len(values)-1))
                
            #add standard deviation to the list
            stdevList.append(stdev)          
            
        #add the phased-average column to the dictionary    
        avgs[col] = averageList
        #make a title for the standard deviation column
        errorcol = str(col) + '_std'
        #add the standard-deviation column to the dictionary
        avgs[errorcol] = stdevList
        
    #make a time sequence
    time = np.array(range(0,p))
    time = time/100.
    avgs['time']=time

    #convert dictionary to a dataframe
    avgs=pd.DataFrame(avgs, columns=avgs.keys())
        
    #save to an excel file
    finalpath = path + '/queen2_results_phaseAvg_fixed_fixedstdev_axisfix.xlsx'
    avgs.to_excel(finalpath)
    
    #confirm done
    print 'Saved file as ' + finalpath.split('/')[-1]
        



#list the test types conducted
tests = ['tail_gap','tail_mid']

#where the test data is stored
directory = 'G:/Nonunif_performance_PIV/Pressure code analysis/'

#for each test,
for t in tests:
    #make the full path for the test subfolder
    path = directory + t
    #find all the folders with test results
    folders = queen2_results_folders(path)
    
    #for each of those folders,
    for f in folders.keys():
        #retrieve the pressure code data and assemble into one table
        data = queen2_results(folders[f])
        
        #phase average the data
        phaseAvg_queen2(folders[f], data)
        