# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 14:05:02 2023

@author: Taylor Gibbons with contributions from Alex Salley
"""

import numpy as np
import os
#import scipy.signal as sc
#This code is a quick way to find the values of the 'weird' variability based on the limits that are user input.

quas_name = 'J123830.71+480003.4' #name of the quasar
quas_direct = quas_name +'/' 
    
quas_direct = os.getcwd() + '/DATA_VARIABILITY/' + quas_direct #folder for the object spectra
    
test = 'n' #prints testing statements

##################################################################################################################################################################################
#These are the limits of finding the values going from left to right
trough_left = 5250 #beginning limit for finding the values
in_trough_left = 1 #The first limit inside of the trough INDEX
min_index = 2
in_trough_right = 3 #the second limit inside of the trough INDEX
trough_right = 5450 #the last limit outside of the trough

    
##################################################################################################################################################################################

#Functions:
    


def findBoundries (flux): #This finds the maximum value of the sides of the trough, sets them as limits and finds the minimum value
                          #then returns the indexes: in_left,in_right, trough_min_index
    
    trough_min_index = np.argmin(flux)
    
    
    in_left = np.argmax(flux[:trough_min_index])
    
    in_right = np.argmax(flux[trough_min_index:]) + trough_min_index
    

    if test == 'y':
        print('in_trough_left= '+str(in_left))
        print('in_trough_right= '+str(in_right))
        print('min_index= ' + str(trough_min_index))
        
    return in_left,in_right, trough_min_index

def fcalc_avg(flux): #This will take in the flux array and return the value of the average from the left side of the array
                     # to the in_trough_right index.
                     
    add_flux = 0
    
    counter = in_trough_left + ((len(flux)-1)- in_trough_right) #number of elements exluding the trough

    for i in range(in_trough_left):
        add_flux += flux[i]
        
    for k in range(in_trough_right, len(flux)):
        add_flux += flux[k]
    
    avg_flux =  add_flux / counter
    
    if test == 'y':
        print('counter: ' + str(counter))
        print('add_flux: ' + str(add_flux))
    return avg_flux


    
#End of the Functions
##################################################################################################################################################################################
#Reads in the dered files and makes an array

dered = []
    
files = os.listdir(quas_direct)
    
for i in range(len(files)):
    if files[i].endswith('dered.dr16') :
                    dered.append(files[i])
                    
##################################################################################################################################################################################      

text = []       #A 2d array with [row][column]
waves = []      #An array with [wavelengths]
flux = []       #An array with [flux]
        
avg_min = 0        

for k in range(len(dered)): 
    text = np.loadtxt(quas_direct + dered[k],) #makes the array filled with the numbers of each line of the text
                                                # for each spectra. This will have two rounds, one for each spec.
                                                 
    waves = text[:,0] #set all of the text wavelength numbers to the waves array
    flux = text[:,1] #set all of the text flux numbers to the flux array
 
    bitmask = np.where((waves > trough_left) & (waves < trough_right)) #placeholder index for the indexes within the left/right limits
    waves = waves[bitmask] #reindexes the waves to only have the bitmask elements
    flux = flux[bitmask]    #reindexes the flux to only have the bitmask elements
    
    in_trough_left,in_trough_right, min_index = findBoundries(flux) #Finds the intrough boundries and index at minimum value
    
    avg_min += flux[min_index]
    
    
    Iok = round(fcalc_avg(flux),2)
    print('I_{o' + str(k+1) + '} = ' + str(Iok))
    
flux_min = round(avg_min / 2,2)   
print('I1 ~~ I2 = ' + str(flux_min))   




