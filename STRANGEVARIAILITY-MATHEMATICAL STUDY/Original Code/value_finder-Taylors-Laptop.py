# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 14:05:02 2023

@author: Taylor Gibbons
"""

import numpy as np
import os
import scipy.signal as sc
#This code is a quick way to find the values of the 'weird' variability based on the limits that are user input.

quas_name = 'J123830.71+480003.4' #name of the quasar
quas_direct = quas_name +'/' 
    
quas_direct = os.getcwd() + '/DATA_VARIABILITY/' + quas_direct #folder for the object spectra
    
    #print(quas_direct) #testing
    
##################################################################################################################################################################################
#These are the limits of finding the values going from left to right
trough_left = 5250 #beginning limit for finding the values
in_trough_left = 1 #The first limit inside of the trough
in_trough_right = 2 #the second limit inside of the trough
trough_right = 5400 #the last limit outside of the trough
    
##################################################################################################################################################################################

#Functions:
    
def fcalc_avg(flux): #This will take in the flux array and return the value of the average.
    
    add_flux = 0
    counter = 0
    for i in range(len(flux)):
        add_flux += add_flux + flux[i]
        counter +=  1
        print (counter)
    avg_flux =  add_flux / counter
    return avg_flux


def findBoundries (flux, k):
    #trough_min = np.min(text[:,1])
    trough_min_index = np.argmin(flux)
    in_trough_left = np.argmax(flux[:trough_min_index])
    print('Spectra #'+str(k+1) + ' trough_left index: ' + str(in_trough_left))
    
    in_trough_right = np.argmax(flux[trough_min_index:])
    print('Spectra #'+str(k+1) + ' trough_right index: ' +str(in_trough_right))
    return 


'''
def reindex(reindex, n):
    tl = trough_left
    itl = in_trough_left_limit
    itr = in_trough_right_limit
    tr = trough_right
    indexed = []
    
    if n == 1: #reindexs for before the trough (zone 1)
        indexed = np.where(tl < reindex < itl)
    elif n == 2: #The case where we are reindexing based on the limits given for the trough (zone 2)
        indexed = np.where(itl < reindex < itr)
    elif n == 3: #reindexs for after the trough (zone 3)
        indexed = np.where(itr < reindex < tr)
    else:
        return reindex
    
    return indexed
'''    
    
#End of the Functions
##################################################################################################################################################################################

dered = []
    
files = os.listdir(quas_direct)
    
for i in range(len(files)):
    if files[i].endswith('dered.dr16') :
                    dered.append(files[i])
                    
        
text = []
waves = []
flux = []

#This for loop fills the text array with the information of EACH line. ie, the first line of numbers is [0] and
#then it will repeat with the next document. So, the calculations might need to be inside of this for loop and
#will be calculating Io1 and Io2 for k + 1. 
#Process (for now):
    #1. Read in the text from the dered files found in the dered[].
    #2. Takes the text and will separate (?) the numbers into more arrays, ie. wave[],flux[],....
    #3. Once seperate, reindex the arrays based on where the wavelengths will start: ie, if the trough_left_limit
        #is somewhere in the middle of the array, it will then become index = 0. and the same with the other lims.
    #4. Reindex the middle of the arrays for the middle limits
    #5. With the newly updated arrays, start with the calculations. Might be able to make a def to have it
        #calculate the values?
        
Ion = []    
for k in range(len(dered)): 
    text = np.loadtxt(quas_direct + dered[k],) #makes the array filled with the numbers of each line of the text
                                                # for each spectra. This will have two rounds, one for each spec.
                                                 


#    n_text = np.where( (trough_left < text[:,0]) & ( text[:,0]< trough_right ))
#    print('n_text[0,0] (' + str(n_text[0,0]))c
    waves = text[:,0]
    flux = text[:,1] #set all of the text flux numbers to the flux array
 
    bitmask = np.where((waves > trough_left) & (waves < trough_right))
    waves = waves[bitmask]
    flux = flux[bitmask]
    
    findBoundries(flux,k)
    Ion[k] = fcalc_avg(flux[trough_left:in_trough_left])
    Ion[k] += fcalc_avg(flux[in_trough_right: trough_right])
    Ion[k] = Ion[k]/2
    print (Ion[k])
    #testing the flux array and making sure that it is all of the numbers we want
    #for i in range(10):
     #   print('In spec #' + str(k+1) +': flux[] ('+ str(text[i,1]) +') = ' + str(flux[i]))

    #Iok = fcalc_avg(flux)
    #print('I_{o' + str(k+1) + '} = ' + str(Iok))
    
    
    
'''

print('number of lines in (n_text): ' + str(len(n_text)))
print('n_text[0]='+ str(n_text[0]))
'''




