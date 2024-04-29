# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 14:42:21 2023

@author: heave
"""
list_09 = [50,100,2950,45,1,0,105,92,76,48]
list_1 = [4,9000,59,100,1,89,105,910,97,1092,92,75,76]
data = [900,4500,9752,61789,4712,1752,972,14,2367,1249,2186,24456,541]

print(list_1[3])
al = [10,20,30,55,78,99,105,42,65] # 1.0 velocity list
nal = [20,78,65] # 0.9 velocity list
matched = [] # matching indexies

for i in range(len(list_1)): # Runs through all the 1.0 values and checks to see if they match a 0.9 value. If they do, save the index
    for j in range(len(list_09)):
        if list_09[j] == list_1[i]:
            matched.append(i)

others = [75,500,27,9000,78,90,45,677,47] # Other 1.0 data (like EW)
sort = []  # New list containing 1.0 data that corresponds to 0.9 detections (like EW)

for i in range(len(matched)):
    sort.append(others[matched[i]])
    # Repeat for each value (EW, depth, velocity, etc)
    
    
