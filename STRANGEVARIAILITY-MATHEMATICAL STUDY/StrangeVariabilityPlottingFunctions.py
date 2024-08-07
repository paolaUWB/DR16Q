# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 14:54:39 2024

@author: Taylor Gibbons
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
#import StrangeVariabilityCalculationFunctions as svc
import time

def visualize(data,xaxis,yaxis, graph_num = 0): #In Progress, and more of a test than an actual function to be used
    """
    
    """
    
    number = 0
    if graph_num == 1 or graph_num == 0:
        plt.figure(number)
        plt.scatter(data[xaxis],data[yaxis])
        plt.title('Minimum Optical Depth by Number of Days Between Observations')
        plt.xlabel('Number of Days between Observations')
        plt.ylabel('Minimum Optical Depth')
        number += 1
    if graph_num == 2 or graph_num == 0:
        plt.figure(number)
        plt.scatter(data['days'],data['new_alpha'])
        plt.title('New Alpha vs Number of Days Between Obervations')
        plt.ylabel('New alpha')
        plt.xlabel('Number of Days Between Observations')
        number += 1
    if graph_num == 3 or graph_num == 0:
        plt.figure(number)
        plt.scatter(data[''])
    if graph_num == 4 or graph_num == 0:
        table = pd.pivot_table(data, values = ['new_alpha','minOpt','minCov1','minCov2'],index = ['objectName'], aggfunc= 'mean')
        xdata = table['new_alpha']
        xbins = np.array([0.25,0.5,0.75,1])

        style = {'facecolor': 'none', 'edgecolor': 'C0','linewidth':3}

        fig,ax = plt.subplots()

        ax.hist(xdata,bins = xbins, **style)
        ax.set_ylabel('Number of Objects with new alpha values')
        ax.set_xlabel('New Alpha')
        
        
        
        
def Cf_tau_grapher(data, alpha_group = 'no'):
    """
    

    Parameters
    ----------
    data : DataFrame
        A dataframe that contains all of the calculated values and initial values
        from the detection code. This requires:
            I01, I02, I1, I2, alpha, minOpt, minCov1 (spectral values 1, sv1), 
            minCov2 (Spectral_value_2, sv2), Spec_Cf1, Spec_Cf2, Spec_Cf2Cf1 
    graph_num : int, optional
        DESCRIPTION. The default is 0. This input changes based on the graph you are
        looking to get

    Returns
    -------
    None.

    """
    
    #Whenever alpha_group is yes, itll take the template dataframe that we pass into it
    # then will plot the basic solutions of the detections, but none of the dotted lines
    #
    #
    #if alpha_group == yes:
    #   then it will prevent the plotting of the dotted lines (minimum values) and start plotting all of the 
    #   lists within the list of the observational data from the template df. ie. 
# =============================================================================
#             if alpha_group == 'yes':
#                 for index in range(len(row.Obs)):
#                     plt.plot(tau, row['Obs'][index])
#             else:
#                 plt.plot(dotted lines and all comes with it)
# =============================================================================
    




    tau = np.arange(0.1,5,0.1)
    ttau = -1 * tau


    template = data #for now

    if alpha_group == 'yes':
        for i in range(len(template)):
            temp_cf1 = template.Obs_Cf1[i]
            for k in range(len(temp_cf1)):
#                plt.plot(tau, temp_cf1[k], color=[0,0.5,1-k])
                plt.plot(tau, template.Obs_Cf1[i][k], color = [0,0.5,1-i] )
        
        for k in range(len(template)):
            plt.plot(tau, template.Spec_Cf1[k], 'r')




                  
        
        
    
    graph_time = []
    total_time_start = time.time()    
    
    
    #alpha, mintau, Sv1, Sv2, Spectral_value_Cf2Cf1, Spectral_Cf1, Spectral_Cf2 =svc.min_value_finder(I01, I02, I1, I2, prev_num)

    
    num = 1
    prev_num = 0
    #Goes through each of the 
    for index, row in data.iterrows():
        time_start = time.time()
        
        if row['I01'] == row['I02']:
            continue
        
        
        print('Plotting: ' + row['objectName'])
        
        #Initializing
        alpha = row['alpha']
        inv_alpha = 1 / alpha
        Spectral_Cf1 = row['Spec_Cf1']
        Spectral_Cf2 = row['Spec_Cf2']
        #Spectral_value_Cf2Cf1 = row['Spec_Cf2Cf1']
        Sv1 = row['minCov1']
        Sv2 = row['minCov2']
        mintau = row['minOpt']  
        
        
        #Setting up plots
        fig = plt.figure(figsize = (12,5))
        
        plt.suptitle(row['objectName'] + r' $\alpha$ = ' + str(alpha))
        
        plt.subplot(1,2,1)
        plt.title(r'$C_{f1}$' + ' vs. ' + r'$\tau$')
        
        plt.subplot(1,2,2)
        plt.title(r'$C_{f2}$' + ' vs. ' + r'$\tau$')
        
        
        if row['I01'] < row['I02']:
            subplt1 = 1
            subplt2 = 2
            alpha_eqn = r"; $I_{o1} = \alpha I_{o2}$"
        else:
            subplt1 = 2
            subplt2 = 1
            alpha_eqn = r"; $I_{o2} = \alpha I_{o1}$"
        
        
        
        #Plotting the Cf1 vs tau graph
        plt.subplot(1,2,1)
        
        plt.ylabel(r'$C_{f1}$')
        plt.xlabel(r'$\tau$')
        
        plt.ylim(0,1)
        plt.xlim(-0.95,5)
        
        plt.text(3,0.05,r"$\alpha =$"+str(alpha)+alpha_eqn,bbox=dict(facecolor='white', alpha=0.7))
        
        #Observational Data/Solution Cf1
        plt.plot(tau,Spectral_Cf1,color = 'r')
        plt.plot([-0.95,5],[Sv1,Sv1],'r--')
        plt.text(2.8,Sv1-0.05,r'$C_{f1}(1-e^{- \tau})=$' + str(Sv1),color='r',fontsize=10)
        
        
        #Minimum Optical Depth
        plt.plot([mintau,mintau],[0,1],'r--')
        plt.text(mintau-1,-0.1,r"mintau= "+str(round(mintau,2)),color = 'r',fontsize=12)
        
        
        #Plotting the Cf2 vs tau graph
        plt.subplot(1,2,2)
        
        plt.ylabel(r'$C_{f2}$')
        plt.xlabel(r'$\tau$')
        
        plt.ylim(0,1)
        plt.xlim(-0.95,5)
        
        plt.text(3,0.05,r"$\alpha =$"+str(alpha)+alpha_eqn,bbox=dict(facecolor='white', alpha=0.7))
        
        #Observation Data/Solution
        plt.plot(tau,Spectral_Cf2,color = 'r')
        plt.plot([-0.95,5],[Sv2,Sv2], 'r--')
        plt.text(2.8,Sv2-0.05,r'$C_{f2}(1-e^{- \tau})=$' + str(Sv2),color='r',fontsize=10)
        
        #Minimum Optical Depth
        plt.plot([mintau,mintau],[0,1],'r--')
        plt.text(mintau-1,-0.1,r"mintau= "+str(round(mintau,2)),color = 'r',fontsize=12)
        
        
        Cf = 0.1
        
        for index in range(1,11): #The numbers in range determine how many solution sets are shown
            
            plt.subplot(1,2,subplt1)
            Cf1 = alpha * Cf + (alpha - 1)/(np.exp(ttau)-1)
            
            plt.plot(tau, Cf1,label = r'$C_{f}$=' + str(Cf), color = [0.05,0.6,index/10.])
            
            plt.subplot(1,2,subplt2)
            Cf2 = inv_alpha * Cf + (inv_alpha - 1)/(np.exp(ttau)-1)
            
            plt.plot(tau, Cf2, color = [0.05,0.6,index/10.])
            
            Cf = round(Cf + 0.1, 1)
        
        
        fig.legend(loc = 'outside upper left')
        
        
        g_time = round(time.time() - time_start,2)
        graph_time.append(g_time)
        print('Plot took ' + str(g_time) + ' seconds')
        
        #plt.show()
        
        
        if index != 0:
            if data.iloc[index - 1]['alpha'] == alpha:
                prev_num += 1
            else:
                prev_num = 0
        
        if row['save_fig'] == 'yes':
            #plt.savefig(os.getcwd() + '/CovOpt_Plots/a=' + str(alpha) + '(' + str(prev_num) + ')' + '.png',dpi = 350)
            plt.savefig(os.getcwd() + 'OUTPUT_FILES/alpha_plots/a=' + str(alpha) + '(' +str(prev_num) + ').png',dpi = 350)            
        elif alpha_group == 'yes':
            plt.savefig(os.getcwd() + 'OUTPUT_FILES/Grouped Alpha/a=' + str(alpha) + '.png', dpi = 350)
        
        plt.clf()
        


        print(str(len(data) - num) + ' detections left')
        num += 1
        
        
    avg_time = np.mean(graph_time)
    print('Average Plotting Time: ' + str(avg_time) + ' seconds\nTotal Plotting Time: ' + 
          str(round(time.time() - total_time_start,2)) + ' seconds')
    
    
    
def alpha_group_grapher(data):
    """
    This function takes in a dataframe of detections with calculations and 
    groups the detections by alpha values, sends it to the Cf_tau_grapher to plot
    all of the observational data on intervaled plots (ie, alpha = 0.25 +- 0.05)

    Parameters
    ----------
    data : DATAFRAME
        A Dataframe of all the detections that has been run through
        StrangeVariabilityCoverageOptical to get all of the calculations.

    Returns
    -------
    None.
    
    """
    
    
    
    #Setting up data
    testing = 'no'
    template_file = 'CSV Files/Alpha_Template.csv'
    data_alpha = data.sort_values(by = ['alpha'])
    template = pd.read_csv(template_file)
    #print(template.iloc[0])
    print('Total Number of Detections: ' + str(len(data_alpha)))
    
    #initiallizing Variables/Lists
     
    total_Cf1_obs = [] #The Total List of lists of Cf1 observational Data
    total_Cf2_obs = [] #the Total List of lists of Cf2 Observational Data
    
    
    #Iterating through the template dataframe rows 
    for index, row in template.iterrows():
        #print(row['alpha'])
        print('alpha interval: [' + str(round(row['alpha'] - 0.05,2)) + ',' + str(round(row['alpha'] + 0.05,2)) + ')')
        cond1 = data_alpha['alpha'] >= round(row['alpha'] - 0.05,2)
        cond2 = data_alpha['alpha'] < round(row['alpha'] + 0.05,2)
        
        data_iv = data_alpha.loc[(cond1) & (cond2)] 
        
        data_iv.reset_index(inplace = True, drop = True) #data frame containing all rows in conditions
        
        print('Number of alpha values in interval: ' + str(len(data_iv)))
        
        temp_obs_Cf1 = data_iv.Spec_Cf1.squeeze()

        temp_obs_Cf2 = data_iv.Spec_Cf2.squeeze()

        #Adds each of the spec_cf1/Cf2 into list 
        total_Cf1_obs.append(temp_obs_Cf1)
        
        total_Cf2_obs.append(temp_obs_Cf2)

        if testing == 'yes':
           print(data_iv['alpha'])
           print('Length of obs_Cf1: ' + str(len(data_iv.Spec_Cf1)))
           print('Length of obs_Cf2: ' + str(len(data_iv.Spec_Cf2)))       
           print('data_iv_Cf1 index 0: ' + str(data_iv.Spec_Cf1.iloc[0]))  
           print('total_Cf1_obs Test [0]:')
           print(total_Cf1_obs[0])          
           
        print('Number of lists in template dataframe (should be ' + str(index + 1) + ' total): ' + str(len(total_Cf1_obs)))
        
        
    template['obs_Cf1'] = total_Cf1_obs
    template['obs_Cf2'] = total_Cf2_obs

    #Cf_tau_grapher(template,'yes')
    #return template,total_Cf1_obs
    
def alpha_Testing(I01, I02, I1, I2, Sv1, Sv2, Spectral_Cf1, Spectral_Cf2, graph):
    """

     I01 : float
         The average value of the flux around the trough of the first observation.
     I02 : float
         The average value of the flux around the trouhg of the second observation.
     I12 : float
         The average value of the flux at the bottom of the trough. At the point of 
         using this code, the values should be found and the quasars that have strange 
         variablity should be found.

    Returns
    -------
    None.

    
    OUTLINE:
       1. Take in values from observations
       2. Makes the two alpha values (at least one of both cases, maybe 2)
       3. graphs the two 

    """
    
    
    a12 = round(I01/I02,2)
    a21 = round(I02/I01,2)
    print('a12 = I01/I02 = ' + str(a12) + '\n'
          + 'a21 = I02/I01 = ' + str(a21))
    
    
    if a12 < 1: #1st Case where the second epoch is larger than the first 
        print('Proceeding with Case #1:')
    elif a21 < 1: #2nd Case where second epoch is less than the first
        print('Proceeding with Case #2')    
    
    if graph == 'yes':
        """
        Plan for this part of the function is to make is so that it ends up being
        a lot more consolidated. Ie, instead of the two graphs being made 
        separately, it will be made in one for loop with the graphs being made
        at generally the same time
        """
        #CASE 1: a12 < 1
        title = 'Alpha12 = ' + str(a12)+ ' Alpha21 = '+str(a21)
        tau_label = r'$\tau$'
        cf1_label = r'$C_{f1}$'
        cf2_label = r'$C_{f2}$'
        
        tau = np.arange(0.1,5.1,0.1)
        ttau=tau*(-1)  #= -tau
        
        plt.figure(figsize= (12,5))
        plt.suptitle(title)

        plt.subplot(1,2,1)
        plt.title(cf1_label + ' using ' + r'$\alpha_{12}$' + '=' + str(a12))
        plt.xlabel(r"$\tau$")
        plt.ylabel(r"$C_{f1}$")
        plt.ylim(0,1)
        plt.xlim(-0.95,5) 
        
        
        plt.subplot(1,2,2)
        plt.title(cf2_label + ' using ' + r'$\alpha_{12}$' + '=' + str(a12))
        plt.xlabel(tau_label)
        plt.ylabel(cf2_label)
        plt.ylim(0,1)
        plt.xlim(-0.95,5)        
        
        Cf=0.1
        
        for index in range(1,11):
    
            plt.subplot(1,2,1)
            Cf1 = a21*Cf + (a21-1)/(np.exp(ttau)-1)
            plt.plot(tau, Cf1 ,label = r'$Cf_2$='+str(Cf),color=[0.05,0.6,index/10.])

            plt.subplot(1,2,2)
            Cf2 = a12*Cf + (a12-1)/(np.exp(ttau)-1)
            plt.plot(tau, Cf2 ,label = cf2_label +'='+str(Cf),color=[0.05,0.6,index/10.])

    
            Cf=round(Cf+0.1,1)
            
        
        
        Cf1=0.1
        ttau=tau*(-1)  #= -tau
        
        #Save figure and clear to start new one
        plt.show()
        plt.savefig(os.getcwd() + '/AlphaChangePlots/' + title + ').png',dpi= 500)

