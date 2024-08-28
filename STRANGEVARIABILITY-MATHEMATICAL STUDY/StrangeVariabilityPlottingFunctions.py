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
    #Starts recording the time it takes to run the program entirely
    graph_time = []
    total_time_start = time.time()      


    #saving files/paths
    file_out = os.getcwd() + '/OutputFiles/CovOpt_Plots/'
    alpha_plot_out = os.getcwd() + '/OutputFiles/CovOpt_Plots/AlphaSorted/a='
    alpha_group_out = os.getcwd() + '/OutputFiles/Grouped Alpha/a='
    
    #Variable initialization
    tau = np.arange(0.1,5,0.1)
    ttau = -1 * tau
    prev_num = 0
    prev_alpha = 0
    num = 1

    #Colors for the grouping 
# =============================================================================
#     pickle = [0,0.2,0]
#     obs_colors = [0,51 - index/10,0]
#     It would be sik af if we can get the color of the observation to be determined based on the closeness of the midpoint of the interval
#       color = [0.8,0,0.8 * obs_alpha/alpha] Red initial
#       
# =============================================================================
    
################################################################################################################################################### 
    #Starts iterating through each of the dataframe rows
    for index, row in data.iterrows():
        time_start = time.time()
        
        #Makes sure that the I0's are not equal/ alpha = 1
        if row['I01'] == row['I02']:
            continue
        
        
        print('Plotting: ' + row['objectName'])
        
        #Inner Initialization
        alpha = row['alpha']
        inv_alpha = 1 / alpha
        Spectral_Cf1 = row['Spec_Cf1']
        Spectral_Cf2 = row['Spec_Cf2']
        #Spectral_value_Cf2Cf1 = row['Spec_Cf2Cf1']
        Sv1 = row['minCov1']
        Sv2 = row['minCov2']
        mintau = row['minOpt']  
        
       
        #Setting up plots for both subplots
        
        #Figure Size
        fig = plt.figure(figsize = (12,5))        
        
        #Sets title of plot
        if alpha_group == 'yes':
            if isinstance(row.obs_Cf1, np.ndarray):
                detect = 1
            else:
                detect = len(row.obs_Cf1)
            
            plt.suptitle(r' $\alpha$ = ' + str(alpha) + ' with ' + str(detect) + ' detections')
        else:
            plt.suptitle(row['objectName'] + r' $\alpha$ = ' + str(alpha))
        
        #Subplot 1 Settings
        plt.subplot(1,2,1)
        plt.title(r'$C_{f1}$' + ' vs. ' + r'$\tau$')
        plt.ylabel(r'$C_{f1}$')
        plt.xlabel(r'$\tau$')
        plt.ylim(0,1)
        plt.xlim(-0.95,5)
        
        #Subplot 2 Settings
        plt.subplot(1,2,2)
        plt.title(r'$C_{f2}$' + ' vs. ' + r'$\tau$')
        plt.ylabel(r'$C_{f2}$')
        plt.xlabel(r'$\tau$')
        plt.ylim(0,1)
        plt.xlim(-0.95,5)      
        
        #Not really necessary but I am too lazy to change the variable names
        subplt1 = 1
        subplt2 = 2
        
        #Sets the equation used to get the alpha value
        if row['I01'] < row['I02']:
            alpha_eqn = r"; $I_{o1} = \alpha I_{o2}$"
            
        else:
            alpha_eqn = r"; $I_{o2} = \alpha I_{o1}$"  
                
        #Starting solution sets of coverage fraction vs tau
        Cf = 0.1 #control coverage fraction value
            
        for index in range(1,11): #The numbers in range determine how many solution sets are shown
                
        
            plt.subplot(1,2,subplt1)
            #solving for Cf1
            Cf1 = alpha * Cf + (alpha - 1)/(np.exp(ttau)-1)
            plt.plot(tau, Cf1,label = r'$C_{f}$=' + str(Cf), color = [0.05,0.6,index/10.])
            

            plt.subplot(1,2,subplt2)
            #Solving for Cf2
            Cf2 = inv_alpha * Cf + (inv_alpha - 1)/(np.exp(ttau)-1)
            plt.plot(tau, Cf2, color = [0.05,0.6,index/10.])
                
            Cf = round(Cf + 0.1, 1)
            
        #LEGEND - DAIRY
        fig.legend(loc = 'outside upper left')
        
        if alpha_group == 'yes':
            #starts iterating through data frame to look into each indexes set of lists
            
            if len(row.obs_Cf1) == 0: #If the interval has no observational data, it will stop this loop and continue
                                 
                continue
                             
                 
            if isinstance(row.obs_Cf1, np.ndarray): #The dataframe saves a single list as an array instead of in a series, so this is one single detection
                    
                    plt.subplot(1,2,1)
                    plt.plot(tau, row.obs_Cf1, color = [0.8,0,0.8])

                    plt.subplot(1,2,2)
                    plt.plot(tau, row.obs_Cf2, color = [0.8,0,0.8]) 
            
            elif isinstance(row.obs_Cf1, pd.Series): #A list of lists is saved as a series in a dataframe. This is anycase where there is more than one obs. data in an interval
                    
                    
                    #iterates through each of the lists inside of the series
                    for k in range(len(row.obs_Cf1)): 

                        #Originally for the use of finding how close an obs data was to the plotted alpha value, but is unnessassary now.
                        if alpha > row.obs_alpha.iloc[k]:
                            alp_ratio = row.obs_alpha.iloc[k]/alpha
                            
                        else:
                            alp_ratio = alpha/row.obs_alpha.iloc[k]
                        
                        if k == 0:
                            print(alp_ratio)
                        
                        blue = 0.8 * np.abs(alpha - row.obs_alpha.iloc[k]) * 20
                        
                        #plots each of the obs data for the interval on Cf1
                        plt.subplot(1,2,1)    
                        plt.plot(tau, row.obs_Cf1.iloc[k], color = [0.8,0,blue])
                        
                        #plots each of the obs data for the interval on Cf2
                        plt.subplot(1,2,2) 
                        plt.plot(tau, row.obs_Cf2.iloc[k], color = [0.8,0,blue])
                       
    
        else: #If alpha_group == 'no'. This is the case for plotting each detection and is defaulted
            #Observational Data/Solution Cf1
            plt.subplot(1,2,1)
            plt.plot(tau,Spectral_Cf1,color = 'r')
            plt.plot([-0.95,5],[Sv1,Sv1],'r--')
            plt.text(2.8,Sv1-0.05,r'$C_{f1}(1-e^{- \tau})=$' + str(Sv1),color='r',fontsize=10)
       
            plt.subplot(1,2,2)
            #Observation Data/Solution Cf2
            plt.plot(tau,Spectral_Cf2,color = 'r')
            plt.plot([-0.95,5],[Sv2,Sv2], 'r--')
            plt.text(2.8,Sv2-0.05,r'$C_{f2}(1-e^{- \tau})=$' + str(Sv2),color='r',fontsize=10)
                    
            
            #Minimum Optical Depth
            plt.plot([mintau,mintau],[0,1],'r--')
            plt.text(mintau-1,-0.1,r"mintau= "+str(round(mintau,2)),color = 'r',fontsize=12)
                   
           
            #Minimum Optical Depth
            plt.plot([mintau,mintau],[0,1],'r--')
            plt.text(mintau-1,-0.1,r"mintau= "+str(round(mintau,2)),color = 'r',fontsize=12)
        
        
        #Plotting the Cf1 vs tau graph
        plt.subplot(1,2,1)
        plt.text(3,0.05,r"$\alpha =$"+str(alpha)+alpha_eqn,bbox=dict(facecolor='white', alpha=0.6))
    
    
    
        #Plotting the Cf2 vs tau graph
        plt.subplot(1,2,2)        
        plt.text(3,0.05,r"$\alpha =$"+str(alpha)+alpha_eqn,bbox=dict(facecolor='white', alpha=0.6))
    
        
        
        #Finds the time it took for each plot
        g_time = round(time.time() - time_start,2)
        graph_time.append(g_time)
        print('Plot took ' + str(g_time) + ' seconds')
        
        #Determines if the alpha value has been seen before (assuming the dataframe is in alpha_numerical order)
        if prev_alpha == alpha:
            prev_num += 1
        else:
            prev_num = 0
        
        
        prev_alpha = alpha
        
        #Saving figures 
        if row['save_fig'] == 'yes':            
            
            if alpha_group == 'yes':
                plt.savefig(alpha_group_out + str(alpha) + '.png', dpi = 350)
                plt.show()
            else:
                plt.savefig(file_out + row.objectName + '.png',dpi = 350)
                plt.savefig(alpha_plot_out + str(alpha) + '(' +str(prev_num) + ').png',dpi = 350)
        
                plt.clf()
        

        #Just to see progress and how many detections are left at a time
        print(str(len(data) - num) + ' detections left')
        num += 1
        
    #Time it takes for the entire plotting function to finish
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
    
    #We might need to add a variable that will let us change the interval size due to detection rate at 
    #specific intervals. IE, if only one in a range then change the interval
    
    
    #Setting up data
    template_file = 'CSV Files/Alpha_Template.pkl'
    
    data_alpha = data.sort_values(by = ['alpha'])
    template = pd.read_pickle(template_file)
    epsilon = 0.05
    
    print('Total Number of Detections: ' + str(len(data_alpha)))
    
    #initiallizing Variables/Lists
    total_Cf1_obs = [] #The Total List of lists of Cf1 observational Data
    total_Cf2_obs = [] #the Total List of lists of Cf2 Observational Data
    obs_alpha = []
    
    
    #Iterating through the template dataframe rows 
    for index, row in template.iterrows():
        
        #Creating an interval and conditions fitting that interval
        #print('alpha interval: [' + str(round(row['alpha'] - 0.05,2)) + ',' + str(round(row['alpha'] + 0.05,2)) + ')')
        cond1 = data_alpha['alpha'] >= round(row['alpha'] - epsilon,2)
        cond2 = data_alpha['alpha'] < round(row['alpha'] + epsilon,2)
        
        #Making a dataframe that contains all of the detections within the conditions/interval
        data_iv = data_alpha.loc[(cond1) & (cond2)]

        #this is where you would maybe start with checking the number of obs_data and start seperating by tens.
        for k in range(len(data_iv.Spec_Cf1)):
            if k % 10 == 0:
                print('something')
            #This is where you would add an index into the data frame and you can copy original/calculated data, then put 10 observational datas into its respective column. 
            #Research:
            #Inserting indexes/rows into a data frame and copying data from the row into it. 
            #Pandas df.squeeze() might be usefull again

        #Adds each of the spec_cf1/Cf2 into list 
        total_Cf1_obs.append(data_iv.Spec_Cf1.squeeze())
        total_Cf2_obs.append(data_iv.Spec_Cf2.squeeze())
        #print(data_iv.alpha.squeeze())
        alpha_sqz = data_iv.alpha.squeeze()
        obs_alpha.append(alpha_sqz)
        
       
          
        
    #Saving list of series into the template dataframe that will go to the graphing function
    #print(obs_alpha)
    template['obs_Cf1'] = total_Cf1_obs
    template['obs_Cf2'] = total_Cf2_obs
    template['obs_alpha'] = obs_alpha

    template.to_pickle(os.getcwd() + '/CSV Files/templateDF.pkl')
    return template

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

#Everything above has been created with respect to ONLY Coverage Fraction vs Optical Depth. Below is the start of Optical Depth vs Coverage Fraction

def tau_Cf_grapher(data, alpha_group = 'no'):



    for index, row in data.iterrows():
        #Inner Initialization
        alpha = row['alpha']
        inv_alpha = 1 / alpha
        Spectral_Cf1 = row['Spec_Cf1']
        Spectral_Cf2 = row['Spec_Cf2']
        #Spectral_value_Cf2Cf1 = row['Spec_Cf2Cf1']
        Sv1 = row['minCov1']
        Sv2 = row['minCov2']
        mintau = row['minOpt']  

        Cf = np.arange(0.1,1.1,0.02)

        for i in range(1,52):
            tau1 =  
        
