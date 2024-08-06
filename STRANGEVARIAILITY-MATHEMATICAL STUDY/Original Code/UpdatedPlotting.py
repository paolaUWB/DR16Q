# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 12:03:11 2024

@author: Taylor Gibbons
"""



# =============================================================================
# Notes for future:
# 7/22/2024:
#   As of the end of the day, I have finished a fairly easy version of the Cf_grapher
#   function, however needs to be detail, so those will need to add. Also
#   There seems to be a problem with some of the spectral values not being even 
#   close to the solution sets. 
# =============================================================================



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import os

def Updated_Cf_grapher(data, graph_num= 0):
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

    graph_time = []
    total_time_start = time.time()    
    
    
    #alpha, mintau, Sv1, Sv2, Spectral_value_Cf2Cf1, Spectral_Cf1, Spectral_Cf2 =svc.min_value_finder(I01, I02, I1, I2, prev_num)
    tau = np.arange(0.1,5,0.1)
    ttau = -1 * tau
    
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
            alpha_eqn = 'I01 = aI02'
        else:
            subplt1 = 2
            subplt2 = 1
            alpha_eqn = 'I02 = aI01'
        
        
        
        #Plotting the Cf1 vs tau graph
        plt.subplot(1,2,1)
        
        plt.ylabel(r'$C_{f1}$')
        plt.xlabel(r'$\tau$')
        
        plt.ylim(0,1)
        plt.xlim(-0.95,5)
        
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
        
        

        
        #plt.show()
        
        
        if index != 0:
            if data.iloc[index - 1]['alpha'] == alpha:
                print('Same Alpha Value Detected')
                prev_num += 1
            else:
                prev_num = 0
        
        if row['save_fig'] == 'yes':
            plt.savefig(os.getcwd() + '/alpha_plots/a=' + str(alpha) + '(' + str(prev_num) + ')' + '.png',dpi = 350)
        #     plt.savefig(os.getcwd() + '/alpha_plots/a=' + str(row['alpha']) + str(num) + '.png',dpi = 350)            
        plt.clf()
            
        print(str(len(data) - num) + ' detections left')
        num += 1
        
        
        g_time = round(time.time() - time_start,2)
        graph_time.append(g_time)
        print('Plot took ' + str(g_time) + ' seconds')        
        
    avg_time = np.mean(graph_time)
    print('Average Plotting Time: ' + str(avg_time) + ' seconds\nTotal Plotting Time: ' + 
          str(round(time.time() - total_time_start,2)) + ' seconds')
    
    data['t_t_plot'] = graph_time
    
    data.to_csv('CovOpt_wTime.csv')
    

        
    