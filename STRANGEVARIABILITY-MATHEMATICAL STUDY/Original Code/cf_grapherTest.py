# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 21:56:00 2024

@author: tayta
"""


import numpy as np
import os
import matplotlib.pyplot as plt
import StrangeVariabilityCalculationFunctions as svc

def Cf_grapher(data, graph_num= 0):
    """
    IN PROGRESS
    
    
    INPUT: O
    
    
    
    """
    prev_num = 0
    
    
    #alpha, mintau, Sv1, Sv2, Spectral_value_Cf2Cf1, Spectral_Cf1, Spectral_Cf2 =svc.min_value_finder(I01, I02, I1, I2, prev_num)
    tau = np.arange(0.1,5,0.1)
    
    alpha_eqn = 'I01 = aI02'
    
    for index, row in data.iterrows():
        alpha = row['alpha']
        Spectral_Cf1 = row['Spectral_Cf1']
        Spectral_Cf2 = row['Spectral_Cf2']
        Spectral_value_Cf2Cf1 = row['Spectral_Cf2Cf1']
        Sv1 = row['Sv1']
        Sv2 = row['Sv2']
        mintau = row['minOpt']
        
        
        
        if graph_num == 1 or graph_num == 0:
            plt.figure(figsize= (18,5))
            plt.suptitle(row['objectName'] + '(' + str(prev_num) + ')')
            plt.subplot(1, 3, 1)
            plt.title(r'$C_{f1}$ vs. $\tau$')
            Cf2=0.1
            ttau=tau*(-1)  #= -tau
            for lala in range(1,11):
        
                Cf1 = (Cf2 - ((alpha-1)/(np.exp(ttau)-1)))/ alpha
                
                plt.plot(tau, Cf1 ,label = r'$Cf_2$='+str(Cf2),color=[0.05,0.6,lala/10.])
        
                plt.xlabel(r"$\tau$")
                
                plt.ylabel(r"$C_{f1}$")
                plt.ylim(0,1)
                plt.xlim(-0.95,5)
                      
                Cf2=round(Cf2+0.1,1)
                
                
            
            plt.plot(tau,Spectral_Cf1,color = 'r')
            plt.text(2.8,Sv1-0.05,r'$C_{f1}(1-e^{- \tau})=$' + str(Sv1),color='r',fontsize=10)
        
            #plt.legend(loc='upper left')
            #plt.text(3,0.19,r"$C_{f2}- \alpha C_{f1} = \frac{\alpha - 1}{e^{- \tau}-1}$",bbox=dict(facecolor='white', alpha=0.7))
            #plt.text(3,0.12,r"assume $\tau_1 = \tau_2$")
            #plt.text(3,0.05,r"$\alpha =$"+str(alpha)+r"; $I_{o1} = \alpha I_{o2}$",bbox=dict(facecolor='white', alpha=0.7))
        
            plt.plot([mintau,mintau],[0,1],'r--',) #RED verticle dotted Line of line (^). Biggest step is understanding  
            plt.text(mintau-1,-0.1,r"mintau= "+str(round(mintau,2)),color = 'r',fontsize=12)
        
        if graph_num == 1 or graph_num == 1221:
            plt.figure(figsize= (12,5))
            plt.suptitle(row['objectName'] + '(' + str(prev_num) + ')')
            plt.subplot(1, 2, 1)
            plt.title(r'$C_{f1}$ vs. $\tau$')
            Cf2=0.1
            ttau=tau*(-1)  #= -tau
            for lala in range(1,11):
        
                Cf1 = (Cf2 - ((alpha-1)/(np.exp(ttau)-1)))/ alpha
                
                plt.plot(tau, Cf1 ,label = r'$Cf_2$='+str(Cf2),color=[0.05,0.6,lala/10.])
        
                plt.xlabel(r"$\tau$")
                
                plt.ylabel(r"$C_{f1}$")
                plt.ylim(0,1)
                plt.xlim(-0.95,5)
                      
                Cf2=round(Cf2+0.1,1)
                
                
            plt.subplot(1,2,1)
            plt.plot(tau,Spectral_Cf1,color = 'r')
            plt.text(2.8,Sv1-0.05,r'$C_{f1}(1-e^{- \tau})=$' + str(Sv1),color='r',fontsize=10)
        
            #plt.legend(loc='upper left')
            #plt.text(3,0.19,r"$C_{f2}- \alpha C_{f1} = \frac{\alpha - 1}{e^{- \tau}-1}$",bbox=dict(facecolor='white', alpha=0.7))
            #plt.text(3,0.12,r"assume $\tau_1 = \tau_2$")
            #plt.text(3,0.05,r"$\alpha =$"+str(alpha)+r"; $I_{o1} = \alpha I_{o2}$",bbox=dict(facecolor='white', alpha=0.7))
        
            plt.plot([mintau,mintau],[0,1],'r--',) #RED verticle dotted Line of line (^). Biggest step is understanding  
            plt.text(mintau-1,-0.1,r"mintau= "+str(round(mintau,2)),color = 'r',fontsize=12)
            
            
            plt.subplot(1, 2, 2)
            plt.title(r'$C_{f2}$ vs. $\tau$')
            Cf1=0.1
            ttau=tau*(-1)  #= -tau
            for lala in range(1,11):
                 
                Cf2 = (Cf1 - ((alpha-1)/(np.exp(ttau)-1)))/ alpha
                
                plt.plot(tau, Cf2 ,label = r'$Cf_2$='+str(Cf2),color=[0.05,0.6,lala/10.])
        
                plt.xlabel(r"$\tau$")
                
                plt.ylabel(r"$C_{f1}$")
                plt.ylim(0,1)
                plt.xlim(-0.95,5)
                      
                Cf1=round(Cf1+0.1,1)
                
                
            
            plt.plot(tau,Spectral_Cf1,color = 'r')
            plt.text(2.8,Sv1-0.05,r'$C_{f1}(1-e^{- \tau})=$' + str(Sv1),color='r',fontsize=10)
        
            #plt.legend(loc='upper left')
            #plt.text(3,0.19,r"$C_{f2}- \alpha C_{f1} = \frac{\alpha - 1}{e^{- \tau}-1}$",bbox=dict(facecolor='white', alpha=0.7))
            #plt.text(3,0.12,r"assume $\tau_1 = \tau_2$")
            #plt.text(3,0.05,r"$\alpha =$"+str(alpha)+r"; $I_{o1} = \alpha I_{o2}$",bbox=dict(facecolor='white', alpha=0.7))
        
            plt.plot([mintau,mintau],[0,1],'r--',) #RED verticle dotted Line of line (^). Biggest step is understanding  
            plt.text(mintau-1,-0.1,r"mintau= "+str(round(mintau,2)),color = 'r',fontsize=12)
    
    ##########################################################################################################################################################
        if graph_num == 2 or graph_num == 0:
            tau= np.arange(0.1,5,0.1)
            ttau=tau*(-1)  #= -tau
        
            
            Cf1=0.1
            plt.subplot(1, 3, 2)
            
            plt.title(r'$C_{f2}$ vs. $\tau$')        
            plt.xlabel(r"$\tau$")
                
            plt.ylabel(r"$C_{f2}$")
            for lala in range(1,11):
        
                Cf2 =(alpha-1.)/(np.exp(ttau)-1.)+alpha*Cf1 
                
                
                plt.plot(tau, Cf2 , label = r'$Cf_1$='+str(Cf1),color=[0.05,0.6,lala/10.])
        
        
                plt.ylim(0,1)
                plt.xlim(-0.95,5)
        
                
                Cf1=round(Cf1+0.1,1)
        
            
        
            Spectral_Cf1 = Sv1/(1-np.exp(-tau)) #this will go specifically into the graph tau vs. Cf1.
            plt.plot(tau,Spectral_Cf2,color='r')
            plt.text(2.8,Sv2-0.05,r"$C_{f2}(1-e^{- \tau})=$"+str(Sv2),color='r',fontsize=10) 
        
            mintau=(-1)*np.log(1.-(Sv2/1.))  # From Cf1(1-exp(-tau))=1-I1/I01; since that is the value of tau for Cf=1 and the observational data
        
            plt.plot([mintau,mintau],[0,1],'r--',) #RED verticle dotted Line of line (^). Biggest step is understanding  
            plt.text(mintau-1,-0.1,r"mintau= "+str(round(mintau,2)),color = 'r',fontsize=12)
        
        
            plt.legend(loc='upper left')
            plt.text(3,0.19,r"$C_{f2}- \alpha C_{f1} = \frac{\alpha - 1}{e^{- \tau}-1}$",bbox=dict(facecolor='white', alpha=0.7))
            plt.text(3,0.12,r"assume $\tau_1 = \tau_2$")
            plt.text(3,0.05,r"$\alpha =$"+str(alpha) + alpha_eqn ,bbox=dict(facecolor='white', alpha=0.7))
    
    #plotting mintau on the first plot.
    
    ##########################################################################################################################################################
        if graph_num == 3 or graph_num == 0:
            Cf1=0.1
            
            plt.subplot(1, 3, 3)  
            plt.title(r'$C_{f2}/C_{f1}$ vs. $\tau$')
                  
            plt.xlabel(r"$\tau$")
            plt.ylabel(r"$C_{f2}/C_{f1}$")
            for lala in range(1,11):
                
                Cf2Cf1=((alpha-1.)/(np.exp(ttau)-1.)+alpha*Cf1)/Cf1  # Cf2/Cf1
                
                
                
                
                plt.plot(tau,Cf2Cf1,label=r'$Cf_1$='+str(Cf1),color=[0.05,0.6,lala/10.])
        
                plt.ylim(0.75 * Spectral_value_Cf2Cf1,1.2 * Spectral_value_Cf2Cf1)
                plt.xlim(-0.95,5)
                
                
                Cf1=round(Cf1+0.1,1)
                
            plt.legend(loc='lower left')
            plt.text(3, 0.84 * Spectral_value_Cf2Cf1,r"$C_{f2}- \alpha C_{f1} = \frac{\alpha - 1}{e^{- \tau}-1}$",bbox=dict(facecolor = 'white', alpha=0.7))
            plt.text(3,0.80 * Spectral_value_Cf2Cf1,r"assume $\tau_1 = \tau_2$")
            plt.text(3,0.76 * Spectral_value_Cf2Cf1,r"$\alpha =$"+str(alpha)+r"; $I_{o1} = \alpha I_{o2}$",bbox=dict(facecolor='white', alpha=0.7))
        
        
            plt.plot([0,6],[Spectral_value_Cf2Cf1,Spectral_value_Cf2Cf1],color='r') #The horizontal red line showing the ratio between spectral_v_2 to spec_v_1
            rr=round(Sv2/Sv1,2)
            plt.text(3.2,Spectral_value_Cf2Cf1 - 0.05,r"$C_{f2}/C_{f1}=$"+str(rr),color='r',fontsize=12)
        
            #plt.plot([mintau,mintau],[0,1.2 * Spectral_value_Cf2Cf1],'r--',)
            #plt.text(mintau-1,0.7 * Spectral_value_Cf2Cf1,r"mintau= "+str(round(mintau,2)),color = 'r',fontsize=12)
            
            
            plt.savefig(os.getcwd() + '/CovOpt_Plots/' + row['objectName'] + ' (' + str(prev_num) + ').png',dpi= 500)
            plt.show()
            plt.clf()