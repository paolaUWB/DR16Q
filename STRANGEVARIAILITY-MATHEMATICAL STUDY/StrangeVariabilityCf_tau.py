# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 22:32:42 2024

@author: tayta
"""


#imports
import numpy as np

import pandas as pd

import matplotlib.pyplot as plt 

import os




#functions

def epsilon_finder(I01, I02, mintau, cf1):
    '''

    Parameters
    ----------
    I01 : float, int
        This is the value for the average flux around the trough on the first observation
    I02 : float, int
        This is the value for the average flux around the trough on the second observation
    mintau : float
        This is the minimum value of the optical depth based on the initial fluxes.

    Returns
    -------
    Epsilon: float
        A value that we can consider to be the most pessimistic value of error in order to determine if the 
        observations are considered useful data. 
    '''
    
    alpha = I01/I02
    epsilon = (alpha*cf1 - (alpha - 1))*(I02(np.exp(-mintau)-1))
    
    return epsilon


def Cf_finder(name,I01,I02,I12): #maybe should just be passed alpha instead of both I01 and I02???
    """
    Cf_finder will calculate the coverage fraction values for both observations.
    Essentially the case in which the optical depths are the same through all observations
    

    Parameters
    ----------
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
    A list of the coverage fraction values for both observations.

    """
    
    plt.figure(figsize= (18,5))
    Spectral_Cf1=round(1.-(I12/I01),4)/(1-np.exp(-tau))
    Spectral_value_2 = round(1-(I12/I02))#  = Cf1(1-exp(-tau)) 
    Spectral_Cf2=(Spectral_value_2)/(1-np.exp(-tau))  #=Cf2(1-exp(-tau))
    alpha = round(I01/I02,1)
    
##############################################################################
###Coverage fraction 1 calculations and plots    

    Cf2 = 0.1
    

    plt.suptitle(name)
    plt.subplot(1,3,1)
    plt.title('Finding Coverage Fraction 1')
    for i in range(1,11):
        #print(index)
        Cf1 = (Cf2 - ((alpha-1)/(np.exp(-tau)-1)))
                
        plt.plot(tau, Cf1 ,label = r'$Cf_2$='+str(Cf2),color=[0.05,0.6,i/10.])
        
        plt.ylabel(r"$C_{f1}$")
        plt.xlabel(r"$\tau$")
        
        plt.ylim(0,1)
        plt.xlim(-0.95,5)

        Cf2 = round(Cf2 + 0.1, 1)
    
    plt.plot(tau,Spectral_Cf1,color= 'r')
    '''
    plt.text(2.8,Spectral_Cf1,r'$C_{f1}(1-e^{- \tau})=$' + str(Spectral_Cf1),color='r',fontsize=10)
    '''
    plt.legend(loc='lower left')
    plt.text(3,0.19,r"$C_{f2}- \alpha C_{f1} = \frac{\alpha - 1}{e^{- \tau}-1}$",bbox=dict(facecolor='white', alpha=0.7))
    plt.text(3,0.12,r"assume $\tau_1 = \tau_2$")
    plt.text(3,0.05,r"$\alpha =$"+str(alpha)+r"; $I_{o1} = \alpha I_{o2}$",bbox=dict(facecolor='white', alpha=0.7))
    
    mintau=(-1)*np.log(1.-(Spectral_value_2/1.))  # From Cf1(1-exp(-tau))=1-I1/I01; since that is the value of tau for Cf=1 and the observational data
    plt.plot([mintau,mintau],[0,1],'r--',) #RED verticle dotted Line of line (^). Biggest step is understanding  
    plt.text(mintau-1,-0.1,r"mintau= "+str(round(mintau,2)),color = 'r',fontsize=12)
    ###########################################################################
    #Coverage Fraction 2 
    
    plt.subplot(1,3,2)
    plt.title('Finding Coverage Fraction 2')
    Cf1 = 0.1
    
    for i in range(1,11):
        
        Cf2 = Cf1 * alpha + ((alpha-1)/(np.exp(-tau)-1))
        
        plt.plot(tau,Cf1, label = r'$Cf_1$=' + str(Cf1), color = [0.05,0.6,i/10.])
        
        plt.xlabel(r'$\tau$')
        plt.ylabel(r'$C_{f2}$')
        plt.ylim(0,1)
        plt.xlim(-0.95,5)
        
        Cf1 = round(Cf1 + 0.1,1)
    
    plt.plot(tau, Spectral_Cf2, color = 'r')
    plt.text(2.8,Spectral_value_2-0.05,r"$C_{f2}(1-e^{- \tau})=$"+str(Spectral_value_2),color='r',fontsize=10)
      
    mintau=(-1)*np.log(1.-(Spectral_value_2/1.))
    
    plt.plot([mintau,mintau],[0,1],'r--',) #RED verticle dotted Line of line (^). Biggest step is understanding  
    plt.text(mintau-1,-0.1,r"mintau= "+str(round(mintau,2)),color = 'r',fontsize=12)
   
    #plt.plot(tau,Cf1)
    
    plt.subplot(1,3,3)
    
    plt.savefig(os.getcwd() + '/FigureSavingTesting/' + name + '.png',dpi= 500)

def tau_finder(name,Io1, Io2, I12):
    plt.subplot(1,3,1)
    Cf = np.arange(0.1,1.1,0.02)
    alpha = round(Io1/Io2,1)
    alpha2 = round(Io2/Io1,1)
    yyy = ((1-alpha)*(1-Cf))/(Cf) 
    xxx = ((1-alpha2)*(1-Cf))/Cf
    tau2 = 0.1
    plt.title(name)
    
    for i in range(1,52):
        tau1 = (-1) * np.log((yyy/alpha) + np.exp(-1*tau2)/alpha)
        
        if round(tau2-0.1,1) & 1 == 0:
            plt.plot(Cf, tau1, lebel = r'$\tau_2$=' + str(tau2),color = [0.05,0.6,i/51])
            
        else:
            plt.plot(Cf,tau1, color = [0.05,0.6,i/51])
            
        plt.ylabel(r'$\tau_1$')
        plt.xlabel(r'$C_{f}$')
        plt.xlim(0,1)
        plt.ylim(0,5)
        
        tau2 = round(tau2 + 0.1, 2)

    tau2 = round(tau2-0.1,2)




    

#Main Programming
'''
If there is going to be anything added, make sure that it is being added here
(unless it is a new function then whatever works bruh)
'''
    

filename = 'StrangeVariability.csv'



#Makes a dataframe of all the data from csv using PANDAS
data = pd.read_csv(filename)
tau = np.arange(0.1,5.1,0.1)
for index,row in data.iterrows():
    print(row['Spec1'])
    
    Cf_finder(row['Spec1'], row['I01'], row['I02'], row['I1'],)





plt.show()

"""
    SV1, SV2 = spectral_finder(tau, Io1, Io2, I1)
plt.figure(1)
plot(tau,SV1,'Spectral Value 1','$\tau$','Spectral Value 1')
plot(tau,SV2,'Spectral Value 2','tau','Spectral Value 2')
"""
    
