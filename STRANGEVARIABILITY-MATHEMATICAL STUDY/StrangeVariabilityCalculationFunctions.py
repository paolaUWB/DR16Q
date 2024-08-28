# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 14:48:40 2024

@author: Taylor Gibbons
"""

import numpy as np
        
def get_spectral_values(I01,I02,I12):
    Spectral_value_Cf2Cf1=round((1.-(I12/I02))/(1.-(I12/I01)),2) #Cf2/Cf1 from writing eqn for each intensity and dividing
    Spectral_value_1=round(1.-(I12/I01),4) #  = Cf1(1-exp(-tau)) 
    Spectral_value_2=round(1.-(I12/I02),4)  #=Cf2(1-exp(-tau))
    
    if Spectral_value_1 < 0:
            print('Negative SV1 Found: ' + str(Spectral_value_1))
            Spectral_value_1 = 0
            
    if Spectral_value_2 < 0:
            print('Negative SV2 Found: ' + str(Spectral_value_2))
            Spectral_value_2 = 0
    
    return Spectral_value_1,Spectral_value_2,Spectral_value_Cf2Cf1


def min_value_finder(I01,I02,I1,I2): #maybe should just be passed alpha instead of both I01 and I02???
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
    prev_num: int
        This variable is the number of times that the object has been seen. This is in
        attempts to prevent the overwriting of plots that have the same object but have
        multiple possibilities of strange variability. If the object has not been 
        seen before, the prev_num = 0. If it has been seen once before, the value
        will be 1, and twice will be 2 and so on. This will then get printed into
        the saved plot and prevent the overwriting. 
    
    Returns
    -------
    A list of the coverage fraction values for both observations.

    """
    
    I12 = round((round(I1,2) + round(I2,2))/2,2)
    I1,I2 = I12, I12
    
    Sv1, Sv2, SvCf2Cf1 = get_spectral_values(I01,I02,I12)
    
    
##############################################################################
###Coverage fraction calculations and plots    

    inv_alpha=round(I01/I02,1) # Io1=alpha * Io2 (Works best when alpha > 1???????)
    
    if I01/I02 < 1:
        alpha = round(I01/I02,2)
    else:
        alpha = round(I02/I01,2)
        
#####################################################################################################################################################
    
    tau= np.arange(0.1,5,0.1)
    
    Spectral_Cf1 = Sv1/(1-np.exp(-tau)) #this will go specifically into the graph tau vs. Cf1.
    Spectral_Cf2 = Sv2/(1.-np.exp(-tau)) 
    
    mintau=(-1)*np.log(1.-(Sv2))  # From Cf1(1-exp(-tau))=1-I1/I01; since that is the value of tau for Cf=1 and the observational data
    
    
    if mintau < 0:
        print('Negative minTau Found: ' + str(mintau))
        mintau = 0

    
            
    return alpha,inv_alpha, round(mintau,2), Sv1, Sv2, Spectral_Cf1, Spectral_Cf2, SvCf2Cf1


def tau_finder(name,I01, I02, I12):#IN PROGRESS
    I1 = I12
    I2 = I12
    
    Spectral_Cf1, Spectral_Cf2, Spectral_Cf2Cf1 = get_spectral_values(I01, I02, I12)
    maxtau1 = (-1)*np.log(1.-Spectral_value_Cf2Cf1)
    mintau1 = (-1)*np.log(1.-(Spectral_value_1/1.))
    Spectral_tau1 = (-1)*np.log(1.-(Spectral_value_1/Cf))
    ###########################################################################

    ## First set of Tau1 vs Cf
#     plt.subplot(1,3,1)
#     Cf = np.arange(0.1,1.1,0.02)
#     alpha = round(Io1/Io2,1)
#     alpha2 = round(Io2/Io1,1)
#     yyy = ((1-alpha)*(1-Cf))/(Cf) 
#     xxx = ((1-alpha2)*(1-Cf))/Cf
#     tau2 = 0.1

#     
#     if plot == 'Yes':
#         plt.title(name)
        
#         for i in range(1,52):
#             tau1 = (-1) * np.log((yyy/alpha) + np.exp(-1*tau2)/alpha)
            
#             if round(tau2-0.1,1) & 1 == 0:
#                 plt.plot(Cf, tau1, lebel = r'$\tau_2$=' + str(tau2),color = [0.05,0.6,i/51])
                
#             else:
#                 plt.plot(Cf,tau1, color = [0.05,0.6,i/51])
                
#             plt.ylabel(r'$\tau_1$')
#             plt.xlabel(r'$C_{f}$')
#             plt.xlim(0,1)
#             plt.ylim(0,5)
            
#             tau2 = round(tau2 + 0.1, 2)
    
#         tau2 = round(tau2-0.1,2)
#         plt.legend(loc = 'upper left')
#         plt.text(0.6,4.6,r"$\frac{1-C_{f}}{C_{f}} = \frac{\alpha e^{- \tau_1} - e^{- \tau_2}}{1 - \alpha}$",fontsize=10,bbox=dict(facecolor='white', alpha=0.7))
#         plt.text(0.6,4.15,r"assume $C_{f1} = C_{f2}$")
#         plt.text(0.6,3.8,r"$\alpha =$"+str(alpha)+r"; $I_{o1} = \alpha I_{o2}$",bbox=dict(facecolor='white', alpha=0.7))
    
        
#         plt.plot(Cf,Spectral_tau1,'r:')
#         plt.plot(Cf[10:],Spectral_tau1[10:],'ro')
#         plt.text(0.23,2.5,r"$C_{f} (1-e^{- \tau_1})=$"+str(SV1),color='r',fontsize=12)

    
#         plt.plot([0,1],[mintau1,mintau1],'r--')
#         plt.plot([0,1],[maxtau1,maxtau1],'r--')
    
#         plt.text(0.05,mintau1-0.25,r'min $\tau_1$ = '+str(round(mintau1,3)),color = 'r')
#         plt.plot([Spectral_value_2,Spectral_value_2],[0,5],'r--') #Verticle Line
    

    

    
    
#     ###########################################################################
#     #Second plot of Cf vs tau2
    
    
#     if plot == 'Yes'.any():
#         plt.subplot(1, 3, 2)
#         Cf=np.arange(0.1,1.1,0.02)
    
#         tau1=0.1
#         plt.title(name)
#         for lala in range(1,52): #The number of lines produced on the first graph
            
#             tau2=(-1)*(np.log(alpha*np.exp(-tau1)-yyy))
    
    
    
#         if round((tau1-0.1),1) % 1 == 0:
#             plt.plot(Cf,tau2,label=r'$\tau_1$='+str(tau1),color=[0.05,0.6,lala/51.])
#         else: 
#             plt.plot(Cf,tau2,color=[0.05,0.6,lala/51.])
#         plt.ylabel(r"$\tau_2$")
#         plt.xlabel(r"$C_{f}$")
#         plt.xlim(0,1)
#         plt.ylim(0,5)
        
#         tau1=round(tau1+0.1,2)
        
        
#         tau1=round(tau1-0.1,2)
    
    
#         plt.legend(loc='upper left')
#         plt.text(0.6,4.6,r"$\frac{1-C_{f}}{C_{f}} = \frac{\alpha e^{- \tau_1} - e^{- \tau_2}}{1 - \alpha}$",fontsize=10,bbox=dict(facecolor='white', alpha=0.7))
#         plt.text(0.6,4.15,r"assume $C_{f1} = C_{f2}$")
#         plt.text(0.6,3.8,r"$\alpha =$"+str(alpha)+r"; $I_{o1} = \alpha I_{o2}$",bbox=dict(facecolor='white', alpha=0.7))
        
#         Spectral_tau2=(-1)*np.log(1.-(Spectral_value_2/Cf))  # From Spectral_value/Cf = 1-exp(-tau)
#         plt.plot(Cf,Spectral_tau2,'r:')
#         plt.plot(Cf[10:],Spectral_tau2[10:],'ro')
#         plt.text(0.23,2.5,r"$C_{f} (1-e^{- \tau_2})=$"+str(Spectral_value_2),color='r',fontsize=12) #RED LINE
        
#         maxtau2 = (-1)*np.log(1.-Spectral_value_Cf2Cf1)
        
#         mintau2=(-1)*np.log(1.-(Spectral_value_2/1.))  #Point where Cf=1 in Cf*(1-exp(-tau2))=Sp2
#         plt.plot([0,1],[mintau2,mintau2],'r--',)
#         plt.plot([0,1],[maxtau2,maxtau2],'r--',)
#         plt.plot([Spectral_value_1,Spectral_value_1],[0,5],'r--') #Verticle Line
    
#         plt.text(0.05,maxtau2-0.25,r'max $\tau_2$ =' + str(round(maxtau2,2)), color = 'r')
#         plt.text(0.05,mintau2-0.25,r'min $\tau_2$ = '+ str(round(mintau2,2)), color = 'r')