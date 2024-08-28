#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 09:39:28 2018

@author: Paola Rodriguez Hidalgo 
with additions by Taylor Gibbons
"""

import numpy as np

import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages



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
    
    alpha = Io1/Io2
    epsilon = (alpha*cf1 - (alpha - 1))*(I02(np.exp(-mintau)-1))
    
    return epsilon


def Cf_finder(tau,I01,I02,I12): #maybe should just be passed alpha instead of both I01 and I02???
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
    alpha = I01/I02
    

def spectral_finder(tau, I01,I02,I1): 
    """
    Calculates the spectral values of each observation. 

    Parameters
    ----------
    tau : Array; dtype= float
        Array of values ranging from 0 to a number sufficiently high to consider
        an outflow saturated. This represents the optical depth of an outflow
    I01 : TYPE
        .
    I02 : TYPE
        DESCRIPTION.
    I1 : TYPE
        DESCRIPTION.

    Returns
    -------
    SV1 : Array;dtype= float
        An array of the spectral values for the first observation.
    SV2 : Float
        An array of the spectral values for the second observation.

    """
    
    SV1 = np.zeros(len(tau), dtype= float)
    SV1 = round(1-I1/I01,4)
    SV1 = SV1/(1 - np.exp(-tau))
    
    SV2 = np.zeros(len(tau), dtype = float)
    SV2 = round(1-I1/I02,4)
    SV2 = SV2/(1 - np.exp(-tau))
    
    return SV1, SV2
    
    
#pp2= PdfPages('Cf_tau_sametau.pdf')
f = plt.figure(figsize=(18,5))
ax = f.add_subplot(121)
ax2 = f.add_subplot(122)

save_format = '.pdf' #Does not work, I think because the code originally uses PDFpages and not just the plot.save stuff
#==========================================
# Input parameters: values from the spectra



case_num = 7

case_name = ''


#I1=I2 is the value at the bottom of the trough, they will be equal
#Io1 is the avged value of the first observation above the trough (not sure if the whole spectra before/after or only a little bit)
#Io2 is the avged value of the second observation above the trough ^
if case_num == 1  :  
    case_name = 'Case spSpec-51959-051-305'
    I1=3
    I2=I1
    Io1=10
    Io2=7

elif case_num == 2:
        case_name = 'Case LBQS'
        I1=8
        I2=I1
        Io1=23
        Io2=18

elif case_num == 3:
        case_name = 'Case 304'
        I1=1.5 
        I2=I1
        Io1=10
        Io2=7.5

elif case_num ==4: 
        case_name = 'Case 901'
        I1=10.5
        I2=I1
        Io1=22
        Io2=16.5
        
elif case_num == 5:
        case_name = 'Case spSpec-51959-051-305'
        I1=3
        I2=I1
        Io1=10
        Io2=7

elif case_num == 6:
        case_name = 'Case LBQS, Absorption complex D'
        I1=8
        I2=I1
        Io1=23
        Io2=18

elif case_num == 7:
        case_name = 'J123830.71+480003.4' #Name of object changes based on.... object
        I1 = 12.7  #average minimum value
        I2 = I1
        Io1 = 20.63 #average values from before and after trough of first observation
        Io2 = 23.66 #average values from before and after trough of second observation
else :
        case_name = 'No inputs'
        I1 = 0
        I2 = I1
        Io1 = 0
        Io2 = 0
        

    

tau = np.arange(0.1,5.1,0.1)
SV1,SV2 = spectral_finder(tau, Io1, Io2, I1)
#Cf1,Cf2 = Cf_finder(Io1,Io2,I2)

plt.plot(tau, SV1, 'r',tau,SV2,'b')
        
        
# ------------------------------------------
#-- First set of plots. 

alpha=round(Io1/Io2,1) # Io1=alpha * Io2 (Works best when alpha > 1???????) 

pp2= PdfPages(case_name + '-tau' + save_format)

Spectral_value_Cf2Cf1=round((1.-(I2/Io2))/(1.-(I1/Io1)),2) #Cf2/Cf1 from writing eqn for each intensity and dividing
Spectral_value_1=round(1.-(I1/Io1),4) #  = Cf1(1-exp(-tau)) 
Spectral_value_2=round(1.-(I2/Io2),4)  #=Cf2(1-exp(-tau))


#-------
#These lines are attempts for Taylor to see if the Spec.Values were wrongly written, however ^ equations are
#Correct (as far as I know) because the lines below become negative, which doesn't work inside of logs.

#Spectral_value_1 = round((I1/Io1)-1.,2) #=Cf1(exp(-tau)-1) Taylors Equation > I1 = (1-Cf1)Io1+Cf1*Io1*exp(-tau)
#Spectral_value_2 = round((I2/Io2)-1.,2) #=Cf2(exp(-tau)-1) ^


##########################################################################################################################################################
#The graph for finding tau vs. Cf1. 

minbufferCf1 = 0 #this will change the location of the mintau by +- how every much. left needs (-), right is +

tau= np.arange(0.1,5,0.1)

ttau=tau*(-1)  #= -tau

Cf2=0.1
#yy= 1-(Cf2*np.exp(ttau))/(alpha+np.exp(ttau)) # right side of eqn
#xx = ((alpha-1)/(np.exp(-tau)-1))
plt.subplot(1, 3, 1)
plt.suptitle(case_name)
for lala in range(1,11):

   # Cf1 = (Cf2/alpha)-(1-(1/alpha))/(np.exp(-tau)-1)
    aCf1 = (Cf2 - ((alpha-1)/(np.exp(ttau)-1)))
    Cf1 = aCf1 / alpha
    plt.plot(tau, Cf1 ,label = r'$Cf_2$='+str(Cf2),color=[0.05,0.6,lala/10.])

    plt.xlabel(r"$\tau$")
    #plt.ylabel(r"$C_{f2}- \alpha C_{f1}$")
    plt.ylabel(r"$C_{f1}$")
    plt.ylim(0,1)
    plt.xlim(-0.95,5)
          
    Cf2=round(Cf2+0.1,1)
    
    
Spectral_Cf1 = Spectral_value_1/(1-np.exp(-tau)) #this will go specifically into the graph tau vs. Cf1.
plt.plot(tau,Spectral_Cf1,color = 'r')
plt.text(2.8,Spectral_value_1-0.05,r'$C_{f1}(1-e^{- \tau})=$' + str(Spectral_value_1),color='r',fontsize=10)

plt.legend(loc='lower left')
plt.text(3,0.19,r"$C_{f2}- \alpha C_{f1} = \frac{\alpha - 1}{e^{- \tau}-1}$",bbox=dict(facecolor='white', alpha=0.7))
plt.text(3,0.12,r"assume $\tau_1 = \tau_2$")
plt.text(3,0.05,r"$\alpha =$"+str(alpha)+r"; $I_{o1} = \alpha I_{o2}$",bbox=dict(facecolor='white', alpha=0.7))

mintau=(-1)*np.log(1.-(Spectral_value_2))  # From Cf1(1-exp(-tau))=1-I1/I01; since that is the value of tau for Cf=1 and the observational data
plt.plot([mintau+minbufferCf1,mintau+minbufferCf1],[0,1],'r--',) #RED verticle dotted Line of line (^). Biggest step is understanding  
plt.text(mintau-1,-0.1,r"mintau= "+str(round(mintau+minbufferCf1,2)),color = 'r',fontsize=12)

#Trying to find the maximum value of tau where there can be values found. What I think could work is by putting the
#spectral value into the eqn (1) from paper. not sure how though. 
##########################################################################################################################################################

minbufferCf2 = 0 #this will change the location of the mintau by +- how every much. left needs (-), right is +

tau= np.arange(0.1,5,0.1)
ttau=tau*(-1)  #= -tau
yy=(alpha-1.)/(np.exp(ttau)-1.)  # right side of eqn

Cf1=0.1
plt.subplot(1, 3, 2)
plt.title(case_name)
for lala in range(1,11):

    Cf2 = yy+alpha*Cf1 #Correct -Taylor
    
    
    plt.plot(tau, Cf2 , label = r'$Cf_1$='+str(Cf1),color=[0.05,0.6,lala/10.])

    plt.xlabel(r"$\tau$")
    #plt.ylabel(r"$C_{f2}- \alpha C_{f1}$")
    plt.ylabel(r"$C_{f2}$")
    plt.ylim(0,1)
    plt.xlim(-0.95,5)

    
    Cf1=round(Cf1+0.1,1)

Spectral_Cf2=Spectral_value_2/(1.-np.exp(-tau)) 
#BEST IDEA YET:
#Spectral_cf2 is the possible values of Cf that can come from one of the observations, thus is why it is the
#line that can show the possibility of the values.

Spectral_Cf1 = Spectral_value_1/(1-np.exp(-tau)) #this will go specifically into the graph tau vs. Cf1.

plt.plot(tau,Spectral_Cf2,color='r')


plt.text(2.8,Spectral_value_2-0.05,r"$C_{f2}(1-e^{- \tau})=$"+str(Spectral_value_2),color='r',fontsize=10) 


mintau=(-1)*np.log(1.-(Spectral_value_2/1.))  # From Cf1(1-exp(-tau))=1-I1/I01; since that is the value of tau for Cf=1 and the observational data

#print('mintau = '+ str(mintau))
plt.plot([mintau+minbufferCf2,mintau+minbufferCf2],[0,1],'r--',) #RED verticle dotted Line of line (^). Biggest step is understanding  
plt.text(mintau-1,-0.1,r"mintau= "+str(round(mintau+minbufferCf2,2)),color = 'r',fontsize=12)

#ee=np.exp(-tau)*Cf2
#sv=min(abs(ee-Spectral_value_2))
#sv_index=np.where(ee == sv+Spectral_value_2)
#sv_index=np.where(ee == -(sv+Spectral_value_2))


plt.legend(loc='upper left')
plt.text(3,0.19,r"$C_{f2}- \alpha C_{f1} = \frac{\alpha - 1}{e^{- \tau}-1}$",bbox=dict(facecolor='white', alpha=0.7))
plt.text(3,0.12,r"assume $\tau_1 = \tau_2$")
plt.text(3,0.05,r"$\alpha =$"+str(alpha)+r"; $I_{o1} = \alpha I_{o2}$",bbox=dict(facecolor='white', alpha=0.7))

#plotting mintau on the first plot.

##########################################################################################################################################################

minbuffer = -0.23 #this will change the location of the mintau by +- how every much. left needs (-), right is +

Cf1=0.1

for lala in range(1,11):
    
    Cf2Cf1=(yy+alpha*Cf1)/Cf1  # Cf2/Cf1
    
    plt.title(case_name)
    
    plt.subplot(1, 3, 3)
    plt.plot(tau,Cf2Cf1,label=r'$Cf_1$='+str(Cf1),color=[0.05,0.6,lala/10.])
    plt.xlabel(r"$\tau$")
    plt.ylabel(r"$C_{f2}/C_{f1}$")
    plt.ylim(0.75 * Spectral_value_Cf2Cf1,1.2 * Spectral_value_Cf2Cf1)
    plt.xlim(-0.95,5)
    
    
    Cf1=round(Cf1+0.1,1)
    
plt.legend(loc='upper right')
plt.text(3, 0.84 * Spectral_value_Cf2Cf1,r"$C_{f2}- \alpha C_{f1} = \frac{\alpha - 1}{e^{- \tau}-1}$",bbox=dict(facecolor='white', alpha=0.7))
plt.text(3,0.80 * Spectral_value_Cf2Cf1,r"assume $\tau_1 = \tau_2$")
plt.text(3,0.76 * Spectral_value_Cf2Cf1,r"$\alpha =$"+str(alpha)+r"; $I_{o1} = \alpha I_{o2}$",bbox=dict(facecolor='white', alpha=0.7))


plt.plot([0,6],[Spectral_value_Cf2Cf1,Spectral_value_Cf2Cf1],color='r') #The horizontal red line showing the ratio between spectral_v_2 to spec_v_1
rr=round(Spectral_value_2/Spectral_value_1,2)
plt.text(3.2,Spectral_value_Cf2Cf1 - 0.05,r"$C_{f2}/C_{f1}=$"+str(rr),color='r',fontsize=12)

plt.plot([mintau+minbuffer,mintau+minbuffer],[0,1.2 * Spectral_value_Cf2Cf1],'r--',)
plt.text(mintau-1,0.7 * Spectral_value_Cf2Cf1,r"mintau= "+str(round(mintau+minbuffer,2)),color = 'r',fontsize=12)

#pp1.savefig()
#pp1.close()
pp2.savefig()
plt.show()
pp2.close()
plt.clf()
##########################################################################################################################################################

# -------------------------- Second set of plots 

# This case seems unreasonable to me:
# If I don't have absorption, e(-tau)=1 --> I1=I01; I2=Io2; I1=I2, but Io1=alpha*Io2
# but if it skips the tau =o region, it mihgt be feasable? (and the red part avoids it...)

pp1= PdfPages(case_name + '-Cf' + save_format)

ff = plt.figure(figsize=(18,5))
ax = f.add_subplot(121)
ax2 = f.add_subplot(122)

##########################################################################################################################################################
#First plot of Cf vs tau1

plt.subplot(1,3,1)
Cf = np.arange(0.1,1.1,0.02)
alpha2 = Io2/Io1
yyy = ((1-alpha)*(1-Cf))/(Cf) 
xxx = ((1-alpha2)*(1-Cf))/Cf
tau2 = 0.1
plt.title(case_name)
for lala in range(1,52):
    
    #tau1 = (-1)*np.log((yyy+np.exp((-1)*tau2))/alpha)
    tau1 = (-1)*np.log((yyy/alpha) + np.exp(-1*tau2)/alpha)
    #taux = (-1)*np.log(alpha2*np.exp(-tau2)-(xxx))
    if round(tau2-0.1,1) % 1 == 0:
        plt.plot(Cf, tau1, label = r'$\tau_2$=' + str(tau2), color = [0.05,0.6,lala/51.])
        #plt.plot(Cf, taux, label = r'$\tau_x$=' + str(tau2), color = 'm')
    else:
        plt.plot(Cf,tau1, color=[0.05,0.6,lala/51.])
        #plt.plot(Cf,taux, color = 'm')
    plt.ylabel(r'$\tau_1$')
    plt.xlabel(r'$C_{f}$')
    plt.xlim(0,1) # 0 < Cf <1
    plt.ylim(0,5)
    
    tau2 = round(tau2+0.1,2)

tau2 = round(tau2-0.1,2)

plt.legend(loc = 'upper left')
plt.text(0.6,4.6,r"$\frac{1-C_{f}}{C_{f}} = \frac{\alpha e^{- \tau_1} - e^{- \tau_2}}{1 - \alpha}$",fontsize=10,bbox=dict(facecolor='white', alpha=0.7))
plt.text(0.6,4.15,r"assume $C_{f1} = C_{f2}$")
plt.text(0.6,3.8,r"$\alpha =$"+str(alpha)+r"; $I_{o1} = \alpha I_{o2}$",bbox=dict(facecolor='white', alpha=0.7))

Spectral_tau1 = (-1)*np.log(1.-(Spectral_value_1/Cf))
plt.plot(Cf,Spectral_tau1,'r:')
plt.plot(Cf[10:],Spectral_tau1[10:],'ro')
plt.text(0.23,2.5,r"$C_{f} (1-e^{- \tau_1})=$"+str(Spectral_value_1),color='r',fontsize=12)
maxtau1 = (-1)*np.log(1.-Spectral_value_Cf2Cf1)
mintau1 = (-1)*np.log(1.-(Spectral_value_1/1.))

plt.plot([0,1],[mintau1,mintau1],'r--')
plt.plot([0,1],[maxtau1,maxtau1],'r--')
#plt.text(0.05,maxtau1-0.25,r'max $\tau_1$ = '+ str(round(maxtau1,2)), color = 'r')
plt.text(0.05,mintau1-0.25,r'min $\tau_1$ = '+str(round(mintau1,3)),color = 'r')
plt.plot([Spectral_value_2,Spectral_value_2],[0,5],'r--') #Verticle Line

##########################################################################################################################################################
#Second plot of Cf vs tau2
plt.subplot(1, 3, 2)
Cf=np.arange(0.1,1.1,0.02)

tau1=0.1

for lala in range(1,52): #The number of lines produced on the first graph
    plt.title(case_name)
    tau2=(-1)*(np.log(alpha*np.exp(-tau1)-yyy)) #Equation is correct, but coming up with an error when ran
    #print('lala = ' + str(lala))
    #print ('tau2 = '+ str(tau2))

    if round((tau1-0.1),1) % 1 == 0:
        plt.plot(Cf,tau2,label=r'$\tau_1$='+str(tau1),color=[0.05,0.6,lala/51.])
    else: 
        plt.plot(Cf,tau2,color=[0.05,0.6,lala/51.])
    plt.ylabel(r"$\tau_2$")
    plt.xlabel(r"$C_{f}$")
    plt.xlim(0,1)
    plt.ylim(0,5)
    
    tau1=round(tau1+0.1,2)
    
    
tau1=round(tau1-0.1,2)
    
plt.legend(loc='upper left')
plt.text(0.6,4.6,r"$\frac{1-C_{f}}{C_{f}} = \frac{\alpha e^{- \tau_1} - e^{- \tau_2}}{1 - \alpha}$",fontsize=10,bbox=dict(facecolor='white', alpha=0.7))
plt.text(0.6,4.15,r"assume $C_{f1} = C_{f2}$")
plt.text(0.6,3.8,r"$\alpha =$"+str(alpha)+r"; $I_{o1} = \alpha I_{o2}$",bbox=dict(facecolor='white', alpha=0.7))


#plt.plot([mintau,mintau],[0,1],'r--',)


Spectral_tau2=(-1)*np.log(1.-(Spectral_value_2/Cf))  # From Spectral_value/Cf = 1-exp(-tau)
plt.plot(Cf,Spectral_tau2,'r:')
plt.plot(Cf[10:],Spectral_tau2[10:],'ro')
plt.text(0.23,2.5,r"$C_{f} (1-e^{- \tau_2})=$"+str(Spectral_value_2),color='r',fontsize=12) #RED LINE

#maxtau2=-(1)*np.log(1.-((Spectral_value_2/Spectral_value_1)*(1.-np.exp(-tau1))))  # Bring Cf=Sp1/(1-exp(-tau1)) into Cf*(1-exp(-tau2))=Sp2
maxtau2=(-1)*np.log(1.-((Spectral_value_2/Spectral_value_1)))  # Bring Cf=Sp1/(1-exp(-tau1)) into Cf*(1-exp(-tau2))=Sp2
#equal to
maxtau2 = (-1)*np.log(1.-Spectral_value_Cf2Cf1)

mintau2=(-1)*np.log(1.-(Spectral_value_2/1.))  #Point where Cf=1 in Cf*(1-exp(-tau2))=Sp2
plt.plot([0,1],[mintau2,mintau2],'r--',)
plt.plot([0,1],[maxtau2,maxtau2],'r--',)
plt.plot([Spectral_value_1,Spectral_value_1],[0,5],'r--') #Verticle Line

plt.text(0.05,maxtau2-0.25,r'max $\tau_2$ =' + str(round(maxtau2,2)), color = 'r')
plt.text(0.05,mintau2-0.25,r'min $\tau_2$ = '+ str(round(mintau2,2)), color = 'r')

##########################################################################################################################################################

tau1=0.1

plt.subplot(1, 3, 3)
for lala in range(1,52):
    
    plt.title(case_name)
    
    tau2=(-1)*np.log(alpha*np.exp(-tau1)-yyy)
    
    tau2tau1=tau2/tau1
    
    
    
    if round((tau1-0.1),1) % 1 == 0:
        plt.plot(Cf,tau2tau1,label=r'$\tau_1$='+str(tau1 ),color=[0.05,0.6,lala/51.])
    else: 
        plt.plot(Cf,tau2tau1,color=[0.05,0.6,lala/51.])
    plt.ylabel(r"$\tau_2/\tau_1$")
    plt.xlabel(r"$C_{f}$")
    plt.xlim(0,1)
    plt.ylim(0,10)
    
    tau1=round(tau1+0.1,2)
    
plt.legend(loc='upper left')
#plt.text(0.55,0.92 * tau2tau1,r"$\frac{1-C_{f}}{C_{f}} = \frac{\alpha e^{- \tau_1} - e^{- \tau_2}}{1 - \alpha}$",fontsize=12,bbox=dict(facecolor='white', alpha=0.7))
#plt.text(0.55,0.83 * tau2tau1,r"assume $C_{f1} = C_{f2}$")
#plt.text(0.55,0.76 * tau2tau1,r"$\alpha =$"+str(alpha)+r"; $I_{o1} = \alpha I_{o2}$",bbox=dict(facecolor='white', alpha=0.7))


pp1.savefig()
plt.show()
pp1.close()