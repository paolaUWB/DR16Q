#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.io import ascii
from astropy.table import unique


# In[2]:


qso = Table.read(ddir + 'alr53/BALs_nonBALs/mastertable_AmyR_extended2_LDR1_210210.fits')
CIVm = np.zeros(len(qso['CIV_EW_eun']), 'bool')
CIVm[np.where((qso['CIV_EW_eun']>10.) & (qso['CIV_blue_eun']>-1200) & (qso['CIV_blue_eun']<7500) & 
              (qso['chisqr']<2.) & (qso['Ngood_CIV']>0))[0]] = 1
CIVm[qso['AlIII_BI'] > 0] = 0
CIVm[qso['FeLoBAL']] = 0

CIVm[np.log10(qso['CIV_EW_eun'])<-0.00023077 * qso['CIV_blue_eun'] + 1.3231] = 0

good = np.zeros(len(qso['BI_BI']), 'bool')
good[np.where((qso['chisqr'] <= 2.) & (CIVm==1))[0]] = 1
good[qso['AlIII_BI'] > 0] = 0

inds = (qso['BI_BI'] > 0) * (qso['catalogue'] == 1) * good * (qso['CIV_EW'] > 10) #* (qso['SNRpaul'] >= 10)


# In[3]:


speeds = ascii.read('DR16EHVO_DR14RankineInfo_wEHVOspeed_corrected.csv', format='commented_header')
speeds


# In[4]:


len(unique(speeds['SDSS_NAME0', 'plate_parentDR16']))


# In[5]:


speedsdr9 = ascii.read('DR9EHVO_DR14RankineInfo_wEHVOspeed_corrected.csv', format='commented_header')
speedsdr9['dr9'] = np.ones(len(speedsdr9), bool)
speedsdr9


# In[6]:


len(unique(speedsdr9['SDSS_NAME0', 'plate_parentDR16']))


# In[7]:


from astropy.table import vstack, unique


# In[8]:


speeds = vstack((speeds, speedsdr9))
speeds['dr9'][speeds['dr9'].mask] = False
speeds


# In[9]:


len(unique(speeds['SDSS_NAME0', 'plate_parentDR16']))


# In[10]:


ehvo = Table.read('EVHO_DR16.fits')
ehvo.info


# In[12]:


plt.scatter(speeds['CIV_blue0'], np.log10(speeds['CIV_EW0']), c=speeds['EHVO_vmax'], # VMAX
            zorder=100, cmap='RdYlBu_r', edgecolor='k', s=50)
xlims = plt.gca().get_xlim()
ylims = plt.gca().get_ylim()


# In[13]:


# remove crazy high EW EHVO
speeds.remove_row(np.where(np.log10(speeds['CIV_EW0']) > 2.1)[0][0])


# In[14]:


colormap='Blues'
import matplotlib.colors as colors

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))
fig, ax = plt.subplots(ncols=2)

cmap = plt.get_cmap(colormap)
new_cmap = truncate_colormap(cmap, 0.1, 1)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap)
plt.show()


# # Vmax
# 

# In[15]:


plt.hist(-speeds['EHVO_vmax']/1e3, cumulative=True)


# In[16]:


minvmax, maxvmax = np.min(-speeds['EHVO_vmax']/1e3), np.max(-speeds['EHVO_vmax']/1e3)
vmax = qso['BI_VMAX'][inds]

vmax[vmax < np.percentile(vmax, 5)] = np.percentile(vmax, 5)
vmax[vmax > np.percentile(vmax, 95)] = np.percentile(vmax, 95)


# In[17]:


min(vmax), max(-speeds['EHVO_vmax'])


# In[18]:


vmax = vmax/10**3


# # Vmin
# 

# In[19]:


plt.hist(-speeds['EHVO_vmin']/1e3, 50, cumulative=True)
plt.axvline(np.percentile(-speeds['EHVO_vmin']/1e3, 10), c='C1')


# In[20]:


np.percentile(-speeds['EHVO_vmin']/1e3, 10)


# In[21]:


vmin = qso['BI_VMIN'][inds] / 10**3
#minvmin, maxvmin = np.percentile(-speeds['EHVO_vmin']/1e3, 1), 50
minvmin, maxvmin = np.min(-speeds['EHVO_vmin']/1e3), np.max(-speeds['EHVO_vmin']/1e3)
minvmin, maxvmin


# In[22]:


vmin[vmin > np.percentile(vmin, 90)] = np.percentile(vmin, 90)


# In[23]:


# BAL or non-BAL in Rankine+20
nb = speeds['nBAL0'] == 'True'
b = speeds['BAL0'] == 'True'
# bal = np.zeros(len(ehvo), 'bool')
# bal[~ehvo['Qinds'].mask][qso['BI_BI'][ehvo['Qinds'][~ehvo['Qinds'].mask]] > 0] = 1
# bal
sum(b), sum(nb)


# In[27]:


fig = plt.figure(figsize=(12, 6))
x0 = 0.05
y0 = 0.05
xw = 0.45
yh = 0.95

ygap = 0.16

cw = 0.03
cgap = 0.02

axmin = fig.add_axes([x0, y0, xw, yh])

cax1min = fig.add_axes([x0 + xw + cgap, y0, cw, yh/2])

cax2min = fig.add_axes([x0 + xw + cgap, y0 + yh/2, cw, yh/2])

axmax = fig.add_axes([x0 + xw + cgap + cw + ygap, y0, xw, yh], sharex=axmin, sharey=axmin)

cax1max = fig.add_axes([x0 + xw + cgap + cw + ygap + xw + cgap, y0, cw, yh/2])

cax2max = fig.add_axes([x0 + xw + cgap + cw + ygap + xw + cgap, y0 + yh/2, cw, yh/2])

axmin.set_xlim(min([-1500, xlims[0]]), max([5000, xlims[1]]))
axmin.set_ylim(min([1, ylims[0]]), max([2.3, ylims[1]]))
axmin.set_xlabel(CIVx)
axmin.set_ylabel(CIVy)
axmax.set_xlabel(CIVx)
axmax.set_ylabel(CIVy)

fig.text(x0 + xw + cgap + cw + 0.06, y0 + yh/2, r'V$_{\text{min}}$ [$10^3$ km\,s$^{-1}$]', 
         rotation=270, fontsize=18, va='center', ha='left')

fig.text(x0 + xw + cgap + cw + ygap + xw + cgap+ cw + 0.06, y0 + yh/2, r'V$_{\text{max}}$ [$10^3$ km\,s$^{-1}$]', 
         rotation=270, fontsize=18, va='center', ha='left')



#### VMIN #####
new_cmap = truncate_colormap(cmap, 0.05, 0.6)
cm=axmin.hexbin(qso['CIV_blue_eun'][inds], np.log10(qso['CIV_EW_eun'][inds]),
                                         gridsize=50,
                                         C=vmin,
                                         cmap=new_cmap,
                                         reduce_C_function=np.median, mincnt=4,
                                         vmin=None, vmax=None,#vmin[i], vmax=vmax[i],
                                         edgecolor=None, zorder=99)

cbar = plt.colorbar(cm, cax=cax1min)
# print(cbar.ax.get_yticklabels())
# cbar.ax.set_yticks(cbar.ax.get_yticks()[1:-2])#(MaxNLocator(prune='both'))
cbar.set_label('Rankine+20 BALs')

cmeh = axmin.scatter(speeds['CIV_blue0'][nb], np.log10(speeds['CIV_EW0'][nb]), c=-speeds['EHVO_vmin'][nb]/1e3, 
                     zorder=100, cmap=truncate_colormap(plt.get_cmap('Purples'), 0.1, 1), 
                     edgecolor='k', s=80,
                     vmin=minvmin, vmax=maxvmin, label='EHVO only')

cmeh = axmin.scatter(speeds['CIV_blue0'][b], np.log10(speeds['CIV_EW0'][b]), c=-speeds['EHVO_vmin'][b]/1e3, 
                     zorder=100, cmap=truncate_colormap(plt.get_cmap('Purples'), 0.1, 1), 
                     edgecolor='k', s=200, marker='*',
                     vmin=minvmin, vmax=maxvmin, label='EHVO+BAL')


cbareh = plt.colorbar(cmeh, cax=cax2min)
cbareh.set_label(r'RH+20, CS+in prep EHVOs')

axmin.legend()

#### VMAX #####
#new_cmap = truncate_colormap(cmap, 0.1, 1)
cm=axmax.hexbin(qso['CIV_blue_eun'][inds], np.log10(qso['CIV_EW_eun'][inds]),
                                         gridsize=50,
                                         C=vmax,
                                         cmap=new_cmap,
                                         reduce_C_function=np.median, mincnt=4,
                                         vmin=None, vmax=None,#vmin[i], vmax=vmax[i],
                                         edgecolor=None, zorder=99)


cbar = plt.colorbar(cm, cax=cax1max)
cbar.set_label('Rankine+20 BALs')

#edgecolorn = ['w' if v/1e3 > 50 else 'k' for v in -speeds['VMAXz'][nb]]
#edgecolorb = ['w' if v/1e3 > 50 else 'k' for v in -speeds['VMAXz'][b]]

cmeh = axmax.scatter(speeds['CIV_blue0'][nb], np.log10(speeds['CIV_EW0'][nb]), c=-speeds['EHVO_vmax'][nb]/1e3, 
                     zorder=100, cmap=truncate_colormap(plt.get_cmap('Purples'), 0.1, 1), 
                     edgecolor='k', 
                     s=80,
                     vmin=minvmax, vmax=maxvmax)

cmeh = axmax.scatter(speeds['CIV_blue0'][b], np.log10(speeds['CIV_EW0'][b]), c=-speeds['EHVO_vmax'][b]/1e3, 
                     zorder=100, cmap=truncate_colormap(plt.get_cmap('Purples'), 0.1, 1), 
                     edgecolor='k', 
                     s=200, marker='*',
                     vmin=minvmax, vmax=maxvmax)

cbareh = plt.colorbar(cmeh, cax=cax2max)
cbareh.set_label(r'RH+20, CS+in prep EHVOs')

plt.savefig('CIV_vmin_vmax_Rankine_crange.pdf')


# In[ ]:





# In[ ]:




