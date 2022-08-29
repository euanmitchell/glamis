############################################
# Script to plot the summary results of
# ground-finding code as plots of beam
# sensitivity vs. bias or RMSE for different
# groundFinder.py updateMethods.
############################################

import numpy as np
import matplotlib.pyplot as plt
import argparse
import warnings

#########################################

# Plotting settings
plt.rcParams['figure.constrained_layout.use']=True
plt.rcParams['figure.figsize']=(6,3.5)
plt.rcParams['legend.fontsize']=5
plt.rcParams['legend.frameon']=False
plt.rcParams['xtick.labelsize']=5
plt.rcParams['xtick.major.size']=2
plt.rcParams['xtick.major.width']=0.4
plt.rcParams['xtick.major.pad']=2
plt.rcParams['ytick.labelsize']=5
plt.rcParams['ytick.major.size']=2
plt.rcParams['ytick.major.width']=0.4
plt.rcParams['ytick.major.pad']=2
plt.rcParams['axes.labelsize']=6
plt.rcParams['axes.linewidth']=0.5
plt.rcParams['axes.labelpad']=3
plt.rcParams['axes.titlesize']=7
plt.rcParams['lines.linewidth']=0.5
plt.rcParams['lines.markersize']=2
#plt.rcParams['savefig.format']='tif'

#########################################

# Get the text files containing the data
inputLS='../data/summary/laselvaSummary.txt'
inputRC='../data/summary/robsonSummary.txt'
inputN='../data/summary/nouraguesSummary.txt'
inputP='../data/summary/paracouSummary.txt'
inputHB='../data/summary/hubbardSummary.txt'
inputWR='../data/summary/windriverSummary.txt'
inputOR='../data/summary/oakridgeSummary.txt'

inputLS10='../data/summary/laselvaSummary10m.txt'
inputRC10='../data/summary/robsonSummary10m.txt'
inputN10='../data/summary/nouraguesSummary10m.txt'

# Text file columns
#0 - beam sensitivity
#1 - uncorrected RMSE
#2 - uncorrected bias
#3 - natural2 RMSE
#4 - natural2 bias
#5 - hybrid-1 RMSE
#6 - hybrid-1 bias
#7 - hybrid-2 RMSE
#8 - hybrid-2 bias
#9 - mean RMSE
#10 - mean bias
#11 - hybrid-3 RMSE
#12 - hybrid-3 bias

# Read in the text files containing the data for each site
laselva=np.loadtxt(inputLS, comments='#', skiprows=1)
robson=np.loadtxt(inputRC, comments='#', skiprows=1)
nouragues=np.loadtxt(inputN, comments='#', skiprows=1)
paracou=np.loadtxt(inputP, comments='#', skiprows=1)
hubbard=np.loadtxt(inputHB, comments='#', skiprows=1)
wind=np.loadtxt(inputWR, comments='#', skiprows=1)
oak=np.loadtxt(inputOR, comments='#', skiprows=1)

laselva10=np.loadtxt(inputLS10, comments='#', skiprows=1)
robson10=np.loadtxt(inputRC10, comments='#', skiprows=1)
nouragues10=np.loadtxt(inputN10, comments='#', skiprows=1)

sites=[laselva,robson,nouragues,paracou,hubbard,wind,oak]
sites10=[laselva10,robson10,nouragues10]
for site in sites:
    site[site==0.00]=np.nan
for site in sites10:
    site[site==0.00]=np.nan
names=['La Selva', 'Robson Creek', 'Nouragues', 'Paracou', 'Hubbard Brook', 'Wind River', 'Oak Ridge']
colours=['tab:blue', 'tab:red', 'tab:green', 'tab:purple', 'tab:orange', 'tab:olive', 'tab:cyan']

########################

# 30m vs 10m Raw Data
fig,axes=plt.subplots(2,3)
subplots=[axes[0,0],axes[0,1],axes[0,2],axes[1,0],axes[1,1],axes[1,2]]

# RMSE 30m
subplots[0].axhline(y=10,color='darkgray',linestyle='--')
for i in range(3):
    site=sites[i]
    symbol='o'
    subplots[0].plot(site[:,0],site[:,1],label=names[i],linestyle='--',color=colours[i],marker=symbol)
subplots[0].set(title='30m Raw Data',ylabel='RMSE (m)')
subplots[0].set_xlim(91,100)
subplots[0].set_ylim(0,30)
subplots[0].set_xticks(np.arange(91,101,1))
subplots[0].set_yticks(np.arange(0,35,5))
handles,labels=subplots[0].get_legend_handles_labels()
handles=[handles[1],handles[2],handles[0]]
labels=[labels[1],labels[2],labels[0]]
subplots[0].legend()

# RMSE 10m
subplots[1].axhline(y=10,color='darkgray',linestyle='--')
for i in range(3):
    site=sites10[i]
    symbol='o'
    subplots[1].plot(site[:,0],site[:,1],label=names[i],linestyle='--',color=colours[i],marker=symbol)
subplots[1].set(title='10m Raw Data')
subplots[1].set_xlim(91,100)
subplots[1].set_ylim(0,30)
subplots[1].set_xticks(np.arange(91,101,1))
subplots[1].set_yticks(np.arange(0,35,5))

# RMSE 10m Corrected
subplots[2].axhline(y=10,color='darkgray',linestyle='--')
for i in range(3):
    site=sites10[i]
    symbol='o'
    subplots[2].plot(site[:,0],site[:,11],label=names[i],linestyle='-',color=colours[i],marker=symbol)
subplots[2].set(title='10m Corrected Data')
subplots[2].set_xlim(91,100)
subplots[2].set_ylim(0,30)
subplots[2].set_xticks(np.arange(91,101,1))
subplots[2].set_yticks(np.arange(0,35,5))

# Bias 30m
subplots[3].axhline(y=3,color='darkgray',linestyle='--')
subplots[3].axhline(y=-3,color='darkgray',linestyle='--')
for i in range(3):
    site=sites[i]
    symbol='o'
    subplots[3].plot(site[:,0],site[:,2],label=names[i],linestyle='--',color=colours[i],marker=symbol)
subplots[3].set(xlabel='Beam Sensitivity (%)',ylabel='Bias (m)')
subplots[3].set_xlim(91,100)
subplots[3].set_ylim(-4,20)
subplots[3].set_xticks(np.arange(91,101,1))
subplots[3].set_yticks(np.arange(-4,21,2))

# Bias 10m
subplots[4].axhline(y=3,color='darkgray',linestyle='--')
subplots[4].axhline(y=-3,color='darkgray',linestyle='--')
for i in range(3):
    site=sites10[i]
    symbol='o'
    subplots[4].plot(site[:,0],site[:,2],label=names[i],linestyle='--',color=colours[i],marker=symbol)
subplots[4].set(xlabel='Beam Sensitivity (%)')
subplots[4].set_xlim(91,100)
subplots[4].set_ylim(-4,20)
subplots[4].set_xticks(np.arange(91,101,1))
subplots[4].set_yticks(np.arange(-4,21,2))

# Bias 10m Corrected
subplots[5].axhline(y=3,color='darkgray',linestyle='--')
subplots[5].axhline(y=-3,color='darkgray',linestyle='--')
for i in range(3):
    site=sites10[i]
    symbol='o'
    subplots[5].plot(site[:,0],site[:,12],label=names[i],linestyle='-',color=colours[i],marker=symbol)
subplots[5].set(xlabel='Beam Sensitivity (%)')
subplots[5].set_xlim(91,100)
subplots[5].set_ylim(-4,20)
subplots[5].set_xticks(np.arange(91,101,1))
subplots[5].set_yticks(np.arange(-4,21,2))

plt.savefig('../data/summary/hybrid10vs30m.eps',dpi=300)
plt.close()
plt.clf()
print('../data/summary/hybrid10vs30m.eps written to disk')
#plt.show()

########################

# Natural2 Method Plots
fig,axes=plt.subplots(2,3)
subplots=[axes[0,0],axes[0,1],axes[0,2],axes[1,0],axes[1,1],axes[1,2]]

# RMSE Raw
subplots[0].axhline(y=10,color='darkgray',linestyle='--')
for index, site in enumerate(sites):
    symbol='o' if index < 4 else 's'
    subplots[0].plot(site[:,0],site[:,1],label=names[index],linestyle='--',color=colours[index],marker=symbol)
subplots[0].set(title='Raw Data',ylabel='RMSE (m)')
subplots[0].set_xlim(91,100)
subplots[0].set_ylim(0,30)
subplots[0].set_xticks(np.arange(91,101,1))
subplots[0].set_yticks(np.arange(0,35,5))
handles,labels=subplots[0].get_legend_handles_labels()
handles=[handles[1],handles[2],handles[3],handles[0],handles[4],handles[5],handles[6]]
labels=[labels[1],labels[2],labels[3],labels[0],labels[4],labels[5],labels[6]]
subplots[0].legend(ncol=2)

# RMSE Natural2 Corrected
subplots[1].axhline(y=10,color='darkgray',linestyle='--')
for index, site in enumerate(sites):
    symbol='o' if index < 4 else 's'
    subplots[1].plot(site[:,0],site[:,3],label=names[index],linestyle='-',color=colours[index],marker=symbol)
subplots[1].set(title='NaturalRV')
subplots[1].set_xlim(91,100)
subplots[1].set_ylim(0,30)
subplots[1].set_xticks(np.arange(91,101,1))
subplots[1].set_yticks(np.arange(0,35,5))

# RMSE Mean Corrected
subplots[2].axhline(y=10,color='darkgray',linestyle='--')
for index, site in enumerate(sites):
    symbol='o' if index < 4 else 's'
    subplots[2].plot(site[:,0],site[:,9],label=names[index],linestyle='-',color=colours[index],marker=symbol)
subplots[2].set(title='MeanWF')
subplots[2].set_xlim(91,100)
subplots[2].set_ylim(0,30)
subplots[2].set_xticks(np.arange(91,101,1))
subplots[2].set_yticks(np.arange(0,35,5))

# Bias Raw
subplots[3].axhline(y=3,color='darkgray',linestyle='--')
subplots[3].axhline(y=-3,color='darkgray',linestyle='--')
for index, site in enumerate(sites):
    symbol='o' if index < 4 else 's'
    subplots[3].plot(site[:,0],site[:,2],label=names[index],linestyle='--',color=colours[index],marker=symbol)
subplots[3].set(xlabel='Beam Sensitivity (%)',ylabel='Bias (m)')
subplots[3].set_xlim(91,100)
subplots[3].set_ylim(-4,20)
subplots[3].set_xticks(np.arange(91,101,1))
subplots[3].set_yticks(np.arange(-4,21,2))

# Bias Natural2 Corrected
subplots[4].axhline(y=3,color='darkgray',linestyle='--')
subplots[4].axhline(y=-3,color='darkgray',linestyle='--')
for index, site in enumerate(sites):
    symbol='o' if index < 4 else 's'
    subplots[4].plot(site[:,0],site[:,4],label=names[index],linestyle='-',color=colours[index],marker=symbol)
subplots[4].set(xlabel='Beam Sensitivity (%)')
subplots[4].set_xlim(91,100)
subplots[4].set_ylim(-4,20)
subplots[4].set_xticks(np.arange(91,101,1))
subplots[4].set_yticks(np.arange(-4,21,2))

# Bias Mean Corrected
subplots[5].axhline(y=3,color='darkgray',linestyle='--')
subplots[5].axhline(y=-3,color='darkgray',linestyle='--')
for index, site in enumerate(sites):
    symbol='o' if index < 4 else 's'
    subplots[5].plot(site[:,0],site[:,10],label=names[index],linestyle='-',color=colours[index],marker=symbol)
subplots[5].set(xlabel='Beam Sensitivity (%)')
subplots[5].set_xlim(91,100)
subplots[5].set_ylim(-4,20)
subplots[5].set_xticks(np.arange(91,101,1))
subplots[5].set_yticks(np.arange(-4,21,2))

plt.savefig('../data/summary/individualCombined.eps',dpi=300)
plt.close()
plt.clf()
print('../data/summary/individualCombined.eps written to disk')
#plt.show()

########################

# Mean Method Plots
fig,axes=plt.subplots(2,2)
subplots=[axes[0,0],axes[0,1],axes[1,0],axes[1,1]]

# RMSE Raw
subplots[0].axhline(y=10,color='darkgray',linestyle='--')
for index, site in enumerate(sites):
    symbol='o' if index < 4 else 's'
    subplots[0].plot(site[:,0],site[:,1],label=names[index],linestyle='--',color=colours[index],marker=symbol)
subplots[0].set(ylabel='RMSE (m)')
subplots[0].set_xlim(91,100)
subplots[0].set_ylim(0,30)
subplots[0].set_xticks(np.arange(91,101,1))
subplots[0].set_yticks(np.arange(0,35,5))
handles,labels=subplots[0].get_legend_handles_labels()
handles=[handles[1],handles[2],handles[3],handles[0],handles[4],handles[5],handles[6]]
labels=[labels[1],labels[2],labels[3],labels[0],labels[4],labels[5],labels[6]]
subplots[0].legend()

# RMSE Corrected
subplots[1].axhline(y=10,color='darkgray',linestyle='--')
for index, site in enumerate(sites):
    symbol='o' if index < 4 else 's'
    subplots[1].plot(site[:,0],site[:,9],label=names[index],linestyle='-',color=colours[index],marker=symbol)
#axes[1].set(xlabel='Beam Sensitivity (%)')
subplots[1].set_xlim(91,100)
subplots[1].set_ylim(0,30)
subplots[1].set_xticks(np.arange(91,101,1))
subplots[1].set_yticks(np.arange(0,35,5))

# Bias Raw
subplots[2].axhline(y=3,color='darkgray',linestyle='--')
subplots[2].axhline(y=-3,color='darkgray',linestyle='--')
for index, site in enumerate(sites):
    symbol='o' if index < 4 else 's'
    subplots[2].plot(site[:,0],site[:,2],label=names[index],linestyle='--',color=colours[index],marker=symbol)
subplots[2].set(xlabel='Beam Sensitivity (%)',ylabel='Bias (m)')
subplots[2].set_xlim(91,100)
subplots[2].set_ylim(-4,20)
subplots[2].set_xticks(np.arange(91,101,1))
subplots[2].set_yticks(np.arange(-4,21,2))

# Bias Corrected
subplots[3].axhline(y=3,color='darkgray',linestyle='--')
subplots[3].axhline(y=-3,color='darkgray',linestyle='--')
for index, site in enumerate(sites):
    symbol='o' if index < 4 else 's'
    subplots[3].plot(site[:,0],site[:,10],label=names[index],linestyle='-',color=colours[index],marker=symbol)
subplots[3].set(xlabel='Beam Sensitivity (%)')
subplots[3].set_xlim(91,100)
subplots[3].set_ylim(-4,20)
subplots[3].set_xticks(np.arange(91,101,1))
subplots[3].set_yticks(np.arange(-4,21,2))

#plt.savefig('../data/summary/meanCombined.png',dpi=300)
#plt.close()
#plt.clf()
#print('../data/summary/meanCombined.png written to disk')
#plt.show()

########################

# Hybrid-3 Method Plots
fig,axes=plt.subplots(2,2,figsize=(6,4))
subplots=[axes[0,0],axes[0,1],axes[1,0],axes[1,1]]

# RMSE Raw
subplots[0].axhline(y=10,color='darkgray',linestyle='--')
for index, site in enumerate(sites):
    symbol='o' if index < 4 else 's'
    subplots[0].plot(site[:,0],site[:,1],label=names[index],linestyle='--',color=colours[index],marker=symbol)
subplots[0].set(title='Raw Data',ylabel='RMSE (m)')
subplots[0].set_xlim(91,100)
subplots[0].set_ylim(0,30)
subplots[0].set_xticks(np.arange(91,101,1))
subplots[0].set_yticks(np.arange(0,35,5))
handles,labels=subplots[0].get_legend_handles_labels()
handles=[handles[1],handles[2],handles[3],handles[0],handles[4],handles[5],handles[6]]
labels=[labels[1],labels[2],labels[3],labels[0],labels[4],labels[5],labels[6]]
subplots[0].legend(ncol=2)

# RMSE Corrected
subplots[1].axhline(y=10,color='darkgray',linestyle='--')
for index, site in enumerate(sites):
    symbol='o' if index < 4 else 's'
    subplots[1].plot(site[:,0],site[:,11],label=names[index],linestyle='-',color=colours[index],marker=symbol)
subplots[1].set(title='Hybrid')
subplots[1].set_xlim(91,100)
subplots[1].set_ylim(0,30)
subplots[1].set_xticks(np.arange(91,101,1))
subplots[1].set_yticks(np.arange(0,35,5))

# Bias Raw
subplots[2].axhline(y=3,color='darkgray',linestyle='--')
subplots[2].axhline(y=-3,color='darkgray',linestyle='--')
for index, site in enumerate(sites):
    symbol='o' if index < 4 else 's'
    subplots[2].plot(site[:,0],site[:,2],label=names[index],linestyle='--',color=colours[index],marker=symbol)
subplots[2].set(xlabel='Beam Sensitivity (%)',ylabel='Bias (m)')
subplots[2].set_xlim(91,100)
subplots[2].set_ylim(-4,20)
subplots[2].set_xticks(np.arange(91,101,1))
subplots[2].set_yticks(np.arange(-4,21,2))

# Bias Corrected
subplots[3].axhline(y=3,color='darkgray',linestyle='--')
subplots[3].axhline(y=-3,color='darkgray',linestyle='--')
for index, site in enumerate(sites):
    symbol='o' if index < 4 else 's'
    subplots[3].plot(site[:,0],site[:,12],label=names[index],linestyle='-',color=colours[index],marker=symbol)
subplots[3].set(xlabel='Beam Sensitivity (%)')
subplots[3].set_xlim(91,100)
subplots[3].set_ylim(-4,20)
subplots[3].set_xticks(np.arange(91,101,1))
subplots[3].set_yticks(np.arange(-4,21,2))

plt.savefig('../data/summary/hybrid3Combined.eps',dpi=300)
plt.close()
plt.clf()
print('../data/summary/hybrid3Combined.eps written to disk')
#plt.show()

########################

# Combined Plots
'''fig,axes=plt.subplots(2,2,figsize=(8,6),sharex=True,sharey=True)
subplots=[axes[0,0],axes[0,1],axes[1,0],axes[1,1]]
# Bias Plot
for index, site in enumerate(sites):
    subplots[0].plot(site[:,0],site[:,2],label=names[index]+' Raw',linestyle='--',color=colours[index],marker='o')
    subplots[0].plot(site[:,0],site[:,4],label=names[index]+' Corrected',linestyle='-',color=colours[index],marker='o')
    subplots[1].plot(site[:,0],site[:,2],label=names[index]+' Raw',linestyle='--',color=colours[index],marker='o')
    subplots[1].plot(site[:,0],site[:,10],label=names[index]+' Corrected',linestyle='-',color=colours[index],marker='o')
    subplots[2].plot(site[:,0],site[:,2],label=names[index]+' Raw',linestyle='--',color=colours[index],marker='o')
    subplots[2].plot(site[:,0],site[:,6],label=names[index]+' Corrected',linestyle='-',color=colours[index],marker='o')
    subplots[3].plot(site[:,0],site[:,2],label=names[index]+' Raw',linestyle='--',color=colours[index],marker='o')
    subplots[3].plot(site[:,0],site[:,12],label=names[index]+' Corrected',linestyle='-',color=colours[index],marker='o')
for axis in axes.flat:
    axis.set(xlabel='Beam Sensitivity (%)',ylabel='Bias (m)')
    axis.set_xlim(91,99)
    axis.set_ylim(-2,20)
    axis.set_yticks(np.arange(-2,22,2))
for axis in axes.flat:
    axis.label_outer()
#plt.tight_layout()
#plt.savefig('../data/summary/biasAllHybrid1&3.png',dpi=300)
#plt.close()
#plt.clf()
#print('../data/summary/biasAllHybrid1&3.png written to disk')

fig,axes=plt.subplots(2,2,figsize=(8,6),sharex=True,sharey=True)
subplots=[axes[0,0],axes[0,1],axes[1,0],axes[1,1]]
# Bias Plot
for index, site in enumerate(sites):
    subplots[0].plot(site[:,0],site[:,1],label=names[index]+' Raw',linestyle='--',color=colours[index],marker='o')
    subplots[0].plot(site[:,0],site[:,3],label=names[index]+' Corrected',linestyle='-',color=colours[index],marker='o')
    subplots[1].plot(site[:,0],site[:,1],label=names[index]+' Raw',linestyle='--',color=colours[index],marker='o')
    subplots[1].plot(site[:,0],site[:,9],label=names[index]+' Corrected',linestyle='-',color=colours[index],marker='o')
    subplots[2].plot(site[:,0],site[:,1],label=names[index]+' Raw',linestyle='--',color=colours[index],marker='o')
    subplots[2].plot(site[:,0],site[:,5],label=names[index]+' Corrected',linestyle='-',color=colours[index],marker='o')
    subplots[3].plot(site[:,0],site[:,1],label=names[index]+' Raw',linestyle='--',color=colours[index],marker='o')
    subplots[3].plot(site[:,0],site[:,11],label=names[index]+' Corrected',linestyle='-',color=colours[index],marker='o')
for axis in axes.flat:
    axis.set(xlabel='Beam Sensitivity (%)',ylabel='RMSE (m)')
    axis.set_xlim(91,99)
    axis.set_ylim(0,30)
    axis.set_yticks(np.arange(0,32,5))
for axis in axes.flat:
    axis.label_outer()
#plt.tight_layout()
#plt.savefig('../data/summary/rmseAllHybrid1&3.png',dpi=300)
#plt.close()
#plt.clf()
#print('../data/summary/rmseAllHybrid1&3.png written to disk')'''

########################

'''# Hybrid-1 Method Plots
fig,axes=plt.subplots(1,3,figsize=(12,4),gridspec_kw={'width_ratios':[3,3,1]})
# Bias Plot
for index, site in enumerate(sites):
    axes[0].plot(site[:,0],site[:,2],label=names[index]+' Raw',linestyle='--',color=colours[index],marker='o')
    axes[0].plot(site[:,0],site[:,6],label=names[index]+' Corrected',linestyle='-',color=colours[index],marker='o')
axes[0].set_title("'Hybrid-1' Beam Sensitivity vs. Bias")
axes[0].set(xlabel='Beam Sensitivity (%)',ylabel='Bias (m)')
axes[0].set_xlim(91,99)
axes[0].set_ylim(-2,20)
axes[0].set_yticks(np.arange(-2,22,2))
#axes[0].legend(bbox_to_anchor=(1,2),loc='upper left')

# RMSE Plot
for index, site in enumerate(sites):
    axes[1].plot(site[:,0],site[:,1],label=names[index]+' Raw',linestyle='--',color=colours[index],marker='o')
    axes[1].plot(site[:,0],site[:,5],label=names[index]+' Corrected',linestyle='-',color=colours[index],marker='o')
axes[1].set_title("'Hybrid-1' Beam Sensitivity vs. RMSE")
axes[1].set(xlabel='Beam Sensitivity (%)',ylabel='RMSE (m)')
axes[1].set_xlim(91,99)
axes[1].set_ylim(0,30)
axes[1].set_yticks(np.arange(0,32,5))
axes[1].legend(bbox_to_anchor=(1,1),loc='upper left')

# Legend
axes[2].set_frame_on(False)
axes[2].get_xaxis().set_visible(False)
axes[2].get_yaxis().set_visible(False)
#plt.savefig('../data/summary/hybrid1.png',dpi=300)
#plt.close()
#plt.clf()
#print('../data/summary/hybrid1.png written to disk')

########################

# Hybrid-2 Method Plots
fig,axes=plt.subplots(1,3,figsize=(12,4),gridspec_kw={'width_ratios':[3,3,1]})
# Bias Plot
for index, site in enumerate(sites):
    axes[0].plot(site[:,0],site[:,2],label=names[index]+' Raw',linestyle='--',color=colours[index],marker='o')
    axes[0].plot(site[:,0],site[:,8],label=names[index]+' Corrected',linestyle='-',color=colours[index],marker='o')
axes[0].set_title("'Hybrid-2' Beam Sensitivity vs. Bias")
axes[0].set(xlabel='Beam Sensitivity (%)',ylabel='Bias (m)')
axes[0].set_xlim(91,99)
axes[0].set_ylim(-2,20)
axes[0].set_yticks(np.arange(-2,22,2))
#axes[0].legend(bbox_to_anchor=(1,2),loc='upper left')

# RMSE Plot
for index, site in enumerate(sites):
    axes[1].plot(site[:,0],site[:,1],label=names[index]+' Raw',linestyle='--',color=colours[index],marker='o')
    axes[1].plot(site[:,0],site[:,7],label=names[index]+' Corrected',linestyle='-',color=colours[index],marker='o')
axes[1].set_title("'Hybrid-2' Beam Sensitivity vs. RMSE")
axes[1].set(xlabel='Beam Sensitivity (%)',ylabel='RMSE (m)')
axes[1].set_xlim(91,99)
axes[1].set_ylim(0,30)
axes[1].set_yticks(np.arange(0,32,5))
axes[1].legend(bbox_to_anchor=(1,1),loc='upper left')

# Legend
axes[2].set_frame_on(False)
axes[2].get_xaxis().set_visible(False)
axes[2].get_yaxis().set_visible(False)
#plt.savefig('../data/summary/hybrid2.png',dpi=300)
#plt.close()
#plt.clf()
#print('../data/summary/hybrid2.png written to disk')

########################'''
