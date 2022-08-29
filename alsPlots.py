############################################
# Quick script to load a single metric file
# per site to generate plots of ALS data.
############################################

import numpy as np
import matplotlib.pyplot as plt
from rasterTools import readFromFile

#########################################

# Read gediMetric output

rootDir='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/'
inputs=[
    rootDir+'laselva/grid30/metric98-1/metric30-98-1Summary.txt',
    rootDir+'robson/grid30/metric98-1/metric30-98-1Summary.txt',
    rootDir+'nouragues/grid30/metric98-1/metric30-98-1Summary.txt',
    rootDir+'paracou/grid30/metric98-1/metric30-98-1Summary.txt',
    rootDir+'hubbard/grid30/metric98-1/metric30-98-1Summary.txt',
    rootDir+'windriver/grid30/metric98-1/metric30-98-1Summary.txt',
    rootDir+'oakridge/grid30/metric98-1/metric30-98-1Summary.txt'
]
inputs10=[
    rootDir+'laselva/grid10/metric98-1/metric10-98-1Summary.txt',
    rootDir+'robson/grid10/metric98-1/metric10-98-1Summary.txt',
    rootDir+'nouragues/grid10/metric98-1/metric10-98-1Summary.txt'
]

lsData=np.loadtxt(inputs[0], comments='#', skiprows=1)              # Columns:
rcData=np.loadtxt(inputs[1], comments='#', skiprows=1)              # 2 = trueGround
noData=np.loadtxt(inputs[2], comments='#', skiprows=1)              # 4 = ALScover
paData=np.loadtxt(inputs[3], comments='#', skiprows=1)              # 5 = rhReal95
huData=np.loadtxt(inputs[4], comments='#', skiprows=1)              #10 = rhReal50
wrData=np.loadtxt(inputs[5], comments='#', skiprows=1)
orData=np.loadtxt(inputs[6], comments='#', skiprows=1)

lsData10=np.loadtxt(inputs10[0], comments='#', skiprows=1)
rcData10=np.loadtxt(inputs10[1], comments='#', skiprows=1)
noData10=np.loadtxt(inputs10[2], comments='#', skiprows=1)

data=[lsData,rcData,noData,paData,huData,wrData,orData]
data10=[lsData10,rcData10,noData10]

#########################################

# Read slope rasters

slopeInputs=[
    '../data/laSelvaGroundError/30-98/ALSGroundSlope.tif',
    '../data/robsonGroundError/30-98/ALSGroundSlope.tif',
    '../data/nouraguesGroundError/30-98/groundErrorALSSlope.tif',
    '../data/paracouGroundError/30-98/groundErrorALSSlope.tif',
    '../data/hubbardGroundError/30-98/ALSGroundSlope.tif',
    '../data/windRiverGroundError/30-98/ALSGroundSlope.tif',
    '../data/oakRidgeGroundError/30-98/ALSGroundSlope.tif'
]

lsSlope,_=readFromFile(slopeInputs[0])
robsonSlope,_=readFromFile(slopeInputs[1])
nouraguesSlope,_=readFromFile(slopeInputs[2])
paracouSlope,_=readFromFile(slopeInputs[3])
hubbardSlope,_=readFromFile(slopeInputs[4])
windSlope,_=readFromFile(slopeInputs[5])
oakSlope,_=readFromFile(slopeInputs[6])

slopes=[lsSlope,robsonSlope,nouraguesSlope,paracouSlope,hubbardSlope,windSlope,oakSlope]

#########################################

# Print stats for each site

names=['La Selva', 'Robson Creek', 'Nouragues', 'Paracou', 'Hubbard Brook', 'Wind River', 'Oak Ridge']
for index, datum in enumerate(data):
    meanCover=np.mean(datum[:,4])
    sdCover=np.std(datum[:,4])
    meanHeight=np.mean(datum[:,5])
    sdHeight=np.std(datum[:,5])
    meanSlope=np.nanmean(slopes[index][slopes[index]>-999])
    sdSlope=np.nanstd(slopes[index][slopes[index]>-999])
    print(names[index])
    print('n Samples:',datum.shape[0])
    print('Canopy cover mean and std. dev.: {0:.2f}, {1:.2f}'.format(meanCover*100,sdCover*100))
    print('Canopy cover > 99%: {0:.2f}%'.format((len(datum[:,4][datum[:,4]>0.99])/datum.shape[0])*100))
    print('Canopy cover > 98%: {0:.2f}%'.format((len(datum[:,4][datum[:,4]>0.98])/datum.shape[0])*100))
    print('Canopy height mean and std. dev.: {0:.2f}, {1:.2f}'.format(meanHeight,sdHeight))
    print('Slope mean and std. dev.: {0:.2f}, {1:.2f}'.format(meanSlope,sdSlope))
    print(' ')
for index, datum in enumerate(data10):
    meanCover=np.mean(datum[:,4])
    sdCover=np.std(datum[:,4])
    meanHeight=np.mean(datum[:,5])
    sdHeight=np.std(datum[:,5])
    print(names[index]+' 10m grid data')
    print('n Samples:',datum.shape[0])
    print('Canopy cover mean and std. dev.: {0:.2f}, {1:.2f}'.format(meanCover*100,sdCover*100))
    print('Canopy cover > 99%: {0:.2f}%'.format((len(datum[:,4][datum[:,4]>0.99])/datum.shape[0])*100))
    print('Canopy cover > 98%: {0:.2f}%'.format((len(datum[:,4][datum[:,4]>0.98])/datum.shape[0])*100))
    print('Canopy height mean and std. dev.: {0:.2f}, {1:.2f}'.format(meanHeight,sdHeight))
    print(' ')

#########################################

# Plotting settings
plt.rcParams['figure.constrained_layout.use']=True
plt.rcParams['figure.figsize']=(4.5,8)
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
plt.rcParams['lines.linewidth']=0.1
plt.rcParams['lines.markersize']=2
#plt.rcParams['savefig.format']='tif'

#########################################

# Canopy Cover Full Range
fig,axes=plt.subplots(4,2,sharex=True)
subplots=[axes[0,0],axes[1,0],axes[2,0],axes[3,0],axes[0,1],axes[1,1],axes[2,1]]
for index, subplot in enumerate(subplots):
    subplot.hist(data[index][:,4],bins=np.arange(0.0,1.0+0.1,0.02),density=True,edgecolor='black',linewidth=0.1,label=names[index])
    subplot.legend()
for axis in axes.flat:
    axis.set(ylabel='Normalised Frequency')
#for axis in axes.flat:
    #axis.label_outer()
axes[3,0].set(xlabel='Canopy Cover')
axes[3,1].set_frame_on(False)
axes[3,1].get_xaxis().set_visible(False)
axes[3,1].get_yaxis().set_visible(False)
plt.savefig('../data/coverHistAll.png',dpi=300)
plt.close()
plt.clf()
print('../data/coverHistAll.png written to disk')
#plt.show()

# Canopy Cover High
fig,axes=plt.subplots(4,2,sharex=True)
subplots=[axes[0,0],axes[1,0],axes[2,0],axes[3,0],axes[0,1],axes[1,1],axes[2,1]]
for index, subplot in enumerate(subplots):
    subplot.hist(data[index][:,4],bins=np.arange(0.6,1.0+0.005,0.005),density=True,edgecolor='black',linewidth=0.1,label=names[index])
    subplot.legend()
for axis in axes.flat:
    axis.set(ylabel='Normalised Frequency')
#for axis in axes.flat:
    #axis.label_outer()
axes[3,0].set(xlabel='Canopy Cover')
axes[3,1].set_frame_on(False)
axes[3,1].get_xaxis().set_visible(False)
axes[3,1].get_yaxis().set_visible(False)
plt.savefig('../data/coverHistHighAll.png',dpi=300)
plt.close()
plt.clf()
print('../data/coverHistHighAll.png written to disk')
#plt.show()

# Canopy Height
fig,axes=plt.subplots(4,2,sharex=True)
subplots=[axes[0,0],axes[1,0],axes[2,0],axes[3,0],axes[0,1],axes[1,1],axes[2,1]]
for index, subplot in enumerate(subplots):
    subplot.hist(data[index][:,5],bins=np.arange(0,70+1,1),density=True,edgecolor='black',linewidth=0.1,label=names[index])
    subplot.legend()
for axis in axes.flat:
    axis.set(ylabel='Normalised Frequency')
#for axis in axes.flat:
    #axis.label_outer()
axes[3,0].set(xlabel='Canopy Height (m)')
axes[3,1].set_frame_on(False)
axes[3,1].get_xaxis().set_visible(False)
axes[3,1].get_yaxis().set_visible(False)
plt.savefig('../data/heightHistAll.png',dpi=300)
plt.close()
plt.clf()
print('../data/heightHistAll.png written to disk')
#plt.show()

# Slope
fig,axes=plt.subplots(4,2,sharex=True)
subplots=[axes[0,0],axes[1,0],axes[2,0],axes[3,0],axes[0,1],axes[1,1],axes[2,1]]
for index, subplot in enumerate(subplots):
    subplot.hist(slopes[index][slopes[index]>-999].flatten(),bins=np.arange(0,60+1,1),density=True,edgecolor='black',linewidth=0.1,label=names[index])
    subplot.legend()
for axis in axes.flat:
    axis.set(ylabel='Normalised Frequency')
#for axis in axes.flat:
    #axis.label_outer()
axes[3,0].set(xlabel='Ground Slope (deg)')
axes[3,1].set_frame_on(False)
axes[3,1].get_xaxis().set_visible(False)
axes[3,1].get_yaxis().set_visible(False)
plt.savefig('../data/slopeHistAll.png',dpi=300)
plt.close()
plt.clf()
print('../data/slopeHistAll.png written to disk')
#plt.show()
