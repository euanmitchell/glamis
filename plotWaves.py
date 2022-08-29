############################################
# Script to read HDF5 output files from
# gediRat/addNoiseHDF and plot waveforms.
############################################

import h5py
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import glob

#########################################

# Defining the command line reading function
def readCommands():
    '''
    Read the arguments passed from the command line
    '''
    p=argparse.ArgumentParser(description=('Specify parameters for program control'),
        formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('--site', dest='site', type=str, default='robson',
        help=('''Site name. Valid values: laselva, robson, nouragues, paracou.'''))
    p.add_argument('--gridSize', dest='gridSize', type=int, default=30,
        help=('The spacing between footprints. Default is 30'))
    cmdargs=p.parse_args()
    return cmdargs

#########################################

class groundComparison(object):
    '''A class to plot waveforms with various ground estimations'''

    def __init__(self):
        '''Initialise the class'''

    def readNoisedWaves(self,input):
        '''Read the HDF5 file containing the raw waveforms with added noise'''

        # Open the HDF5 file and extract data to arrays
        print('Reading HDF5 file {} ...'.format(input))
        with h5py.File(input,'r') as f:
            self.bins=int(np.array(f['NBINS']))
            self.lats=np.array(f['LAT0'])
            self.lons=np.array(f['LON0'])
            self.z0=np.array(f['Z0'])
            self.zN=np.array(f['ZN'])
            self.zG=np.array(f['ZG'])
            self.waves=np.array(f['RXWAVECOUNT'])

    def readMetricFiles(self,rootDir):
        '''Read the gediMetric output files for the various algorithm options'''

        print('Reading gediMetric output files ...')
        metric1=rootDir+'metric98-1/metric{0}-98-1Summary.txt'.format(args.gridSize)
        metric2=rootDir+'metric98-2/metric{0}-98-2Summary.txt'.format(args.gridSize)
        metric3=rootDir+'metric98-3/metric{0}-98-3Summary.txt'.format(args.gridSize)
        metric4=rootDir+'metric98-4/metric{0}-98-4Summary.txt'.format(args.gridSize)
        metric5=rootDir+'metric98-5/metric{0}-98-5Summary.txt'.format(args.gridSize)
        metric6=rootDir+'metric98-6/metric{0}-98-6Summary.txt'.format(args.gridSize)
        metric7=rootDir+'metric98-7/metric{0}-98-7Summary.txt'.format(args.gridSize)
        metric8=rootDir+'metric98-8/metric{0}-98-8Summary.txt'.format(args.gridSize)
        metric9=rootDir+'metric98-9/metric{0}-98-9Summary.txt'.format(args.gridSize)
        metric10=rootDir+'metric98-10/metric{0}-98-10Summary.txt'.format(args.gridSize)
        metric150=rootDir+'metric98-150/metric{0}-98-150Summary.txt'.format(args.gridSize)

        data1=np.loadtxt(metric1, comments='#', skiprows=1)
        data2=np.loadtxt(metric2, comments='#', skiprows=1)
        data3=np.loadtxt(metric3, comments='#', skiprows=1)
        data4=np.loadtxt(metric4, comments='#', skiprows=1)
        data5=np.loadtxt(metric5, comments='#', skiprows=1)
        data6=np.loadtxt(metric6, comments='#', skiprows=1)
        data7=np.loadtxt(metric7, comments='#', skiprows=1)
        data8=np.loadtxt(metric8, comments='#', skiprows=1)
        data9=np.loadtxt(metric9, comments='#', skiprows=1)
        data10=np.loadtxt(metric10, comments='#', skiprows=1)
        data150=np.loadtxt(metric150, comments='#', skiprows=1)
        self.metricArrays=[data1,data2,data3,data4,data5,data6,data7,data8,data9,data150,data10]

    def plotWaveform(self):
        '''Plot the noised waves and the corresponding inferred ground elevations'''

        # Plotting settings
        plt.rcParams['figure.figsize']=(3,2.4)
        plt.rcParams['legend.fontsize']=5
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
        plt.rcParams['axes.titlesize']=8
        plt.rcParams['lines.linewidth']=0.5
        #plt.rcParams['savefig.format']='tif'

        counter=0
        print('Writing matching waveforms ...')
        for i in range(self.lats.shape[0]):
            if self.lats[i]==8106540 and self.lons[i]==351010:
                # Get the lat/lon of the next wave
                lat=self.lats[i]
                lon=self.lons[i]

                # Read in a text file containing lat/lon of outlier footprints
                inFile='../data/robsonGroundError/30-98/naturalSD1F1OutlierCoords.txt'
                with open(inFile,'r') as textfile:
                    outliers=textfile.readlines()[1:]
                    for outlier in outliers:
                        outLon=float(outlier.split()[0])
                        outLat=float(outlier.split()[1])
                        outNG=float(outlier.split()[2])
                        if (lat==outLat and lon==outLon):

                            # Get the nine ground options for this wave
                            groundOptions=[]
                            for array in self.metricArrays:
                                useInd=np.where((array[:,0]==lon) & (array[:,1]==lat))
                                elevation=array[useInd,6]
                                groundOptions.append(elevation)
                            print('Number of signal processing options:',len(groundOptions))

                            # Plot the waveform and options
                            z=np.linspace(self.z0[i],self.zN[i],num=self.bins)
                            wave=self.waves[i]
                            plt.plot(wave,z,label='Noised Waveform')
                            for index, ground in enumerate(groundOptions):
                                try:
                                    if index==0:
                                        plt.axhline(y=ground,color='gray',linestyle='--',label='Ground Estimates')
                                    else:
                                        plt.axhline(y=ground,color='gray',linestyle='--')
                                except:
                                    continue
                            plt.axhline(y=self.zG[i],color='red',label='ALS Ground')
                            #plt.axhline(y=outNG,color='r',label='newGround',linewidth=1)
                            plt.title('Footprint: {0:.0f}, {1:.0f}'.format(lon,lat))
                            plt.xlabel('DN')
                            plt.ylabel('Elevation (m)')
                            plt.legend()
                            #plt.tight_layout()
                            plt.savefig('../data/exampleWaveform.eps',dpi=300)
                            plt.close()
                            plt.clf()
                            counter+=1
                            #plt.show()
        print('{} waveforms written to disk'.format(counter))


#########################################
# The main block
if __name__ == '__main__':
    args=readCommands()
    plots=groundComparison()

    file='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/{0}/grid{1}/addNoise/350500.8105400.h5'.format(args.site,args.gridSize)
    waveforms='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/{0}/grid{1}/addNoise/'.format(args.site,args.gridSize)
    metrics='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/{0}/grid{1}/'.format(args.site,args.gridSize)

    inputs='*.h5'
    path=os.path.join(waveforms,inputs)
    fileList=glob.glob(path)

    #for file in fileList:
    plots.readNoisedWaves(file)
    plots.readMetricFiles(metrics)
    plots.plotWaveform()
