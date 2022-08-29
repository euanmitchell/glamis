#######################################
# Script to read gediMetric text files
# and investigate spatial distribution
# of canopy cover statistics.
#######################################

#import h5py
import numpy as np
import matplotlib.pyplot as plt
from pyproj import Proj, transform
from osgeo import gdal
from osgeo import osr
from math import sqrt
from rasterTools import writeToFile, defineRegions
import argparse

#######################################

# Defining the command line reading function
def readCommands():
    '''
    Read the arguments passed from the command line
    '''
    p=argparse.ArgumentParser(description=('Specify parameters for program control'),
        formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('--res', dest='res', type=int, default=30,
        help=('The grid size of input data and resolution for output GeoTiffs. Default is 30.'))
    p.add_argument('--epsg', dest='epsg', type=int, default=4326,
        help=('The EPSG code for input data and output GeoTiffs. Default is 4326 (WGS84).'))
    p.add_argument('--stats', dest='stats', action='store_true', default=False,
        help=('Call loadData() to generate data size and stats print out. Requires --res to read correct input files.'))
    p.add_argument('--noImage', dest='noImage', action='store_true', default=False,
        help=('Skip writing the raw image array (i.e., canopy cover raster) to file.'))
    p.add_argument('--binary', dest='binary', action='store_true', default=False,
        help=('Generate binary raster based on threshold cover value specified by "threshold".'))
    p.add_argument('--regions', dest='regions', action='store_true', default=False,
        help=('Call computeRegions() to generate region raster based on threshold cover value specified by "threshold".'))
    p.add_argument('--threshold', dest='threshold', type=float, default=0.98,
        help=('The threshold value for generating binary rasters. Output will be zero for cells less than threshold and one for cells greater. Default is 0.98.'))
    cmdargs=p.parse_args()
    return cmdargs

#######################################

class coverStats(object):
    '''A class to handle reading and plotting ALS and GEDI ground estimates for multiple sites'''

    def __init__(self):
        '''Initialise the class'''

    def loadData(self,dir):
        '''Read in data from gediMetric output text files'''

        # Columns:
        # 0 = x
        # 1 = y
        # 2 = trueGround
        # 3 = beam sensitivity
        # 4 = ALScover
        # 5 = rhReal95

        if args.res==30:
            grid30Input=dir+'grid30/metric/grid30MetricSummary.txt'
            self.grid30Data=np.loadtxt(grid30Input,comments='#',skiprows=1)
            self.useInd30=np.where((self.grid30Data[:,4] >= args.threshold))
            print('Grid30 shape:',self.grid30Data.shape)
            print('Grid30 threshold shape:',self.grid30Data[self.useInd30].shape)
            print('{0:.1f}% of 30m grid waveforms are greater than {1}% cover'.format(((self.grid30Data[self.useInd30].shape[0]/self.grid30Data.shape[0])*100),args.threshold*100))

        if args.res==20:
            grid20Input=dir+'grid20/metric/grid20MetricSummary.txt'
            self.grid20Data=np.loadtxt(grid20Input,comments='#',skiprows=1)
            self.useInd20=np.where((self.grid20Data[:,4] >= args.threshold))
            print('Grid20 shape:',self.grid20Data.shape)
            print('Grid20 threshold shape:',self.grid20Data[self.useInd20].shape)
            print('{0:.1f}% of 20m grid waveforms are greater than {1}% cover'.format(((self.grid20Data[self.useInd20].shape[0]/self.grid20Data.shape[0])*100),args.threshold*100))

        if args.res==10:
            grid10Input=dir+'grid10/metric/grid10MetricSummary.txt'
            self.grid10Data=np.loadtxt(grid10Input,comments='#',skiprows=1)
            self.useInd10=np.where((self.grid10Data[:,4] >= args.threshold))
            print('Grid10 shape:',self.grid10Data.shape)
            print('Grid10 threshold shape:',self.grid10Data[self.useInd10].shape)
            print('{0:.1f}% of 10m grid waveforms are greater than {1}% cover'.format(((self.grid10Data[self.useInd10].shape[0]/self.grid10Data.shape[0])*100),args.threshold*100))

        if args.res==5:
            grid5Input=dir+'2019/grid5/metric/grid5_2019MetricSummary.txt'
            self.grid5Data=np.loadtxt(grid5Input,comments='#',skiprows=1)
            self.useInd5=np.where((self.grid5Data[:,4] >= args.threshold))
            print('Grid5 shape:',self.grid5Data.shape)
            print('Grid5 threshold shape:',self.grid5Data[self.useInd5].shape)
            print('{0:.1f}% of 5m grid waveforms are greater than {1}% cover'.format(((self.grid5Data[self.useInd5].shape[0]/self.grid5Data.shape[0])*100),args.threshold*100))

    def plotCover(self):
        plt.figure(figsize=(10,8))
        plt.scatter(self.grid10Data[:,0],self.grid10Data[:,1],s=0.1,c=self.grid10Data[:,4],cmap='Greens')
        plt.colorbar()
        plt.scatter(self.grid10Data[self.useInd10,0],self.grid10Data[self.useInd10,1],s=0.1,c='r',label='Cover > 98%')
        plt.title('10m Grid by ALS Cover')
        plt.xlabel('Easting (m)')
        plt.ylabel('Northing (m)')
        plt.legend()
        plt.show()

    def writeGeoTiff(self,input,res,epsg,outImage,outBinary,outRegion):
        '''Write canopy cover values to GeoTIFF of chosen resolution'''

        # Read in specified input data file
        data=np.loadtxt(input, comments='#', skiprows=1)

        # Determine data bounds and image dimensions
        lons=data[:,0]
        lats=data[:,1]
        covers=data[:,4]
        self.minX=np.min(lons)
        self.maxX=np.max(lons)
        self.nX=int((self.maxX-self.minX)/res)+1
        self.minY=np.min(lats)
        self.maxY=np.max(lats)
        self.nY=int((self.maxY-self.minY)/res)+1

        # Make array of no data values
        self.imageArray=np.full((self.nY,self.nX),-999.0)

        # Determine pixel index of each latitude and longitude
        xInds=((lons[:]-(self.minX-(args.res/2)))//res).astype(np.int)
        yInds=(((self.maxY+(args.res/2))-lats[:])//res).astype(np.int)

        # Fill the image array with the average data value at each pixel
        counter=np.full((self.nY,self.nX),0)
        totals=np.full((self.nY,self.nX),0.0)
        for i in range(covers.shape[0]):
            y=yInds[i]
            x=xInds[i]
            value=covers[i]
            if self.imageArray[y,x]==-999.0:
                self.imageArray[y,x]=value
                totals[y,x]=value
                counter[y,x]+=1
            else:
                counter[y,x]+=1
                totals[y,x]=totals[y,x]+value
                value=totals[y,x]/counter[y,x]
                self.imageArray[y,x]=value

        # Get size of the image portion of the image array (i.e., not NoData value cells)
        self.imageData=np.sum(np.where(self.imageArray==-999.0,0,1))

        # Generate binary raster (0 = cover < 0.98; 1 = cover >= 0.98)
        self.binaryArray=np.full((self.nY,self.nX),-999.0)
        self.binaryArray[self.imageArray<args.threshold]=0.0
        self.binaryArray[self.imageArray>=args.threshold]=1.0

        # Write to file
        if args.binary:
            writeToFile(self.binaryArray,self.minX-(args.res/2),self.maxY+(args.res/2),self.nX,self.nY,res,epsg,outBinary)
        if args.noImage:
            pass
        else:
            writeToFile(self.imageArray,self.minX-(args.res/2),self.maxY+(args.res/2),self.nX,self.nY,res,epsg,outImage)

        # Define regions with imported method and write to file
        if args.regions:
            regions=defineRegions(self.imageArray,self.binaryArray,res,epsg)
            writeToFile(regions,self.minX-(args.res/2),self.maxY+(args.res/2),self.nX,self.nY,res,epsg,outRegion)

    def summaryStatsPlot(self):
        '''Plot summary statistics from region analysis'''

        gridSize=[10,20,30]

        # Percent coverage of 2x2 squares at 95% cover threshold for 10/20/30 m grids
        thresh95_4=[35.7372,29.8603,24.7648]    #[] - RC
        thresh96_4=[25.1564,15.8594,9.8375]     #[]
        thresh97_4=[13.9035,4.3740,1.6920]      #[]
        thresh98_4=[4.8447,0.5517,0.0779]       #[52.0304,45.4753,41.6976]
        thresh99_4=[0.6448,0.0198,0]            #[18.6500,10.0857,6.0008]

        thresh95_9=[16.2388,11.7544,7.6696]     #[]
        thresh96_9=[7.4632,3.3053,1.4582]       #[]
        thresh97_9=[2.1778,0.3124,0.1057]       #[]
        thresh98_9=[0.3060,0.0174,0.0]          #[32.9454,25.9165,22.9980]
        thresh99_9=[0.0056,0.0,0.0]             #[6.2778,2.1035,0.9176]

        thresh95_16=[4.4269,2.6643,1.3803]     #[]
        thresh96_16=[1.1077,0.3025,0.0891]     #[]
        thresh97_16=[0.1247,0.0198,0.0]        #[]
        thresh98_16=[0.0,0.0,0.0]              #[17.6086,11.7709,10.3354]
        thresh99_16=[0.0,0.0,0.0]              #[1.5485,0.3913,0.1449]

        thresh95_25=[0.9314,0.2914,0.1391]     #[]
        thresh96_25=[0.1235,0.0310,0.0]        #[]
        thresh97_25=[0.0078,0.0,0.0]           #[]
        thresh98_25=[0.0,0.0,0.0]              #[8.0665,4.5942,3.6101]
        thresh99_25=[0.0,0.0,0.0]              #[0.2646,0.0810,0.0]

        thresh95_36=[0.0844,0.0,0.0]           #[]
        thresh96_36=[0.0,0.0,0.0]              #[]
        thresh98_36=[0.0,0.0,0.0]              #[3.0406,1.9955,1.4338]
        thresh99_36=[0.0,0.0,0.0]              #[0.0381,0,0.0]

        thresh96_49=[0.0,0.0,0.0]              #[]
        thresh98_49=[0.0,0.0,0.0]              #[1.0355,0.7475,0.2113]
        thresh99_49=[0.0,0.0,0.0]              #[0.0216,0,0.0]

        thresh96_64=[0.0,0.0,0.0]              #[]
        thresh98_64=[0.0,0.0,0.0]              #[0.3048,0.1727,0.0]
        thresh99_64=[0.0,0.0,0.0]              #[0,0,0.0]

        thresh96_81=[0.0,0.0,0.0]              #[]
        thresh98_81=[0.0,0.0,0.0]              #[0.1144,0.0,0.0]
        thresh99_81=[0.0,0.0,0.0]              #[0,0,0.0]

        thresh95_total=[57.4187,44.5704,33.9538]        #[]
        thresh96_total=[33.8508,19.4982,11.3848]        #[]
        thresh97_total=[16.2138,4.7062,1.7977]          #[]
        thresh98_total=[5.1507,0.5691,0.0779]           #[]
        thresh99_total=[0.6504,0.0198,0.0]              #[]

        plt.plot(gridSize,thresh95_4,'-',color='purple',label='95% cover threshold')
        plt.plot(gridSize,thresh96_4,'-',color='green',label='96% cover threshold')
        plt.plot(gridSize,thresh97_4,'-',color='blue',label='97% cover threshold')
        plt.plot(gridSize,thresh98_4,'-',color='orange',label='98% cover threshold')
        plt.plot(gridSize,thresh99_4,'-',color='red',label='99% cover threshold')
        plt.plot(gridSize,thresh95_total,'--',color='purple')#,label='Cumulative 95% coverage')
        plt.plot(gridSize,thresh96_total,'--',color='green')#,label='Cumulative 96% coverage')
        plt.plot(gridSize,thresh97_total,'--',color='blue')#,label='Cumulative 97% coverage')
        plt.plot(gridSize,thresh98_total,'--',color='orange')#,label='Cumulative 98% coverage')
        plt.plot(gridSize,thresh99_total,'--',color='red')#,label='Cumulative 99% coverage')
        plt.title('Raster Coverage by 2x2 Squares')
        plt.xlabel('Grid Size (m)')
        plt.xticks([10,20,30])
        plt.ylabel('Area of Raster Above Threshold (%)')
        plt.legend()
        plt.show()

#######################################

# The main block
if __name__ == '__main__':
    args=readCommands()

    loadDir='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/grids/nouragues/'
    inFilename='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/grids/nouragues/grid30/metric/grid30MetricSummary.txt'
    outFilename='../data/tiffs/Nouragues{}.tif'.format(args.res)
    outBinaryName='../data/tiffs/binary_regions/NouraguesBinary{0}-{1:.0f}.tif'.format(args.res,args.threshold*100)
    outRegionName='../data/tiffs/binary_regions/NouraguesRegions{0}-{1:.0f}.tif'.format(args.res,args.threshold*100)

    data=coverStats()
    if args.stats:
        data.loadData(loadDir)
    #data.plotCover()
    data.writeGeoTiff(inFilename,args.res,args.epsg,outFilename,outBinaryName,outRegionName)
    #data.summaryStatsPlot()
