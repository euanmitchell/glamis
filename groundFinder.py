############################################
# Script to load and combine results of
# denoising data for a range of parameters
# and determine the best ground estimate for
# each waveform in a grid of footprints.
############################################

import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
from osgeo import gdal
from osgeo import osr
from math import sqrt,tan,pi
from scipy.stats import norm
from statistics import mean
from rasterTools import writeToFile, defineRegions, readFromFile
import argparse
import sys
import csv
import warnings

#########################################

# Defining the command line reading function
def readCommands():
    '''
    Read the arguments passed from the command line
    '''
    p=argparse.ArgumentParser(description=('Specify parameters for program control'),
        formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('--suppress', dest='suppress', action='store_true', default=False,
        help=('Suppress "RuntimeWarnings" in code.'))
    p.add_argument('--site', dest='site', type=str, default='',
        help=('''Shortcut to set --rootDir, --epsg, and --plotRoot for a whole field site. Valid values: laselva, robson, nouragues, paracou, hubbard, windriver, oakridge. Default is blank.'''))
    p.add_argument('--name', dest='name', type=str, default='',
        help=('''The name of the site, set by --site. Default is blank.'''))
    p.add_argument('--rootDir', dest='rootDir', type=str,
        default='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/laselva/grid30/',
        help=('''The root directory to look in for sub-directories with individual gediMetric output text files. Sub directories should follow the naming convention: "metric{covThresh}-n" and text files should follow the naming convention: "metric{gridSize}-{covThresh}-nSummary.txt" where n is an integer identifier.'''))
    p.add_argument('--epsg', dest='epsg', type=int, default=32616,
        help=('EPSG code of the data. Default is 32616'))
    p.add_argument('--gridSize', dest='gridSize', type=int, default=30,
        help=('The spacing between footprints. Default is 30'))
    p.add_argument('--covThresh', dest='covThresh', type=int, default=98,
        help=('The canopy cover threshold for input datasets (i.e., the -linkNoise beam sensitivity). Default is 98'))
    p.add_argument('--varScale', dest='varScale', action='store_true', default=False,
        help=('Use the extra data with varScale = 1 and varScale = 1.5.'))
    p.add_argument('--outlierThresh', dest='outlierThresh', type=int, default=10,
        help=('The threshold distance from the mean for calculating the outlier proportion beyond. Default is 10'))
    p.add_argument('--fThresh', dest='fThresh', type=int, default=3,
        help=('The threshold number of low std dev cells within "searchSize" window before correcting target high std dev cell. Default is 3'))
    p.add_argument('--prefAlgm', dest='prefAlgm', type=int, default=1,
        help=('The preferred algorithm to use to fill the newGround array at low std dev footprints. Default is 1'))
    p.add_argument('--plotRoot', dest='plotRoot', type=str, default='../data/laSelvaGroundError/30-98/',
        help=('The root of the directory to save plots and tiffs to. Default is "../data/laSelvaGroundError/30-98/"'))
    p.add_argument('--plots', dest='plots', action='store_true', default=False,
        help=('Run plotGroundRange() to make plots for different algorithm settings'))
    p.add_argument('--sdThresh', dest='sdThresh', type=float, default=3.0,
        help=('The standard deviation threshold for generating binary and region rasters. Default is 3.0'))
    p.add_argument('--sdScatter', dest='sdScatter', action='store_true', default=False,
        help=('Generate standard deviation vs ground error scatter plots for different algorithm settings'))
    p.add_argument('--sdHist', dest='sdHist', action='store_true', default=False,
        help=('Generate standard deviation histogram plots for the nine algorithms'))
    p.add_argument('--groundHist', dest='groundHist', action='store_true', default=False,
        help=('Generate ground elevation histogram plots for different algorithm settings'))
    p.add_argument('--coverHist', dest='coverHist', action='store_true', default=False,
        help=('Generate ALS cover histogram plots for the ground truth data'))
    p.add_argument('--heightHist', dest='heightHist', action='store_true', default=False,
        help=('Generate ALS tree height (RH95) histogram plots for the ground truth data'))
    p.add_argument('--rh50Hist', dest='rh50Hist', action='store_true', default=False,
        help=('Generate ALS RH50 height histogram plots for the ground truth data'))
    p.add_argument('--errorHist', dest='errorHist', action='store_true', default=False,
        help=('Generate ground error histogram plots for different algorithm settings'))
    p.add_argument('--fhdHist', dest='fhdHist', action='store_true', default=False,
        help=('Generate FHD histogram plots for different algorithm settings'))
    p.add_argument('--sdBox', dest='sdBox', action='store_true', default=False,
        help=('Generate standard deviation vs ground error box plots for different algorithm settings.'))
    p.add_argument('--sdRaster', dest='sdRaster', action='store_true', default=False,
        help=('Call stdDevRasters() to generate rasters from standard deviation array.'))
    p.add_argument('--update', dest='update', action='store_true', default=False,
        help=('Call updateGround() to correct high standard deviation cells from standard deviation array.'))
    p.add_argument('--updateMethod', dest='updateMethod', type=str, default='mean',
        help=('''Choose the method for estimating ground elevation at high std dev footprints. Valid values: hybrid, ideal, idealAll, mean, nearest, natural, linear, quadratic. Default is mean'''))
    p.add_argument('--forceToWave', dest='forceToWave', action='store_true', default=False,
        help=('Force the chosen update method to return an identified peak in the waveform when running updateGround().'))
    p.add_argument('--noUpdatePlots', dest='noUpdatePlots', action='store_true', default=False,
        help=('Skip saving plots and writing tiffs when running updateGround().'))
    p.add_argument('--writeCoords', dest='writeCoords', action='store_true', default=False,
        help=('Write the coordinates of large newGroundError footprints to file when running updateGround().'))
    cmdargs=p.parse_args()
    return cmdargs

#########################################

class groundFinder(object):
    '''A class to estimate ground elevations from a regular grid of footprints'''

    def __init__(self):
        '''Initialise the class'''

    def loadData(self,rootDir,step,threshold):
        '''Load the different algorithm datasets into arrays and compute various properties'''

        # Get the files with the individual algorithm setting outputs
        algorithm1=rootDir+'metric'+str(args.covThresh)+'-1/metric'+str(args.gridSize)+'-'+str(args.covThresh)+'-1Summary.txt'
        algorithm2=rootDir+'metric'+str(args.covThresh)+'-2/metric'+str(args.gridSize)+'-'+str(args.covThresh)+'-2Summary.txt'
        algorithm3=rootDir+'metric'+str(args.covThresh)+'-3/metric'+str(args.gridSize)+'-'+str(args.covThresh)+'-3Summary.txt'
        algorithm4=rootDir+'metric'+str(args.covThresh)+'-4/metric'+str(args.gridSize)+'-'+str(args.covThresh)+'-4Summary.txt'
        algorithm5=rootDir+'metric'+str(args.covThresh)+'-5/metric'+str(args.gridSize)+'-'+str(args.covThresh)+'-5Summary.txt'
        algorithm6=rootDir+'metric'+str(args.covThresh)+'-6/metric'+str(args.gridSize)+'-'+str(args.covThresh)+'-6Summary.txt'
        algorithm7=rootDir+'metric'+str(args.covThresh)+'-7/metric'+str(args.gridSize)+'-'+str(args.covThresh)+'-7Summary.txt'
        algorithm8=rootDir+'metric'+str(args.covThresh)+'-8/metric'+str(args.gridSize)+'-'+str(args.covThresh)+'-8Summary.txt'
        algorithm9=rootDir+'metric'+str(args.covThresh)+'-9/metric'+str(args.gridSize)+'-'+str(args.covThresh)+'-9Summary.txt'
        algorithm10=rootDir+'metric'+str(args.covThresh)+'-10/metric'+str(args.gridSize)+'-'+str(args.covThresh)+'-10Summary.txt'
        algorithm150=rootDir+'metric'+str(args.covThresh)+'-150/metric'+str(args.gridSize)+'-'+str(args.covThresh)+'-150Summary.txt'

        # Read each input file into its own array                               # Columns:
        print('Reading raw data into arrays ...')                               # 0 = x
        data1=np.loadtxt(algorithm1, comments='#', skiprows=1)                  # 1 = y
        data2=np.loadtxt(algorithm2, comments='#', skiprows=1)                  # 2 = trueGround
        data3=np.loadtxt(algorithm3, comments='#', skiprows=1)                  # 3 = beam sensitivity
        data4=np.loadtxt(algorithm4, comments='#', skiprows=1)                  # 4 = ALScover
        data5=np.loadtxt(algorithm5, comments='#', skiprows=1)                  # 5 = rhReal95
        data6=np.loadtxt(algorithm6, comments='#', skiprows=1)                  # 6 = groundHeight
        data7=np.loadtxt(algorithm7, comments='#', skiprows=1)                  # 7 = cover
        data8=np.loadtxt(algorithm8, comments='#', skiprows=1)                  # 8 = rh95
        data9=np.loadtxt(algorithm9, comments='#', skiprows=1)                  # 9 = FHD
        data10=np.loadtxt(algorithm10, comments='#', skiprows=1)                #10 = rhReal50
        data150=np.loadtxt(algorithm150, comments='#', skiprows=1)              #11 = gSlope; 12 = waveEnergy; 13 = meanNoise; 14 = noiseStdDev

        if args.varScale:
            dataArrays=[data1,data2,data3,data4,data5,data6,data7,data8,data9,data150,data10]
        else:
            dataArrays=[data1,data2,data3,data4,data5,data6,data7,data8,data9]
        print('Raw data1 array shape:',data1.shape)
        self.numRecords=data1.shape[0]
        for index, array in enumerate(dataArrays):
            if array.shape != dataArrays[0].shape:
                print('WARNING: Mismatched array shapes!')
                print('data{} array shape:'.format(index+1),array.shape)

        # Make re-shaped 3D arrays for each of the six inputs
        self.minX=(np.min(data1[:,0]))
        self.maxX=(np.max(data1[:,0]))
        binsX=((self.maxX-self.minX)/args.gridSize)+1       # +1 so that each footprint is centred in a raster pixel later
        self.minY=(np.min(data1[:,1]))
        self.maxY=(np.max(data1[:,1]))
        binsY=((self.maxY-self.minY)/args.gridSize)+1
        dataCube1=np.full((int(binsY),int(binsX),data1.shape[1]),np.nan)
        dataCube2=np.full((int(binsY),int(binsX),data2.shape[1]),np.nan)
        dataCube3=np.full((int(binsY),int(binsX),data3.shape[1]),np.nan)
        dataCube4=np.full((int(binsY),int(binsX),data4.shape[1]),np.nan)
        dataCube5=np.full((int(binsY),int(binsX),data5.shape[1]),np.nan)
        dataCube6=np.full((int(binsY),int(binsX),data6.shape[1]),np.nan)
        dataCube7=np.full((int(binsY),int(binsX),data7.shape[1]),np.nan)
        dataCube8=np.full((int(binsY),int(binsX),data8.shape[1]),np.nan)
        dataCube9=np.full((int(binsY),int(binsX),data9.shape[1]),np.nan)
        dataCube10=np.full((int(binsY),int(binsX),data10.shape[1]),np.nan)
        dataCube150=np.full((int(binsY),int(binsX),data150.shape[1]),np.nan)
        if args.varScale:
            dataCubes=[dataCube1,dataCube2,dataCube3,dataCube4,dataCube5,dataCube6,dataCube7,dataCube8,dataCube9,dataCube150,dataCube10]
        else:
            dataCubes=[dataCube1,dataCube2,dataCube3,dataCube4,dataCube5,dataCube6,dataCube7,dataCube8,dataCube9]

        # Fill the six data cubes with the re-shaped original data
        def fillDataCube(input,output):
            for i in range(0,input.shape[0]):
                row=input[i,:]
                xDist=row[0]-(self.minX-(args.gridSize/2))      # Start counting from half a raster cell beyond data edge
                yDist=(self.maxY+(args.gridSize/2))-row[1]
                xBin=(int(xDist//step))
                yBin=(int(yDist//step))
                output[yBin,xBin,:]=row
            return output

        print('Reshaping raw data arrays into data cubes ...')
        for index, array in enumerate(dataArrays):
            dataCubes[index]=fillDataCube(array,dataCubes[index])
        print('Data cube shape:',dataCube1.shape)

        # Build a 3D array of just ground elevation for the nine options
        if args.varScale:
            self.ground=np.stack((
                dataCube1[:,:,6],dataCube2[:,:,6],dataCube3[:,:,6],
                dataCube4[:,:,6],dataCube5[:,:,6],dataCube6[:,:,6],
                dataCube7[:,:,6],dataCube8[:,:,6],dataCube9[:,:,6],
                dataCube150[:,:,6],dataCube10[:,:,6]),
                axis=2)
            self.canopyHeight=np.stack((
                dataCube1[:,:,8],dataCube2[:,:,8],dataCube3[:,:,8],
                dataCube4[:,:,8],dataCube5[:,:,8],dataCube6[:,:,8],
                dataCube7[:,:,8],dataCube8[:,:,8],dataCube9[:,:,8],
                dataCube150[:,:,8],dataCube10[:,:,8]),
                axis=2)
            self.fhd=np.stack((
                dataCube1[:,:,9],dataCube2[:,:,9],dataCube3[:,:,9],
                dataCube4[:,:,9],dataCube5[:,:,9],dataCube6[:,:,9],
                dataCube7[:,:,9],dataCube8[:,:,9],dataCube9[:,:,9],
                dataCube150[:,:,9],dataCube10[:,:,9]),
                axis=2)
            self.bsMinusCover=np.stack((
                dataCube1[:,:,3]-dataCube1[:,:,4],dataCube2[:,:,3]-dataCube2[:,:,4],dataCube3[:,:,3]-dataCube3[:,:,4],
                dataCube4[:,:,3]-dataCube4[:,:,4],dataCube5[:,:,3]-dataCube5[:,:,4],dataCube6[:,:,3]-dataCube6[:,:,4],
                dataCube7[:,:,3]-dataCube7[:,:,4],dataCube8[:,:,3]-dataCube8[:,:,4],dataCube9[:,:,3]-dataCube9[:,:,4],
                dataCube150[:,:,3]-dataCube150[:,:,4],dataCube10[:,:,3]-dataCube10[:,:,4]),
                axis=2)
            self.cover=np.stack((
                dataCube1[:,:,7],dataCube2[:,:,7],dataCube3[:,:,7],
                dataCube4[:,:,7],dataCube5[:,:,7],dataCube6[:,:,7],
                dataCube7[:,:,7],dataCube8[:,:,7],dataCube9[:,:,7],
                dataCube150[:,:,7],dataCube10[:,:,7]),
                axis=2)
            self.gSlope=np.stack((
                dataCube1[:,:,11],dataCube2[:,:,11],dataCube3[:,:,11],
                dataCube4[:,:,11],dataCube5[:,:,11],dataCube6[:,:,11],
                dataCube7[:,:,11],dataCube8[:,:,11],dataCube9[:,:,11],
                dataCube150[:,:,11],dataCube10[:,:,11]),
                axis=2)
            self.waveEnergy=np.stack((
                dataCube1[:,:,12],dataCube2[:,:,12],dataCube3[:,:,12],
                dataCube4[:,:,12],dataCube5[:,:,12],dataCube6[:,:,12],
                dataCube7[:,:,12],dataCube8[:,:,12],dataCube9[:,:,12],
                dataCube150[:,:,12],dataCube10[:,:,12]),
                axis=2)
            self.meanNoise=np.stack((
                dataCube1[:,:,13],dataCube2[:,:,13],dataCube3[:,:,13],
                dataCube4[:,:,13],dataCube5[:,:,13],dataCube6[:,:,13],
                dataCube7[:,:,13],dataCube8[:,:,13],dataCube9[:,:,13],
                dataCube150[:,:,13],dataCube10[:,:,13]),
                axis=2)
            self.noiseStdDev=np.stack((
                dataCube1[:,:,14],dataCube2[:,:,14],dataCube3[:,:,14],
                dataCube4[:,:,14],dataCube5[:,:,14],dataCube6[:,:,14],
                dataCube7[:,:,14],dataCube8[:,:,14],dataCube9[:,:,14],
                dataCube150[:,:,14],dataCube10[:,:,14]),
                axis=2)
        else:
            self.ground=np.stack((
                dataCube1[:,:,6],dataCube2[:,:,6],dataCube3[:,:,6],
                dataCube4[:,:,6],dataCube5[:,:,6],dataCube6[:,:,6],
                dataCube7[:,:,6],dataCube8[:,:,6],dataCube9[:,:,6]),
                axis=2)
            self.canopyHeight=np.stack((
                dataCube1[:,:,8],dataCube2[:,:,8],dataCube3[:,:,8],
                dataCube4[:,:,8],dataCube5[:,:,8],dataCube6[:,:,8],
                dataCube7[:,:,8],dataCube8[:,:,8],dataCube9[:,:,8]),
                axis=2)
            self.fhd=np.stack((
                dataCube1[:,:,9],dataCube2[:,:,9],dataCube3[:,:,9],
                dataCube4[:,:,9],dataCube5[:,:,9],dataCube6[:,:,9],
                dataCube7[:,:,9],dataCube8[:,:,9],dataCube9[:,:,9]),
                axis=2)
            self.bsMinusCover=np.stack((
                dataCube1[:,:,3]-dataCube1[:,:,4],dataCube2[:,:,3]-dataCube2[:,:,4],dataCube3[:,:,3]-dataCube3[:,:,4],
                dataCube4[:,:,3]-dataCube4[:,:,4],dataCube5[:,:,3]-dataCube5[:,:,4],dataCube6[:,:,3]-dataCube6[:,:,4],
                dataCube7[:,:,3]-dataCube7[:,:,4],dataCube8[:,:,3]-dataCube8[:,:,4],dataCube9[:,:,3]-dataCube9[:,:,4]),
                axis=2)
            self.cover=np.stack((
                dataCube1[:,:,7],dataCube2[:,:,7],dataCube3[:,:,7],
                dataCube4[:,:,7],dataCube5[:,:,7],dataCube6[:,:,7],
                dataCube7[:,:,7],dataCube8[:,:,7],dataCube9[:,:,7]),
                axis=2)
            self.gSlope=np.stack((
                dataCube1[:,:,11],dataCube2[:,:,11],dataCube3[:,:,11],
                dataCube4[:,:,11],dataCube5[:,:,11],dataCube6[:,:,11],
                dataCube7[:,:,11],dataCube8[:,:,11],dataCube9[:,:,11]),
                axis=2)
            self.waveEnergy=np.stack((
                dataCube1[:,:,12],dataCube2[:,:,12],dataCube3[:,:,12],
                dataCube4[:,:,12],dataCube5[:,:,12],dataCube6[:,:,12],
                dataCube7[:,:,12],dataCube8[:,:,12],dataCube9[:,:,12]),
                axis=2)
            self.meanNoise=np.stack((
                dataCube1[:,:,13],dataCube2[:,:,13],dataCube3[:,:,13],
                dataCube4[:,:,13],dataCube5[:,:,13],dataCube6[:,:,13],
                dataCube7[:,:,13],dataCube8[:,:,13],dataCube9[:,:,13]),
                axis=2)
            self.noiseStdDev=np.stack((
                dataCube1[:,:,14],dataCube2[:,:,14],dataCube3[:,:,14],
                dataCube4[:,:,14],dataCube5[:,:,14],dataCube6[:,:,14],
                dataCube7[:,:,14],dataCube8[:,:,14],dataCube9[:,:,14]),
                axis=2)

        # Build 2D arrays of lat and long
        self.eastings=dataCube1[:,:,0]
        self.northings=dataCube1[:,:,1]

        # Build 2D arrays of true ground, cover, RH95, and RH50 from ALS - should be the same for all algorithms
        self.groundTruth=dataCube1[:,:,2]
        self.alsCover=dataCube1[:,:,4]
        self.alsCanopyHeight=dataCube1[:,:,5]
        self.rh50=dataCube1[:,:,10]

        # Build 3D array of ground errors - i.e., residuals
        self.groundError=np.full((self.ground.shape),0.0)
        for i in range(self.ground.shape[2]):
            self.groundError[:,:,i]=self.ground[:,:,i]-self.groundTruth

        # Get arrays with RMSE, Bias, and outlier fraction for each of the algorithms
        means=np.nanmean(self.ground[:,:,0:9],axis=2)
        self.rmse=np.full((self.ground.shape[2]),0.0)
        self.bias=np.full((self.ground.shape[2]),0.0)
        self.outFrac=np.full((self.ground.shape[2]),0.0)
        for i in range(self.ground.shape[2]):
            self.rmse[i]=sqrt(np.nansum(self.groundError[:,:,i]**2)/np.sum(~np.isnan(self.groundError[:,:,i])))
            self.bias[i]=(np.nansum(self.groundError[:,:,i])/np.sum(~np.isnan(self.groundError[:,:,i])))
            self.outFrac[i]=(np.nansum(np.abs((self.ground[:,:,i]-means))>args.outlierThresh)/np.sum(~np.isnan(self.groundTruth)))*100
        with np.printoptions(precision=2,floatmode='fixed'):
            print('RMSE array:',self.rmse)
            print('Bias array:',self.bias)
            print('Outlier ({} m) proportion array:'.format(args.outlierThresh),self.outFrac)

        # Build 2D array of standard deviations of ground estimates at each location - ignore the last two datasets (varScale = 1 & 1.5)
        self.groundStdDev=np.nanstd(self.ground[:,:,0:9],axis=2)

        # Generate a binary standard deviation raster based on threshold (0 < threshold; 1 > threshold)
        self.binaryStdDev=np.full((self.groundStdDev.shape),-999,dtype=int)
        self.binaryStdDev[self.groundStdDev<threshold]=0
        self.binaryStdDev[self.groundStdDev>=threshold]=1

    def plotGroundRange(self,threshold):
        '''Make plots of ground error vs standard deviation for each algorithm'''

        # Figure formatting for publications
        plt.rcParams['figure.figsize']=(6,4.8)
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
        plt.rcParams['boxplot.boxprops.linewidth']=0.5
        plt.rcParams['boxplot.capprops.linewidth']=0.5
        plt.rcParams['boxplot.medianprops.linewidth']=0.5
        plt.rcParams['boxplot.whiskerprops.linewidth']=0.5
        plt.rcParams['boxplot.whiskers']=0.5
        plt.rcParams['axes.labelpad']=3
        plt.rcParams['axes.titlesize']=8
        #plt.rcParams['savefig.format']='tif'

        # Low standard deviation mask
        lowStdDev=np.where((self.groundStdDev<args.sdThresh))

        # Write csv file for CloudCompare
        '''allZ=self.groundTruth.flatten()
        allX=self.eastings.flatten()
        allY=self.northings.flatten()
        alsOutput=np.stack((allX,allY,allZ),axis=-1)
        header=['X','Y','Z']
        print(alsOutput.shape)
        f=open('../data/NouraguesALS.txt','w')
        writer=csv.writer(f)
        writer.writerow(header)
        writer.writerows(alsOutput)
        f.close()
        subZ=self.ground[:,:,0][lowStdDev].flatten()
        subX=self.eastings[lowStdDev].flatten()
        subY=self.northings[lowStdDev].flatten()
        waveOutput=np.stack((subX,subY,subZ),axis=-1)
        print(waveOutput.shape)
        f=open('../data/NouraguesWave.txt','w')
        writer=csv.writer(f)
        writer.writerow(header)
        writer.writerows(waveOutput)
        f.close()'''

        # Generate a 3x3 figure with histograms of ground elevation for each algorithm
        if args.groundHist:
            fig,axes=plt.subplots(3,3,figsize=(15,12),sharex=True,sharey=True)
            subplots=[axes[0,0],axes[0,1],axes[0,2],axes[1,0],axes[1,1],axes[1,2],axes[2,0],axes[2,1],axes[2,2]]
            for i in range(self.ground.shape[2]):
                subplots[i].hist(self.ground[:,:,i][~np.isnan(self.ground[:,:,i])].flatten(),bins=50,label='Algorithm {0}'.format(i+1))
                subplots[i].axvline(x=np.nanmin(self.ground[:,:,i]),color='g',label='Min: {0:.1f}'.format(np.nanmin(self.ground[:,:,i])),linestyle='--')
                subplots[i].axvline(x=np.nanmean(self.ground[:,:,i]),color='r',label='Mean: {0:.1f}'.format(np.nanmean(self.ground[:,:,i])),linestyle='--')
                subplots[i].axvline(x=np.nanmean(self.groundTruth),color='grey',label='ALS Mean: {0:.1f}'.format(np.nanmean(self.groundTruth)),linestyle='--')
                subplots[i].axvline(x=np.nanmax(self.ground[:,:,i]),color='b',label='Max: {0:.1f}'.format(np.nanmax(self.ground[:,:,i])),linestyle='--')
                handles,labels=subplots[i].get_legend_handles_labels()
                handles=[handles[4],handles[0],handles[1],handles[2],handles[3]]
                labels=[labels[4],labels[0],labels[1],labels[2],labels[3]]
                subplots[i].legend(handles,labels,framealpha=1)
            for axis in axes.flat:
                axis.set(xlabel='Ground Elevation (m)',ylabel='Frequency')
            for axis in axes.flat:
                axis.label_outer()
            plt.savefig(args.plotRoot+'GroundHist_All_Scaled.png',dpi=300)
            plt.close()
            plt.clf()
            print('Ground elevation histogram plots written to '+args.plotRoot+'GroundHist_All_Scaled.png')
            #plt.show()

        # Generate a ground standard deviation histogram
        if args.sdHist:
            plt.hist(self.groundStdDev[~np.isnan(self.groundStdDev)].flatten(),bins=60)
            plt.title('{} Ground Standard Deviation'.format(args.name))
            plt.xlabel('Ground Standard Deviation (m)')
            plt.ylabel('Frequency')
            plt.savefig(args.plotRoot+'StdDevHist.png',dpi=300)
            plt.close()
            plt.clf()
            print('Standard deviation histogram plots written to '+args.plotRoot+'StdDevHist.png')
            #plt.show()

        # Generate a 3x3 figure with histograms of ground elevation error for each algorithm
        if args.errorHist:
            fig,axes=plt.subplots(3,3,figsize=(15,12))
            subplots=[axes[0,0],axes[0,1],axes[0,2],axes[1,0],axes[1,1],axes[1,2],axes[2,0],axes[2,1],axes[2,2]]
            for i in range(self.ground.shape[2]):
                _, bins, _=subplots[i].hist(self.groundError[:,:,i][~np.isnan(self.groundError[:,:,i])].flatten(),bins=50,label='Alg. {0}'.format(i+1))
                subplots[i].hist(self.groundError[:,:,i][self.groundStdDev>=threshold].flatten(),bins=bins,label='High SD')
                #subplots[i].axvline(x=np.nanmin(self.ground[:,:,i]),color='g',label='Min: {0:.1f}'.format(np.nanmin(self.ground[:,:,i])),linestyle='--')
                subplots[i].axvline(x=0,color='r',linestyle='--')
                #subplots[i].axvline(x=np.nanmean(self.groundTruth),color='grey',label='ALS Mean: {0:.1f}'.format(np.nanmean(self.groundTruth)),linestyle='--')
                #subplots[i].axvline(x=np.nanmax(self.ground[:,:,i]),color='b',label='Max: {0:.1f}'.format(np.nanmax(self.ground[:,:,i])),linestyle='--')
                #handles,labels=subplots[i].get_legend_handles_labels()
                #handles=[handles[4],handles[0],handles[1],handles[2],handles[3]]
                #labels=[labels[4],labels[0],labels[1],labels[2],labels[3]]
                subplots[i].legend()
            for axis in axes.flat:
                axis.set(xlabel='Ground Error (m)',ylabel='Frequency')
            #for axis in axes.flat:
                #axis.label_outer()
            plt.savefig(args.plotRoot+'ErrorHistSD{0:.0f}_All_Unscaled.png'.format(args.sdThresh),dpi=300)
            plt.close()
            plt.clf()
            print('Ground error histogram plots written to '+args.plotRoot+'ErrorHistSD{0:.0f}_All_Unscaled.png'.format(args.sdThresh))
            #plt.show()

        # Generate a 3x3 figure of standard deviation vs ground error scatter plots for each algorithm
        if args.sdScatter:
            fig,axes=plt.subplots(3,3,figsize=(15,12),sharey=True)      # Include "sharey=True" for same y-axis on all subplots
            subplots=[axes[0,0],axes[0,1],axes[0,2],axes[1,0],axes[1,1],axes[1,2],axes[2,0],axes[2,1],axes[2,2]]
            for i in range(9):
                subplots[i].scatter(self.groundStdDev,self.groundError[:,:,i],s=1,label='Algorithm {0}\nBias: {1:.2f}\nRMSE: {2:.2f}\nOutlier%: {3:.2f}'.format(i+1,self.bias[i],self.rmse[i],self.outFrac[i]))
                subplots[i].axhline(y=0.0,color='r',linestyle='--')
                subplots[i].legend(loc='lower left')
                #subplots[i].set_title('Algorithm {}'.format(i+1))
            for axis in axes.flat:
                axis.set(xlabel='Ground Standard Deviation (m)',ylabel='Ground Error (m)')
            for axis in axes.flat:
                axis.label_outer()        # If sharey=True
            plt.savefig(args.plotRoot+'Scatter_All_Scaled.png',dpi=300)
            plt.close()
            plt.clf()
            print('Standard deviation scatter plots written to '+args.plotRoot+'Scatter_All_Scaled.png')
            #plt.show()

        # Generate a 3x3 figure of standard deviation vs ground error box plots for each algorithm
        if args.sdBox:
            fig,axes=plt.subplots(3,3,figsize=(6,4.8),sharey=True)      # Include "sharey=True" for same y-axis on all subplots
            subplots=[axes[0,0],axes[0,1],axes[0,2],axes[1,0],axes[1,1],axes[1,2],axes[2,0],axes[2,1],axes[2,2]]
            for i in range(9):
                data=[]
                for k in range(15):
                    data.append(self.groundError[:,:,i][np.where((self.groundStdDev < k+1) & (self.groundStdDev > k) & (~np.isnan(self.groundError[:,:,i])))])
                    #x1=self.groundError[:,:,i][np.where((self.groundStdDev < args.sdThresh) & (~np.isnan(self.groundError[:,:,i])))]
                    #x2=self.groundError[:,:,i][np.where((self.groundStdDev >= args.sdThresh) & (self.groundStdDev < 6.0) & (~np.isnan(self.groundError[:,:,i])))]
                    #x3=self.groundError[:,:,i][np.where((self.groundStdDev >= 6.0) & (self.groundStdDev < 9.0) & (~np.isnan(self.groundError[:,:,i])))]
                    #x4=self.groundError[:,:,i][np.where((self.groundStdDev >= 9.0) & (self.groundStdDev < 12.0) & (~np.isnan(self.groundError[:,:,i])))]
                    #x5=self.groundError[:,:,i][np.where((self.groundStdDev >= 12.0) & (self.groundStdDev < 15.0) & (~np.isnan(self.groundError[:,:,i])))]
                    x6=self.groundError[:,:,i][np.where((self.groundStdDev >= 15.0) & (~np.isnan(self.groundError[:,:,i])))]
                    #data=[x1, x2, x3, x4, x5]

                subplots[i].boxplot(data,sym='',)
                subplots[i].axhline(y=0.0,color='r',linestyle='--',linewidth=0.5)
            for axis in axes.flat:
                axis.set(xlabel='Ground Standard Deviation (m)',
                    ylabel='Ground Error (m)',
                    xticks=np.arange(1,16,step=2),
                    xticklabels=['1','3','5','7','9','11','13','15'])
            for axis in axes.flat:
                axis.label_outer()
            plt.tight_layout()
            plt.savefig(args.plotRoot+'Box_All_Scaled.png',dpi=300)
            plt.close()
            plt.clf()
            print('Standard deviation box plots written to '+args.plotRoot+'Box_All_Scaled.png')
            print('Number of footprints in the > 15 m standard deviation range is {0} ({1:.2f}%)'.format(len(x6),len(x6)/self.numRecords*100))
            #plt.show()

        # Generate a 3x3 figure with histograms of FHD for each algorithm
        if args.fhdHist:
            fig,axes=plt.subplots(3,3,figsize=(15,12),sharex=True,sharey=True)
            subplots=[axes[0,0],axes[0,1],axes[0,2],axes[1,0],axes[1,1],axes[1,2],axes[2,0],axes[2,1],axes[2,2]]
            for i in range(self.ground.shape[2]-2):
                subplots[i].hist(self.fhd[:,:,i][~np.isnan(self.fhd[:,:,i])].flatten(),bins=50,label='Alg. {0}'.format(i+1))
                subplots[i].legend(loc='upper left')
            for axis in axes.flat:
                axis.set(xlabel='Foliage Height Diversity',ylabel='Frequency')
            for axis in axes.flat:
                axis.label_outer()
            plt.savefig(args.plotRoot+'FHDHist_All_Scaled.png',dpi=300)
            plt.close()
            plt.clf()
            print('FHD histogram plots written to '+args.plotRoot+'FHDHist_All_Scaled.png')
            #plt.show()

        # Generate a canopy cover histogram from ALS data
        if args.coverHist:
            _, bins, _=plt.hist(self.alsCover[self.alsCover>=0.8].flatten(),bins=40,edgecolor='black',linewidth=0.5,
                label='All Data\n  Median: {0:.3f}'.format(np.nanmedian(self.alsCover[self.alsCover>=0.8])))
            plt.hist(self.alsCover[(self.alsCover>=0.8) & (self.groundStdDev<args.sdThresh)].flatten(),bins=bins,edgecolor='black',linewidth=0.5,
                label='Std. Dev. < {0:.0f} Data\n  Median: {1:.3f}'.format(args.sdThresh,np.nanmedian(self.alsCover[(self.alsCover>=0.8) & (self.groundStdDev<args.sdThresh)])))
            plt.hist(self.alsCover[(self.alsCover>=0.8) & (self.groundStdDev<1.0)].flatten(),bins=bins,edgecolor='black',linewidth=0.5,
                label='Std. Dev. < 1 Data\n  Median: {0:.3f}'.format(np.nanmedian(self.alsCover[(self.alsCover>=0.8) & (self.groundStdDev<1.0)])))
            plt.title('{} ALS Canopy Cover'.format(args.name))
            plt.xlabel('ALS Cover')
            plt.ylabel('Frequency')
            plt.legend()
            plt.savefig(args.plotRoot+'canopyCoverHighHistComp.png',dpi=300)
            plt.close()
            plt.clf()
            print('Standard deviation histogram plots written to '+args.plotRoot+'canopyCoverHighHistComp.png')
            #plt.show()

            _, bins, _=plt.hist(self.alsCover[self.alsCover>=0.0].flatten(),bins=100,edgecolor='black',linewidth=0.5,
                label='All Data\n  Median: {0:.3f}'.format(np.nanmedian(self.alsCover[self.alsCover>=0.0])))
            plt.hist(self.alsCover[(self.alsCover>=0.0) & (self.groundStdDev<args.sdThresh)].flatten(),bins=bins,edgecolor='black',linewidth=0.5,
                label='Std. Dev. < {0:.0f} Data\n  Median: {1:.3f}'.format(args.sdThresh,np.nanmedian(self.alsCover[(self.alsCover>=0.0) & (self.groundStdDev<args.sdThresh)])))
            plt.hist(self.alsCover[(self.alsCover>=0.0) & (self.groundStdDev<1.0)].flatten(),bins=bins,edgecolor='black',linewidth=0.5,
                label='Std. Dev. < 1 Data\n  Median: {0:.3f}'.format(np.nanmedian(self.alsCover[(self.alsCover>=0.0) & (self.groundStdDev<1.0)])))
            plt.title('{} ALS Canopy Cover'.format(args.name))
            plt.xlabel('ALS Cover')
            plt.ylabel('Frequency')
            plt.legend()
            plt.savefig(args.plotRoot+'canopyCoverHistComp.png',dpi=300)
            plt.close()
            plt.clf()
            print('Standard deviation histogram plots written to '+args.plotRoot+'canopyCoverHistComp.png')
            #plt.show()

        # Generate a canopy height histogram from ALS data
        if args.heightHist:
            _, bins, _=plt.hist(self.alsCanopyHeight[~np.isnan(self.alsCanopyHeight)].flatten(),bins=50,edgecolor='black',linewidth=0.5,
                label='All Data\n  Median: {0:.1f}'.format(np.nanmedian(self.alsCanopyHeight[~np.isnan(self.alsCanopyHeight)])))
            plt.hist(self.alsCanopyHeight[(~np.isnan(self.alsCanopyHeight)) & (self.groundStdDev<args.sdThresh)].flatten(),bins=bins,
                edgecolor='black',linewidth=0.5,
                label='Std. Dev. < {0:.0f} Data\n  Median: {1:.1f}'.format(args.sdThresh,np.nanmedian(self.alsCanopyHeight[(~np.isnan(self.alsCanopyHeight)) & (self.groundStdDev<args.sdThresh)])))
            plt.hist(self.alsCanopyHeight[(~np.isnan(self.alsCanopyHeight)) & (self.groundStdDev<1.0)].flatten(),bins=bins,
                edgecolor='black',linewidth=0.5,
                label='Std. Dev. < 1 Data\n  Median: {0:.1f}'.format(np.nanmedian(self.alsCanopyHeight[(~np.isnan(self.alsCanopyHeight)) & (self.groundStdDev<1.0)])))
            plt.title('{} ALS Canopy Height'.format(args.name))
            plt.xlabel('ALS RH95 (m)')
            plt.ylabel('Frequency')
            plt.legend()
            plt.savefig(args.plotRoot+'alsHeightHistComp.png',dpi=300)
            plt.close()
            plt.clf()
            print('Standard deviation histogram plots written to '+args.plotRoot+'alsHeightHistComp.png')
            #plt.show()

        # Generate an RH50 histogram from ALS data
        if args.rh50Hist:
            _, bins, _=plt.hist(self.rh50[~np.isnan(self.rh50)].flatten(),bins=50,edgecolor='black',linewidth=0.5,
                label='All Data\n  Median: {0:.1f}'.format(np.nanmedian(self.rh50[~np.isnan(self.rh50)])))
            plt.hist(self.rh50[(~np.isnan(self.rh50)) & (self.groundStdDev<args.sdThresh)].flatten(),bins=bins,
                edgecolor='black',linewidth=0.5,
                label='Std. Dev. < {0:.0f} Data\n  Median: {1:.1f}'.format(args.sdThresh,np.nanmedian(self.rh50[(~np.isnan(self.rh50)) & (self.groundStdDev<args.sdThresh)])))
            plt.hist(self.rh50[(~np.isnan(self.rh50)) & (self.groundStdDev<1.0)].flatten(),bins=bins,
                edgecolor='black',linewidth=0.5,
                label='Std. Dev. < 1 Data\n  Median: {0:.1f}'.format(np.nanmedian(self.rh50[(~np.isnan(self.rh50)) & (self.groundStdDev<1.0)])))
            plt.title('{} ALS RH50 Height'.format(args.name))
            plt.xlabel('ALS RH50 (m)')
            plt.ylabel('Frequency')
            plt.legend()
            plt.savefig(args.plotRoot+'alsRH50HistComp.png',dpi=300)
            plt.close()
            plt.clf()
            print('Standard deviation histogram plots written to '+args.plotRoot+'alsRH50HistComp.png')
            #plt.show()

        # Plot of beam sensitivity minus cover vs. ground std. dev. for all algorithms
        '''fig,axes=plt.subplots(3,3,figsize=(15,12),sharex=True,sharey=True)      # Include "sharey=True" for same y-axis on all subplots
        subplots=[axes[0,0],axes[0,1],axes[0,2],axes[1,0],axes[1,1],axes[1,2],axes[2,0],axes[2,1],axes[2,2]]
        for i in range(self.ground.shape[2]):
            subplots[i].scatter(self.bsMinusCover[:,:,i],self.groundStdDev,s=1,label='Algorithm {0}'.format(i+1))
            subplots[i].legend(loc='upper right')
        for axis in axes.flat:
            axis.set(xlabel='Beam Sensitivity - ALS Cover',ylabel='Ground Standard Deviation (m)')
        for axis in axes.flat:
            axis.label_outer()        # If sharey=True
        plt.savefig(args.plotRoot+'BSmC_All.png',dpi=300)
        plt.close()
        plt.clf()
        print('Beam sensitivity scatter plots written to '+args.plotRoot+'BSmC_All.png')
        #plt.show()'''

        # Beam sens minus cover vs proportion of low SD footprints
        '''xData=[]    # beamSens - ALS cover
        yData=[]    # frac. low SD footprints
        yStep=int(np.ceil(self.groundStdDev.shape[0]/10))
        xStep=int(np.ceil(self.groundStdDev.shape[1]/10))
        for i in np.arange(0,yStep*10,yStep):
            for j in np.arange(0,xStep*10,xStep):
                maxY=min(self.groundStdDev.shape[0],i+yStep)
                maxX=min(self.groundStdDev.shape[1],j+xStep)
                x=np.nanmean(self.bsMinusCover[i:maxY,j:maxX,args.prefAlgm-1])
                numFootprints=(self.binaryStdDev[i:maxY,j:maxX]!=-999).sum()
                y=(self.binaryStdDev[i:maxY,j:maxX]==0).sum()/numFootprints
                xData.append(x)
                yData.append(y)
        plt.scatter(xData,yData)
        plt.title('{0} {1}% Beam Sensitivity'.format(args.name,args.covThresh))
        plt.xlabel('Avg. Beam Sensitivity - ALS Cover')
        plt.ylabel('Fraction of Low Std. Dev. Footprints')
        plt.savefig(args.plotRoot+'BSmCvsFracLowSD{0:.0f}BS{1}.png'.format(args.sdThresh,args.covThresh),dpi=300)
        plt.close()
        plt.clf()
        print('Beam sensitivity vs. proportion of low SD footprints plot written to '+args.plotRoot+'BSmCvsFracLowSD{0:.0f}BS{1}.png'.format(args.sdThresh,args.covThresh))
        #plt.show()'''

        '''plt.scatter(self.groundStdDev,minError,marker='o',s=1,c='blue',label='Minimum Error')
        plt.scatter(self.groundStdDev,maxError,marker='o',s=1,c='red',label='Maximum Error')
        plt.scatter(maxError-minError,maxError,marker='o',s=1,label='Error Max-Min Diff.')
        plt.xlabel('Ground Error Difference (Max - Min) (m)')
        plt.ylabel('Maximum Ground Error (m)')
        plt.legend()
        plt.savefig(args.plotRoot+'maxErrorDiff.png',dpi=300)
        plt.close()
        plt.clf()'''

    def stdDevRasters(self):
        '''Generate rasters of raw standard deviation plus binary and region arrays'''

        # Write the raw standard deviation raster to file
        # Adjust raster starting position by half a gridSize to match fillDataCube() method above
        writeToFile(self.groundStdDev,self.minX-(args.gridSize/2),self.maxY+(args.gridSize/2),self.groundStdDev.shape[1],self.groundStdDev.shape[0],args.gridSize,args.epsg,args.plotRoot+'StdDev.tif',noData=-999.0)

        # Write the ALS cover raster to file
        writeToFile(self.alsCover,self.minX-(args.gridSize/2),self.maxY+(args.gridSize/2),self.groundStdDev.shape[1],self.groundStdDev.shape[0],args.gridSize,args.epsg,args.plotRoot+'ALSCover.tif',noData=-999.0)

        # Write the binary raster to file
        writeToFile(self.binaryStdDev,self.minX-(args.gridSize/2),self.maxY+(args.gridSize/2),self.groundStdDev.shape[1],self.groundStdDev.shape[0],args.gridSize,args.epsg,args.plotRoot+'binaryStdDevSD{0:.0f}BS{1:.0f}.tif'.format(args.sdThresh,args.covThresh),noData=-999.0)

        # Generate region array from this and write to file
        regions=defineRegions(self.groundStdDev,self.binaryStdDev,args.gridSize,args.epsg)
        writeToFile(regions,self.minX-(args.gridSize/2),self.maxY+(args.gridSize/2),self.groundStdDev.shape[1],self.groundStdDev.shape[0],args.gridSize,args.epsg,args.plotRoot+'HighStdRegions.tif',noData=-999.0)

        # Write the ALS true ground raster to file
        writeToFile(self.groundTruth,self.minX-(args.gridSize/2),self.maxY+(args.gridSize/2),self.groundStdDev.shape[1],self.groundStdDev.shape[0],args.gridSize,args.epsg,args.plotRoot+'ALSGround.tif',noData=-999.0)

    def updateGround(self,medianCanopyHeight):
        '''Calculate new ground elevation array'''

        # Plotting setttings for publications
        plt.rcParams['figure.figsize']=(6.3,2.3)
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
        plt.rcParams['axes.titlesize']=8
        plt.rcParams['lines.linewidth']=0.75
        #plt.rcParams['savefig.format']='tif'

        # Low and high standard deviation masks
        highStdDev=np.where((self.groundStdDev>=args.sdThresh))
        lowStdDev=np.where((self.groundStdDev<args.sdThresh))

        # Get the statistics of the high std dev footprints alone before correction
        means=np.nanmean(self.ground[:,:,0:9],axis=2)
        medians=np.nanmedian(self.ground[:,:,0:9],axis=2)
        self.rmseHSD=np.full((self.ground.shape[2]),0.0)
        self.biasHSD=np.full((self.ground.shape[2]),0.0)
        self.outFracHSD=np.full((self.ground.shape[2]),0.0)
        for i in range(self.ground.shape[2]):
            self.rmseHSD[i]=sqrt(np.nansum(self.groundError[:,:,i][highStdDev]**2)/np.sum(~np.isnan(self.groundError[:,:,i][highStdDev])))
            self.biasHSD[i]=(np.nansum(self.groundError[:,:,i][highStdDev])/np.sum(~np.isnan(self.groundError[:,:,i][highStdDev])))
            self.outFracHSD[i]=(np.nansum(np.abs((self.ground[:,:,i][highStdDev]-means[highStdDev]))>args.outlierThresh)/np.sum(~np.isnan(self.groundTruth[highStdDev])))*100
        with np.printoptions(precision=2,floatmode='fixed'):
            print('High std. dev. RMSE array:',self.rmseHSD)
            print('High std. dev. Bias array:',self.biasHSD)
            print('High std. dev. Outlier ({} m) proportion array:'.format(args.outlierThresh),self.outFracHSD)

        # Get the statistics of ALS RH95 and low std dev footprint waveform RH95
        '''useInd=np.where((self.groundStdDev<args.sdThresh)&(self.canopyHeight[:,:,args.prefAlgm-1]>2.0))
        useIndHi=np.where((self.groundStdDev>args.sdThresh)&(self.canopyHeight[:,:,args.prefAlgm-1]>2.0))
        meanCanopyHeight=np.nanmean(self.canopyHeight[:,:,args.prefAlgm-1][useInd])
        stdCanopyHeight=np.nanstd(self.canopyHeight[:,:,args.prefAlgm-1][useInd])
        print('Dataset ALS RH95 mean: {0:.2f}'.format(np.nanmean(self.alsCanopyHeight[self.canopyHeight[:,:,args.prefAlgm-1]>2.0])))
        print('Dataset ALS RH95 standard deviation: {0:.2f}'.format(np.nanstd(self.alsCanopyHeight[self.canopyHeight[:,:,args.prefAlgm-1]>2.0])))
        print('Dataset low std. dev. footprint waveform RH95 mean: {0:.2f}'.format(meanCanopyHeight))
        print('Dataset low std. dev. footprint waveform RH95 standard deviation: {0:.2f}'.format(stdCanopyHeight))
        print('Dataset high std. dev. footprint waveform RH95 mean: {0:.2f}'.format(np.nanmean(self.canopyHeight[:,:,args.prefAlgm-1][useIndHi])))
        print('Dataset high std. dev. footprint waveform RH95 standard deviation: {0:.2f}'.format(np.nanstd(self.canopyHeight[:,:,args.prefAlgm-1][useIndHi])))'''

        # Create null value array to record when newGround was populated with final value
        updateRecord=np.full((self.groundStdDev.shape),-999,dtype=int)
        # Fill with 1s where std dev is low - no need to update these further
        # Fill with 0s where std dev is high - these need to be corrected
        updateRecord[lowStdDev]=1
        updateRecord[highStdDev]=0

        # Create the newGround array and fill with preferred values at low std dev footprints
        self.newGround=np.full((self.groundStdDev.shape),np.nan)
        self.newGround[lowStdDev]=self.ground[:,:,args.prefAlgm-1][lowStdDev]
        # Write the initial newGround raster to file
        #writeToFile(self.newGround,self.minX-(args.gridSize/2),self.maxY+(args.gridSize/2),self.groundStdDev.shape[1],self.groundStdDev.shape[0],args.gridSize,args.epsg,args.plotRoot+'newGroundInitial.tif',noData=-999.0)

        # Get initial stats on binary array
        zeros,ones,noData=(self.binaryStdDev==0).sum(),(self.binaryStdDev==1).sum(),(self.binaryStdDev==-999).sum()
        print('Total zeros in binary std dev array:',zeros)
        print('Total ones in binary std dev array: {0} ({1:.2f}% of footprints)'.format(ones,(ones/(ones+zeros))*100))
        print('Total noData values in binary std dev array:',noData)

        # Implement natural neighbour interpolation on the entire dataset
        if args.updateMethod=='natural':
            lowSDx=self.eastings[lowStdDev].flatten()
            lowSDy=self.northings[lowStdDev].flatten()
            lowSDz=self.ground[:,:,args.prefAlgm-1][lowStdDev].flatten()
            # Predict at every footprint, both high and low standard deviation?
            pred=griddata((lowSDx,lowSDy),lowSDz,(self.eastings.flatten(),self.northings.flatten()),method='cubic',rescale=True).reshape(self.binaryStdDev.shape)
            if args.forceToWave:
                # Get the nearest algorithm value for every footprint
                for i in range(updateRecord.shape[0]):
                    for j in range(updateRecord.shape[1]):
                        if updateRecord[i,j]==0:
                            nearest=self.ground[i,j,np.abs(self.ground[i,j,:]-pred[i,j]).argmin()]
                            self.newGround[i,j]=nearest
                        else:
                            continue
            else:
                # Set newGround array at high SD footprints to interpolated value at same
                self.newGround[highStdDev]=pred[highStdDev]
            # Write the newGround raster to file for visual inspection
            #writeToFile(self.newGround,self.minX-(args.gridSize/2),self.maxY+(args.gridSize/2),self.groundStdDev.shape[1],self.groundStdDev.shape[0],args.gridSize,args.epsg,args.plotRoot+'newGroundNatural2SD{0:.0f}.tif'.format(args.sdThresh),noData=-999.0)

        # Implement a hybrid method combining mean and natural2 methods
        elif args.updateMethod=='hybrid':
            # Create temporary new newGround arrays for each update method
            newGroundMean=np.full((self.groundStdDev.shape),np.nan)
            newGroundMean[lowStdDev]=self.ground[:,:,args.prefAlgm-1][lowStdDev]
            newGroundNatural=np.full((self.groundStdDev.shape),np.nan)
            newGroundNatural[lowStdDev]=self.ground[:,:,args.prefAlgm-1][lowStdDev]
            # Create arrays for the min and max ground estimate from canopy height statistics
            groundMin=np.full((self.groundStdDev.shape),np.nan)
            groundMax=np.full((self.groundStdDev.shape),np.nan)

            # Implement the "natural2" update
            print('Implementing first update method ...')
            lowSDx=self.eastings[lowStdDev].flatten()
            lowSDy=self.northings[lowStdDev].flatten()
            lowSDz=self.ground[:,:,args.prefAlgm-1][lowStdDev].flatten()
            # Predict at every footprint, both high and low standard deviation?
            pred=griddata((lowSDx,lowSDy),lowSDz,(self.eastings.flatten(),self.northings.flatten()),method='cubic',rescale=True).reshape(self.binaryStdDev.shape)
            newGroundNatural[highStdDev]=pred[highStdDev]

            # Implement the "mean" update
            print('Implementing second update method ...')
            statsSummary={}
            statsSummary['1']={'Remaining':ones,'Added':0}
            print('Starting first loop, {} footprints to be corrected'.format(ones))
            searchSize=3
            while True:
                lowStep=searchSize//2
                hiStep=searchSize//2+1
                counter=0
                for i in range(updateRecord.shape[0]):
                    for j in range(updateRecord.shape[1]):
                        if updateRecord[i,j]==0:
                            lowx=max(0,j-lowStep)
                            lowy=max(0,i-lowStep)
                            hix=min(j+hiStep,updateRecord.shape[1])
                            hiy=min(i+hiStep,updateRecord.shape[0])
                            window=self.binaryStdDev[lowy:hiy,lowx:hix]
                            if (window==0).sum()<args.fThresh:
                                continue
                            else:
                                mean=np.mean(self.ground[lowy:hiy,lowx:hix,args.prefAlgm-1][window==0])
                                # How many of the 11 algorithms to use in this method - make this a cli option?
                                nearest=self.ground[i,j,np.abs(self.ground[i,j,0:11]-mean).argmin()]
                                updateRecord[i,j]=searchSize
                                newGroundMean[i,j]=nearest
                                counter+=1
                        else:
                            continue
                remaining=statsSummary[str(searchSize-2)]['Remaining']-counter
                print('Loop with step size of {} complete: {} footprints updated, {} remaining.'.format(searchSize,counter,remaining))
                statsSummary[str(searchSize)]={'Remaining':remaining,'Added':counter}
                if remaining > 0:
                    searchSize+=2
                else:
                    break
            print(statsSummary)

            # Combine results from each method
            print('Combining results from each update method ...')
            fromMean=0
            fromNatural=0
            fromMedian=0
            for i in range(updateRecord.shape[0]):
                for j in range(updateRecord.shape[1]):
                    if updateRecord[i,j]==1 or updateRecord[i,j]==-999:
                        continue
                    else:
                        # Base choice between methods on the range defined by the 9/10/11 algorithm estimates or by canopy height mean and std dev
                        # Make this a cli option?
                        algMin=np.min(self.ground[i,j,0:11])
                        algMax=np.max(self.ground[i,j,0:11])
                        if newGroundNatural[i,j]>algMax or newGroundNatural[i,j]<algMin:
                        #if newGroundNatural[i,j]>groundMax[i,j] or newGroundNatural[i,j]<groundMin[i,j]:
                            self.newGround[i,j]=newGroundMean[i,j]
                            fromMean+=1
                        else:
                            self.newGround[i,j]=newGroundNatural[i,j]
                            fromNatural+=1
            print('Hybrid update complete')
            print('{0} footprints updated by mean method, {1} footprints updated by natural method, {2} footprints updated by median method'.format(fromMean,fromNatural,fromMedian))

        # Otherwise, run one of the loop update methods
        else:
            # Initialise dictionary to record statistics about the number of newGround values added in each loop
            statsSummary={}
            # First entry is the number of high std dev footprints that need to be corrected
            statsSummary['1']={'Remaining':ones,'Added':0}
            print('Starting first loop, {} footprints to be corrected'.format(ones))

            probs=[]
            errors=[]

            # Set searchSize and window bounds
            searchSize=3
            while True:
                lowStep=searchSize//2
                hiStep=searchSize//2+1
                counter=0
                skipped=0
                pSig=0.7
                fSig=7.5
                # Loop through update record array looking for 0s - i.e., high std dev footprints
                for i in range(updateRecord.shape[0]):
                    for j in range(updateRecord.shape[1]):
                        # Get an array of the closest algorithm estimate to true ground for every footprint
                        if args.updateMethod=='idealAll':
                            if updateRecord[i,j]>-999:
                                nearest=self.ground[i,j,np.abs(self.ground[i,j,:]-self.groundTruth[i,j]).argmin()]
                                updateRecord[i,j]=searchSize
                                self.newGround[i,j]=nearest
                                counter+=1
                        elif updateRecord[i,j]==0:
                            lowx=max(0,j-lowStep)
                            lowy=max(0,i-lowStep)
                            hix=min(j+hiStep,updateRecord.shape[1])
                            hiy=min(i+hiStep,updateRecord.shape[0])
                            # Look for 0s in binary raster window - if none continue
                            window=self.binaryStdDev[lowy:hiy,lowx:hix]

                            # Get an array of the closest algorithm estimate to true ground for each high std dev footprint
                            if args.updateMethod=='ideal':
                                nearest=self.ground[i,j,np.abs(self.ground[i,j,:]-self.groundTruth[i,j]).argmin()]
                                updateRecord[i,j]=searchSize
                                self.newGround[i,j]=nearest
                                counter+=1

                            # Implement ground finding through Gaussian amplitude probability
                            if args.updateMethod=='amplitude':
                                # Calculate the Gaussian amplitude and probability for each algorithm ground feature
                                cdf=np.full((self.ground.shape[2]),0.0)
                                for k in range(self.ground.shape[2]):
                                    gInt=self.waveEnergy[i,j,k]*(1-self.cover[i,j,k])*(0.57/0.4)
                                    gWidth=sqrt(pSig**2+fSig**2*tan(self.gSlope[i,j,k])**2)
                                    gAmp=gInt/(gWidth*sqrt(2*pi))
                                    cdf[k]=norm.cdf(gAmp+self.meanNoise[i,j,k],loc=self.meanNoise[i,j,k],scale=self.noiseStdDev[i,j,k])
                                    probs.append(cdf[k])
                                    errors.append(self.ground[i,j,k]-self.groundTruth[i,j])
                                # Get the lowest of the high probability ground values
                                try:
                                    lowest=np.min(self.ground[i,j,:][cdf>=0.95])
                                    self.newGround[i,j]=lowest
                                except:
                                    skipped+=1
                                    print('#{0} - Max cdf: {1:.3f}'.format(skipped,np.max(cdf)))
                                counter+=1
                                plt.scatter(cdf,self.ground[i,j,:]-self.groundTruth[i,j],s=1)

                            # Implement ground finding by comparing to mean of neighbouring low std dev footprints
                            if args.updateMethod=='mean':
                                if (window==0).sum()<args.fThresh:
                                    continue
                                else:
                                    # We have at least one good original ground value in the window
                                    mean=np.mean(self.ground[lowy:hiy,lowx:hix,args.prefAlgm-1][window==0])
                                    nearest=self.ground[i,j,np.abs(self.ground[i,j,:]-mean).argmin()]
                                    updateRecord[i,j]=searchSize
                                    if args.forceToWave:
                                        self.newGround[i,j]=nearest
                                    else:
                                        self.newGround[i,j]=mean
                                    counter+=1

                            # Implement ground finding by comparing to mean of neighbouring low std dev footprints
                            if args.updateMethod=='height':
                                if (window==0).sum()<args.fThresh:
                                    continue
                                else:
                                    # We have at least one good original ground value in the window
                                    # Get the closest algorithm estimate as used in the 'mean' method
                                    mean=np.mean(self.ground[lowy:hiy,lowx:hix,args.prefAlgm-1][window==0])
                                    meanNearest=self.ground[i,j,np.abs(self.ground[i,j,:]-mean).argmin()]
                                    # Get the closest algorithm estimate from canopy top and tree height estimate
                                    heightDiff=(self.canopyHeight[i,j,args.prefAlgm-1])-meanCanopyHeight
                                    groundFromHeightDiff=self.ground[i,j,args.prefAlgm-1]+heightDiff
                                    groundMin=groundFromHeightDiff-stdCanopyHeight
                                    groundMax=groundFromHeightDiff+stdCanopyHeight
                                    heightNearest=self.ground[i,j,np.abs(self.ground[i,j,:]-groundFromHeightDiff).argmin()]
                                    updateRecord[i,j]=searchSize
                                    if args.forceToWave:
                                        self.newGround[i,j]=(heightNearest+meanNearest)/2
                                    else:
                                        self.newGround[i,j]=(mean+groundFromHeightDiff)/2
                                    counter+=1

                            # Implement ground finding by taking value/mean of nearest neighbour low std dev footprint(s)
                            if args.updateMethod=='nearest':
                                # Do we have at least one good footprint in the window?
                                if (window==0).sum()<args.fThresh:
                                    continue
                                else:
                                    # Make array to hold the distance to each surrounding point in the window
                                    dists=np.full((self.groundStdDev.shape),np.nan)
                                    for k in range(window.shape[0]):
                                        for l in range(window.shape[1]):
                                            dists[lowy+k,lowx+l]=sqrt(((self.eastings[lowy+k,lowx+l]-self.eastings[i,j])**2)+((self.northings[lowy+k,lowx+l]-self.northings[i,j])**2))
                                    # Get the distance to the nearest good footprint
                                    low=np.nanmin(dists[lowy:hiy,lowx:hix][window==0])
                                    # Take the mean of all good footprints in the window at the lowest distance
                                    mean=np.mean(self.ground[lowy:hiy,lowx:hix,args.prefAlgm-1][(window==0) & (dists[lowy:hiy,lowx:hix]==low)])
                                    nearest=self.ground[i,j,np.abs(self.ground[i,j,:]-mean).argmin()]
                                    updateRecord[i,j]=searchSize
                                    if args.forceToWave:
                                        self.newGround[i,j]=nearest
                                    else:
                                        self.newGround[i,j]=mean
                                    counter+=1
                                    dists[:]=np.nan

                            # Implement ground finding by fitting linear plane (1st-order polynomial) to min. of 3 neighbouring low std dev footprints
                            if args.updateMethod=='linear':
                                # Do we have the required number of good footprints in the window (min. 3 here?)?
                                if (window==0).sum()<args.fThresh:
                                    continue
                                else:
                                    # Get x-y-z in window where low std dev
                                    x=self.eastings[lowy:hiy,lowx:hix][window==0]
                                    y=self.northings[lowy:hiy,lowx:hix][window==0]
                                    z=self.ground[lowy:hiy,lowx:hix,args.prefAlgm-1][window==0]
                                    # Fit linear surface
                                    data=np.c_[x.flatten(),y.flatten(),z.flatten()]
                                    A=np.c_[data[:,0],data[:,1],np.ones(data.shape[0])]
                                    coeffs,residues,rank,s=scipy.linalg.lstsq(A,data[:,2])
                                    # Predict at chosen location
                                    pred=coeffs[0]*self.eastings[i,j] + coeffs[1]*self.northings[i,j] + coeffs[2]
                                    # Pick nearest option to prediction and update
                                    nearest=self.ground[i,j,np.abs(self.ground[i,j,:]-pred).argmin()]
                                    updateRecord[i,j]=searchSize
                                    if args.forceToWave:
                                        self.newGround[i,j]=nearest
                                    else:
                                        self.newGround[i,j]=pred
                                    counter+=1

                            # Implement ground finding by fitting quadratic surface (2nd-order polynomial) to min. of 3 neighbouring low std dev footprints
                            # Updated to follow Simon Mudd's suggestion of fitting curves to each of the nearest neighbours and taking a mean
                            if args.updateMethod=='quadratic':
                                # Do we have the required number of good footprints in the window (min. 3 here?)?
                                if (window==0).sum()<args.fThresh:
                                    continue
                                else:
                                    # Make empty list to hold the result of each individual prediction
                                    predictedValues=[]
                                    # Make array to hold the distance to each surrounding point in the window
                                    dists=np.full((self.groundStdDev.shape),np.nan)
                                    for k in range(window.shape[0]):
                                        for l in range(window.shape[1]):
                                            if window[k,l]==0:
                                                dists[lowy+k,lowx+l]=sqrt(((self.eastings[lowy+k,lowx+l]-self.eastings[i,j])**2)+((self.northings[lowy+k,lowx+l]-self.northings[i,j])**2))
                                    # For each of the n lowest values in 'dists', get their nearest neighbours
                                    nNeighbours=np.argsort(dists.flatten()) #The indices of dists, flattened and sorted
                                    for nN in range(args.fThresh):
                                        # Get the x and y of the first neighbour
                                        nNx=self.eastings.flatten()[nNeighbours[nN]]
                                        nNy=self.northings.flatten()[nNeighbours[nN]]

                                        # Next need the (minimum three) closest points to this neighbour to build prediction from
                                        nWindowSize=3
                                        while True:
                                            njIndex=int(j-((self.eastings[i,j]-nNx)/args.gridSize))
                                            niIndex=int(i+((self.northings[i,j]-nNy)/args.gridSize))
                                            lowNx=max(0,njIndex-(nWindowSize//2))
                                            lowNy=max(0,niIndex-(nWindowSize//2))
                                            hiNx=min(njIndex+(nWindowSize//2+1),updateRecord.shape[1])
                                            hiNy=min(niIndex+(nWindowSize//2+1),updateRecord.shape[0])
                                            nnWindow=self.binaryStdDev[lowNy:hiNy,lowNx:hiNx]
                                            nnDists=np.full((self.groundStdDev.shape),np.nan)
                                            for kk in range(nnWindow.shape[0]):
                                                for ll in range(nnWindow.shape[1]):
                                                    if nnWindow[kk,ll]==0:
                                                        nnDists[lowNy+kk,lowNx+ll]=sqrt(((self.eastings[lowNy+kk,lowNx+ll]-nNx)**2)+((self.northings[lowNy+kk,lowNx+ll]-nNy)**2))
                                            if (nnDists[lowNy:hiNy,lowNx:hiNx]>0).sum() < args.fThresh:
                                                nWindowSize+=1
                                            else:
                                                break
                                        # I now have at least three neighbours to this first original neighbour - sort these?
                                        nnNeighbours=np.argsort(nnDists.flatten())
                                        # Get the x,y,z of these three neighbours
                                        x=self.eastings.flatten()[nnNeighbours[1:args.fThresh+1]]
                                        y=self.northings.flatten()[nnNeighbours[1:args.fThresh+1]]
                                        z=self.ground[:,:,args.prefAlgm-1].flatten()[nnNeighbours[1:args.fThresh+1]]
                                        # Predict from these back to the original footprint
                                        # Fit polynomial
                                        data=np.c_[x,y,z]
                                        mn=np.min(data,axis=0)
                                        mx=np.max(data,axis=0)
                                        A=np.c_[np.ones(data.shape[0]),data[:,:2],np.prod(data[:,:2],axis=1),data[:,:2]**2]
                                        coeffs,residues,rank,s=scipy.linalg.lstsq(A,data[:,2])
                                        # Predict at chosen location
                                        pred=np.dot(np.c_[1,self.eastings[i,j],self.northings[i,j],self.eastings[i,j]*self.northings[i,j],self.eastings[i,j]**2,self.northings[i,j]**2],coeffs)
                                        # Add to predictedValues list
                                        predictedValues.append(pred[0])
                                    predictions=np.array(predictedValues)
                                    meanPrediction=np.mean(predictions)
                                    # Pick nearest option to prediction and update
                                    nearest=self.ground[i,j,np.abs(self.ground[i,j,:]-meanPrediction).argmin()]
                                    updateRecord[i,j]=searchSize
                                    if args.forceToWave:
                                        self.newGround[i,j]=nearest
                                    else:
                                        self.newGround[i,j]=meanPrediction
                                    counter+=1

                        else:
                            continue
                remaining=statsSummary[str(searchSize-2)]['Remaining']-counter
                print('Loop with step size of {} complete: {} footprints updated, {} remaining.'.format(searchSize,counter,remaining))
                statsSummary[str(searchSize)]={'Remaining':remaining,'Added':counter}

                if remaining > 0:
                    searchSize+=2
                else:
                    break

            print(statsSummary)

            if args.updateMethod=='amplitude':
                print('Total uncorrected waveforms:',skipped)

                plt.axhline(y=0.0,color='r',linestyle='--')
                plt.axvline(x=0.95,color='orange',linestyle='--')
                plt.xlabel('Probability')
                plt.ylabel('Ground Error (m)')
                plt.title('{0}: {1:.0f}% Beam Sensitivity'.format(args.name,args.covThresh))
                plt.savefig(args.plotRoot+'ProbVsGroundError.png',dpi=300)
                plt.close()
                plt.clf()
                print('Probability vs error plot written to '+args.plotRoot+'ProbVsGroundError.png')
                #plt.show()

                probArray=np.asarray(probs)
                errorArray=np.asarray(errors)
                x1=errorArray[np.where((probArray<=0.55))]
                x2=errorArray[np.where((probArray>0.55)&(probArray<=0.60))]
                x3=errorArray[np.where((probArray>0.60)&(probArray<=0.65))]
                x4=errorArray[np.where((probArray>0.65)&(probArray<=0.70))]
                x5=errorArray[np.where((probArray>0.70)&(probArray<=0.75))]
                x6=errorArray[np.where((probArray>0.75)&(probArray<=0.80))]
                x7=errorArray[np.where((probArray>0.80)&(probArray<=0.85))]
                x8=errorArray[np.where((probArray>0.85)&(probArray<=0.90))]
                x9=errorArray[np.where((probArray>0.90)&(probArray<=0.95))]
                x10=errorArray[np.where((probArray>0.95))]
                data=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]

                plt.boxplot(data,sym='',)
                plt.axhline(y=0.0,color='r',linestyle='--')
                plt.xticks(np.arange(1,11,step=1),['0.525','0.575','0.625','0.675','0.725','0.775','0.825','0.875','0.925','0.975'])
                plt.title('{0}: {1:.0f}% Beam Sensitivity'.format(args.name,args.covThresh))
                plt.xlabel('Probability')
                plt.ylabel('Ground Error (m)')
                plt.savefig(args.plotRoot+'ProbVsErrorBox.png',dpi=300)
                plt.close()
                plt.clf()
                print('Probability vs error box plot written to '+args.plotRoot+'ProbVsErrorBox.png')
                #plt.show()

        ### END OF UPDATE METHOD CONDITIONALS ###

        # Get stats on the newGround array and the high std dev footprints alone
        newGroundError=self.newGround-self.groundTruth
        print('newGround RMSE: {0:.2f}'.format(sqrt(np.nansum((newGroundError)**2)/np.sum(~np.isnan(newGroundError)))))
        print('newGround Bias: {0:.2f}'.format(np.nansum(newGroundError)/np.sum(~np.isnan(newGroundError))))
        print('newGround ({0} m) outlierProportion: {1:.2f}'.format(args.outlierThresh,np.nansum(np.abs((self.newGround-means))>args.outlierThresh)/np.sum(~np.isnan(self.groundTruth))*100))

        print('newGround high std. dev. RMSE: {0:.2f}'.format(sqrt(np.nansum((newGroundError[highStdDev])**2)/np.sum(~np.isnan(newGroundError[highStdDev])))))
        print('newGround high std. dev. Bias: {0:.2f}'.format(np.nansum(newGroundError[highStdDev])/np.sum(~np.isnan(newGroundError[highStdDev]))))
        print('newGround high std. dev. ({0} m) outlierProportion: {1:.2f}'.format(args.outlierThresh,np.nansum(np.abs((self.newGround[highStdDev]-means[highStdDev]))>args.outlierThresh)/np.sum(~np.isnan(self.groundTruth[highStdDev]))*100))

        # Write list of coordinates for outliers
        if args.writeCoords:
            outFile=args.plotRoot+'{0}SD{1:.0f}F{2:.0f}OutlierCoordsGt10.txt'.format(args.updateMethod,args.sdThresh,args.fThresh)
            file=open(outFile,'w')
            file.write('Eastings Northings newGround\n')
            file.close()
            with open(outFile,'a') as file:
                for i in range(newGroundError.shape[0]):
                    for j in range(newGroundError.shape[1]):
                        if np.abs(newGroundError[i,j])<=20 and np.abs(newGroundError[i,j])>10:
                            coords=str(self.eastings[i,j])+' '+str(self.northings[i,j])+' '+str(self.newGround[i,j])+'\n'
                            file.write(coords)

        if args.noUpdatePlots:
            pass
        else:
            if args.forceToWave or args.updateMethod=='hybrid':
                rastername=args.updateMethod.title()
            else:
                rastername=args.updateMethod.title()+'RV'
            # Write the updateRecord raster to file - args.fThresh in filename
            writeToFile(updateRecord,self.minX-(args.gridSize/2),self.maxY+(args.gridSize/2),self.groundStdDev.shape[1],self.groundStdDev.shape[0],args.gridSize,args.epsg,args.plotRoot+'updateRecord{0}F{1:.0f}SD{2:.0f}BS{3:.0f}.tif'.format(rastername,args.fThresh,args.sdThresh,args.covThresh),noData=-999.0)
            # Write the newGround raster to file
            writeToFile(self.newGround,self.minX-(args.gridSize/2),self.maxY+(args.gridSize/2),self.groundStdDev.shape[1],self.groundStdDev.shape[0],args.gridSize,args.epsg,args.plotRoot+'newGround{0}F{1:.0f}SD{2:.0f}BS{3:.0f}.tif'.format(rastername,args.fThresh,args.sdThresh,args.covThresh),noData=-999.0)
            # Write the newGroundError raster to file
            writeToFile(newGroundError,self.minX-(args.gridSize/2),self.maxY+(args.gridSize/2),self.groundStdDev.shape[1],self.groundStdDev.shape[0],args.gridSize,args.epsg,args.plotRoot+'newGroundError{0}F{1:.0f}SD{2:.0f}BS{3:.0f}.tif'.format(rastername,args.fThresh,args.sdThresh,args.covThresh),noData=-999.0)

        # Map out where the very high and low error footprints occur
        '''hiError=np.full((self.groundStdDev.shape),-999.0)
        hiError[newGroundError>10.0]=newGroundError[newGroundError>10.0]
        lowError=np.full((self.groundStdDev.shape),-999.0)
        lowError[newGroundError<-10.0]=newGroundError[newGroundError<-10.0]
        writeToFile(lowError,
            self.minX-(args.gridSize/2),self.maxY+(args.gridSize/2),
            self.groundStdDev.shape[1],self.groundStdDev.shape[0],
            args.gridSize,args.epsg,
            args.plotRoot+'newGroundErrorLow{}.tif'.format(args.fThresh),noData=-999.0)
        writeToFile(hiError,
            self.minX-(args.gridSize/2),self.maxY+(args.gridSize/2),
            self.groundStdDev.shape[1],self.groundStdDev.shape[0],
            args.gridSize,args.epsg,
            args.plotRoot+'newGroundErrorHigh{}.tif'.format(args.fThresh),noData=-999.0)'''

        if args.noUpdatePlots:
            pass
        else:
            minVal=min(np.nanmin(self.groundError[:,:,args.prefAlgm-1]),np.nanmin(newGroundError))
            maxVal=max(np.nanmax(self.groundError[:,:,args.prefAlgm-1]),np.nanmax(newGroundError))
            valRange=(minVal,maxVal)
            numBins=50 if (maxVal-minVal) <=50 else 100     #"Ternary operator"
            # Plot histograms of ground error before and after correcting - args.fThresh in filename
            fig,axes=plt.subplots(1,2)
            _, bins, _=axes[0].hist(self.groundError[:,:,args.prefAlgm-1][~np.isnan(self.groundError[:,:,args.prefAlgm-1])].flatten(),bins=numBins,range=valRange,
                label='Alg. {0}\n Mean: {1:.2f}\n Std. Dev: {2:.2f}'.format(args.prefAlgm,np.nanmean(self.groundError[:,:,args.prefAlgm-1]),np.nanstd(self.groundError[:,:,args.prefAlgm-1])))
            axes[0].hist(newGroundError[~np.isnan(newGroundError)],bins=bins,range=valRange,alpha=0.5,
                label='newGround\n Mean: {0:.2f}\n Std. Dev: {1:.2f}'.format(np.nanmean(newGroundError),np.nanstd(newGroundError)))
            axes[0].axvline(x=0,color='r',linestyle='--')
            axes[0].set(title='All Footprints',ylabel='Frequency')
            axes[0].legend()

            _, bins, _=axes[1].hist(self.groundError[:,:,args.prefAlgm-1][highStdDev].flatten(),bins=numBins,range=valRange,
                label='Alg. {0}\n Mean: {1:.2f}\n Std. Dev: {2:.2f}'.format(args.prefAlgm,np.nanmean(self.groundError[:,:,args.prefAlgm-1][highStdDev]),np.nanstd(self.groundError[:,:,args.prefAlgm-1][highStdDev])))
            axes[1].hist(newGroundError[highStdDev].flatten(),bins=bins,range=valRange,alpha=0.5,
                label='newGround\n Mean: {0:.2f}\n Std. Dev: {1:.2f}'.format(np.nanmean(newGroundError[highStdDev]),np.nanstd(newGroundError[highStdDev])))
            axes[1].axvline(x=0,color='r',linestyle='--')
            axes[1].set(title='High Std. Dev. Footprints')
            axes[1].legend()
            for axis in axes.flat:
                axis.set(xlabel='Ground Error (m)')
                axis.set_xlim(-40,40)
                axis.set_ylim(0,10000)
            #for axis in axes.flat:
                #axis.label_outer()
            if args.forceToWave or args.updateMethod=='hybrid':
                plotname=args.updateMethod.title()
            else:
                plotname=args.updateMethod.title()+'RV'
            plt.savefig(args.plotRoot+'ErrorHist{0}SD{1:.0f}BS{2}.pdf'.format(plotname,args.sdThresh,args.covThresh),dpi=300)
            plt.close()
            plt.clf()
            print('Ground error histogram plots written to '+args.plotRoot+'ErrorHist{0}SD{1:.0f}BS{2}.pdf'.format(plotname,args.sdThresh,args.covThresh))
            #plt.show()

    def slopeFit(self):
        '''Experiment with linear slope fitting over a region'''

        # Ground should be broadly linear over short length scales so try plane fitting over small regions?
        # Generate a small test dataset
        testX=self.eastings[101:104,101:104][self.binaryStdDev[101:104,101:104]==0]
        print('X:')
        print(testX)
        testY=self.northings[101:104,101:104][self.binaryStdDev[101:104,101:104]==0]
        print('Y:')
        print(testY)
        testZ=self.ground[101:104,101:104,6][self.binaryStdDev[101:104,101:104]==0]
        print('Z:')
        print(testZ)
        testTrueZ=self.groundTruth[101:104,101:104][self.binaryStdDev[101:104,101:104]==0]
        print('True Z:')
        print(testTrueZ)
        testStd=self.binaryStdDev[101:104,101:104]
        print('Binary Std. Dev.:')
        print(testStd)

        print('Test point Z:',self.ground[102,102,6])
        print('Test point true Z:',self.groundTruth[102,102])

        # Fit linear surface
        data=np.c_[testX.flatten(),testY.flatten(),testZ.flatten()]
        mn=np.min(data,axis=0)
        mx=np.max(data,axis=0)
        X,Y=np.meshgrid(np.linspace(mn[0],mx[0],100),np.linspace(mn[1],mx[1],100))
        XX=X.flatten()
        YY=Y.flatten()

        # Linear plane (1st-order)
        A=np.c_[data[:,0],data[:,1],np.ones(data.shape[0])]
        coeffs,residues,rank,s=scipy.linalg.lstsq(A,data[:,2])
        # Evaluate on grid
        Z = coeffs[0]*X + coeffs[1]*Y + coeffs[2]
        # Predict at chosen location
        pred=coeffs[0]*824790 + coeffs[1]*1153410 + coeffs[2]
        print('Predicted height at 824790,1153410:',pred)

        # Quadratic curve (2nd-order)
        #A=np.c_[np.ones(data.shape[0]),data[:,:2],np.prod(data[:,:2],axis=1),data[:,:2]**2]
        #coeffs,residues,rank,s=scipy.linalg.lstsq(A,data[:,2])
        # Evaluate on grid
        #Z = np.dot(np.c_[np.ones(XX.shape),XX,YY,XX*YY,XX**2,YY**2],coeffs).reshape(X.shape)
        # Predict at chosen location
        #pred = np.dot(np.c_[1,824790,1153410,824790*1153410,824790**2,1153410**2],coeffs)
        #print('Predicted height at 824790,1153410:',pred)

        print('Least-squares coefficients:',coeffs[0],coeffs[1],coeffs[2])
        print('Residues:',residues)
        print('Rank:',rank)
        print('s:',s)

        # Plot test data in 3D
        fig=plt.figure(figsize=(10,10))
        ax=fig.gca(projection='3d')
        ax.plot_surface(X,Y,Z,rstride=1,alpha=0.2)
        ax.scatter(testX,testY,testTrueZ)
        plt.xlabel('Easting (m)')
        plt.ylabel('Northing (m)')
        ax.set_zlabel('Elevation (m)')
        plt.show()

    def ampTest(self):
        '''Experiment with extracting the amplitude of each Gaussian and calculating probability'''

        x=351010
        y=8106570
        useInd=np.where((self.eastings==x) & (self.northings==y))
        with np.printoptions(precision=3,floatmode='fixed'):
            print('')
            print('')
            print('Ground options:',self.ground[useInd])
            print('Cover options:',self.cover[useInd])
            print('GroundSlope options:',self.gSlope[useInd])
            print('WaveEnergy options:',self.waveEnergy[useInd])
            print('MeanNoise options:',self.meanNoise[useInd])
            print('NoiseStdDev options:',self.noiseStdDev[useInd])
            print('')

        pSig=0.7
        fSig=7.5
        ground=self.ground[useInd]
        groundTruth=self.groundTruth[useInd]
        cover=self.cover[useInd]
        gSlope=self.gSlope[useInd]
        waveEnergy=self.waveEnergy[useInd]
        meanNoise=self.meanNoise[useInd]
        noiseStdDev=self.noiseStdDev[useInd]
        #print('cover, gSlope, waveEnergy, meanNoise, noiseStdDev')
        #print(cover, gSlope, waveEnergy, meanNoise, noiseStdDev)

        cdf1=np.full((11),0.0)
        cdf2=np.full((11),0.0)
        for i in range(11):
            gInt=waveEnergy[0,i]*(1-cover[0,i])*(0.57/0.4)
            gWidth=sqrt(pSig**2+fSig**2*tan(gSlope[0,i])**2)
            gAmp=gInt/(gWidth*sqrt(2*pi))
            cdf1[i]=norm.cdf(gAmp+meanNoise[0,i],loc=meanNoise[0,i],scale=noiseStdDev[0,i])
            cdf2[i]=norm.cdf(gAmp,loc=meanNoise[0,i],scale=noiseStdDev[0,i])
            print('cdf1, cdf2, amplitude, ground, cover, gSlope, waveEnergy, meanNoise, noiseStdDev, groundTruth')
            print(cdf1[i], cdf2[i], gAmp, ground[0,i], cover[0,i], gSlope[0,i], waveEnergy[0,i], meanNoise[0,i], noiseStdDev[0,i], groundTruth[0])
            print('')
        #print('cdf2:',cdf2)

        print('Sorted ground options:',np.sort(ground[0]))
        print('High probability ground options:',ground[0,cdf1>=0.95])
        print('Low probability ground options:',ground[0,cdf1<0.95])
        print('Lowest high probability option:',np.min(ground[0,cdf1>=0.95]))

        plt.scatter(cdf1,ground-groundTruth,s=1)
        plt.axhline(y=0.0,color='r',linestyle='--')
        plt.axvline(x=0.95,color='orange',linestyle='--')
        plt.xlabel('Probability')
        plt.ylabel('Ground Error (m)')
        plt.title('Footprint: {0:.0f}, {1:.0f}'.format(x,y))
        plt.show()

#########################################

args=readCommands()
# The main block
if __name__ == '__main__':
    #args=readCommands()

    if args.suppress:
        warnings.filterwarnings('ignore',category=RuntimeWarning)

    if args.site=='laselva':
        args.rootDir='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/laselva/grid{}/'.format(args.gridSize)
        args.epsg=32616
        args.plotRoot='../data/laSelvaGroundError/{}-{}/'.format(args.gridSize,args.covThresh)
        args.name='La Selva'
        siteCanopyHeight=29
    elif args.site=='robson':
        args.rootDir='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid{}/'.format(args.gridSize)
        args.epsg=32755
        args.plotRoot='../data/robsonGroundError/{}-{}/'.format(args.gridSize,args.covThresh)
        args.name='Robson Creek'
        siteCanopyHeight=33
    elif args.site=='nouragues':
        args.rootDir='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/nouragues/grid{}/'.format(args.gridSize)
        args.epsg=32622
        args.plotRoot='../data/nouraguesGroundError/{}-{}/'.format(args.gridSize,args.covThresh)
        args.name='Nouragues'
        siteCanopyHeight=38
    elif args.site=='paracou':
        args.rootDir='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/paracou/grid{}/'.format(args.gridSize)
        args.epsg=32622
        args.plotRoot='../data/paracouGroundError/{}-{}/'.format(args.gridSize,args.covThresh)
        args.name='Paracou'
        siteCanopyHeight=34
    elif args.site=='hubbard':
        args.rootDir='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/hubbard/grid{}/'.format(args.gridSize)
        args.epsg=32619
        args.plotRoot='../data/hubbardGroundError/{}-{}/'.format(args.gridSize,args.covThresh)
        args.name='Hubbard Brook'
        siteCanopyHeight=''
    elif args.site=='windriver':
        args.rootDir='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/windriver/grid{}/'.format(args.gridSize)
        args.epsg=32610
        args.plotRoot='../data/windRiverGroundError/{}-{}/'.format(args.gridSize,args.covThresh)
        args.name='Wind River'
        siteCanopyHeight=''
    elif args.site=='oakridge':
        args.rootDir='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/oakridge/grid{}/'.format(args.gridSize)
        args.epsg=32616
        args.plotRoot='../data/oakRidgeGroundError/{}-{}/'.format(args.gridSize,args.covThresh)
        args.name='Oak Ridge'
        siteCanopyHeight=''
    else:
        print('WARNING: No valid site name given, exiting programme ...')
        sys.exit()

    ground=groundFinder()
    ground.loadData(args.rootDir,args.gridSize,args.sdThresh)
    if args.plots:
        ground.plotGroundRange(args.sdThresh)
    if args.sdRaster:
        ground.stdDevRasters()
    if args.update:
        ground.updateGround(siteCanopyHeight)
    #ground.slopeFit()
    #ground.ampTest()
