#######################################
# Library of tools for
# working with raster data.
#######################################

import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
from osgeo import osr

#######################################

def writeToFile(data,minX,maxY,nX,nY,res,epsg,outFile,noData=-999.0):
    '''Write 2-D numpy array to GeoTiff raster image'''

    # Set geolocation and load data to geotiff object
    geotransform=(minX, res, 0, maxY, 0, -res)

    dst_ds=gdal.GetDriverByName('GTiff').Create(outFile, nX, nY, 1, gdal.GDT_Float32)
    dst_ds.SetGeoTransform(geotransform)
    srs=osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    dst_ds.SetProjection(srs.ExportToWkt())
    dst_ds.GetRasterBand(1).WriteArray(data)
    dst_ds.GetRasterBand(1).SetNoDataValue(noData)
    dst_ds.FlushCache()
    dst_ds=None

    print('Image written to',outFile)

def readFromFile(input):
    '''Read a GeoTiff into a 2D numpy array'''

    print('Reading {} into array ...'.format(input))
    ds=gdal.Open(input)
    proj=osr.SpatialReference(wkt=ds.GetProjection())
    epsg=int(proj.GetAttrValue('AUTHORITY',1))
    nX=ds.RasterXSize
    nY=ds.RasterYSize
    transform_ds=ds.GetGeoTransform()
    xOrigin=transform_ds[0]
    yOrigin=transform_ds[3]
    pixelWidth=transform_ds[1]
    pixelHeight=transform_ds[5]
    dataArray=ds.GetRasterBand(1).ReadAsArray(0,0,nX,nY)
    print('Array shape:',dataArray.shape)
    print('Data min/max:',np.amin(dataArray[dataArray>-999]),np.amax(dataArray[dataArray>-999]))

    return dataArray,epsg

def defineRegions(imageArray,binaryArray,res,epsg):
    '''Generate new GeoTiff with regions of contiguous cells greater than threshold value.
    Requires the original image array and the binary array generated from that as inputs.'''

    # Get size of the image portion of the image array (i.e., not NoData value cells)
    imageData=np.sum(np.where(imageArray==-999.0,0,1))
    print('Data cells in the image array is',imageData)

    regionArray=binaryArray.copy()
    numCells=np.sum(np.where(binaryArray==0.0,0,1))
    print('Number of nonzero cells in binary array is',numCells)
    numAllCells=np.size(binaryArray)
    print('Total number of cells in binary array is',numAllCells)
    print('Nonzero cells are {:.2f}% of the image'.format(numCells/imageData*100))
    rows=regionArray.shape[0]
    cols=regionArray.shape[1]

    summaryStats={}
    increments=[]
    imagePerc=[]
    nonzeroPerc=[]

    counter=0
    increment=2
    while True:
        for i in range(rows):
            for j in range(cols):
                if binaryArray[i,j]==1.0:
                    maxx=min(j+increment,cols)
                    maxy=min(i+increment,rows)
                    sum=np.sum(binaryArray[i:maxy,j:maxx])
                    if sum==float(increment**2):
                        regionArray[i:maxy,j:maxx]=float(increment**2)
                        counter+=1
        print('There were {} regions of {}x{} contiguous cells identified'.format(counter,increment,increment))
        imageArea=(np.sum(np.where(regionArray==float(increment**2),1,0)))/imageData*100
        nonzeroArea=(np.sum(np.where(regionArray==float(increment**2),1,0)))/numCells*100
        print('These regions account for {:.2f}% of the data raster'.format(imageArea))
        print('These regions account for {:.2f}% of the nonzero raster'.format(nonzeroArea))
        summaryStats[increment]=counter
        increments.append(increment)
        imagePerc.append(imageArea)
        nonzeroPerc.append(nonzeroArea)
        if counter==0:
            break
        else:
            counter=0
            increment+=1

    print(summaryStats)
    plt.bar(increments,nonzeroPerc,label='Nonzero Cells')
    plt.bar(increments,imagePerc,label='Entire Image')
    plt.xlabel('Square Side Length')
    plt.ylabel('Area of Raster (%)')
    plt.legend()
    plt.show()

    return regionArray
