import numpy as np

minX=821700
maxX=828900
minY=1149300
maxY=1156500

step=1200
buffer=5

inList='~/glamis/code/groundBounds.nasa_laselva2009_newground.txt'
outRoot='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/laselva/grid10/'
outFile='shell/gediRatLaSelvaNoise10.sh'

with open(outFile,'a') as outf:
    outf.write('#!/usr/bin/env bash\n\n')

counter=1
for x in np.arange(minX,maxX,step):
    for y in np.arange(minY,maxY,step):
        output=outRoot+str(x)+'.'+str(y)+'.h5'
        with open(outFile,'a') as outf:
            outf.write('echo "Processing subset ' + str(counter) + ' of 36"\n')
            outf.write('gediRat -inList {0} -output {1} -ground -hdf -gridBound {2} {3} {4} {5} -gridStep 10 -pSigma 0.7 -fSigma 2.5 -checkCover -pBuff 2 -countOnly\n\n'.format(inList,output,x,x+step-buffer,y,y+step-buffer))
        counter+=1

with open(outFile,'a') as outf:
    outf.write('echo "All done"')
