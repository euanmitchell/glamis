############################################
# Script to generate a shell script to
# concatenate and gawk the output of
# gediMetric runs to generate the text
# files needed by groundFinder.py.
############################################

import argparse

# Defining the command line reading function
def readCommands():
    '''
    Read the arguments passed from the command line
    '''
    p=argparse.ArgumentParser(description=('Specify parameters for program control'),
        formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('--cat', dest='cat', action='store_true', default=False,
        help=('Concatenate gediMetric output text files if not already done.'))
    p.add_argument('--removeOld', dest='removeOld', action='store_true', default=False,
        help=('Remove any existing gediMetric summary files.'))
    p.add_argument('--gridSize', dest='gridSize', type=int, default=30,
        help=('''The footprint spacing. Default is 30'''))
    p.add_argument('--scriptName', dest='scriptName', type=str, default='gawk.sh',
        help=('''The name (and file path) of the output shell script. Default is "gawk.sh".'''))
    cmdargs=p.parse_args()
    return cmdargs

args=readCommands()

#########################################

#sites=['laselva','robson','nouragues','paracou','hubbard','windriver','oakridge']
#thresholds=['99','98','97','96','95','94','93','92']
sites=[]
thresholds=[]
algorithms=['1','10','150','2','3','4','5','6','7','8','9']

rootDir='/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/'
outShell=args.scriptName

with open(outShell,'w') as outf:
    outf.write('#!/usr/bin/env bash\n\n')

if args.removeOld:
    # Remove existing summary files
    for site in sites:
        with open(outShell,'a') as outf:
            outf.write('echo "Deleting old {0} files ...."\n\n'.format(site))
        for threshold in thresholds:
            with open(outShell,'a') as outf:
                outf.write('echo "Threshold {0} ..."\n\n'.format(threshold))
            for algorithm in algorithms:
                file1='{0}{1}/grid{2}/metric{3}-{4}/metric{2}-{3}-{4}All.txt'.format(rootDir,site,args.gridSize,threshold,algorithm)
                file2='{0}{1}/grid{2}/metric{3}-{4}/metric{2}-{3}-{4}Summary.txt'.format(rootDir,site,args.gridSize,threshold,algorithm)
                with open(outShell,'a') as outf:
                    outf.write('rm {0}\n\n'.format(file1))
                    outf.write('rm {0}\n\n'.format(file2))

if args.cat:
    # Concatenate raw metric files
    for site in sites:
        with open(outShell,'a') as outf:
            outf.write('echo "Concatenating {0} data ...."\n\n'.format(site))
        for threshold in thresholds:
            with open(outShell,'a') as outf:
                outf.write('echo "Threshold {0} ..."\n\n'.format(threshold))
            for algorithm in algorithms:
                inFile='{0}{1}/grid{2}/metric{3}-{4}/*.txt'.format(rootDir,site,args.gridSize,threshold,algorithm)
                outFile='{0}{1}/grid{2}/metric{3}-{4}/metric{2}-{3}-{4}All.txt'.format(rootDir,site,args.gridSize,threshold,algorithm)
                with open(outShell,'a') as outf:
                    outf.write('cat {0} > {1}\n\n'.format(inFile,outFile))

# Gawk the concatenated output
for site in sites:
    with open(outShell,'a') as outf:
        outf.write('echo "Processing {0} data ...."\n\n'.format(site))
    for threshold in thresholds:
        with open(outShell,'a') as outf:
            outf.write('echo "Threshold {0} ..."\n\n'.format(threshold))
        for algorithm in algorithms:
            inFile='{0}{1}/grid{2}/metric{3}-{4}/metric{2}-{3}-{4}All.txt'.format(rootDir,site,args.gridSize,threshold,algorithm)
            outFile='{0}{1}/grid{2}/metric{3}-{4}/metric{2}-{3}-{4}Summary.txt'.format(rootDir,site,args.gridSize,threshold,algorithm)
            with open(outShell,'a') as outf:
                outf.write('''gawk 'BEGIN{print "x,y,trueGround,beamSens,ALScover,rhReal95,gHeight,cover,rhGauss95,FHD,rhReal50,gSlope,waveEnergy,meanNoise,noiseStdDev"}($1!="#"){if($2>-1000000)print $107,$108,$2,$113,$5,$96,$6,$11,$33,$117,$87,$148,$112,$120,$121}' < '''+inFile+''' > '''+outFile+'''\n\n''')

with open(outShell,'a') as outf:
    outf.write('echo "All Done"')

print('Shell file written to',outShell)
