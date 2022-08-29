#!/usr/bin/env bash

# Run gediMetric over the output of gediRat grid sims.
# Each loop over files generates one of the 11 algorithms used.

files=/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/*.h5

for file in $files
do
  echo "Processing file $file"
  root="/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/metric98-1/$(basename -- $file .h5)"
  gediMetric -input $file -readHDFgedi -ground -linkNoise 0 0.98 -linkFsig 2.5 -linkPsig 0.7 -minGsig 0.7 -sWidth 0.5 -varScale 2 -minWidth 3 -outRoot $root
done

###############

files=/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/*.h5

for file in $files
do
  echo "Processing file $file"
  root="/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/metric98-10/$(basename -- $file .h5)"
  gediMetric -input $file -readHDFgedi -ground -linkNoise 0 0.98 -linkFsig 2.5 -linkPsig 0.7 -minGsig 0.7 -sWidth 0.5 -varScale 1 -minWidth 3 -outRoot $root
done

###############

files=/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/*.h5

for file in $files
do
  echo "Processing file $file"
  root="/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/metric98-150/$(basename -- $file .h5)"
  gediMetric -input $file -readHDFgedi -ground -linkNoise 0 0.98 -linkFsig 2.5 -linkPsig 0.7 -minGsig 0.7 -sWidth 0.5 -varScale 1.5 -minWidth 3 -outRoot $root
done

##############

files=/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/*.h5

for file in $files
do
  echo "Processing file $file"
  root="/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/metric98-2/$(basename -- $file .h5)"
  gediMetric -input $file -readHDFgedi -ground -linkNoise 0 0.98 -linkFsig 2.5 -linkPsig 0.7 -minGsig 0.7 -sWidth 0.5 -varScale 3 -minWidth 3 -outRoot $root
done

##############

files=/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/*.h5

for file in $files
do
  echo "Processing file $file"
  root="/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/metric98-3/$(basename -- $file .h5)"
  gediMetric -input $file -readHDFgedi -ground -linkNoise 0 0.98 -linkFsig 2.5 -linkPsig 0.7 -minGsig 0.7 -sWidth 0.5 -varScale 6 -minWidth 3 -outRoot $root
done

#############

files=/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/*.h5

for file in $files
do
  echo "Processing file $file"
  root="/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/metric98-4/$(basename -- $file .h5)"
  gediMetric -input $file -readHDFgedi -ground -linkNoise 0 0.98 -linkFsig 2.5 -linkPsig 0.7 -minGsig 0.7 -sWidth 1.0 -varScale 2 -minWidth 3 -outRoot $root
done

##############

files=/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/*.h5

for file in $files
do
  echo "Processing file $file"
  root="/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/metric98-5/$(basename -- $file .h5)"
  gediMetric -input $file -readHDFgedi -ground -linkNoise 0 0.98 -linkFsig 2.5 -linkPsig 0.7 -minGsig 0.7 -sWidth 1.0 -varScale 3 -minWidth 3 -outRoot $root
done

#############

files=/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/*.h5

for file in $files
do
  echo "Processing file $file"
  root="/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/metric98-6/$(basename -- $file .h5)"
  gediMetric -input $file -readHDFgedi -ground -linkNoise 0 0.98 -linkFsig 2.5 -linkPsig 0.7 -minGsig 0.7 -sWidth 1.0 -varScale 6 -minWidth 3 -outRoot $root
done

###############

files=/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/*.h5

for file in $files
do
  echo "Processing file $file"
  root="/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/metric98-7/$(basename -- $file .h5)"
  gediMetric -input $file -readHDFgedi -ground -linkNoise 0 0.98 -linkFsig 2.5 -linkPsig 0.7 -minGsig 0.7 -sWidth 0.5 -varScale 2 -minWidth 5 -outRoot $root
done

#############

files=/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/*.h5

for file in $files
do
  echo "Processing file $file"
  root="/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/metric98-8/$(basename -- $file .h5)"
  gediMetric -input $file -readHDFgedi -ground -linkNoise 0 0.98 -linkFsig 2.5 -linkPsig 0.7 -minGsig 0.7 -sWidth 0.5 -varScale 3 -minWidth 5 -outRoot $root
done

###############

files=/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/*.h5

for file in $files
do
  echo "Processing file $file"
  root="/exports/csce/datastore/geos/groups/3d_env/data/glamis/emitchell/noisedGrids/robson/grid10/metric98-9/$(basename -- $file .h5)"
  gediMetric -input $file -readHDFgedi -ground -linkNoise 0 0.98 -linkFsig 2.5 -linkPsig 0.7 -minGsig 0.7 -sWidth 0.5 -varScale 6 -minWidth 5 -outRoot $root
done

echo "All done"

