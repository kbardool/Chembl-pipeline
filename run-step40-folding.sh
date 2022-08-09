#!/bin/bash

LEADERFOLLOWER="python ../leader_follower/cluster.py"
VERSION=chembl_29
# cd ..


echo Creating folds...
python mtcv.py  --y         output/${VERSION}/${VERSION}_thresholds.csv \
                --clusters  output/${VERSION}/${VERSION}_clustering.npy \
                --folding   output/${VERSION}/${VERSION}_folding.npy \
                --compounds output/${VERSION}/${VERSION}_X_cmpds.csv

echo "run.sh - Step 40 finished "