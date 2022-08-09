#!/bin/bash

VERSION=chembl_29_dev
LEADERFOLLOWER="../leader_follower/cluster.py"
# LEADERFOLLOWER_ORIGINAL="../leader_follower_original/"
# LEADERFOLLOWER_PROGRAM="cluster.py"

echo Step 30 - Clustering compounds...
# original line : $LEADERFOLLOWER --x output/${VERSION}_X.npy --out output/clustering.npy --dist 0.5 0.7 
# line expands to the following:
#    python ../leader_follower/cluster.py --x output/${VERSION}/${VERSION}_X.npy --out output/${VERSION}_clustering.npy --dist 0.5 0.7 
python $LEADERFOLLOWER \
    --x   output/${VERSION}/${VERSION}_X.npy \
    --out output/${VERSION}/${VERSION}_clustering.npy \
    --dist 0.5 0.7 
echo "run.sh - Step 30 finished "
 

# python $LEADERFOLLOWER \
    #  --x   output/${VERSION}/adam_${VERSION}_X.npy \
    #  --out output/${VERSION}/adam_clustering.npy \
    #  --dist 0.5 0.7 

