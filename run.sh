#!/bin/bash

LEADERFOLLOWER="python ../leader_follower/cluster.py"

if [ -d "input" ] 
then
    echo "Directory ./input exists." 
else
    echo "Directory ./input does not exists."
    mkdir input
fi

cd input

if [ -n "${VERSION}" ]; then
    echo "Requested version: ${VERSION}"
else
    echo "Default: chembl_29"
    VERSION=chembl_29
fi

if [ -f "${VERSION}.sdf.gz" ] ; then
    echo "Chemical structure data exists."
else
    echo "Downloading Chemical structure version ${VERSION}"
    wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/${VERSION}/${VERSION}.sdf.gz
    echo Gunziping the sdf...
    gunzip -k ${Version}.sdf.gz
fi

if [ -f "${VERSION}_sqlite.tar.gz" ]; then
    echo "Database exists."
else
    echo "Downloading sqlite database version ${VERSION}"
    wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/${VERSION}/${VERSION}_sqlite.tar.gz
    echo Gunziping the database...
    tar -xf ${VERSION}_sqlite.tar.gz
fi

cd ..

# echo Step 10 - Create Dataset ...
# python 10_create_dataset.py \
#         --sqlite input/${VERSION}_sqlite/${VERSION}.db \
#         --datadir Step10
#         --prefix ${VERSION}
## Writes the followig files:
##        output/${VERSION}/${VERSION}_thresh
##        output/${VERSION}/${VERSION}_compounds
##        output/${VERSION}/${VERSION}_targets
# echo "run.sh - Step 10 finished "


# echo Step 20 - SDF to ECFP ...
# python 20_sdfToECFP.py -s input/${VERSION}.sdf \
#                        -o output/${VERSION}_X.csv \
#                        -c output/${VERSION}/${VERSION}_compounds.csv

# mv fp_32000.npy    output/${VERSION}/${VERSION}_folding.npy
# mv cmpd_list_X.csv output/${VERSION}/${VERSION}_X_cmpds.csv
# echo "run.sh - Step 20 finished "

# echo Step 30 - Clustering compounds...
# $LEADERFOLLOWER \
#     --x   output/${VERSION}/${VERSION}_X.npy \
#     --out output/${VERSION}/${VERSION}_clustering.npy \
#     --dist 0.5 0.7 
# ####   above line expands to the following:
# ####   python ../leader_follower/cluster.py --x output/${VERSION}/${VERSION}_X.npy --out output/${VERSION}_clustering.npy --dist 0.5 0.7 
# echo "run.sh - Step 30 finished "

# echo Creating folds...
# python mtcv.py  --y         output/${VERSION}/${VERSION}_thresh.csv \
#                 --clusters  output/${VERSION}/${VERSION}_clustering.npy \
#                 --folding   output/${VERSION}/${VERSION}_folding.npy \
#                 --compounds output/${VERSION}/${VERSION}_X_cmpds.csv

echo "run.sh finished "