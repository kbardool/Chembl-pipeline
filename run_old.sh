#!/bin/bash

LEADERFOLLOWER="python ../leader_follower/cluster.py"

# if [ -d "input" ] 
# then
#     echo "Directory ./input exists." 
# else
#     echo "Directory ./input does not exists."
#     mkdir input
# fi

cd input

if [ -n "${VERSION}" ]; then
    echo "Requested version: ${VERSION}"
else
    echo "Default: chembl_29"
    VERSION=chembl_29
fi

# if [ -f "${VERSION}.sdf.gz" ] ; then
#     echo "Chemical structure data exists."
# else
#     echo "Downloading Chemical structure version ${VERSION}"
#     wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/${VERSION}.sdf.gz
#     echo Gunziping the sdf...
#     gunzip -k ${Version}.sdf.gz
# fi

# if [ -f "${VERSION}_sqlite.tar.gz" ]; then
#     echo "Database exists."
# else
#     echo "Downloading sqlite database version ${VERSION}"
#     wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/${VERSION}_sqlite.tar.gz
#     echo Gunziping the database...
#     tar -xf ${VERSION}_sqlite.tar.gz
# fi

cd ..
# echo Creating Dataset
# python 10_create_dataset.py --sqlite input/${VERSION}_sqlite/${VERSION}.db --prefix ${VERSION}

# echo sdfT0ECFP ....
# python 20_sdfToECFP.py -s input/${VERSION}.sdf -o output/${VERSION}_X.csv -c output/${VERSION}/${VERSION}_compounds.csv
# mv fp_32000.npy    output/${VERSION}_X.npy
# mv cmpd_list_X.csv output/${VERSION}_X_cmpds.csv

echo Clustering compounds...
# original line : $LEADERFOLLOWER --x output/${VERSION}_X.npy --out output/clustering.npy --dist 0.5 0.7 
python  ../leader_follower/cluster.py   \
        --x    output/${VERSION}_X.npy  \
        --out  output/clustering.npy    \
        --dist 0.5 0.7 

echo Creating folds...
python 40_mtcv.py --y output/${VERSION}/${VERSION}_thresh.csv\
               --clusters  output/clustering.npy \
               --folding   output/folding.npy \
               --compounds output/${VERSION}_X_cmpds.csv

echo "run.sh finished "