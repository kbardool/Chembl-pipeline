import numpy as np
import pandas as pd
import scipy.sparse

def mtcv(mask, nfolds=None, pfolds=None, seed=None):
    """
    Return a vector of folds

    Args:
    mask     binary mask matrix of [compounds x targets]
    nfolds   number of folds (integer)
    pfolds   array, specifying fold sizes (probability)
             If specified nfolds is ignored.
             Must sum to 1.
    """
    if nfolds is None and pfolds is None:
        raise ValueError("nfolds or pfolds must be specified.")
    if pfolds is not None and np.abs(np.sum(pfolds) - 1.0) > 1e-5:
        raise ValueError("pfolds must sum to 1.0.")

    if pfolds is None:
        pfolds = np.ones(nfolds) / nfolds
    else:
        nfolds = len(pfolds)

    target_sizes = np.array(mask.sum(0)).flatten()
    comp_sizes   = np.array(mask.sum(1)).flatten()

    if seed is not None:
        np.random.seed(seed)

    df = pd.DataFrame({"row": np.arange(mask.shape[0]), "size": comp_sizes})
    ## Randomize the dataframe
    df1 = df.sample(frac=1)
    
    ## Sort in descending compound count order 
    df1.sort_values("size", inplace=True, ascending=False)

    fold_sizes = np.zeros((nfolds, mask.shape[1]), dtype=np.int8)
    #max_fold_sizes = np.ceil(target_sizes / nfolds).astype(np.int)
    max_fold_sizes = np.ceil(np.outer(target_sizes, pfolds))

    ## output vector (initialized to -1)
    folds = np.zeros(mask.shape[0], dtype=np.int8) - 1
    farray = np.arange(nfolds)    

    for i in range(mask.shape[0]):
        idx   = df1.index[i]

        #for j in np.random.permutation(nfolds):
        # make a random list of the folds, and go through the list 
        for j in np.random.choice(farray, size=nfolds, replace=False, p=pfolds):
            if mask[idx,:].dot((fold_sizes[j,:] + mask[idx,:] > max_fold_sizes[:,j]).transpose()) == 0:
                ## compound fits into fold j
                folds[idx] = j
                fold_sizes[j,:] += mask[idx,:]
                break
        if folds[idx] == -1:
            folds[idx] = np.random.randint(0, nfolds)
            fold_sizes[folds[idx],:] += mask[idx,:]

    return folds

def mtcv_clustered(Y, clusters, nfolds=None, pfolds=None, seed=None):
    """splits rows of Y, based on clusters into either 
    a) equally into nfolds 
    b) or according to the ratios defined by pfolds
    """

    assert clusters.shape[0] == Y.shape[0]
    cl_uniq = clusters.unique()
    cl2id   = pd.Series(np.arange(cl_uniq.shape[0]), index=cl_uniq)
    cid     = cl2id[clusters]
    print(" Number of unique clusters: {cl_uniq}")
    
    
    ## creating cluster2compound matrix
    ##
    ## Number of rows: cid.values
    C    = scipy.sparse.csr_matrix(
            (np.ones(cid.shape[0], dtype=np.int8), (cid.values, np.arange(cid.shape[0])))
            )


    Ybin      = Y.copy()
    Ybin.data = np.ones(Ybin.data.shape[0], dtype=np.int8)
    cl_counts = C.dot(Ybin)

    ## get cluster folds
    folds = mtcv(mask=cl_counts, nfolds=nfolds, pfolds=pfolds, seed=seed)

    return folds[cid]

import argparse

parser = argparse.ArgumentParser(description="Multi-task clustered crossvalidation")
parser.add_argument('--y', type=str, required=True)
parser.add_argument('--compounds', type=str, required=True)
parser.add_argument('--clusters', type=str, required=True)
parser.add_argument('--folding', type=str, default="output/folding.npy")
parser.add_argument('--nfolds', type=int, default=5)
conf = parser.parse_args()
print(conf)

#Load y matrix and create the mask
df = pd.read_csv(conf.y)
print(df)

#Load the clustering
clusters = np.load(conf.clusters, allow_pickle=True)
print(clusters)

#Load compound list 
cmpd_list = pd.read_csv(conf.compounds, header=None)
cmpd_list["cid"] = range(len(cmpd_list))
print(cmpd_list)

#Create mask
assert len(clusters) == len(cmpd_list[0]), "Number of compound must agree with the size of clustering vector"
Ncmpd = len(set(cmpd_list[0]))
Ncol  = len(df["target_id"].unique()) * len(df["variable"].unique())

# Create dataframe of unique target_ids
target_list = pd.DataFrame(df["target_id"].unique())
target_list["tid"] = range(len(target_list))

# create data frame of unique thresholds
variable_list = pd.DataFrame(df["variable"].unique())
variable_list["vid"] = range(len(variable_list))
Nvar = len(variable_list)

## Add unique target_id identifier to df of chembl_29_thresholds
print("joining tables...")
tmp1 = pd.merge(df,target_list, left_on="target_id", right_on=0)

## Add unique threshold identifier to df of chembl_29 thresholds
tmp2 = pd.merge(tmp1, cmpd_list, left_on="cmpd_id", right_on=1), 

join = pd.merge(tmp2[0], variable_list, left_on="variable", right_on=0) #WHY [0]



#-----------------------------------------------------------------------------------------------------
# Build COOrdiante formatted sparse matrices for Ymask and  Y
#
# Ymask =  coo_matrix((data, (i, j)), [shape=(M, N)])
# to construct from three arrays:
# data[:] the entries of the matrix, in any order
# 
# i[:] the row indices of the matrix entries
# j[:] the column indices of the matrix entries
# 
# Where A[i[k], j[k]] = data[k]. When shape is not specified, it is inferred from the index arrays
#
# Number of rows    : Ncmpd - Unique compound ids 
# Number of columns : Ncol  - Number of unique target ids x number of unqiue thresholds 
#
#
#-----------------------------------------------------------------------------------------------------
I = join["cid"].to_numpy()
J = (Nvar * join["tid"] + join["vid"]).to_numpy()
V = np.ones(len(I))

Ymask = scipy.sparse.coo_matrix((V,(I,J)),(Ncmpd,Ncol))

Y = scipy.sparse.coo_matrix((join["value"].to_numpy(),(I,J)),(Ncmpd, Ncol))
print(f"write Y file to {conf.y[:-3]}")
np.save(conf.y[:-3]+"npy", Y)


print("Calculate folding... nfolds: {conf.nfolds}")
folds = mtcv_clustered(Ymask, pd.Series(clusters), conf.nfolds)

#TODO: save folding
np.save(conf.folding, folds)