# Dataset filters:
# 1) Select all:
#      assay.assay_organism="Homo sapiens", 
#      target_dictionary.targe_type="SINGLE PROTEIN" 
#      proteins that have at least N standard_type="IC50" measurements (N=100 or 200)
#TODO: check: Rattus norvegicus, Mus musculus
#TODO: check these: ED50, AC50, EC50, MIC, Ki, Activity, Inhibition, GI50, Potency (later the higher number)
# 2) Filter according to the standard_units="nM" #TODO: ug.mL-1
# 3) Pick the minimum IC50 for all cells (take care of missing, NA, ...)
# 4) Filter credible values: 10^9 > IC50 >= 10^-5  #TODO: tighten #3<->4
# 5) pIC50 = 9 - log10(IC50)
# 6) Refilter proteins that have at least N compounds
# 7) Remove empty rows
# OUT: matrix market data, compound name, protein name lists

import configargparse
import sqlite3
import pandas as pd
import numpy as np
import scipy.io
import os
import logging


p = configargparse.ArgParser(default_config_files=["default.ini"])
p.add('-c', '--config', required=False, is_config_file=True, help='Config file path')
p.add('--sqlite'      , required=True, type=str, help="ChEMBL sqlite database")
#p.add("--organism"   , required=True, help="Organisms for protein filtering" )
#p.add("--targettype" , required=True, help="Target type for protein filtering")
p.add('--mincmpdcount', required=True, help='Minimal number of compounds required for an assays', type=int)
p.add('--thresholds'  , required=True, help="Thresholds for classification", type=float, action="append")
p.add('--datadir'     , required=True, help="Data directory to write to (append prefix)", type=str)
p.add('--prefix'      , required=True, help="Prefix for the current dataset", type=str)
options = p.parse_args()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

outdir = options.datadir + "/" + options.prefix

if not os.path.exists(outdir):
    logging.info("Building output dir: '%s'" % outdir)
    os.makedirs(outdir)

logging.info("Querying sqlite database '%s'" % options.sqlite)

conn = sqlite3.connect(options.sqlite)
# Homo sapiens, Single Protein, IC50, nM
df = pd.read_sql_query("""SELECT target_dictionary.chembl_id as target_id, 
                                 molecule_dictionary.chembl_id as cmpd_id, 
                                 activities.standard_units as stu,
                                 CASE activities.standard_units
                                    WHEN 'nM' THEN activities.standard_value
                                    WHEN 'ug.mL-1' THEN activities.standard_value / compound_properties.full_mwt * 1E6
                                 END ic50,
                                 CASE activities.standard_relation 
                                    WHEN '<'  THEN '<'
                                    WHEN '<=' THEN '<'
                                    WHEN '='  THEN '='
                                    WHEN '>'  THEN '>'
                                    WHEN '>='  THEN '>' 
                                    ELSE 'drop' 
                                 END relation
                                 FROM molecule_dictionary 
                            JOIN activities on activities.molregno == molecule_dictionary.molregno 
                            JOIN assays on assays.assay_id == activities.assay_id 
                            JOIN target_dictionary on target_dictionary.tid == assays.tid
                            JOIN compound_properties on compound_properties.molregno = molecule_dictionary.molregno
                           WHERE 
                                target_dictionary.organism='Homo sapiens' AND 
                                target_dictionary.target_type='SINGLE PROTEIN' AND
                                activities.standard_type = 'IC50' AND 
                                activities.standard_units IN  ('nM','ug.mL-1') AND
                                activities.standard_relation IN ('<', '<=', '=','>', '>=')  AND
                                ic50 < 10e9 AND 
                                ic50 >= 10e-5 
                                ORDER BY target_id, cmpd_id """, conn)
conn.close()
df.to_csv('Step10/1_sql_output.csv')

#import ipdb; ipdb.set_trace()
logging.info("Step10/ Filtering and thresholding activity data")

# Pick the minimum
logging.info("Step10/ Grouping by target_id [cmpd_id].min() ")
df = df.groupby(["target_id","cmpd_id"]).min().reset_index()
df.to_csv('2_groupby_min.csv')

# at least N compounds per assay
logging.info("Step10/Grouping by target_id [cmpd_id].nunique() ")
c  = df.groupby("target_id")["cmpd_id"].nunique()
c.to_csv('Step10/3_groupby_nunique.csv')

# extract indexes for compounds with N compounds per assay
logging.info("Step10/ Filtering by min compound count : '%f'" % options.mincmpdcount)
i  = c[c >= options.mincmpdcount].index
df = df[df.target_id.isin(i)]
df.to_csv("Step10/4_target_id_is_in.csv")

# 
df["pic50"] = 9 - np.log10(df["ic50"])
df.to_csv('Step10/5_temporary_continuous.csv')

# df = pd.read_csv("temporary_continuous.csv")


#Thresholding
value_vars = []
for thr in options.thresholds:
    logging.info("Step10/ Processing threshold : '%1.1f'" % thr)
    value_vars.append("%1.1f" % thr)
    thr_str = "%1.1f" % thr
    ## using +1 and -1 for actives and inactives
    df[thr_str] = (df["pic50"] >= thr) * 2.0 - 1.0
    df[thr_str] = np.where(np.logical_and((df["relation"] == '<'), (df['pic50'] < thr)), np.nan, df[thr_str]) 
    df[thr_str] = np.where(np.logical_and((df["relation"] == '>'), (df['pic50'] > thr)), np.nan, df[thr_str]) 

df.to_csv("Step10/6_df.csv")

logging.info("Step10/ Unpivot target_id and cmpd_id ")
# Unpivot df from wide to long format
melted = pd.melt(df, id_vars=['target_id','cmpd_id'], value_vars=value_vars).dropna()

# Write threshold file
melted.to_csv('%s/%s_thresh.csv' % (outdir, options.prefix), index = False)

# Write unique compound IDs
np.savetxt("%s/%s_compounds.csv" % (outdir, options.prefix), melted["cmpd_id"].unique(), fmt="%s")

# Write unique target  IDs
np.savetxt("%s/%s_targets.csv"   % (outdir, options.prefix), melted["target_id"].unique(), fmt="%s")


