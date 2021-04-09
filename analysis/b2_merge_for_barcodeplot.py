# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess
sys.path.append('/home/unix/maxwshen/')
import numpy as np
from collections import defaultdict
from mylib import util
import pandas as pd

# Default params
inp_dir = _config.OUT_PLACE + 'b_normalize/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'master_exp_design.csv')
conds = exp_design['Name']

##
# Functions
##
def run_tests(target):
  id_cols = [
    'Peptide.number.within.gene',
    'Protein',
    'Name',
    'Design category',
    'Gene',
  ]

  conds = exp_design[exp_design['Target'] == target]['Name']

  tmdf = None
  cmdf = None
  for cond in conds:
    tdf = pd.read_csv(inp_dir + f'{cond}.csv', index_col = 0)
    cdf = pd.read_csv(inp_dir + f'{cond}_control.csv', index_col = 0)

    stat_cols = [col for col in tdf.columns if 'logit fraction' in col]
    stat_cols += [col.replace(' logit', '') for col in stat_cols]

    tdf = tdf[id_cols + stat_cols]
    tdf = tdf.rename(columns = {col : f'{col}, {cond}' for col in stat_cols})
    if tmdf is None:
      tmdf = tdf
    else:
      tmdf = tmdf.merge(tdf, on = id_cols, how = 'outer')

    cdf = cdf[id_cols + stat_cols]
    cdf = cdf.rename(columns = {col : f'{col}, {cond}' for col in stat_cols})
    if cmdf is None:
      cmdf = cdf
    else:
      cmdf = cmdf.merge(cdf, on = id_cols, how = 'outer')

  cmdf.to_csv(out_dir + f'{target}_control_merged.csv')
  tmdf.to_csv(out_dir + f'{target}_merged.csv')

  return


##
# Main
##
@util.time_dec
def main():
  print(NAME)
  
  targets = [
    'LibBwithssODN',
    'LibC',
  ]

  # Function calls
  for target in targets:
    print(target)
    run_tests(target)

  return


if __name__ == '__main__':
  main()