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
inp_dir = _config.OUT_PLACE + 'a2_gather/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'master_exp_design.csv')
conds = exp_design['Name']

##
# Functions
##
def normalize(cond):
  df = pd.read_csv(inp_dir + f'{cond}.csv', index_col = 0)
  df['Name'] = df.index

  lib_design = _config.get_lib_design(cond)
  df = df.merge(lib_design, how = 'outer', on = 'Name')

  stat_cols = _config.stat_cols
  df['Total'] = df[stat_cols + ['Wild-type']].apply(np.nansum, axis = 'columns')

  from scipy.special import logit
  for stat_col in stat_cols:
    df[f'{stat_col} fraction'] = df[stat_col] / df['Total']
    df[f'{stat_col} logit fraction'] = logit(df[stat_col] / df['Total'])

  df = df.replace([np.inf, -np.inf], np.nan)

  # Form control df
  ctrl_df = df[df['Design category'] == '1_wt_control_peptides']
  ctrl_df = ctrl_df[ctrl_df['Total'] >= 20]
  print('Num. control rows:', len(ctrl_df))

  for stat_col in stat_cols:
    sclf = f'{stat_col} logit fraction'
    mean = np.mean(ctrl_df[sclf])
    ctrl_df[sclf] -= mean
    df[sclf] -= mean
    print(stat_col, mean)

  nc_df = df[df['Design category'] != '1_wt_control_peptides']
  nc_df = nc_df[nc_df['Total'] >= 20]

  nc_df.to_csv(out_dir + f'{cond}.csv')
  ctrl_df.to_csv(out_dir + f'{cond}_control.csv')

  return


##
# Main
##
@util.time_dec
def main():
  print(NAME)
  
  # Function calls
  for cond in conds:
    print(cond)
    normalize(cond)

  return


if __name__ == '__main__':
  main()