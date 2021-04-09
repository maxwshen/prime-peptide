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
inp_dir = _config.OUT_PLACE + 'mle_d_solve/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'master_exp_design.csv')
conds = exp_design['Name']

##
# Support
##


##
# Main
##
@util.time_dec
def main():
  print(NAME)
  
  targets = _config.targets
  exp_subsets = _config.exp_subsets

  # Function calls
  for target in targets:
    stat_cols = _config.stat_cols

    for stat_col in stat_cols:
      for exp_subset in exp_subsets:

        mdf = None
        for start_idx in range(0, 12000, 500):
          end_idx = start_idx + 500

          try:
            df = pd.read_csv(inp_dir + f'{target}_{stat_col}_{exp_subset}_{start_idx}_{end_idx}.csv', index_col = 0)
          except:
            print(f'Missing {start_idx}')

          if mdf is None:
            mdf = df
          else:
            mdf = mdf.append(df, ignore_index = True)

        mdf.to_csv(out_dir + f'{target}_{stat_col}_{exp_subset}.csv')
        print(f'{target}_{stat_col}_{exp_subset}.csv')

  return


if __name__ == '__main__':
  # if len(sys.argv) > 2:
  #   main(sys.argv[1:])
  # else:
  #   gen_qsubs()
  main()