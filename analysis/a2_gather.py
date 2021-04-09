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
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

master_exp_design = pd.read_csv(_config.DATA_DIR + 'master_exp_design.csv')

##
# Main
##
@util.time_dec
def main():
  print(NAME)
  
  '''
    dfs are in same order
  '''

  # Function calls
  for idx, row in master_exp_design.iterrows():
    cond = row['Name']
    dir_nm = row['Directory']
    local_nm = row['Local name']

    print(cond)

    # Simple copy
    command = f'cp {dir_nm}/{local_nm}.csv {out_dir}/{cond}.csv'
    subprocess.check_output(command, shell = True)

  return


if __name__ == '__main__':
  main()