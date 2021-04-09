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
inp_dir_1 = '/ahg/regevdata/projects/CRISPR-libraries/prj/rich-peptide-pe-2003xx/out/d2_combine/'
# inp_dir_0409 = '/ahg/regevdata/projects/CRISPR-libraries/prj/rich-peptide-200409/out/f_summary_stats/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')
conds = exp_design['Name']

##
# Main
##
@util.time_dec
def main():
  print(NAME)
  
  '''
    dfs are in same order
  '''

  stat_cols = [
    'prime_edit',
    'Total',
  ]

  # Function calls
  for cond in conds:
    print(cond)
    df1 = pd.read_csv(inp_dir_0311 + f'data_{cond}.csv', index_col = 0)
    df2 = pd.read_csv(inp_dir_0319 + f'data_{cond}.csv', index_col = 0)

    if 'pre-gRNA' in cond:
      cond = cond.replace('pre-gRNA', 'pregRNA')
    df3 = pd.read_csv(inp_dir_0409 + f'{cond}.csv', index_col = 0)

    assert sum(df1['Peptide name'] == df2['Peptide name']) == len(df1), 'dfs are not aligned by peptide name'
    assert sum(df1['Peptide name'] == df3['Peptide name']) == len(df3), 'dfs are not aligned by peptide name'

    mdf = df1[stat_cols] + df2[stat_cols] + df3[stat_cols]
    id_cols = [col for col in df1.columns if col not in stat_cols]
    dfs = pd.concat([mdf, df1[id_cols]], axis = 1, sort = False)

    dfs['Edited fraction'] = dfs['Edited total'] / dfs['Total reads']

    dfs.to_csv(out_dir + f'{cond}.csv')


  return


if __name__ == '__main__':
  main()