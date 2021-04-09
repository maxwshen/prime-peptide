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
inp_dir_control = _config.OUT_PLACE + 'mle_d_solve_control/'
inp_dir_treat = _config.OUT_PLACE + 'mle_d2_combine_solve/'

NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'master_exp_design.csv')
conds = exp_design['Name']

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
        for inp_dir in [inp_dir_control, inp_dir_treat]:
          df = pd.read_csv(inp_dir + f'{target}_{stat_col}_{exp_subset}.csv')

          '''
            Remove peptides that are seen in too few conditions
          '''
          scs = _config.get_subset(df.columns, exp_subset)
          df['Num replicates with data'] = df[scs].count(axis = 'columns')
          df = df[df['Num replicates with data'] > 2]

          '''
            Unify MLE alpha/beta columns from sympy and fast approx.

            Use sympy if verified and in bounds and not too different from fast approx
            Otherwise, use fast approx
          '''

          keep_cols = [col for col in df.columns if 'MLE' not in col]

          inf_cols = [
            'alpha',
            'beta',
            'alpha+beta',
            'mean',
            'std',
            'var',
          ]
          sources = []
          dd = defaultdict(list)
          for idx, row in df.iterrows():
            source = 'MLE fast approx'
            if 'MLE Sympy verified' not in row.index or row['MLE Sympy verified'] == False:
              source = 'MLE fast approx'
            else:
              if row['MLE alpha'] > 0 and row['MLE beta'] > 0 and row['MLE mean'] <= 1:
                if abs(row['MLE mean'] - row['MLE fast approx mean']) < 0.05:
                  source = 'MLE'

            dd['Inference source'].append(source)
            for inf_col in inf_cols:
              dd[f'Posterior {inf_col}'].append(row[f'{source} {inf_col}'])

          df = df[keep_cols]
          for col in dd:
            df[col] = dd[col]

          '''
            Handle out of bounds variables: When variance is negative, set it to a low number
          '''
          crit = (df['Posterior var'] < 0)
          df.loc[crit, 'Posterior alpha+beta'] = 5000
          df.loc[crit, 'Posterior alpha'] = 5000 * df[crit]['Posterior mean']
          df.loc[crit, 'Posterior beta'] = 5000 * (1 - df[crit]['Posterior mean'])
          var_numer = (df[crit]['Posterior alpha'] * df[crit]['Posterior beta'])
          var_denom1 = (df[crit]['Posterior alpha'] + df[crit]['Posterior beta'])**2
          var_denom2 = df[crit]['Posterior alpha'] + df[crit]['Posterior beta'] + 1
          df.loc[crit, 'Posterior var'] = var_numer / (var_denom1 * var_denom2)
          df.loc[crit, 'Posterior std'] = np.sqrt(df[crit]['Posterior var'])

          # Save
          print(f'{target}_{stat_col}_{exp_subset}')
          if 'control' in inp_dir:
            df.to_csv(out_dir + f'{target}_{stat_col}_{exp_subset}_control.csv')
          else:
            df.to_csv(out_dir + f'{target}_{stat_col}_{exp_subset}.csv')

  return


if __name__ == '__main__':
  # if len(sys.argv) > 2:
  #   main(sys.argv[1:])
  # else:
  #   gen_qsubs()
  main()