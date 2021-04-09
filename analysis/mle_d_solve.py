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
# Support
##
def combine_dfs(exp_conds, stat_col):
  id_cols = [
    'Peptide.number.within.gene',
    'Protein',
    'Name',
    'Design category',
    'Gene',
  ]

  mdf = None
  for cond in exp_conds:
    df = pd.read_csv(inp_dir + f'{cond}.csv', index_col = 0)

    df[f'{stat_col} fraction {cond}'] = df[f'{stat_col} fraction']
    df[f'Total {cond}'] = df['Total']

    keep_cols = [
      f'{stat_col} fraction {cond}',
      f'Total {cond}',
    ]
    df = df[id_cols + keep_cols]

    if mdf is None:
      mdf = df
    else:
      mdf = mdf.merge(df, on = id_cols, how = 'outer')

  return mdf


'''
  Timeout
  https://stackoverflow.com/questions/2281850/timeout-function-if-it-takes-too-long-to-finish
'''
from functools import wraps
import errno
import os
import signal

class TimeoutError(Exception):
  pass


def timeout(seconds = 10, error_message = os.strerror(errno.ETIME)):
  def decorator(func):
    def _handle_timeout(signum, frame):
      raise TimeoutError(error_message)

    def wrapper(*args, **kwargs):
      signal.signal(signal.SIGALRM, _handle_timeout)
      signal.alarm(seconds)
      try:
        result = func(*args, **kwargs)
      finally:
        # print(f'Timed out!')
        signal.alarm(0)
      return result

    return wraps(func)(wrapper)

  return decorator


##
# Statistics
##
@timeout(seconds = 1)
def mle_betabinom(xs, ns):
  num_reps = len(xs)

  from sympy.solvers import solve, nsolve
  from sympy import Symbol, ln, digamma
  alpha = Symbol('alpha')
  beta = Symbol('beta')
  
  d_ll_wrt_alpha = -num_reps * (digamma(alpha) - digamma(alpha + beta)) + \
                   sum([(digamma(x_i + alpha) - digamma(n_i + alpha + beta)) for x_i, n_i in zip(xs, ns)])
  d_ll_wrt_beta = -num_reps * (digamma(beta) - digamma(alpha + beta)) + \
                   sum([(digamma(n_i - x_i + beta) - digamma(n_i + alpha + beta)) for x_i, n_i in zip(xs, ns)])

  try:
    ans = nsolve(
        (d_ll_wrt_alpha, d_ll_wrt_beta),
        (alpha, beta),
        (5, 5),
    )
    verified = True
  except ValueError:
    print(f'Failed, retrying with no verification')
    print(xs, ns)
    ans = nsolve(
        (d_ll_wrt_alpha, d_ll_wrt_beta),
        (alpha, beta),
        (5, 5),
        verify = False,
    )
    verified = False
  return ans[0], ans[1], {'Sympy verified': verified}


def fast_approx_mle_estimate(xs, ns):
  fs = np.array(xs) / np.array(ns)

  obs_var = np.var(fs)
  exp_var = np.mean(fs) * (1- np.mean(fs)) / np.mean(ns)

  var_factor = obs_var / exp_var
  approx_beta_sum = (np.mean(ns) - 1) / (var_factor - 1) - 1

  approx_mean = np.mean(fs)

  approx_alpha = approx_mean * approx_beta_sum
  approx_beta = (1 - approx_mean) * approx_beta_sum

  return approx_alpha, approx_beta, {}


##
# Primary
##
def solve_mle_beta_parameters(target, stat_col, exp_subset, start_idx, end_idx):
  '''
    Note: Approach still doesn't adjust for read count.

    Annotate with:
    - read counts -> confidence score.
    - fractions -> effect size.

    Raivo Kolde, Sven Laur, Priit Adler, Jaak Vilo, Robust rank aggregation for gene list integration and meta-analysis, Bioinformatics, Volume 28, Issue 4, 15 February 2012, Pages 573–580, https://doi.org/10.1093/bioinformatics/btr709
  '''

  id_cols = [
    'Peptide.number.within.gene',
    'Protein',
    'Name',
    'Design category',
    'Gene',
  ]

  exp_conds = _config.get_subset_conds(exp_subset)
  exp_conds = [s for s in exp_conds if target in s]
  df = combine_dfs(exp_conds, stat_col)

  df = df.iloc[int(start_idx) : int(end_idx)]

  n_cols = sorted([col for col in df.columns if 'Total' in col])
  n_cols = _config.get_subset(n_cols, exp_subset)

  f_cols = sorted([col for col in df.columns if f'{stat_col} fraction' in col])
  f_cols = _config.get_subset(f_cols, exp_subset)

  dd = defaultdict(list)
  timer = util.Timer(total = len(df))
  for idx, row in df.iterrows():
    
    ns, xs = [], []
    for f, n in zip(row[f_cols], row[n_cols]):
      if np.isnan(f):
        continue
      xs.append(round(f * n))
      ns.append(n)

    funcs = {
      # 'MLE': mle_betabinom,
      'MLE fast approx': fast_approx_mle_estimate,
    }

    for func_nm in funcs:
      func = funcs[func_nm]

      try:
        alpha, beta, solution_stats = func(xs, ns)
      except TimeoutError:
        alpha = np.nan
        beta = np.nan
        if func_nm == 'MLE':
          solution_stats = {'Sympy verified': False}
        else:
          solution_stats = {}

      dd[f'{func_nm} alpha'].append(alpha)
      dd[f'{func_nm} beta'].append(beta)
      dd[f'{func_nm} alpha+beta'].append(beta + beta)

      dd[f'{func_nm} mean'].append(alpha / (alpha + beta))
      var_numer = alpha * beta
      var_denom = ((alpha + beta)**2) * (alpha + beta + 1)
      mle_var = var_numer / var_denom
      dd[f'{func_nm} var'].append(mle_var)
      dd[f'{func_nm} std'].append(mle_var**(1/2))

      for sol_stat in solution_stats:
        dd[f'{func_nm} {sol_stat}'].append(solution_stats[sol_stat])

    timer.update()

  for stat in dd:
    df[stat] = dd[stat]

  df.to_csv(out_dir + f'{target}_{stat_col}_{exp_subset}_{start_idx}_{end_idx}.csv')

  return

##
# qsub
##
def gen_qsubs():
  # Generate qsub shell scripts and commands for easy parallelization
  print('Generating qsub scripts...')
  qsubs_dir = _config.QSUBS_DIR + NAME + '/'
  util.ensure_dir_exists(qsubs_dir)
  qsub_commands = []

  num_scripts = 0
  
  targets = _config.targets
  exp_subsets = _config.exp_subsets

  # Function calls
  for target in targets:
    stat_cols = _config.stat_cols

    for stat_col in stat_cols:
      for exp_subset in exp_subsets:

        for start_idx in range(0, 12000, 500):
          end_idx = start_idx + 500

          command = 'python %s.py %s %s %s %s %s' % (NAME, target, stat_col, exp_subset, start_idx, end_idx)
          script_id = NAME.split('_')[0]

          # Write shell scripts
          sh_fn = qsubs_dir + 'q_%s_%s_%s_%s_%s.sh' % (script_id, target, stat_col, exp_subset, start_idx)
          with open(sh_fn, 'w') as f:
            f.write('#!/bin/bash\n%s\n' % (command))
          num_scripts += 1

          # Write qsub commands
          qsub_commands.append('qsub -P regevlab -V -l h_rt=8:00:00 -wd %s %s &' % (_config.SRC_DIR, sh_fn))

  # Save commands
  commands_fn = qsubs_dir + '_commands.sh'
  with open(commands_fn, 'w') as f:
    f.write('\n'.join(qsub_commands))

  subprocess.check_output('chmod +x %s' % (commands_fn), shell = True)

  print('Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir))
  return


##
# Main
##
@util.time_dec
def main(argv):
  print(NAME)
  
  [target, stat_col, exp_subset, start_idx, end_idx] = argv
  print(target, stat_col)
  solve_mle_beta_parameters(target, stat_col, exp_subset, start_idx, end_idx)

  return


if __name__ == '__main__':
  if len(sys.argv) > 2:
    main(sys.argv[1:])
  else:
    gen_qsubs()
