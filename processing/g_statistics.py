# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, copy, pickle
import numpy as np
from collections import defaultdict
sys.path.append('/home/unix/maxwshen/')
from mylib import util, compbio
import pandas as pd

# Default params
inp_dir = _config.OUT_PLACE + 'f_summary_stats/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')
target_design = pd.read_csv(_config.DATA_DIR + 'lib_targets_design.csv')
lib_design = pd.read_csv(_config.DATA_DIR + 'lib_12k_peptides_design.csv')

##
# Primary
##
def stats_tests(cond):

  df = pd.read_csv(inp_dir + f'{cond}.csv', index_col = 0)

  # Add columns
  ins_cols = ['ins_1nt_precise', 'ins_1nt_imprecise', 'ins_2ntplus']
  df['ins'] = df[ins_cols].apply(np.nansum, axis = 'columns')

  if 'ssODN' in cond:
    stat_cols = ['hdr', 'del_MH', 'del_MHless', 'ins']
  else:
    stat_cols = ['del_MH', 'del_MHless', 'ins']

  for stat_col in stat_cols:
    df[f'{stat_col} fraction'] = df[stat_col] / df['Edited total']
    df[f'{stat_col} efficiency'] = df[stat_col] / df['Total reads']

  # Form control df
  ctrl_df = df[df['Design category'] == '1_wt_control_peptides']
  ctrl_df = ctrl_df[~pd.isna(ctrl_df['Edited fraction'])]
  ctrl_df = ctrl_df[ctrl_df['Edited total'] >= 100]
  print('Num. control rows:', len(ctrl_df))

  '''
    Form control editing frequencies for chi-squared.
  '''
  mean_edit_count = np.mean(ctrl_df['Edited total'])

  dd = defaultdict(lambda: 0)
  for idx, row in ctrl_df.iterrows():
    for stat_col in stat_cols:
      dd[stat_col] += row[f'{stat_col} fraction'] * mean_edit_count
  tot = sum(dd.values())
  for stat in dd:
    dd[stat] = dd[stat] / tot
  p_exp = np.array([dd[s] for s in stat_cols])

  # Perform chi-squared
  from scipy.stats import chisquare
  print('Running chi-squared tests...')
  nc_df = df[df['Design category'] != '1_wt_control_peptides']

  dd = defaultdict(list)
  timer = util.Timer(total = len(nc_df))
  for idx, row in nc_df.iterrows():
    f_obs = [row[f'{stat_col}'] for stat_col in stat_cols]
    f_exp = p_exp * row['Edited total']
    chisq, pval = chisquare(f_obs, f_exp = f_exp)
    dd['Chi-squared statistic'].append(chisq)
    dd['Chi-squared pval'].append(pval)
    timer.update()
  for col in dd:
    nc_df[col] = dd[col]

  # Chi-squared pval to FDR
  print('Calculating FDR ...')
  nc_df = nc_df.sort_values(by = 'Chi-squared pval')
  num_tests = len(nc_df['Chi-squared pval'].dropna())
  assert num_tests > 0, 'NO GOOD NON-CONTROL ROWS'

  dd = defaultdict(list)
  for fdr_threshold in [0.01, 0.05, 0.1]:
    fdr_accept = True
    for idx, row in nc_df.iterrows():
      pval = row['Chi-squared pval']
      k_over_m = idx / num_tests
      if pval > k_over_m * fdr_threshold:
        fdr_accept = False
      dd[f'Chi-squared FDR {fdr_threshold}'].append(fdr_accept)
    print(fdr_threshold, sum(dd[f'Chi-squared FDR {fdr_threshold}']))
  for col in dd:
    nc_df[col] = dd[col]

  '''
    Efficiency.

  '''
  from scipy.stats import binom_test
  print('Running binomial tests...')
  ctrl_effs = {stat: np.mean(ctrl_df[f'{stat} efficiency']) for stat in stat_cols}

  dd = defaultdict(list)
  timer = util.Timer(total = len(nc_df))
  for idx, row in nc_df.iterrows():
    for stat in stat_cols:
      try:
        pval = binom_test(row[stat], row['Total reads'], p = ctrl_effs[stat])
        dd[f'Binomial pval {stat}'].append(pval)
      except ValueError:
        dd[f'Binomial pval {stat}'].append(np.nan) 
    timer.update()
  for col in dd:
    nc_df[col] = dd[col]

  # Binomial pval to FDR
  for stat_col in stat_cols:
    dd = defaultdict(list)
    pval_col = f'Binomial pval {stat_col}'
    nc_df = nc_df.sort_values(by = pval_col)
    num_tests = len(nc_df[pval_col].dropna())
    for fdr_threshold in [0.01, 0.05, 0.1]:
      fdr_accept = True
      for idx, row in nc_df.iterrows():
        pval = row[pval_col]
        k_over_m = idx / num_tests
        if pval > k_over_m * fdr_threshold:
          fdr_accept = False
        dd[f'Binomial pval {stat_col} FDR {fdr_threshold}'].append(fdr_accept)
      print(stat_col, fdr_threshold, sum(dd[f'Binomial pval {stat_col} FDR {fdr_threshold}']))
    for col in dd:
      nc_df[col] = dd[col]

  nc_df = nc_df.sort_values(by = 'Unnamed: 0.1')
  nc_df.to_csv(out_dir + f'{cond}.csv')

  ctrl_df.to_csv(out_dir + f'{cond}_control.csv')
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
  for condition in exp_design['Name']:
    exp_row = exp_design[exp_design['Name'] == condition].iloc[0]
    lib_nm = exp_row['Library']

    command = f'python {NAME}.py {condition}'
    script_id = NAME.split('_')[0]

    # Write shell scripts
    sh_fn = qsubs_dir + f'q_{script_id}_{condition}.sh'
    with open(sh_fn, 'w') as f:
      f.write('#!/bin/bash\n%s\n' % (command))
    num_scripts += 1

    # Write qsub commands
    qsub_commands.append('qsub -P regevlab -V -l h_rt=4:00:00,h_vmem=2G -wd %s %s &' % (_config.SRC_DIR, sh_fn))

  # Save commands
  commands_fn = qsubs_dir + '_commands.sh'
  with open(commands_fn, 'w') as f:
    f.write('\n'.join(qsub_commands))

  subprocess.check_output('chmod +x %s' % (commands_fn), shell = True)

  print('Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir))
  return

@util.time_dec
def main(argv):

  cond = argv[0]

  stats_tests(cond)

  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1:])
  else:
    gen_qsubs()