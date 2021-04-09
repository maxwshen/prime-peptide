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
inp_dir = _config.OUT_PLACE + 'mle_e_unify_inferences/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'master_exp_design.csv')
conds = exp_design['Name']

hparams = {
  # 'n': 50000,
  'n': 1000,
}

##
# Statistics
##
def get_mc_control_effects(cdf):
  n_per_row = int(hparams['n'] // len(cdf)) + 1

  samples = []
  for idx, row in cdf.iterrows():
    alpha_col = 'Posterior alpha' 
    beta_col = 'Posterior beta'
        
    try:
      samples += list(np.random.beta(row[alpha_col], row[beta_col], size = n_per_row))
    except:
      continue

  samples = np.array(samples)
  return samples[:hparams['n']]


##
# Primary
##
def effect_tests(target, stat_col, exp_subset):
  tdf = pd.read_csv(inp_dir + f'{target}_{stat_col}_{exp_subset}.csv', index_col = 0)
  cdf = pd.read_csv(inp_dir + f'{target}_{stat_col}_{exp_subset}_control.csv', index_col = 0)

  mc_control_effects = get_mc_control_effects(cdf)
  mean_control = np.mean(mc_control_effects)
  std_control = np.std(mc_control_effects)

  tdf = tdf.dropna(subset = ['Posterior mean'])

  dd = defaultdict(list)
  timer = util.Timer(total = len(tdf))
  for idx, row in tdf.iterrows():
    alpha, beta = row['Posterior alpha'], row['Posterior beta']

    samples = np.random.beta(alpha, beta, size = hparams['n'])

    dd['p(effect > control)'].append(sum(samples > mc_control_effects) / len(samples))
    dd['p(effect < control)'].append(sum(samples < mc_control_effects) / len(samples))

    timer.update()

  for col in dd:
    tdf[col] = dd[col]

  tdf.to_csv(out_dir + f'{target}_{stat_col}_{exp_subset}.csv')

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

        command = 'python %s.py %s %s %s' % (NAME, target, stat_col, exp_subset)
        script_id = NAME.split('_')[0]

        # Write shell scripts
        sh_fn = qsubs_dir + 'q_%s_%s_%s_%s.sh' % (script_id, target, stat_col, exp_subset)
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
  
  [target, stat_col, exp_subset] = argv
  print(target, stat_col)
  effect_tests(target, stat_col, exp_subset)

  return


if __name__ == '__main__':
  if len(sys.argv) > 2:
    main(sys.argv[1:])
  else:
    gen_qsubs()
