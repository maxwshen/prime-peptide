# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess
sys.path.append('/home/unix/maxwshen/')
import numpy as np
from collections import defaultdict
import pandas as pd
from mylib import util
import pickle

# Default params
inp_dir = _config.OUT_PLACE + 'cx_overall_quant/'
NAME = util.get_fn(__file__)
out_place = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_place)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')
target_design = pd.read_csv(_config.DATA_DIR + 'lib_targets_design.csv')
lib_design = None
target = None


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

    for jdx in range(0, 12000, 500):
      start_jdx = jdx
      end_jdx = jdx + 500
      command = 'python %s.py %s %s %s' % (NAME, condition, start_jdx, end_jdx)
      script_id = NAME.split('_')[0]

      out_fn = out_place + f'stats_{condition}_{start_jdx}_{end_jdx}.csv'
      if os.path.isfile(out_fn):
        continue

      # Write shell scripts
      sh_fn = qsubs_dir + 'q_%s_%s_%s.sh' % (script_id, condition, start_jdx)
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
  # Parse condition-specific settings

  for nm in exp_design['Name']:
    stats_mdf = None
    data_mdf = None
    print(nm)
    for start in range(0, 12000, 500):
      print(start)
      end = start + 500

      try:
        stats_df = pd.read_csv(inp_dir + f'stats_{nm}_{start}_{end}.csv', index_col = 0)
        data_df = pd.read_csv(inp_dir + f'data_{nm}_{start}_{end}.csv', index_col = 0)
      except:
        print(f'Did not find files for {nm}, {start}, {end}')
        continue

      if stats_mdf is None:
        stats_mdf = stats_df
      else:
        stats_mdf = stats_mdf.merge(stats_df, how = 'outer', left_index = True, right_index = True)

      if data_mdf is None:
        data_mdf = data_df
      else:
        data_mdf = data_mdf.merge(data_df, how = 'outer', left_index = True, right_index = True)

    # import code; code.interact(local=dict(globals(), **locals()))
    stats_mdf['Count'] = stats_mdf.apply(np.nansum, axis = 'columns')
    stats_mdf = stats_mdf[['Count']]
    stats_mdf.to_csv(out_place + f'stats_{nm}.csv')

    data_mdf['Count'] = data_mdf.apply(np.nansum, axis = 'columns')
    data_mdf = data_mdf[['Count']]
    data_mdf.to_csv(out_place + f'data_{nm}.csv')

    data_mdfs = data_mdf[data_mdf['Count'] >= 5]
    data_mdfs.to_csv(out_place + f'datafilt_{nm}.csv')

  return


if __name__ == '__main__':
  # if len(sys.argv) > 1:
  #   main(sys.argv[1:])
  # else:
  #   gen_qsubs()
  main([])