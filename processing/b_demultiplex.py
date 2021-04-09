# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, imp
sys.path.append('/home/unix/maxwshen/')
import numpy as np
from collections import defaultdict
from mylib import util, compbio
import pandas as pd

# Default params
inp_dir = _config.OUT_PLACE + 'a_split/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')


##
# Functions
##
def match(index, nms_to_index1, nms_to_index2):
  for nm in nms_to_index1:
    index1 = index.split('+')[0]
    index2 = index.split('+')[1]
    index1_score = sum([1 for i in range(len(index1)) if index1[i] != nms_to_index1[nm][i]])
    index2_score = sum([1 for i in range(len(index2)) if index2[i] != nms_to_index2[nm][i]])
    if index1_score <= 1 and index2_score <= 1:
      # print('Expected')
      return nm
  return 'other'

##
# primary
##
def demultiplex(parent_fn, read_type, split):
  r1_fn = inp_dir + f'{parent_fn}_{read_type}_{split}.fq'

  eds = exp_design[exp_design['Parent file'] == parent_fn]
  nms = eds['Name']
  index1s = [s.upper() for s in eds['Index 1']]
  index2s = [s.upper() for s in eds['Index 2']]
  nms_to_index1 = {nm: index for nm, index in zip(nms, index1s)}
  nms_to_index2 = {nm: index for nm, index in zip(nms, index2s)}

  for name in list(eds['Name']) + ['other']:
    util.ensure_dir_exists(out_dir + name)
    util.exists_empty_fn(out_dir + name + f'/{read_type}_{split}.fq')

  lc = util.line_count(r1_fn)
  num_bad_q, num_tot = 0, 0
  timer = util.Timer(total = lc)
  with open(r1_fn) as f:
    for i, line in enumerate(f):
      if i % 4 == 0:
        header = line.strip()
        index = header.split(':')[-1]
        pass
      if i % 4 == 1:
        read = line.strip()
      if i % 4 == 3:
        qs = line.strip()

        num_tot += 1

        demultiplex_id = match(index, nms_to_index1, nms_to_index2)
        
        new_header = header[1:]

        out_fn = out_dir +  f'{demultiplex_id}/{read_type}_{split}.fq'
        with open(out_fn, 'a') as f:
          f.write('@%s\n%s\n+\n%s\n' % (new_header, read, qs))

      timer.update()

  print('Rejected %s fraction of reads' % (num_bad_q / num_tot))

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

  parent_fns = set(exp_design['Parent file'])

  for parent_fn in parent_fns:
    for read_type in ['R1', 'R2']:
      for split_num in range(0, 60):
        command = f'python {NAME}.py {parent_fn} {read_type} {split_num}'
        script_id = NAME.split('_')[0]

        # Write shell scripts
        sh_fn = qsubs_dir + f'q_{script_id}_{parent_fn}_{read_type}_{split_num}.sh'
        with open(sh_fn, 'w') as f:
          f.write('#!/bin/bash\n%s\n' % (command))
        num_scripts += 1

        # Write qsub commands
        qsub_commands.append('qsub -P regevlab -V -l h_rt=20:00:00,h_vmem=2G -wd %s %s &' % (_config.SRC_DIR, sh_fn))

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

  [parent_fn, read_type, split] = argv

  demultiplex(parent_fn, read_type, split)
  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1:])
  else:
    gen_qsubs()
