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
inp_dir = _config.OUT_PLACE + 'b_demultiplex/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')

##
# Functions
##
def demultiplex_stats(nm):
  num_lines = 0
  for fn in os.listdir(inp_dir + nm + '/'):
    if 'R1' not in fn: continue
    lc = util.line_count(inp_dir + nm + '/' + fn)
    if lc % 2 == 1:
      print('Error: fq num lines is odd')
      # import code; code.interact(local=dict(globals(), **locals()))
    num_lines += lc

  # divide by 4 for fastq
  num_reads = num_lines / 4
  print(f'{nm}: {num_reads} reads')
  return

##
# qsub
##
# def gen_qsubs():
#   # Generate qsub shell scripts and commands for easy parallelization
#   print('Generating qsub scripts...')
#   qsubs_dir = _config.QSUBS_DIR + NAME + '/'
#   util.ensure_dir_exists(qsubs_dir)
#   qsub_commands = []

#   num_scripts = 0
#   for idx in range(0, 60):
#     command = 'python %s.py %s' % (NAME, idx)
#     script_id = NAME.split('_')[0]

#     # Write shell scripts
#     sh_fn = qsubs_dir + 'q_%s_%s.sh' % (script_id, idx)
#     with open(sh_fn, 'w') as f:
#       f.write('#!/bin/bash\n%s\n' % (command))
#     num_scripts += 1

#     # Write qsub commands
#     qsub_commands.append('qsub -V -wd %s %s' % (_config.SRC_DIR, sh_fn))

#   # Save commands
#   with open(qsubs_dir + '_commands.txt', 'w') as f:
#     f.write('\n'.join(qsub_commands))

#   print('Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir))
#   return

##
# Main
##
@util.time_dec
def main():
  print(NAME)

  for nm in exp_design['Name']:
    demultiplex_stats(nm)

  demultiplex_stats('other')


  return out_dir


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(split = sys.argv[1])
  else:
    main()