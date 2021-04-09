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
inp_place = _config.OUT_PLACE + 'c_alignment/'
NAME = util.get_fn(__file__)
out_place = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_place)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')
target_design = pd.read_csv(_config.DATA_DIR + 'lib_targets_design.csv')
lib_design = None
target = None
target_row = None

##
# Basic sequence operations 
##
def q_to_score(q):
  return ord(q) - 33

def calc_matches(s1, s2):
  assert len(s1) == len(s2)
  nm = 0
  for c1, c2 in zip(s1, s2):
    if c1 == c2 and c2 != '-' and c1 != '-':
      nm += 1
  return nm, nm / len(s1)


def check_mismatch_quality_support(read, ref, qs):
  threshold = 30
  for c, r, q in zip(read, ref, qs):
    if c != r and r != '-' and c != '-':
      if q_to_score(q) < threshold:
        return False
  return True

def check_indel_quality_support(read, ref, qs):
  threshold = 30
  radii = 3
  events = []
  for c, r, q in zip(read, ref, qs):
    if r == '-' or c == '-':
      events.append('indel')
    else:
      events.append(q_to_score(q))

  for idx, e in enumerate(events):
    left = max(0, idx - radii)
    right = min(len(events), idx + radii + 1)
    if 'indel' in events[left : right]:
      if e != 'indel' and e < threshold:
        return False
  return True


def count_indels(s1, s2):
  num_dels, num_ins = 0, 0
  ts1, ts2 = trim_start_end_dashes(s1), trim_start_end_dashes(s2)
  num_dels = len(ts1.replace('-', ' ').split()) - 1
  num_ins = len(ts2.replace('-', ' ').split()) - 1
  return num_dels, num_ins


def trim_start_end_dashes(seq):
  alphabet = set(list(seq))
  if '-' in alphabet:
    alphabet.remove('-')
  alphabet = list(alphabet)
  start = min([seq.index(s) for s in alphabet])
  end = min([seq[::-1].index(s) for s in alphabet])
  if end > 0:
    return seq[start:-end]
  else:
    return seq[start:]

##
# Classification
##
def categorize_alignment(read, ref, qs):
  '''

  '''
  outcomes = {
    'Prime edit': target_row['PE edit search term'],
    'Wild-type': target_row['Wild-type search term'],
    'Hairpin RT': target_row['Hairpin RT search term'],
  }
  for nm in outcomes:
    seq = outcomes[nm]
    if seq in read:
      return nm

  cutsite = target_row['Cutsite index']
  pidx = 0
  trunc_read = ''
  for idx, (c, r) in enumerate(zip(read, ref)):
    if r != '-':
      pidx += 1
    if pidx == cutsite:
      break
  trunc_read = read[idx:]
  trunc_ref = ref[idx:]
  nd, ni = count_indels(trunc_read, trunc_ref)
  if nd + ni > 0:
    return 'indel'
  # import code; code.interact(local=dict(globals(), **locals()))
  return 'other'


###
# Main function
###
def count_data(inp_dir, peptide_nm, stats):
  '''
    Count alignment for a single peptide
  '''
  data = defaultdict(lambda: 0)
  target_prefix = target[:5]
  match_rate_threshold = 0.70
  num_match_threshold = 20

  for split in os.listdir(inp_dir):
    if split == 'aligns': continue
    inp_fn = inp_dir + '%s/%s.txt' % (split, peptide_nm)
    with open(inp_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 0:
          header = line.strip()
        if i % 4 == 1:
          read = line.strip()
        if i % 4 == 2:
          ref = line.strip()
        if i % 4 == 3:
          qs = line.strip()

          # Process alignment
          stats['Total'] += 1

          # Trim random initial section
          if target_prefix not in read[:14]:
            stats['1. Target prefix not found early in read'] += 1
            continue
          tidx = read[:14].index(target_prefix)
          read = read[tidx:]
          ref = ref[tidx:]
          qs = qs[tidx:]

          # Filter low match rate
          num_matches, match_rate = calc_matches(read, ref)
          if num_matches < num_match_threshold:
            stats[f'2a. Num. matches below {num_match_threshold}'] += 1
            continue
          if match_rate < match_rate_threshold:
            stats[f'2b. Match rate below {match_rate_threshold}'] += 1
            continue

          # Check quality support for mismatches and indels
          if not check_mismatch_quality_support(read, ref, qs):
            stats[f'3a. Mismatches had insufficient quality support'] += 1
            continue

          if not check_indel_quality_support(read, ref, qs):
            stats[f'3b. Indels had insufficient quality support'] += 1
            continue

          # Accept
          stats['Accepted'] += 1

          # Categorize type of alignment
          cat = categorize_alignment(read, ref, qs)
          data[cat] += 1


  data_df = pd.DataFrame(data, index = [peptide_nm])
  return data_df


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

    step = 100
    for jdx in range(0, 12000, step):
      start_jdx = jdx
      end_jdx = jdx + step
      command = 'python %s.py %s %s %s' % (NAME, condition, start_jdx, end_jdx)
      script_id = NAME.split('_')[0]

      # out_fn = out_place + f'stats_{condition}_{start_jdx}_{end_jdx}.csv'
      # if os.path.isfile(out_fn):
      #   continue

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
  nm = argv[0]
  start = int(argv[1])
  end = int(argv[2])
  inp_dir = inp_place + nm + '/'

  print(NAME)
  print(nm)

  # Parse condition-specific settings
  exp_row = exp_design[exp_design['Name'] == nm].iloc[0]
  lib_nm = exp_row['Library']
  target_nm = exp_row['Target']

  # Library design
  global lib_design
  lib_design = pd.read_csv(_config.DATA_DIR + f'lib_{lib_nm}_design.csv')
  peptide_nms = list(lib_design['Name'])
  peptide_nms = peptide_nms[start : end]

  # Target
  global target_row
  global target
  target_row = target_design[target_design['Target'] == target_nm].iloc[0]
  target = target_row['Sequence']

  stats = defaultdict(lambda: 0)
  timer = util.Timer(total = len(peptide_nms))
  mdf = pd.DataFrame()
  for peptide_nm in peptide_nms:
    df = count_data(inp_dir, peptide_nm, stats)
    mdf = mdf.append(df)
    timer.update()

  mdf = mdf.fillna(value = 0)
  mdf.to_csv(out_place + f'data_{nm}_{start}_{end}.csv')

  # Save, using start/end
  stats_df = pd.DataFrame(stats, index = ['Value']).T
  stats_df.to_csv(out_place + f'stats_{nm}_{start}_{end}.csv')

  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1:])
  else:
    gen_qsubs()