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
inp_dir = _config.OUT_PLACE + 'a_split/'
NAME = util.get_fn(__file__)
out_place = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_place)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')
target_design = pd.read_csv(_config.DATA_DIR + 'lib_targets_design.csv')
lib_design = pd.read_csv(_config.DATA_DIR + f'lib_12k_peptides_design.csv')

match_len = 8
prefixes = [s[:match_len] for s in lib_design['Sequence']]
peptide_nms = list(lib_design['Name'])
prefix_to_peptide = {prefix: nm for prefix, nm in zip(prefixes, peptide_nms)}
prefixes_set = set(prefixes)

##
# Support
##
def find_peptide_nm(peptide_read):
  '''
    First 8-mer of Sequence in lib_262_peptides.csv is unique
  '''
  # constant_region = 'CCACCGGTACCATGGCTAGC'
  constant_region = 'ACCATGGCTAGC'
  # constant_region = 'ACAAGTGCTGGT'
  if constant_region not in peptide_read:
    return 'missing_constant'

  peptide_seq_start_idx = peptide_read.index(constant_region) + len(constant_region)

  obs_prefix = peptide_read[peptide_seq_start_idx : peptide_seq_start_idx + match_len]

  if obs_prefix in prefixes_set:
    return prefix_to_peptide[obs_prefix]
  else:
    return 'unmatched'

##
# Alignments
##
def alignment(read, target):
  seq_align_tool = '/ahg/regevdata/projects/CRISPR-libraries/tools/seq-align/bin/needleman_wunsch'
  
  align = subprocess.check_output(seq_align_tool + ' --match 1 --mismatch -1 --gapopen -5 --gapextend -0 --freestartgap ' + read + ' ' + target, shell = True)

  align = align.decode('utf-8')
  if align[-2:] == '\n\n':
    align = align[:-2]

  return align

##
# IO
##
def store_alignment(alignment_buffer, peptide_nm, align_header, align, read_q):
  # Place quality scores in same configuration as read alignment
  read_align = align.split()[0]
  rq = ''
  curr_idx = 0
  for ch in read_align:
    if ch == '-':
      rq += ch
    else:
      rq += read_q[curr_idx]
      curr_idx += 1

  align_string = f'{align_header}\n{align}\n{rq}\n'
  if align_string.count('\n') != 4:
    import code; code.interact(local=dict(globals(), **locals()))
  alignment_buffer[peptide_nm].append(align_string)
  return

def init_alignment_buffer():
  alignment_buffer = defaultdict(list)
  return alignment_buffer

def flush_alignments(alignment_buffer, out_dir):
  print(f'Flushing... \n{datetime.datetime.now()}')
  for peptide_nm in alignment_buffer:
    with open(out_dir + f'{peptide_nm}.txt', 'a') as f:
      for align in alignment_buffer[peptide_nm]:
        f.write(align)
  print(f'Done flushing.\n{datetime.datetime.now()}')
  return

def prepare_outfns(out_dir, peptide_nms):
  for peptide_nm in list(peptide_nms) + ['unmatched', 'missing_constant']:
    out_fn = out_dir + f'{peptide_nm}.txt'
    util.exists_empty_fn(out_fn)
  return

##
# Primary
##
def matchmaker(nm, split):
  print(split)
  stdout_fn = _config.SRC_DIR + f'nh_c_{nm}_{split}.out'
  util.exists_empty_fn(stdout_fn)
  out_dir = f'{out_place}{nm}/{split}/'
  util.ensure_dir_exists(out_dir)

  # Parse condition-specific settings
  exp_row = exp_design[exp_design['Name'] == nm].iloc[0]
  parent_fn = exp_row['Parent file']
  lib_nm = exp_row['Library']
  target_nm = exp_row['Target']

  # Library design
  global lib_design

  # Target 
  target_row = target_design[target_design['Target'] == target_nm].iloc[0]
  target = target_row['Sequence']
  target_strand = target_row['gRNA orientation']

  read1_fn = inp_dir + f'{nm}_R1_{split}.fq'
  read2_fn = inp_dir + f'{nm}_R2_{split}.fq'

  alignment_buffer = init_alignment_buffer()
  prepare_outfns(out_dir, peptide_nms)

  quality_pass = 0
  tot_reads = 0

  tot_lines = util.line_count(read1_fn)
  timer = util.Timer(total = tot_lines)
  with open(read1_fn) as f1, open(read2_fn) as f2:
    for i, (line1, line2) in enumerate(zip(f1, f2)):
      if i % 4 == 0:
        h1 = line1.strip()
        h2 = line2.strip()
      if i % 4 == 1:
        peptide_read = line1.strip()
        peptide_nm = find_peptide_nm(peptide_read)

        target_read = line2.strip()
        if target_strand == '-':
          target_read = compbio.reverse_complement(target_read)

      if i % 4 == 3:
        tot_reads += 1
        q1, q2 = line1.strip(), line2.strip()
        peptide_q = q1
        read_q = q2

        qs = [ord(s)-33 for s in q1 + q2]
        if np.mean(qs) >= 30:
          quality_pass += 1

          # Run alignment and store in buffer
          align_header = f'>1_{peptide_read}'
          align = alignment(target_read, target)
          store_alignment(alignment_buffer, peptide_nm, align_header, align, read_q)

      if i % int(tot_lines / 100) == 1 and i > 1:
        # Flush alignment buffer
        flush_alignments(alignment_buffer, out_dir)
        alignment_buffer = init_alignment_buffer()

        # Stats for the curious
        with open(stdout_fn, 'a') as outf:
          outf.write(f'Time: {datetime.datetime.now()}\n')
          outf.write(f'Progress: {i / int(tot_lines / 100)}\n')
          outf.write(f'Frac. quality passing: {(quality_pass / tot_reads)}\n')

      timer.update()
  
  # Final flush
  flush_alignments(alignment_buffer, out_dir)

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
  
  for nm in exp_design['Name']:
    for idx in range(0, _config.num_splits):

      command = 'python %s.py %s %s' % (NAME, nm, idx)
      script_id = NAME.split('_')[0]

      # Write shell scripts
      sh_fn = qsubs_dir + 'q_%s_%s_%s.sh' % (script_id, nm, idx)
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

@util.time_dec
def main(argv):
  print(NAME)

  nm = argv[0]
  split = int(argv[1])

  # Function calls
  matchmaker(nm, split) 
  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1:])
  else:
    gen_qsubs()