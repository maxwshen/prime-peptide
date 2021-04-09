import sys, os

SRC_DIR = os.path.dirname(os.path.realpath(__file__)) + '/'
PRJ_DIR = '/'.join(SRC_DIR.split('/')[:-2]) + '/'

#######################################################
# Note: Directories should end in / always
#######################################################
DATA_DIR = PRJ_DIR + 'data/'
OUT_PLACE = PRJ_DIR + 'out/'
RESULTS_PLACE = PRJ_DIR + 'results/'
QSUBS_DIR = PRJ_DIR + 'qsubs/'
#######################################################
#######################################################

# which data are we using? import that data's parameters
# DATA_FOLD = 'rename_me/'
DATA_FOLD = ''

sys.path.insert(0, DATA_DIR + DATA_FOLD)
print('Using data folder:\n', DATA_DIR + DATA_FOLD)
DATA_DIR += DATA_FOLD
OUT_PLACE += DATA_FOLD
RESULTS_PLACE += DATA_FOLD
QSUBS_DIR += DATA_FOLD

#######################################################
# Project-specific parameters
#######################################################

import pandas as pd
stat_cols = ['prime_edit']
exp_design = pd.read_csv(DATA_DIR + 'master_exp_design.csv')
targets = set(exp_design['Target'])
exp_subsets = ['all']

num_permutations = 100

def get_lib_design(cond):
  lib_nm = exp_design[exp_design['Name'] == cond]['Library'].iloc[0]
  df = pd.read_csv(DATA_DIR + f'lib_{lib_nm}_design.csv')
  return df

def get_subset_conds(exp_subset):
  conds = list(exp_design['Name'])
  if exp_subset == 'all':
    return conds
  if exp_subset == 'noB':
    return [c for c in conds if c[-1] != 'B']
  assert(False, 'Not implemented')

def get_subset(cols, exp_subset):
  scs = get_subset_conds(exp_subset)
  ok_cols = []
  for col in cols:
    for sc in scs:
      if sc in col:
        ok_cols.append(col)
        break
  return ok_cols

