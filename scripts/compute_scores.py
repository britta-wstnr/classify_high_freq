"""Compute classifier outcomes.

Compute the mean accuracy across folds, the confusion matrix, and the average
variable importances. Store output as mat files to be able to visualize them
with MATLAB/FieldTrip.

For SVM analysis: change file names.
Note: no variable importance for SVM classification.

AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
LICENCE: BSD 3-clause
"""


import os.path as op
import numpy as np
import scipy.io as io

# import project settings - data is read from and written to results_dir
from py_project_settings import (results_dir, subjects)

###############################################################################
# ACCURACY
acc_fname = 'ACCURACY_rf.npy'
load_path = op.join(results_dir, acc_fname)
test_scores = (np.load(load_path))

# print mean and range:
print(np.mean(test_scores) * 100,
      np.min(test_scores) * 100, np.max(test_scores) * 100)

# save as .mat file
acc_savename = op.join(results_dir, 'accuracies_rf.mat')
accsave = {}
accsave['accuracies'] = test_scores
io.savemat(acc_savename, accsave)

###############################################################################
# CONFMAT
confm_fname = 'CONFM_rf_{}.npy'

confmats = []
confmatssum = np.zeros([2, 2])
for subj in subjects:
    load_path = op.join(results_dir, confm_fname.format(subj))
    dummy = (np.load(load_path))
    confmats.append(dummy)
    confmatssum = np.add(confmatssum, dummy)

###############################################################################
# For Random Forests only:
# VARIABLE IMPORTANCES
vi_fname = 'VI_rf_{}.npy'

vimean = np.zeros([1, 20780])
for subj in subjects:
    load_path = op.join(results_dir, vi_fname.format(subj))
    dummy = (np.load(load_path))
    vimean = np.add(vimean, dummy)

vimean = vimean / (len(subjects))

# save for plotting in Matlab:
vi_savename = op.join(results_dir, 'vi_mean_rf.mat')
visave = {}
visave['vimean'] = vimean
io.savemat(vi_savename, visave)
