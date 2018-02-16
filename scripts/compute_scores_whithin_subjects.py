"""Compute classifier outcome for within subjects classification.

Compute the mean accuracy across subjects, and make a (preliminary) plot of the
single subjects' accuracies. Save accuracies as .mat file to visualize them in
MATLAB.

Note: for supplementary material.

AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
LICENCE: BSD 3-clause
"""


import numpy as np
import os.path as op
import matplotlib.pyplot as plt
import scipy.io as io

# import project settings - data is read from and written to results_dir
from py_project_settings import (results_dir, subjects)

accuracy = []
for subj in subjects:
    base_fname_acc = 'singlesub_acc_{}.npy'
    load_fname_acc = base_fname_acc.format(str(subj))
    load_path_acc = op.join(results_dir, load_fname_acc)
    accuracy.append(np.load(load_path_acc))


mean_acc = np.mean(accuracy) * 100

savepath_acc = op.join(results_dir, 'accuracies_single_subjs.mat')
acc_save = {}
acc_save['test'] = accuracy
io.savemat(savepath_acc, acc_save)


###############################################################################
# Plotting accuracies
plt.bar(range(len(accuracy)), accuracy, color="cornflowerblue")
plt.title("Single subject classification, mean accuracy %.2f %%" % mean_acc)
plt.xlabel("subjects")
plt.ylabel("accuracy")
plt.xticks(range(0, len(accuracy), 4), subjects[::4])
plt.show()
