"""Using Random Forest to classify high frequency activity across subjects.
"""

# imports
import os
import random
import datetime
import numpy as np
import scipy.io as io
import os.path as op
from sklearn import svm, metrics
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

# project settings:
from py_project_settings import (results_dir, source_dir, datafile,
                                 subjects)


# file names for saving:
save_name_acc = op.join(results_dir, 'ACCURACY_linSVM')

# if random state should be controlled:
rseed = None  # rseed = 510

data_list = []
for sub in subjects:

    # defs and data load
    inputpath = op.join(source_dir, sub)
    os.chdir(inputpath)

    data = io.loadmat(datafile)
    data = data['data_rf']

    class_a = np.where(data[:, data.shape[1] - 1] == 1)
    class_b = np.where(data[:, data.shape[1] - 1] == 0)

    len_class_a = len(class_a[0])
    len_class_b = len(class_b[0])

    # downsample smaller class
    if len_class_a > len_class_b:

        # random.seed(rseed)
        indx = random.sample(class_a[0], len_class_b)
        indx.extend(class_b[0])
        indx = sorted(indx)
        data = data[indx, :]

    elif len_class_b > len_class_a:

        # random.seed(rseed)
        indx = random.sample(class_b[0], len_class_a)
        indx.extend(class_a[0])
        indx = sorted(indx)
        data = data[indx, :]

    data_list.append(data)


# SVM decoding
test_score = []
for sub_n, sub in enumerate(subjects):
    # time keeper
    t1 = datetime.datetime.now()
    print('starting with subject %i at %i:%i.' % (sub_n, t1.hour, t1.minute))

    # set everything up to save data
    save_name_confm = op.join(results_dir, 'CONFM_linSVM_' + sub)

    # cross val:
    train_cross = np.concatenate(data_list[:sub_n] + data_list[sub_n + 1:])
    test_cross = data_list[sub_n]

    train_features = train_cross[:, :train_cross.shape[1] - 1]
    test_features = test_cross[:, :test_cross.shape[1] - 1]

    # responses:
    train_response = train_cross[:, train_cross.shape[1] - 1]
    train_response = train_response.astype(int)
    test_response = test_cross[:, test_cross.shape[1] - 1]
    test_response = test_response.astype(int)

    # SVM:
    # default is RBF kernel, specify "linear" for linear kernel
    clf = make_pipeline(StandardScaler(), svm.SVC(kernel="linear",
                                                  probability=True,
                                                  class_weight='balanced'))
    probs_predicted = (clf.fit(train_features,
                       train_response).predict_proba(test_features))

    test_score.append(clf.score(test_features, test_response))

    # confusion matrix
    test_predicted = (clf.predict(test_features))
    confm = metrics.confusion_matrix(test_response,
                                     test_predicted)  # true, predicted

    print('Saving subject %i.' % sub_n)

    # save to disk
    np.save(save_name_confm, confm)

# save accuracy list of all subjects:
np.save(save_name_acc, test_score)
