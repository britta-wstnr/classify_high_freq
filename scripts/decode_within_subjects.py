"""Using Random Forest to classify high frequency activity within subjects.
 """

# imports
import numpy as np
import scipy.io as io
import os.path as op
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier

# path setup:
from py_project_settings import (results_dir, source_dir, datafile_single,
                                 subjects)

# if random state should be controlled:
rseed = None  # rseed = 510

# if random state should be controlled
rseed = None  # rseed = 510
data_list = []
trees = 25000
p_sub = 'sqrt'

# loop over subjs and load the two folds. then fit 2 rfs to those folds-loop.
# then mean over accuracies
for sub in subjects:
    print('Computing Random Forest for subject %s' % sub)

    # defs and data load
    input_path = op.join(source_dir, sub)

    # output:
    savename_acc = op.join(results_dir, 'singlesub_acc_{}')
    savename_vi = op.join(results_dir, 'singlesub_vi_{}')
    savename_confm = op.join(results_dir, 'singlesub_confm_{}')

    data_folds = []
    response_folds = []

    for ii in range(1, 3):

        # load folds
        datafile_single = op.join(input_path, datafile_single.format(str(ii)))
        fold = io.loadmat(datafile_single)
        fold = fold['data_rf']
        response_folds.append(fold[:, -1])
        data_folds.append(fold[:, :-1])

    forest = RandomForestClassifier(n_estimators=trees, max_features=p_sub,
                                    n_jobs=3, oob_score=True,
                                    class_weight='balanced',
                                    random_state=rseed)

    test_score = []
    conf_m = []
    vi = []

    for ii, (train_data, train_response) in enumerate(zip(data_folds,
                                                          response_folds)):
        test_idx = range(2)
        test_idx.pop(ii)
        test_data = data_folds[test_idx[0]]
        test_response = response_folds[test_idx[0]]

        # Probabilities
        probs_pred = (forest.fit(train_data, train_response).predict_proba(
                                                                    test_data))
        # test score
        test_score.append(forest.score(test_data, test_response))

        # confusion matrix
        test_pred = forest.predict(test_data)
        conf_m.append(metrics.confusion_matrix(test_response, test_pred))

        # variable importances
        vi.append(forest.feature_importances_)

    # mean over folds
    test_acc = sum(test_score) / len(test_score)
    vi_perc = sum(vi) / len(vi)
    confmat = sum(conf_m) / len(conf_m)

    print("Test score: %.2f" % test_acc)
    print("Saving subject %s" % sub)

    # save to disk
    np.save(savename_vi.format(sub), vi_perc)
    np.save(savename_acc.format(sub), test_acc)
    np.save(savename_confm.format(sub), confmat)
