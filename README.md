# Across-subjects classification of stimulus modality from human MEG high frequency activity

Britta U. Westner, Sarang S. Dalal, Simon Hanslmayr, & Tobias Staudigl

## Abstract

Single-trial analyses have the potential to uncover meaningful brain dynamics that are
obscured when averaging across trials. However, low signal-to-noise ratio (SNR) can
impede the use of single-trial analyses and decoding methods. In this study, we
investigate the applicability of a single-trial approach to decode stimulus modality from
magnetoencephalographic (MEG) high frequency activity. In order to classify the
auditory versus visual presentation of words, we combine beamformer source
reconstruction with the random forest classification method. To enable group level
inference, the classification is embedded in an across-subjects framework.
We show that single-trial gamma SNR allows for good classification performance
(accuracy across subjects: 66.44 %). This implies that the characteristics of high
frequency activity have a high consistency across trials and subjects. The random forest
classifier assigned informational value to activity in both auditory and visual cortex
with high spatial specificity. Across time, gamma power was most informative during
stimulus presentation. Among all frequency bands, the 75 Hz to 95 Hz band was the
most informative frequency band in visual as well as in auditory areas. Especially in
visual areas, a broad range of gamma frequencies (55 Hz to 125 Hz) contributed to the
successful classification.
Thus, we demonstrate the feasibility of single-trial approaches for decoding the
stimulus modality across subjects from high frequency activity and describe the
discriminative gamma activity in time, frequency, and space.

## Repository

This repository contains the data analysis scripts for [Westner et al., 2017 (preprint, not peer-reviewed)](https://www.biorxiv.org/content/early/2017/10/12/202424), currently submitted.

The data analysis scripts are organized as follows:

### Configuration files
* `project_settings.m`  settings for all MATLAB files
* `py_project_settings.py`  settings for all PYTHON files

### Data processing
Data processing is completely done in MATLAB, using [FieldTrip](https://github.com/fieldtrip/fieldtrip).
* `get_source_power.m`   compute single trial gamma power on source level
* `get_source_power_singlesubj.m`   compute single trial gamma power on source level for supplementary within subject analysis

### Classification
Decoding analysis is completely done in Python, using [scikit-learn](https://github.com/scikit-learn/scikit-learn).
* `decode_high_freq_RF.py`   classification using random forests
* `decode_high_freq_SVM.py`   classification using SVMs for comparison
* `decode_within_subjects.py`  supplementary analysis within subjects

### Report, evaluation, and plotting
Evaluation of the results is done in Python, plotting is done in MATLAB.
* `compute_scores.py`  classification scores
* `compute_scores_within_subjects.py`   classification scores for supplementary analysis
* `use_fisher.py`  Fisher's exact test on classifier outputs
* `plot_decoding_results.m`  plot results from random forest
* `plot_svm_results.m`   plot results from SVM
* `plot_singlesubj_acc.m`   supplementary figure
* `plot_underlying_activity.m`  plot gamma power in source space


## Data

The data and further information will soon be available at [Open Science Framework](https://osf.io/).

## Dependencies
* FieldTrip
* numpy
* scipy
* scikit-learn
