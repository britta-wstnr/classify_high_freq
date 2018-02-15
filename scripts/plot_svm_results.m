% Plot results from the SVM classification
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot comparison of test accuracies of across-subjects classification between
% random forest classification and SVM models (FIGURE 4).
%
% AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
% LICENCE: BSD 3-clause
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

project_settings

% file names
svm_lin_results = 'accuracies_linSVM.mat';
svm_rbf_results = 'accuracies_rbfSVM.mat';

%% load data

% svm results
load(fullfile(results_dir, svm_rbf_results)); % or svm_lin_results
testgamma_svm = test * 100;
clear test

% rf results
load(fullfile(results_dir, 'accuracies_rf.mat'));
testgamma_rf = test * 100;
clear test


%% plot SVM and RF accuracies together

% compute chance level
trial_num = 4270;
alpha = 0.05;
classes = 2;
ch_lev = binoinv(1 - alpha, trial_num, 1 / classes) * (100 / trial_num);

% set x-axis spacing
x_rf = 1:4:20 * 4;
x_svm = 3:4:20 * 4;

h = figure; hold all
b = bar(x_rf, testgamma_rf);
c = bar(x_svm, testgamma_svm);
set(b, 'facecolor', [33, 113, 181] / 255, 'edgecolor', [33, 113, 181] / 255);
set(c, 'facecolor', [239, 138, 98] / 255, 'edgecolor', [239, 138, 98] / 255);
set(gca, 'ylim', [45 87])
set(gca, 'xlim', [0 20 * 4 + 30])

line(get(gca, 'xlim'), [mean(testgamma_rf), mean(testgamma_rf)], ...
     'color', [5, 48, 97]/255, 'linewidth', 3); hold on
line(get(gca, 'xlim'), [mean(testgamma_svm), mean(testgamma_svm)], ...
     'color', [244, 109, 67]/255, 'linewidth', 3);
line(get(gca, 'xlim'),[ch_lev, ch_lev], 'color', [0 0 0]/255, ...
     'linewidth', 2.5, 'linestyle', '--');
set(gca, 'YTick', [50 60 70 80 90]);
set(gca, 'XTick', [5 10 15 20] * 4);
set(gca, 'XTickLabel', {'5' '10' '15' '20'})

xlabel('Folds/subjects')
ylabel('Accuracy')

set(h, 'Position', [100, 100, 1049, 895]);
a = gca;
set(a,'box','off','color','none')
b = axes('Position', get(a,'Position'), 'box', 'on', 'xtick', [], 'ytick', []);
axes(a)
linkaxes([a b])
set(h, 'color', [1 1 1])
set(gca,'FontSize', 25)
print(h,'-dpdf', '-bestfit',fullfile(fig_dir, 'gammaaccuracies_svm_rbf.pdf'))
