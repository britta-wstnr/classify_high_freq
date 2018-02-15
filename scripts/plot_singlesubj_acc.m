% Plot results from within subject classification
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot individual test accuracies for within subject classification (split-
% half cross-validation) and their respective corrected chance levels
% (SUPPLEMENTARY FIGURE 1).
%
% AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
% LICENCE: BSD 3-clause
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


project_settings

%% load data

load(fullfile(results_dir, 'accuracies_single_subjs.mat'));
testgamma = test * 100;

load(fullfile(results_dir, 'trial_num.mat'));


%% chance levels

bonf_alpha = 0.05/20;

for ii=1:length(testgamma)
    chance_levs(ii) = binoinv(1-bonf_alpha, trial_num(ii), 1/2) ...
        * 100/trial_num(ii);
end


%% Plot individual test accuracies

h = figure;
b = bar(testgamma); hold all
set(b, 'facecolor', [33, 113, 181] / 255, 'edgecolor', [33, 113, 181] / 255)
set(gca, 'ylim', [40 93])
set(gca, 'xlim', [0 30])

% plot individual chance levels
for ii=1:length(testgamma)
    line([ii - 0.5, ii + 0.5], [chance_levs(ii), chance_levs(ii)], ...
         'color', [0 0 0]/255, 'linewidth', 3); hold on
end

set(gca, 'YTick', [50 60 70 80 90]);
set(gca, 'XTick', [5 10 15 20]);
xlabel('Folds/subjects')
ylabel('Accuracy')

set(h, 'Position', [100, 100, 1049, 895]);

a = gca;
set(a, 'box', 'off', 'color', 'none')
b = axes('Position', get(a, 'Position'), 'box', 'on', 'xtick', [], ...
         'ytick', []);
axes(a)
linkaxes([a b])

set(h, 'color', [1 1 1])
set(gca,'FontSize', 25)

print(h, '-dpdf', '-bestfit', ...
      fullfile(fig_dir, 'gammaaccuracies_singlesubj.pdf'))
