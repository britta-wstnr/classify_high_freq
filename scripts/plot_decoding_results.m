% Plot results from the random forest classification
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot confusion matrices and accuracies across folds (FIGURE 1),  variable
% importances overlayed on brain (FIGURE 2), and time-frequency plots of
% specific voxels (FIGURE 3).
%
% AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
% LICENCE: BSD 3-clause
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load project settings from file
project_settings

% load template grid
load(template_grid_file)

% make dummy structures for plotting
make_plotting_dummies


%% load data

vifile = 'vi_mean_rf.mat';
load(fullfile(results_dir, vifile));
atlas = ft_read_atlas(atlas_path);
mri = ft_read_mri(mri_path);


%% Make own colormaps for plotting
% color spaces from colorbrewer2.org

% for plotting of variable importances:
T = [
    255, 255, 217
    237, 248, 177
    199, 233, 180
    127, 205, 187
    65,  182, 196
    29,  145, 192
    34,   94, 168
    ] / 255;
x = [20:42:255, 255];
yellow_blue_map =  interp1(fliplr(x) / 255, T, linspace(1, 0.08, 250));
yellow_blue_map = flipud(yellow_blue_map);

% for plotting of confusion matrix
T = [
    247, 251, 255
    222, 235, 247
    198, 219, 239
    158, 202, 225
    107, 174, 214
    66,  146, 198
    33,  113, 181
    8,    81, 156
    8,    48, 107
    ] / 255;
x = [0:32:255, 255];
blue_map =  interp1(fliplr(x) / 255, T, linspace(1, 0.08, 250));
blue_map = flipud(yellow_blue_map);


%% reshape vimean

freqs = [35 65 85 115 135];
time_points = 4;
data_resh = reshape(vimean, 1039, length(freqs), time_points);


%% get cuttoff values and maskind

cutoff_idx = round(0.98 * length(vimean));
masking_cutoff = sort(vimean);
masking_cutoff = masking_cutoff(cutoff_idx);

freq_ind = 1:length(freqs);
time_cell = {'0-250', '250-500', '500-750', '750-1000'};


%% plot variable importances on brain

for ii = 1:length(freqs)

    for tt = 1:time_points

        source_dummy.avg.pow(source_dummy.inside,:) = data_resh(:, ii, tt);
        mri.coordsys = 'spm';

        %% interpolate data onto MRI
        cfg = [];
        cfg.voxelcoord ='no';
        cfg.interpmethod = 'nearest'; % no smoothing
        cfg.parameter = 'pow';
        source_intp = ft_sourceinterpolate(cfg, source_dummy, mri);


        %% plot it

        % masking
        source_intp.mask = source_intp.pow > masking_cutoff;

        % only plot if something to plot:
        if(sum(sum(sum(source_intp.mask)))>=1)

            cfg               = [];
            cfg.axes          = 'off'
            cfg.method        = 'ortho';
            cfg.funparameter  = 'pow';
            cfg.maskparameter = 'mask';
            cfg.funcolorlim   = [0 max(vimean)];
            cfg.funcolormap   = yellow_blue_map;
            cfg.atlas         = atlas;
            cfg.crosshair     = 'yes';
            ft_sourceplot(cfg, source_intp);

            % title
            suptitle(sprintf('Frequency: %d Hz, %s ms', freqs(ii), ...
                             time_cell{tt}));
            % enlarge figure
            set(gcf, 'Position', [100, 100, 1049, 895]);

        end
    end
end


%% plot TFR in specific voxels
% TODO: all this could be refactored, shared code

% reload dummy
make_plotting_dummies
data_dummy.powspctrm = data_resh;

pos = [
    72 -18 12  % audi gamma
    -4 -100 12   % visual gamma
    ];
pos_it = dsearchn(template_grid.pos(template_grid.inside,:), pos(1,:));

% get data of voxel
cfg = [];
cfg.channel = num2str(pos_it);
data_dummy = ft_selectdata(cfg, data_dummy);

% prepare data structure
freqbands = [25 45; 55 75; 75 95; 105 125; 125 145];
data_dummy.freq = mean(freqbands, 2)';
data_dummy.time = [0.125 0.375 0.625 0.875];
data_dummy.freqbands = freqbands;

dat = squeeze(data_dummy.powspctrm);

% find the gaps in the frequency bands (50 Hz and 100 Hz omits)
for ii = 2:size(data_dummy.freqbands, 1)
    differ(ii-1) = data_dummy.freqbands(ii, 1) - ...
                   data_dummy.freqbands(ii - 1, 2);
end
idx = find(differ~=0);

% take care of gaps in frequency dimension
if(size(idx, 2) ~= 0)      % if there are gaps in data
    idx = fliplr(idx);     % flip the index to avoid shifting of rows
    for ii = 1:size(idx, 2)
        % fill the gap in the freqbands definition
        gapfreq = [data_dummy.freqbands(idx(ii), 2), ...
                   data_dummy.freqbands(idx(ii) + 1, 1)];
        data_dummy.freqbands = [data_dummy.freqbands(1:idx(ii), :); ...
                                gapfreq; ...
                                data_dummy.freqbands(idx(ii) + 1:end, :)];

        % fill in NaNs at appropriate positions
        timedim = size(dat, 2);
        nandummy = NaN(1, timedim);
        dat = [dat(1:idx(ii), :); nandummy; dat(idx(ii) + 1:end, :)];
    end
end

data_dummy.freq = mean(data_dummy.freqbands, 2)';
data_dummy.powspctrm = [data_dummy.powspctrm, data_dummy.powspctrm(1,1:2,:)];
data_dummy.powspctrm(1, :, :) = dat;

data_dummy.mask =  data_dummy.powspctrm >  masking_cutoff;
data_dummy.mask = logical(data_dummy.mask);

% plotting with Fieldtrip is possible but does not yield desired result:
% data_dummy.freq(1) = 40;
% h=figure;
% cfg = [];
% cfg.colorbar = 'no';
% % cfg.maskstyle = 'outline';
% % cfg.maskparameter = 'mask';
% % cfg.maskalpha = 0.6;
% cfg.zlim = [0 max(vimean)];
% cfg.colormap = yellow_blue_map; %colormap(summer);
% ft_singleplotTFR(cfg, data_dummy);

% plot home-brew version instead:
% double the rows with real data for prettier result
idx = [1 1 2 3 3 4 4 5 6 6 7 7];
dat = data_dummy.powspctrm(:, idx, :);
dat = squeeze(dat);

% add one row and one column to make pcolor behave
dat = [dat; dat(end, :)];
dat = [dat, dat(:, end)];

h = figure
pcolor(dat)
colormap(yellow_blue_map)
shading flat
caxis([0 max(vimean)])
set(gca, 'YTick', [1 3 4 6 8 9 11 13]);
set(gca, 'YTickLabel', {[25 45 55 75 95 105 125 145]});
set(gca, 'XTick', [1:5]);
set(gca, 'XTickLabel', {[0 250 500 750 1000]});
set(gcf, 'color', [1 1 1])
set(gca,'FontSize',18)

ylabel('Frequency (Hz)')
xlabel('Time (ms)')

colorbar

set(h, 'Position', [100, 100, 1049, 895]);
print(h,'-dpdf', '-bestfit', fullfile(fig_dir, 'gamma_vis.pdf'))


%% load accuracy estimates

load(fullfile(results_dir, 'accuracies_rf.mat'));
testgamma = test_scores * 100; % convert to percent


%% plot accuracy estimates

% compute chance level
trial_num = 4270;
alpha = 0.05;
classes = 2;
ch_lev = binoinv(1 - alpha, trial_num, 1 / classes) * (100 / trial_num);

% plot folds
h = figure;
b = bar(testgamma); hold all
set(b, 'facecolor', [33, 113, 181] / 255, 'edgecolor', [33, 113, 181] / 255)
set(gca, 'ylim', [45 87])
set(gca, 'xlim', [0 30])

line(get(gca, 'xlim'), [mean(testgamma), mean(testgamma)], ...
     'color', [8, 48, 107] / 255, 'linewidth', 3); hold on
line(get(gca, 'xlim'), [ch_lev, ch_lev], 'color', [0 0 0] / 255, ...
     'linewidth', 2.5, 'linestyle', '--');
set(gca, 'YTick', [50 60 70 80 90]);
set(gca, 'XTick', [5 10 15 20]);

xlabel('Folds/subjects')
ylabel('Accuracy')

% make figure bearable
set(h, 'Position', [100, 100, 1049, 895]);
a = gca;
set(a, 'box', 'off', 'color', 'none')
b = axes('Position', get(a,'Position'), 'box', 'on', 'xtick', [], 'ytick', []);
axes(a)
linkaxes([a b])
set(h, 'color', [1 1 1])
set(gca, 'FontSize', 25)
print(h,'-dpdf', '-bestfit',fullfile(fig_dir, 'gammaaccuracies.pdf'))


%% confusion matrix plot

% trial numbers gamma
gamma = [1486, 649; 786, 1349];

subplot(1, 2, 2)
b = imagesc(gamma / (trial_num / 2));
colormap(blue_map);
caxis([0 1]);
colormap;
set(gcf, 'color', [1 1 1]); colorbar;
