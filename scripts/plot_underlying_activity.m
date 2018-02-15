%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot underlying high frequency activity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load project settings from file
project_settings

% load template_grid
load(template_grid_file);

% make dummy structures for plotting
make_plotting_dummies


%% load data

vifile = 'vi_mean_rf.mat'
load(fullfile(results_dir, vifile));
atlas = ft_read_atlas(atlas_path);
mri = ft_read_mri(mri_path);


%% Make own colormap for plotting
% color space from colorbrewer2.org

colors = [103,   0,  31
          178,  24,  43
          214,  96,  77
          244, 165, 130
          253, 219, 199
          247, 247, 247
          209, 229, 240
          146, 197, 222
           67, 147, 195
           33, 102, 172
            5,  48,  97] / 255;
colors = colors(fliplr(1:11), :);
x = [0:27:250, 255];
redmap=  interp1(fliplr(x)/255, colors, linspace(1,0,250));


%% load data of Vp


pow_file = 'pow_gamma.mat'

peak_pos = [68 -20 10    % auditory gamma
            -4 -100 12]  % visual gamma
timepoints = 4;
freqs = [40 60 70 80 90]

% prepare dummy structure
data_dummy.freq = freqs;
dummy_time = [0.125 0.375 0.625 0.875];
data_dummy.dimord = 'rpt_chan_freq_time';

% preallocation:
data_audi_avg = cell(1, length(Vp));
data_vis_avg  = cell(1, length(Vp));
for sub = 1:length(Vp)

    load(fullfile(source_dir, Vp{sub}, pow_file))

    % get auditory and visual trials
    data_audi = data_rf(find(data_rf(:, end) == 0), :);
    data_vis  = data_rf(find(data_rf(:, end) == 1), :);

    % reshape data
    data_audi_resh = reshape(data_audi(:, 1:end - 1), ...
                             length(find(data_rf(:, end) == 0)), 1039, ...
                             length(freqs), timepoints);
    data_vis_resh  = reshape(data_vis(:, 1:end - 1), ...
                             length(find(data_rf(:, end) == 1)), 1039, ...
                             length(freqs), timepoints);

    % auditory trials - get average for subject
    data_dummy.powspctrm = data_audi_resh;
    data_dummy.time = repmat({dummy_time}, size(data_rf,1),1);

    cfg = [];
    cfg.keeptrials = 'no';
    data_audi_avg{sub} = ft_freqdescriptives(cfg, data_dummy);
    data_audi_avg{sub}.time = dummy_time;

    % visual trials - get average for subject
    data_dummy.powspctrm = data_vis_resh;
    data_dummy.time = repmat({dummy_time}, size(data_rf,1),1);

    cfg = [];
    cfg.keeptrials = 'no';
    data_vis_avg{sub} = ft_freqdescriptives(cfg, data_dummy);
    data_vis_avg{sub}.time = dummy_time;
end


%% average across subjects

cfg = [];
data_audi_gamma = ft_freqgrandaverage(cfg, data_audi_avg{:});
data_vis_gamma  = ft_freqgrandaverage(cfg, data_vis_avg{:});

% get difference
data_plot = data_audi_gamma;
data_plot.powspctrm = data_vis_gamma.powspctrm - data_audi_gamma.powspctrm;


%% plotting brain data

% plot 4th frequency band and 2nd time point (or loop, see below)
freq_is = 4; time_is = 2;

for ii = freq_is      % loop over frequencies
    for tt = time_is  % loop over timepoints

        source_dummy.avg.pow(source_dummy.inside,:) = data_plot.powspctrm(:, ii, tt);
        source_dummy.pos = template_grid.pos;

        %% interpolate source onto MRI
        mri.coordsys = 'spm';

        cfg              = [];
        cfg.voxelcoord   ='no';
        cfg.interpmethod = 'nearest';
        cfg.parameter    = 'pow';
        source_interp    = ft_sourceinterpolate(cfg, source_dummy, mri);

        %% plot orthogonal view

        cfg               = [];
        cfg.axes          = 'off'
        cfg.method        = 'ortho';
        cfg.atlas         = atlas;
        cfg.funparameter  = 'pow';
        cfg.maskparameter = 'mask';
        cfg.funcolorlim   = [-0.35 .35];
        cfg.funcolormap   = redmap;
        cfg.opacitylim    = [-.02 .02];
        cfg.location      = peak_pos(1,:);  % change for 2
        cfg.crosshair     = 'yes';
        ft_sourceplot(cfg, source_interp);

        set(gcf, 'Position', [100, 100, 1049, 895]);

    end
end


%% plot TFR
% TODO: all this could be refactored, shared code

% reload dummy
make_plotting_dummies
data_dummy.powspctrm = data_plot.powspctrm;

% find positions in grid
pos_it = dsearchn(template_grid.pos(template_grid.inside,:), ...
                  peak_pos(1,:));

% select plotting data
cfg         = [];
cfg.channel = num2str(pos_it);
data_dummy  = ft_selectdata(cfg, data_dummy);

% custom color limits
my_own_clim = max(abs([min(min(min(data_dummy.powspctrm))), ...
                    max(max(max(data_dummy.powspctrm)))]));

% frequency band limits
freq_bands = [25 45; 55 75; 75 95; 105 125; 125 145];

% manipulate data_dummy
data_dummy.freq = mean(freq_bands, 2);
data_dummy.time = [0.125 0.375 0.625 0.875];
data_dummy.freq_bands = freq_bands;
dat = squeeze(data_dummy.powspctrm);

% find the gaps in the frequency bands (50 Hz and 100 Hz omits)
for ii = 2:size(data_dummy.freq_bands, 1)
    differ(ii-1) = data_dummy.freq_bands(ii, 1) - ...
                   data_dummy.freq_bands(ii - 1, 2);
end
idx = find(differ~=0);

% take care of gaps in frequency dimension
if(size(idx,2)~=0)      % if there are gaps in data
    idx = fliplr(idx);  % flip the index to avoid shifting of rows
    for ii = 1:size(idx, 2)
        % fill the gap in the freq_bands definition
        gapfreq = [data_dummy.freq_bands(idx(ii), 2), ...
                   data_dummy.freq_bands(idx(ii) + 1, 1)];
        data_dummy.freqbands = [data_dummy.freq_bands(1:idx(ii), :); ...
                                gapfreq; ...
                                data_dummy.freq_bands(idx(ii) + 1:end, :)];

        % fill in NaNs  at appropriate position
        timedim = size(dat, 2);
        nandummy = NaN(1, timedim);
        dat = [dat(1:idx(ii), :); nandummy; dat(idx(ii) + 1:end,:)];
    end
end

data_dummy.freq = mean(data_dummy.freq_bands, 2);
data_dummy.powspctrm = [data_dummy.powspctrm, ...
                        data_dummy.powspctrm(1, 1:2, :)];
data_dummy.powspctrm(1, :, :) = dat;

% double the rows with real data
idx = [1 1 2 3 3 4 4 5 6 6 7 7];
dat = data_dummy.powspctrm(:,idx,:);
dat = squeeze(dat);

% add one row and one column to make pcolor behave
dat = [dat; dat(end, :)];
dat = [dat, dat(:, end)];

h=figure
pcolor(dat)
colormap(redmap)
shading flat
caxis([-my_own_clim my_own_clim])
set(gca, 'YTick', [1 3 4 6 8 9 11 13]);
set(gca, 'YTickLabel', {[25 45 55 75 95 105 125 145]});
set(gca, 'XTick', [1:5]);
set(gca, 'XTickLabel', {[0 250 500 750 1000]});
set(gcf, 'color', [1 1 1])
set(gca,'FontSize',18)

ylabel('Frequency (Hz)')
xlabel('Time (ms)')

caxis([-0.35 0.35])
colorbar

set(h, 'Position', [100, 100, 1049, 895]);

print(h,'-dpdf', '-bestfit', fullfile(fig_dir, 'underlyingvisual.pdf'))
