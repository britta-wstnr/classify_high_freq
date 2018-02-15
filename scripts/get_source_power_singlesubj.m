%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute gamma power within subjects in source space %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

project_settings

% load the custom template grid for source analysis
load(template_grid_file)

% ultimate data length
trial_num = zeros(20, 8);


%% LOOP over subjects

for nn = 1:length(Vp)

    keep Vp nn template_grid toimin toimax trial_num

    % reload project settings
    project_settings

    % data lenghts
    covlat = [-0.5 1.0];

    virtlat = [-.5 1.5];   % this is incl. zero padding
    virtlatbl = [-1 0.5];

    zeropadlats = [0 1];
    zeropadlatsbl = [-0.5 0];


    % TFR specs
    timewindows = {[0 0.25], [0.25 0.5], [0.5 0.75], [0.75 1]};

    headmodelpath = fullfile(headmodel_dir, Vp{nn});
    loadpath = fullfile(raw_dir, Vp{nn});
    datasavepath = fullfile(source_dir, Vp{nn});

    % ensure that path to save exists
    if(exist(fullfile(datasavepath), 'dir') == 0)
        mkdir(datasavepath)
    end


    %% load and concatenate data

    % auditory stimulation
    files={'matchica_audi_shifted'; 'mismatchica_audi_shifted'}

    load(fullfile(loadpath, files{1}))
    load(fullfile(loadpath, files{2}))

    data_audi=ft_appenddata([], matchica_audi, mismatchica_audi)

    % visual stimulation
    files={'matchica_vis'; 'mismatchica_vis'}

    load(fullfile(loadpath, files{1}))
    load(fullfile(loadpath, files{2}))

    data_vis=ft_appenddata([], matchica_vis, mismatchica_vis)

    % correct for stimulation delay in data
    % 23 ms delay for vis stims
    cfg=[];
    cfg.offset = length(nearest(data_vis.time{1}, 0): ...
                        nearest(data_vis.time{1}, 0.023)).*(-1)
    cfg.minlength = 'maxperlen';
    data_vis = ft_redefinetrial(cfg, data_vis)

    % append both, audi and vis
    data = ft_appenddata([], data_audi, data_vis);

    % cut data to equal trial length
    cfg = [];
    cfg.latency = [-1.45 3.95];
    data = ft_selectdata(cfg, data)

    % keep workspace clean
    clear matchica_vis mismatchica_vis
    clear matchica_audi mismatchica_audi


    %% audi vs vis
    trials_audi = find(data.trialinfo(:, 5) == 0);
    trials_vis = find(data.trialinfo(:, 5) == 1);


    %% prepare the headmodel

    if(exist(fullfile(headmodelpath, 'hdm.mat'), 'file') == 0)

        % load stuff for headmodel
        load(fullfile(headmodelpath, 'segmentedmri.mat'));

        cfg        = [];
        cfg.method = 'singleshell';
        hdm        = ft_prepare_headmodel(cfg, segmentedmri);
        hdm = ft_convert_units(hdm, 'mm');

        save(fullfile(headmodelpath, 'hdm.mat'), 'hdm')

    else

        load(fullfile(headmodelpath, 'hdm.mat'))
        % make sure it's in the  right units
        hdm = ft_convert_units(hdm, 'mm');
    end


    %%  prepare sourcemodel

    if(exist(fullfile(headmodelpath, 'sourcemodel_warped.mat'), 'file') == 0)

        load(fullfile(headmodelpath, 'mriRhs.mat'));

        cfg                = [];
        cfg.grid.warpmni   = 'yes';
        cfg.grid.template  = template_grid;
        cfg.grid.nonlinear = 'yes'; % use non-linear normalization
        cfg.mri            = mriRhs;
        sourcemodel        = ft_prepare_sourcemodel(cfg);   % in cm

        save sourcemodel_warped sourcemodel

    else
        load(fullfile(headmodelpath, 'sourcemodel_warped.mat'))
        sourcemodel = ft_convert_units(sourcemodel, 'mm');
    end


    %% preprocessing data

    cfg = [];
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 35;
    cfg.hpfilttype = 'fir';
    datahp = ft_preprocessing(cfg, data);

    %% divide manually in 2 folds here before cov mat and BF

    % downsample the more prevalent class

    if length(trials_audi) > length(trials_vis)
        % if more auditory than visual trials: downsample auditory
        trials_audi = randsample(trials_audi, length(trials_vis), false)
        trials_audi = sort(trials_audi);
    elseif length(trials_vis) > length(trials_audi)
        % if more visual than auditory trials: downsample visual
        trials_vis = randsample(trials_vis, length(trials_audi), false)
        trials_vis = sort(trials_vis);
    end

    % create two folds:
    all_trials = sort([trials_audi; trials_vis]);
    trial_num(nn) = length(all_trials);

    fold1 = sort(randsample(all_trials, length(all_trials)/2, false));
    fold2 = setdiff(all_trials, fold1);
    folds = {fold1, fold2};

    trialinfo = datahp.trialinfo(:, 5);
    trialinfo = [[1:length(trialinfo)]', trialinfo];

    for ff = 1:2

        foldtrials = folds{ff};
        trialinfo_thisfold = trialinfo(foldtrials,:);

        % keep track of visual vs auditory stimulation
        vis_idx = find(trialinfo_thisfold(:, 2) == 1);
        audi_idx = find(trialinfo_thisfold(:, 2) == 0);

        % get the trial numbers
        trials_vis_fold = trialinfo_thisfold(vis_idx,1);
        trials_audi_fold = trialinfo_thisfold(audi_idx,1);


        %% covariance matrices

        cfg = [];
        cfg.trials = foldtrials;
        cfg.covariance = 'yes';
        cfg.covariancewindow = covlat;
        cov = ft_timelockanalysis(cfg, datahp);

        if(0) % optional saving
            covsave = cov.cov;

            covname = ['covmat_singlesubj_fold', num2str(ff)];
            save(fullfile(datasavepath, covname), 'covsave');
        end


        %% beamformer


        cfg                           = [];
        cfg.method                    = 'lcmv';
        cfg.grid                      = sourcemodel;
        cfg.vol                       = hdm;
        cfg.reducerank                = 'no';
        cfg.(cfg.method).keepfilter   = 'yes';
        cfg.(cfg.method).fixedori     = 'yes';
        cfg.(cfg.method).projectnoise = 'yes';
        cfg.(cfg.method).lambda       = '5%';
        cfg.(cfg.method).weightnorm   = 'nai';
        source                        = ft_sourceanalysis(cfg, cov);


        %% cut data for virtual electrodes

        cfg             = [];
        cfg.trials      = foldtrials;
        cfg.latency     = virtlat;
        data_act        = ft_selectdata(cfg, data);

        cfg              = [];
        cfg.trials       = foldtrials;
        cfg.latency      = virtlatbl;
        data_bl          = ft_selectdata(cfg, data);


        %% get indices for zero-padding

        [tp_act] = dsearchn(data_act.time{1}', zeropadlats');
        [tp_bl]  = dsearchn(data_bl.time{1}', zeropadlatsbl');


        %% virtual electrodes: create structure for further analyses

        spatialfilter = cat(1, source.avg.filter{:});

        virtsens = [];
        virtsensbl = [];
        for ii=1:length(foldtrials)

            virtsens.trial{ii}   = spatialfilter * data_act.trial{ii};
            virtsensbl.trial{ii} = spatialfilter * data_bl.trial{ii};

            % zero padding
            virtsens.trial{ii}(:, 1:tp_act(1)) = 0;
            virtsens.trial{ii}(:, tp_act(2):end) = 0;
            virtsensbl.trial{ii}(:, 1:tp_bl(1)) = 0;
            virtsensbl.trial{ii}(:, tp_bl(2):end) = 0;

        end

        virtsens.time = data_act.time;
        virtsens.fsample = data_act.fsample;
        virtsens.trialinfo = data_act.trialinfo;

        virtsensbl.time = data_bl.time;
        virtsensbl.fsample = data_bl.fsample;
        virtsensbl.trialinfo = data_bl.trialinfo;

        % fake labels:
        for i = 1:size(spatialfilter, 1)
            virtsens.label{i} = [num2str(i)];
        end
        virtsensbl.label = virtsens.label;


        %%   GAMMA POWER

        foilims_gamma = [35 65 85 115 135];
        freqsmooth = 10;

        for ii = 1:length(timewindows)
            cfg = [];
            cfg.latency =  timewindows{ii};
            datacut = ft_selectdata(cfg, virtsens);

            cfg = [];
            cfg.method = 'mtmfft';
            cfg.foi = foilims_gamma;
            cfg.taper = 'dpss';
            cfg.tapsmofrq = freqsmooth;
            cfg.keeptrials = 'yes';
            data_tf_g = ft_freqanalysis(cfg, datacut);
            pow_g{ii} = data_tf_g.powspctrm;

        end

        baselinetimewindows = {[-0.5 -0.25], [-0.25 -0.005]}
        for ii = 1:2
            cfg = [];
            cfg.latency =  baselinetimewindows{ii};
            datacut = ft_selectdata(cfg, virtsensbl);

            cfg = [];
            cfg.method = 'mtmfft';
            cfg.foi = foilims_gamma;
            cfg.taper = 'dpss';
            cfg.tapsmofrq = freqsmooth
            cfg.keeptrials = 'yes';
            data_tf_bl_g = ft_freqanalysis(cfg, datacut);
            pow_bl{ii} = data_tf_bl_g.powspctrm;
        end

        %%  put all that together

        dummy = zeros(size(data_tf_g.trialinfo,1), 1039, ...
                      length(data_tf_g.freq), ...
                      [length(baselinetimewindows) + length(timewindows)]);

        % baseline:
        dummy(:,:,:,1)   = pow_bl{1};    % BL1
        dummy(:,:,:,2)   = pow_bl{2};    % BL2

        % activation:
        dummy(:,:,:,3)   = pow_g{1};
        dummy(:,:,:,4)   = pow_g{2};
        dummy(:,:,:,5)   = pow_g{3};
        dummy(:,:,:,6)   = pow_g{4};

        %make structure
        data_tf_g.powspctrm = dummy;
        data_tf_g.freq = [data_tf_g.freq];
        data_tf_g.time = [-.375 -.125 .125 .375 .625 .875];
        data_tf_g.dimord = 'rpt_chan_freq_time';

        clear dummy


        %% Baseline correction

        cfg = [];
        cfg.baseline = [-0.5 -0.01];
        cfg.baselinetype = 'relchange';
        cfg.parameter = 'powspctrm';
        data_corr_bl = ft_freqbaseline(cfg, data_tf_g);


        %% reshape to matrix

        data_rf = reshape(data_corr_bl.powspctrm(:, :, :, 3:end), ...
                          size(data_tf_g.trialinfo, 1), ...
                          [length(data_corr_bl.label) * ...
                           length(data_corr_bl.freq) * ...
                           (length(data_corr_bl.time) - 2)]);


        %% disentangle trials:

        if(intersect(trials_audi_fold, trials_vis_fold));
            error('Something went wrong with setting up the folds.')
        end

        response = data_tf_g.trialinfo(:,5);

        data_rf = [data_rf, response];

        display(sprintf('Saving subject %i, fold %i', nn, ii));
        save_fname = ['pow_gamma_singlesubj_fold', num2str(ff)];
        save(fullfile(datasavepath, save_fname), 'data_rf');
    end
end

save(fullfile(results_dir, 'trial_num.mat'), 'trial_num')
