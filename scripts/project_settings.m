% Project settings for MATLAB files
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Paths, file names, and variables shared between MATLAB scripts. Also takes
% care of adding FieldTrip to the MATLAB path.
%
% AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
% LICENCE: BSD 3-clause
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Toolbox paths
fieldtrip_path = '/path/to/fieldtrip/';

% Data-related paths
raw_dir = '/path/to/sensor/data/';
results_dir = '/path/to/results/';
headmodel_dir = '/path/to/headmodels/';
source_dir = '/path/to/sourcespace/data/';
fig_dir = '/path/for/saving/figures';

% path to custom template grid
template_grid_file = fullfile(raw_dir, 'templategrid_15mm.mat');

% path to Fieldtrip templates
atlas_path = fullfile(fieldtrip_path, 'template/atlas/aal/ROI_MNI_V4.nii');
mri_path = fullfile(fieldtrip_path, 'template/anatomy/single_subj_T1.nii');

% subjects list
Vp = {'S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', ...
      'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20'}


%% Load Fieldtrip if not yet loaded
try
    ft_defaults
catch
    warning('Fieldtrip is not on your path yet, adding it.');
    addpath(fieldtrip_path)
    ft_defaults
end

[ft_ver, ft_path] = ft_version;
display(sprintf('You are using Fieldtrip on path %s', ft_path));
