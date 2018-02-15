% Project settings for MATLAB files
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Paths, file names, and variables shared between MATLAB scripts. Also takes
% care of loading FieldTrip.
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
Vp = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05', 'VP06', 'VP08', 'VP09', ...
      'VP10', 'VP11', 'VP12', 'VP13', 'VP18', 'VP20', 'VP22', 'VP23', ...
      'VP25', 'VP26', 'VP115', 'VP117'}


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
