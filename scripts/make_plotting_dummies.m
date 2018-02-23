% Make dummy structures
% %%%%%%%%%%%%%%%%%%%%%
%
% Make dummy FieldTrip structures for frequency and source space data to use
% for subsequent plotting of classifier outputs.
%
% AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
% LICENCE: BSD 3-clause
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a fake structure for TFR data
% a lot of this will be overwritten before plotting, this is just to provide a
% shape to fill things in.
num_gridpoints = sum(template_grid.inside);  %  number of grid points

data_dummy = [];
for ii = 1:num_gridpoints
    data_dummy.label{ii} = num2str(ii);  % fake label
end
data_dummy.time = [1:4];  % time points
data_dummy.freq = [1:5];  % freqs
data_dummy.powspctrm = NaN(num_gridpoints, length(data_dummy.freq), ...
                           length(data_dummy.time));
data_dummy.dimord = 'chan_freq_time';


% make fake source data structure
% crucial part will be overwritten with meaningful data
% this structure is based on the template grid (MNI pos)
source_dummy.freq =  NaN;
source_dummy.dim = template_grid.dim;
source_dummy.inside = template_grid.inside;
source_dummy.pos = template_grid.pos;
source_dummy.methode = 'average';
source_dummy.avg.pow = NaN(size(template_grid.pos, 1), 1);
source_dummy.filterdimord = '{pos}_ori_chan';
