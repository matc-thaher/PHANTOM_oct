function result = colossus_query(quantity, input_array, cosmo_id, model, z, python_exe, bridge_path, extra)
% colossus_query  Call Colossus from MATLAB via Python bridge.
%
% result = colossus_query('power_spectrum', k, 'planck15', 'eisenstein98', 0, python_exe, bridge_path)
% result = colossus_query('variance',       R, 'planck15', 'eisenstein98', 0, python_exe, bridge_path)
% result = colossus_query('correlation_function', r, 'planck15', 'eisenstein98', 0, python_exe, bridge_path)
% result = colossus_query('concentration', M_vec, 'planck15', 'ishiyama21',  z, python_exe, bridge_path);

if nargin < 4 || isempty(model),       model       = 'eisenstein98'; end
if nargin < 5 || isempty(z),           z           = 0.0;            end
if nargin < 6 || isempty(python_exe),  python_exe  = 'python';       end
if nargin < 7 || isempty(bridge_path)
    bridge_path = fullfile(fileparts(mfilename('fullpath')), 'colossus_bridge.py');
end
if nargin < 8,  extra = struct();  end

% Write task JSON
task_file   = [tempname, '.json'];
output_file = [tempname, '.txt'];

% Build task struct
task.cosmo_id  = cosmo_id;
task.quantity  = quantity;
task.model     = model;
task.z         = z;


switch quantity
    case 'power_spectrum'
        task.k = input_array(:)';
    case 'variance'
        task.R = input_array(:)';
    case 'correlation_function'
        task.r = input_array(:)';
    case 'concentration'
        task.M    = input_array(:)';
        task.mdef = '200c';
    case {'peakheight', 'neff'}
        task.M    = input_array(:)';
    case 'profile_batch'
        task.r    = input_array(:)';
        task.M    = extra.M;
        task.c    = extra.c;
        task.mdef = extra.mdef;
    case 'hmf'
        task.M    = input_array(:)';
        % mdef is resolved model-aware on the Python side.
        % Pass explicitly only if you need to override the default.
        % e.g. colossus_query('hmf', M, cosmo_id, 'tinker08', z, ..., struct('mdef','500c'))
        if isfield(extra, 'mdef')
            task.mdef = extra.mdef;
        end
        % model is already set from the model argument above
    case 'delta_c'
        task.z = z;   % z already set above; nothing else needed

    case 'sigma_fof'
    task.M = input_array(:)';
end

% Write JSON
fid = fopen(task_file, 'w');
fprintf(fid, '%s', jsonencode(task));
fclose(fid);

% Call Python
cmd = sprintf('"%s" "%s" "%s" "%s"', python_exe, bridge_path, task_file, output_file);
[status, msg] = system(cmd);
if status ~= 0
    error('colossus_query: Python call failed.\n%s', msg);
end

% Read result
data = dlmread(output_file, '', 1, 0);   % skip 1 header line
% result = data(:, 2);   % second column is the quantity
if size(data, 2) > 2
    result = data(:, 2:end);   % profile_batch: return all quantity columns
else
    result = data(:, 2);       % all other quantities: single column
end

% Clean up
delete(task_file);
delete(output_file);
end