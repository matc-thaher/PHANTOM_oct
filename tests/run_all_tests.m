%% run_all_tests.m
% Master script to run all PHANTOM test files.
%
% Place BOTH run_all_tests.m AND run_test_isolated.m in F:\PHANTOM\tests\
% Then run:
%   >> cd F:\PHANTOM\tests
%   >> run_all_tests

clear; clc; close all;

TESTS_DIR = fileparts(mfilename('fullpath'));

test_files = {
    'test_Bullock'
    'test_Diemer15'
    'test_Diemer19'
    'test_Duffy'
    'test_Dutton'
    'test_Prada'
    'test_cCDM'
    'test_cFDM'
    'test_child18'
    'test_ishiyama21'
    'test_klypin11'
    'test_klypin16'
    'test_ludlow16'
    'test_profiles'
    'test_chmr'
    'test_halo_bias'
};

n         = numel(test_files);
results   = cell(n, 1);
err_msgs  = cell(n, 1);
durations = zeros(n, 1);

addpath(TESTS_DIR);

fprintf('\n========================================================\n');
fprintf('  PHANTOM — Running %d tests\n', n);
fprintf('========================================================\n\n');

for i = 1:n
    tname    = test_files{i};
    testpath = fullfile(TESTS_DIR, [tname '.m']);

    fprintf('[%2d/%2d]  %-22s ... ', i, n, tname);

    % run_test_isolated() wraps run() inside a FUNCTION scope.
    % This means "clear" inside any test script only clears that
    % function's local workspace — it can NEVER touch results{i},
    % durations(i), t0, or any other variable in this loop.
    [status, elapsed, msg] = run_test_isolated(testpath);

    results{i}   = status;
    durations(i) = elapsed;
    err_msgs{i}  = msg;

    if strcmp(status, 'PASS')
        fprintf('PASS  (%.2f s)\n', elapsed);
    elseif strcmp(status, 'SKIP')
        fprintf('SKIP  (file not found)\n');
    else
        fprintf('FAIL  (%.2f s)\n', elapsed);
        fprintf('         >> %s\n', msg);
    end

    close all;
end

%% Summary
n_pass = sum(strcmp(results, 'PASS'));
n_fail = sum(strcmp(results, 'FAIL'));
n_skip = sum(strcmp(results, 'SKIP'));

fprintf('\n========================================================\n');
fprintf('  SUMMARY:  %d passed  |  %d failed  |  %d skipped\n', ...
        n_pass, n_fail, n_skip);
fprintf('========================================================')





%% supporting function

function [status, elapsed, msg] = run_test_isolated(testpath)
% RUN_TEST_ISOLATED  Runs a test script inside a function workspace.
%
% KEY REASON THIS WORKS:
%   Each test script starts with "clear; clc; close all;"
%   When run() is called from inside a FUNCTION, that "clear" only
%   clears this function's local variables — it cannot touch anything
%   in the calling script (run_all_tests.m). So t0, results, durations
%   are all completely safe.
%
% Inputs:
%   testpath — full path to the .m test script
%
% Outputs:
%   status  — 'PASS', 'FAIL', or 'SKIP'
%   elapsed — elapsed time in seconds
%   msg     — error message string (empty on PASS)

status  = '';
elapsed = 0;
msg     = '';

if ~isfile(testpath)
    status = 'SKIP';
    return;
end

t0 = tic;
try
    run(testpath);
    elapsed = toc(t0);
    status  = 'PASS';
    msg     = '';
catch ME
    elapsed = toc(t0);
    status  = 'FAIL';
    msg     = ME.message;
end

end