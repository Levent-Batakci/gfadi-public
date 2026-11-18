function setup_paths()
%SETUP_PATHS  Put the code folders on the MATLAB path and create the output folders.
%   SETUP_PATHS() adds the repository root and its code folders to the path and creates
%   results/, figs/ and logs/ in the root if they do not exist.
%   Run once per MATLAB session, from anywhere:  run('path/to/repository/setup_paths')
root = fileparts(mfilename('fullpath'));
addpath(root);
for d = {'zolo', 'gfadi', 'shifts', 'solvers', 'problems', 'experiments', 'figures', 'testbed', 'examples', 'teaser'}
    addpath(fullfile(root, d{1}));
end
for d = {'results', 'figs', 'logs'}
    if ~exist(fullfile(root, d{1}), 'dir'), mkdir(fullfile(root, d{1})); end
end
end
