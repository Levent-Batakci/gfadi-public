function r = app_root()
%APP_ROOT  The root folder of the repository.
%   R = APP_ROOT() returns the folder that holds this file.  Runs write to R/results, figure
%   scripts to R/figs, and drivers' logs to R/logs.
r = fileparts(mfilename('fullpath'));
end
