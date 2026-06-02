function package_dir = init_code_openacc(script_dir)
%INIT_CODE_OPENACC Configure paths for the released code_openacc package.

    if nargin < 1 || isempty(script_dir)
        script_dir = fileparts(mfilename('fullpath'));
    end

    package_dir = script_dir;
    while ~isfolder(fullfile(package_dir, 'Truncation'))
        parent_dir = fileparts(package_dir);
        if strcmp(parent_dir, package_dir)
            error('Cannot locate the code_openacc package root from %s.', script_dir);
        end
        package_dir = parent_dir;
    end

    addpath(genpath(package_dir));

    if exist('mosekopt', 'file') ~= 3 && exist('mosekopt', 'file') ~= 2
        warning(['MOSEK was not found on the MATLAB path. ', ...
            'Install MOSEK and run its setup script before executing optimization examples.']);
    end
end
