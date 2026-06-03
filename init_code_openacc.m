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
        error(['MOSEK was not found on the MATLAB path. Install MOSEK and add ', ...
            'its MATLAB toolbox directory before running the examples.']);
    end

    try
        [rcode, ~] = mosekopt('symbcon echo(0)');
    catch ME
        error(['MOSEK was found but could not be loaded. On Windows, make sure ', ...
            'the MOSEK platform bin directory is on PATH before calling mosekopt. ', ...
            'For example, if MOSEK is installed in C:\Program Files\Mosek\11.0, ', ...
            'run: setenv(''PATH'', [getenv(''PATH'') '';C:\Program Files\Mosek\11.0\tools\platform\win64x86\bin'']); ', ...
            'addpath(''C:\Program Files\Mosek\11.0\toolbox\r2017a''). ', ...
            'If you use MOSEK 10, replace 11.0 with 10.0. ', ...
            'Original error: %s'], ...
            ME.message);
    end

    if rcode ~= 0
        error('MOSEK returned a nonzero status during initialization: %d.', rcode);
    end
end
