% INITPATH Initialize paths for ASPIRE toolbox
%
% Usage
%    initpath();

% Find where the package installed.
[pathstr, ~, ~] = fileparts(mfilename('fullpath'));

addpath(pathstr);

addpath(genpath(fullfile(pathstr,'abinitio')));
addpath(genpath(fullfile(pathstr,'basis')));
addpath(genpath(fullfile(pathstr,'common')));
addpath(genpath(fullfile(pathstr,'examples')));
addpath(genpath(fullfile(pathstr,'fourier')));
addpath(genpath(fullfile(pathstr,'install')));
addpath(genpath(fullfile(pathstr,'io')))
addpath(genpath(fullfile(pathstr,'projections')))
addpath(genpath(fullfile(pathstr,'sinograms')))
addpath(genpath(fullfile(pathstr,'reconstruction')))
addpath(genpath(fullfile(pathstr,'refinement')))
addpath(genpath(fullfile(pathstr,'workflow')))
addpath(genpath(fullfile(pathstr,'cwf_denoise')))
addpath(genpath(fullfile(pathstr,'kam_cryo')))

if exist(fullfile(pathstr, 'extern', 'nufftall-1.33'))
    addpath(fullfile(pathstr, 'extern', 'nufftall-1.33'));
end

if exist(fullfile(pathstr, 'extern', 'nfft'))
    addpath(fullfile(pathstr, 'extern', 'nfft', 'lib'));
    addpath(fullfile(pathstr, 'extern', 'nfft', 'share', 'nfft', ...
        'matlab', 'nfft'));
end

if exist(fullfile(pathstr, 'extern', 'finufft'))
    addpath(fullfile(pathstr, 'extern', 'finufft', 'matlab'));
end

if exist(fullfile(pathstr, 'extern', 'SDPLR-1.03-beta'))
    addpath(fullfile(pathstr,'extern','SDPLR-1.03-beta'));
end

%%%%%%%%
% MANOPT
% this supposes we are under the ASPIRE root directory
% and manopt is up at the same level as ASPIRE root directory
%cd  ../manopt;
%addpath(pwd);
% Recursively add Manopt directories to the Matlab path.
%cd manopt;
%addpath(genpath(pwd));
%cd ..;
%cd ../aspire;
addpath(genpath(fullfile(pathstr,'../manopt')));
%%%%%%%%

%%%%%%%%
% SHT - Spherical Harmonic Transform TOolbox
%cd ../SHT;
%addpath(genpath(pwd));
%cd ../aspire;
addpath(genpath(fullfile(pathstr,'../SHT')));
%

clear pathstr;
