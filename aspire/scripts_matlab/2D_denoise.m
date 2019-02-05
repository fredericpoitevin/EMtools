clear all, close all

totTime = tic;


addpath '/home/users/fpoitevi/Toolkit/EMtools/aspire/aspire/';
initpath

%% PATH INFO TO DATA
%%star_path='/scratch/users/fpoitevi/data/f20_oct2016/Extract/job015/particles_subset.star';
star_path='/scratch/users/fpoitevi/data/f20_oct2016/Extract/job015/particles.star';
prefix='/scratch/users/fpoitevi/data/f20_oct2016/';


% PARAMETERS
num_pool = num_cores()  ;       % number of workers
n_im = 1000;              % number of images
n_noise = 100;        % number of images for noise estimation
downsampleddim = 120; % image size
%c = 0.25;             % bandlimit (?)
%R = 150;              % (?)
L0 = downsampleddim;  % (?) 
%n_r = ceil(4*c*R);    % (?)

% READ STAR AND INIT SOME ARRAYS
[CTFdata] = readSTAR(star_path); % Note: These are in real space. Convert to Fourier before denoising
CTFdata.data = CTFdata.data(1:n_im);
ctfid = zeros(numel(CTFdata.data),1); % Index for ctf group
projs = zeros( downsampleddim , downsampleddim , numel(CTFdata.data) ); 
% open_log('test_log.txt') this guy does not seem to be used...

lastStackProcessed = '';

for k = 1:numel(CTFdata.data)
    imageID = CTFdata.data{k}.rlnImageName;
    imparts = strsplit(imageID,'@');
    imageidx = str2double(imparts{1});
    stackname = imparts{2};
    if ~strcmp(stackname,lastStackProcessed) % Save only one CTF per stack. each stack has the same parameters.
        mrcstack = ReadMRC((fullfile(prefix,stackname)));
        lastStackProcessed = stackname;
    end
    projs(:,:,k) = mrcstack(:,:,imageidx); % PROJ K IS STORED IN projs(:,:,k)
end


projs = cryo_normalize_background(projs,round(downsampleddim/2)-10); %we would better pre-normalize according to the instructions
% Normalize images  ********* SKIP BACKGROUND NORMALIZATION ??? ******* ==> % log_message('Normalize background'); % projs = cryo_normalize_background(projs,round(n/2)-10); 
% Estimate Power Spectral Density of images 
psd = cryo_noise_estimation(projs(:,:,1:n_noise)); %
%disp(psd(n,:));
% Apply whitening filter to images (inverse square root of PSD)

[ prewhitened_projs , whiten_filter , nzidx ] = cryo_prewhiten( projs , psd ); %prewhitened_projs is in real space


energy_thresh = 0.99;
[c, R] = choose_support_v6(cfft2(prewhitened_projs), energy_thresh);%estimate band limit and compact support size
c = c*(0.5/floor(downsampleddim/2)); %rescaling between 0 and 0.5
n_r = ceil(4*c*R);


% PRECOMPUTE Fourier-Bessel BASIS (given c and R)
tic_basis = tic;
[ basis, sample_points ] = precomp_fb( n_r, R, c );
timing.basis = toc(tic_basis)

% PARSE MRCS FILES AND EXTRACT RELEVANT INFO
tic_parse = tic;
lastStackProcessed = '';
def_grp=1;
for k = 1:numel(CTFdata.data)
    imageID = CTFdata.data{k}.rlnImageName;
    imparts = strsplit(imageID,'@');
    imageidx = str2double(imparts{1});
    stackname = imparts{2};
    if ~strcmp(stackname,lastStackProcessed) % Save only one CTF per stack. each stack has the same parameters.
        mrcstack = ReadMRC((fullfile(prefix,stackname)));
        lastStackProcessed = stackname;
        [ voltage , DefocusU , DefocusV , DefocusAngle , Cs , pixA , A ] = cryo_parse_Relion_CTF_struct( CTFdata.data{k} );
        pixelscaling = size(projs,1)/downsampleddim;
        pixAdownsampled = pixA*pixelscaling;
        ctfs(:,:,def_grp) = cryo_CTF_Relion( downsampleddim , voltage , DefocusU , DefocusV , DefocusAngle , Cs , pixAdownsampled , A );
        ctfs_rad(:,def_grp)=cryo_CTF_Relion_radial( downsampleddim , voltage , DefocusU , DefocusV , DefocusAngle , Cs , pixAdownsampled , A , sample_points.r ); % Note: no scaling for sample_points.r
        sprintf('Read stack %d : %s',def_grp, stackname)
        def_grp = def_grp+1; 
    end
    ctfid(k) = def_grp-1; % CTF GROUP OF PROJ K
%    projs(:,:,k) = mrcstack(:,:,imageidx); % PROJ K IS STORED IN projs(:,:,k)
end
clear mrcstack
timing.parse=toc(tic_parse)
disp('Done reading')

%% PREWHITENING
disp('Downsampling')
tic_whit = tic;


%[ projs , whiten_filter , nzidx ] = Prewhiten_image2d_tejal( projs , psd ); % [ projs , whiten_filter , nzidx ] = cryo_prewhiten( projs , psd );
[ noise_v_r ] = estimate_noise_real(prewhitened_projs(:,:,1:n_noise)); % Variance of noise in prewhitened projections
% write 'mean' to mrc
mean_img = mean(projs, 3); % mean prewhitened image
WriteMRC(mean_img,pixA,'prewhitened_mean.mrc');
% I believe the next 3 lines fill missing part of the filter with ones
idx = [1:numel(whiten_filter)];
zidx = setdiff( idx , nzidx );
whiten_filter(zidx) = ones(size(whiten_filter(zidx)));
% seems here that w_f is the white filter properly sized and everything...
projs = cfft2(prewhitened_projs); % Images in fourier space
w_f = cryo_downsample( whiten_filter , size(projs,1) ).^(-1); % Image is divided by the filter elementwise, so -1 % This still has some large pixels values
timing.whit = toc(tic_whit)

%% ESTIMATION OF MEAN (mu)
ndef = def_grp - 1;                   % number of defocus groups
w_CTF = ctfs .* repmat(w_f,1,1,ndef); % this must be the 'whitened' CTF filter
regu = 1;                             % this must be the regularization parameter lambda
tic_mean = tic;
mean_image_f = mean_LS( ctfs , ctfid , projs , regu ); % Solve better conditioned system to get W\mu
mean_image_f = mean_image_f ./ w_f;                    % then get \mu
timing.mean = toc(tic_mean)
mean_image_f = double(mean_image_f);
projs = double(projs);
% write to mrc file after inverting back to real space.
mean_image_f_real = real(icfft2(mean_image_f));
WriteMRC(mean_image_f_real,pixA,'mean.mrc');

%% DEMEANING (centering on mu)
tic_demean = tic;
[projs] = demean_y_v6( projs , w_CTF , mean_image_f , ctfid );
projs = real(icfft2(projs)); % at that stage, projs are the CTF-corrected mean-subtracted images, in real space
WriteMRC(projs,pixA,'CTF_corrected_and_mean_substracted_projs.mrcs');
timing.demean = toc(tic_demean)
% coeff_ymu: coeff of (CTF-corrected and mean-subtracted) projs on basis



tic_coeffymu = tic;
[ coeff_ymu ] = coeff_demean( projs , R , basis , sample_points , num_pool); %now projs is in real space
timing.coeffymu = toc(tic_coeffymu)
% coeff_mean: coeff of real-space mean (mu) on same basis
[coeff_mean] = coeff_demean( mean_image_f_real , R , basis , sample_points , num_pool );

%% CTF in new basis: numerical integration
if mod(L0,2)==1
    w_f_rad = interp1( [0:floor(L0/2)] , w_f(floor(L0/2)+1 , floor(L0/2)+1:end) , sample_points.r*((floor(L0/2))/0.5) , 'linear' );
else
    w_f_rad = interp1( [0:floor(L0/2)] , w_f(floor(L0/2) , floor(L0/2):end) , sample_points.r*((floor(L0/2))/0.5) , 'linear' );
end

%% CCWF - Covariance Wiener Filtering
tic_ccwf = tic;
[ denoised_coeff_ccwf , ~ , ~ , num_eigs ] = jobscript_CCWF_cgshrink_jsb( ctfid , w_f_rad , ctfs_rad , basis , sample_points , coeff_mean , coeff_ymu , noise_v_r ); % was: %[ denoised_coeff_ccwf , ~ , ~ , num_eigs, C_FB ] = jobscript_CCWF_cgshrink_jsb( ctfid , w_f_rad , ctfs_rad , basis , sample_points , coeff_mean , coeff_ymu , noise_v_r );
C_FB{1} = denoised_coeff_ccwf{1}*denoised_coeff_ccwf{1}'/n_im;
timing.ccwf = toc(tic_ccwf)

%% Save data for paper figures
tic_recon_imgs = tic;
[recon] = recon_images_FB( c , R , L0 , denoised_coeff_ccwf , 1 , n_im );  %n is the number of images to be denoised
recon = -recon; % images were already contrast-inverted
WriteMRC(recon,pixA, 'denoised.mrcs'); 
timing.recon = toc(tic_recon_imgs);

totTime = toc(totTime);
clearvars -except C_FB recon c R mean_img
