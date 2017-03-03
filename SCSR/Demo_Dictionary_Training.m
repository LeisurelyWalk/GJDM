% ========================================================================
% Demo codes for dictionary training by joint sparse coding
% 
% Reference
%   J. Yang et al. Image super-resolution as sparse representation of raw
%   image patches. CVPR 2008.
%   J. Yang et al. Image super-resolution via sparse representation. IEEE 
%   Transactions on Image Processing, Vol 19, Issue 11, pp2861-2873, 2010
%
% Jianchao Yang
% ECE Department, University of Illinois at Urbana-Champaign
% For any questions, send email to jyang29@uiuc.edu
% =========================================================================

clear all; clc; close all;
addpath(genpath('RegularizedSC'));

%TR_IMG_PATH = 'Data/Training';
TR_IMG_PATH = 'Data/newtrain';
dict_size   = 256;          % dictionary size
lambda      = 0.15;         % sparsity regularization
patch_size  = 5;            % image patch size
nSmp        = 100000;       % number of patches to sample
upscale     = 2;            % upscaling factor

% randomly sample image patches
[Xh, Xl] = rnd_smp_patch(TR_IMG_PATH, '*.jpg', patch_size, nSmp, upscale);

% prune patches with small variances, threshould chosen based on the
% training data
[Xh, Xl] = patch_pruning(Xh, Xl, 10);

% joint sparse coding 
%256
[Dh, Dl] = train_coupled_dict(Xh, Xl, dict_size, lambda, upscale);
dict_path = ['Dictionary/new1D_' num2str(dict_size) '_' num2str(lambda) '_' num2str(patch_size) '_s' num2str(upscale) '.mat' ];
save(dict_path, 'Dh', 'Dl');

dict_size   = 512; 
[Dh, Dl] = train_coupled_dict(Xh, Xl, dict_size, lambda, upscale);
dict_path = ['Dictionary/new1D_' num2str(dict_size) '_' num2str(lambda) '_' num2str(patch_size) '_s' num2str(upscale) '.mat' ];
save(dict_path, 'Dh', 'Dl');

dict_size   = 1024; 
[Dh, Dl] = train_coupled_dict(Xh, Xl, dict_size, lambda, upscale);
dict_path = ['Dictionary/new1D_' num2str(dict_size) '_' num2str(lambda) '_' num2str(patch_size) '_s' num2str(upscale) '.mat' ];
save(dict_path, 'Dh', 'Dl');

dict_size   = 2048; 
[Dh, Dl] = train_coupled_dict(Xh, Xl, dict_size, lambda, upscale);
dict_path = ['Dictionary/new1D_' num2str(dict_size) '_' num2str(lambda) '_' num2str(patch_size) '_s' num2str(upscale) '.mat' ];
save(dict_path, 'Dh', 'Dl');

upscale     = 3; 
% randomly sample image patches
[Xh, Xl] = rnd_smp_patch(TR_IMG_PATH, '*.jpg', patch_size, nSmp, upscale);

% prune patches with small variances, threshould chosen based on the
% training data
[Xh, Xl] = patch_pruning(Xh, Xl, 10);

% joint sparse coding 
dict_size   = 1024; 
[Dh, Dl] = train_coupled_dict(Xh, Xl, dict_size, lambda, upscale);
dict_path = ['Dictionary/new1D_' num2str(dict_size) '_' num2str(lambda) '_' num2str(patch_size) '_s' num2str(upscale) '.mat' ];
save(dict_path, 'Dh', 'Dl');