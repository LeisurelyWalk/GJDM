% =========================================================================
% Simple demo codes for image super-resolution via sparse representation
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

clear all; clc;

% read test image
newpath = 'fig';
mysize = 1024;
imname = 'new215_2.tiff';
im = imread(fullfile('Data/Testing',imname));
im_l = imresize(im,1/3);
lrname = ['LR' imname];
imwrite(im_l, fullfile('fig',lrname));
imwrite(im,fullfile('fig',imname));
% im_l = imresize(im_l,1/3);%Àı–°≤‚ ‘£°
% set parameters
lambda = 0.2;                   % sparsity regularization
overlap = 4;                    % the more overlap the better (patch size 5x5)
up_scale = 3;                   % scaling factor, depending on the trained dictionary
maxIter = 20;                   % if 0, do not use backprojection

% load dictionary
load('Dictionary/D_1024_0.15_5_s3.mat');

% change color space, work on illuminance only
im_l_ycbcr = rgb2ycbcr(im_l);
im_l_y = im_l_ycbcr(:, :, 1);
im_l_cb = im_l_ycbcr(:, :, 2);
im_l_cr = im_l_ycbcr(:, :, 3);

% image super-resolution based on sparse representation
[im_h_y] = ScSR(im_l_y, up_scale, Dh, Dl, lambda, overlap);
im_h_unB = im_h_y;
[im_h_y] = backprojection(im_h_y, im_l_y, maxIter);
% image resuper-resolution base on NL models
% im_h_y = presetimage(im_h_y);
% upscale the chrominance simply by "bicubic" 
[nrow, ncol] = size(im_h_y);
im_h_cb = imresize(im_l_cb, [nrow, ncol], 'bicubic');
im_h_cr = imresize(im_l_cr, [nrow, ncol], 'bicubic');

im_h_ycbcr = zeros([nrow, ncol, 3]);
im_h_ycbcr(:, :, 1) = im_h_y;
im_h_ycbcr(:, :, 2) = im_h_cb;
im_h_ycbcr(:, :, 3) = im_h_cr;
im_h = ycbcr2rgb(uint8(im_h_ycbcr));

im_h_unBc = zeros([nrow,ncol,3]);
im_h_unBc(:,:,1) = im_h_unB;
im_h_unBc(:,:,2) = im_h_cb;
im_h_unBc(:,:,3) = im_h_cr;
im_h_unB = ycbcr2rgb(uint8(im_h_unBc));
resname = ['F_' num2str(up_scale) '_' num2str(mysize) 'S' imname];
imwrite(im_h_unB,fullfile('Results',resname));
imwrite(im_h_unB,fullfile('fig',resname));

% bicubic interpolation for reference
im_b = imresize(im_l, [nrow, ncol], 'bicubic');

% read ground truth image
% im = imread('Data/Testing/gnd.bmp');
im = imresize(im,[nrow, ncol]);
% compute PSNR for the illuminance channel
bb_rmse = compute_rmse(im, im_b);
sp_rmse = compute_rmse(im, im_h);

bb_psnr = 20*log10(255/bb_rmse);
sp_psnr = 20*log10(255/sp_rmse);

fprintf('PSNR for Bicubic Interpolation: %f dB\n', bb_psnr);
fprintf('PSNR for Sparse Representation Recovery: %f dB\n', sp_psnr);
im_b = uint8(im_b);
bname = ['F_' num2str(up_scale) '_' num2str(mysize) 'B' imname];
imwrite(im_b,fullfile('Results',bname));
imwrite(im_b,fullfile(newpath,bname));

% show the images
figure, imshow(im_h);
title('Sparse Recovery');
figure, imshow(im_b);
title('Bicubic Interpolation');