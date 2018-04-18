clear all; clc;

% read test image
imname = 'new3325.tiff';
% im = imread(fullfile('Data/Testing',imname));
im = imread(imname);
im_l = imresize(im,1/2);
% im_l = imresize(im_l,1/3);%Àı–°≤‚ ‘£°
% set parameters
lambda = 0.2;                   % sparsity regularization
overlap = 4;                    % the more overlap the better (patch size 5x5)
up_scale = 2;                   % scaling factor, depending on the trained dictionary
maxIter = 20;                   % if 0, do not use backprojection

% load dictionary
load('Dictionary/D_1024_0.15_5_s2.mat');
im = double(im);
im_l = double(im_l);
% change color space, work on illuminance only
% im_l_ycbcr = rgb2ycbcr(im_l);
% im_l_y = im_l_ycbcr(:, :, 1);
% im_l_cb = im_l_ycbcr(:, :, 2);
% im_l_cr = im_l_ycbcr(:, :, 3);

% image super-resolution based on sparse representation
[im_h] = ScSR(im_l, 2, Dh, Dl, lambda, overlap);
im_h_unB = im_h;
[im_h] = backprojection(im_h, im_l, maxIter);
% image resuper-resolution base on NL models
% im_h_y = presetimage(im_h_y);
% upscale the chrominance simply by "bicubic" 
 [nrow, ncol] = size(im_h);
% im_h_cb = imresize(im_l_cb, [nrow, ncol], 'bicubic');
% im_h_cr = imresize(im_l_cr, [nrow, ncol], 'bicubic');

% im_h_ycbcr = zeros([nrow, ncol, 3]);
% im_h_ycbcr(:, :, 1) = im_h_y;
% im_h_ycbcr(:, :, 2) = im_h_cb;
% im_h_ycbcr(:, :, 3) = im_h_cr;
% im_h = ycbcr2rgb(uint8(im_h_ycbcr));
% 
% im_h_unBc = zeros([nrow,ncol,3]);
% im_h_unBc(:,:,1) = im_h_unB;
% im_h_unBc(:,:,2) = im_h_cb;
% im_h_unBc(:,:,3) = im_h_cr;
% im_h_unB = ycbcr2rgb(uint8(im_h_unBc));
im_h_unB = uint8(im_h_unB);
resname = ['TDSS' imname];
imwrite(im_h_unB,fullfile('Results',resname));
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
im_h = uint8(im_h);
im_b = uint8(im_b);
imwrite(im_b,fullfile('Results','TDBB3325.tiff'));

% show the images
figure, imshow(im_h);
title('Sparse Recovery');
figure, imshow(im_b);
title('Bicubic Interpolation');