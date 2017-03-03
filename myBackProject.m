close all;
addpath('Data/Testing');
addpath('Results');
imname = 'TDSSnew3325.tiff';
im_lname = 'input.bmp';
im_hname = 'new3325.tiff';
im_h = imread(imname);
im_h2 = im_h;
% im_l = imread(im_lname);
% im = imread(fullfile('Data/Testing',im_hname));
im = imread(im_hname);
im_l = imresize(im,1/2);
imshow(im_h);
% change color space, work on illuminance only
% im_l_ycbcr = rgb2ycbcr(im_l);
% im_l_y = im_l_ycbcr(:, :, 1);
% im_l_cb = im_l_ycbcr(:, :, 2);
% im_l_cr = im_l_ycbcr(:, :, 3);
% 
% im_h_ycbcr = rgb2ycbcr(im_h);
% im_h_y = im_h_ycbcr(:, :, 1);
% im_h_cb = im_h_ycbcr(:, :, 2);
% im_h_cr = im_h_ycbcr(:, :, 3);
alfa = 1.5;
lambda = 0.005;
im_h = GualAsent(im_h,alfa,im_l,lambda);

% image resuper-resolution base on NL models
% im_h_y = presetimage(im_h_y);
% upscale the chrominance simply by "bicubic" 
[nrow, ncol] = size(im_h);
% im_h_ycbcr = zeros([nrow, ncol, 3]);
% im_h_ycbcr(:, :, 1) = im_h_y;
% im_h_ycbcr(:, :, 2) = im_h_cb;
% im_h_ycbcr(:, :, 3) = im_h_cr;
% im_h = ycbcr2rgb(uint8(im_h_ycbcr));
im_h = uint8(im_h);
figure;
imshow(im_h);
resname = ['TDG' imname];
path = ['Data/results/' resname];
imwrite(im_h,path);
% im = imread(fullfile('Data/Testing',im_hname));
im = imresize(im,[nrow, ncol]);
% compute PSNR for the illuminance channel
sp_rmse = compute_rmse(im, im_h);
sb_rmse = compute_rmse(im, im_h2);
sb_psnr = 20*log10(255/sb_rmse);
sp_psnr = 20*log10(255/sp_rmse);
fprintf('PSNR for Sparse Representation Recovery(alfa %4f lambda %4f): %f dB\n',alfa, lambda,sp_psnr);
fprintf('PSNR not for MyBackProject: %f dB\n', sb_psnr);