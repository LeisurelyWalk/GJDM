function res = Jmd(imname)
% read test image
newpath = 'newfig';
mysize = 1024;
% imname = 'new2102.tiff';
im = imread(fullfile('testTD',imname));
im_l = imresize(im,1/2);
%lrname = ['LR' imname];
% im_l = imresize(im_l,1/3);%缩小测试！
% set parameters
lambda = 0.2;                   % sparsity regularization
overlap = 4;                    % the more overlap the better (patch size 5x5)
up_scale = 2;                   % scaling factor, depending on the trained dictionary
maxIter = 20;                   % if 0, do not use backprojection

% load dictionary
load('Dictionary/D_1024_0.15_5_2.mat');%1024


% change color space, work on illuminance only
simage = size(im);
if length(simage) == 3

    im_l_ycbcr = rgb2ycbcr(im_l);
    im_l_y = im_l_ycbcr(:, :, 1);
    im_l_cb = im_l_ycbcr(:, :, 2);
    im_l_cr = im_l_ycbcr(:, :, 3);
else
    im_l_y = im_l;
end

% image super-resolution based on sparse representation
[im_h_y] = ScSR_JMD(im_l_y, up_scale, Dh, Dl, lambda, overlap);

% im_h_y = presetimage(im_h_y);
% upscale the chrominance simply by "bicubic" \

[nrow, ncol] = size(im_h_y);
%HR 其他两个分量
if length(simage) == 3
    im_h_cb = imresize(im_l_cb, [nrow, ncol], 'bicubic');
    im_h_cr = imresize(im_l_cr, [nrow, ncol], 'bicubic');
    % G ScSR HR
    im_h_ycbcr = zeros([nrow, ncol, 3]);
    im_h_ycbcr(:, :, 1) = im_h_y;
    im_h_ycbcr(:, :, 2) = im_h_cb;
    im_h_ycbcr(:, :, 3) = im_h_cr;
    im_h = ycbcr2rgb(uint8(im_h_ycbcr));
else
    im_h = im_h_y;
end
im = imresize(im,[nrow, ncol]);
% bicubic interpolation for reference


%处理S
resname = ['F_' num2str(up_scale) '_' num2str(mysize) 'JMd' imname];
%imwrite(im_h_unB,fullfile('Results',resname));
ss_rmse = compute_rmse(im, im_h);
ss_psnr = 20*log10(255/ss_rmse);
im_h = uint8(im_h);
imwrite(im_h,fullfile(newpath,resname));
res = ss_psnr;


% read ground truth image
% im = imread('Data/Testing/gnd.bmp');

% compute PSNR for the illuminance channel
myres = 'Done';
fprintf(myres);
fprintf('\n');

