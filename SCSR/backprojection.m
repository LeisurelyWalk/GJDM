function [im_h] = backprojection(im_h, im_l, maxIter)

[row_l, col_l] = size(im_l);
[row_h, col_h] = size(im_h);

p = fspecial('gaussian', 5, 1);
p = p.^2;
p = p./sum(p(:));

im_l = double(im_l);
im_h = double(im_h);
% B = Compute_NLM_Matrix(im_h,5);

for ii = 1:maxIter,

    
    im_l_s = imresize(im_h, [row_l, col_l], 'bicubic');
    im_diff = im_l - im_l_s;
    
    im_diff = imresize(im_diff, [row_h, col_h], 'bicubic');
%     if (mod(ii,10) == 0)
%         B = Compute_NLM_Matrix(im_h,5);
%     end
%     im_h = im_h(:);
%     im_h = im_h + B*im_h;
%     im_h = reshape(im_h,row_h,col_h);
    im_h = im_h + conv2(im_diff, p, 'same');
end
    