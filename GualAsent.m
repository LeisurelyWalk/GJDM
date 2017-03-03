function [x] = GualAsent(im,alfa,im_l,lambda)
im = double(im);
im_l = double(im_l);
[ml,nl] = size(im_l);
[m,n] = size(im);
p = fspecial('gaussian', 5, 1);
p = p.^2;
p = p./sum(p(:));

B = Compute_NLM_Matrix(im,5);
BTB = B'*B;
% I = eye(size(BTB));
% x = zeros(m,n);
z = zeros(ml,nl);
% INim = inv(I + BTB);
x = im;
for ii = 1:100
%     if (mod(ii,10) == 0)
%         B = Compute_NLM_Matrix(im_h,5);
%         BTB = B'*B;
%     end
for j = 1:10
    c = conv2(imresize(z,[m,n]),p,'same');
    G = 2*BTB*x(:) + 2*x(:) - 2*im(:) + c(:);
    x = x(:) - lambda*G;
end
%     x = im - 1/2*conv2(imresize(z,[m,n]),p,'same');
%     x = INim * x(:);
    x = reshape(x,[m,n]);
    z = z + alfa*(imresize(x,[ml,nl]) - im_l);
end