img_path = 'testTD';
type = '*.tiff';
img_dir = dir(fullfile(img_path, type));%?��???��??bmp

img_num = length(img_dir);%?��?�??
% nper_img = zeros(1, img_num);
TD = 1;
for ii = 1:length(img_dir),
    imname = img_dir(ii).name;
    my = Demo_mySR(TD,imname);
%     im = imread(fullfile(img_path, img_dir(ii).name));
%     nper_img(ii) = prod(size(im));%计�??��???????�?��
end
TD = 2;
for ii = 1:length(img_dir),
    imname = img_dir(ii).name;
    my = Demo_mySR(TD,imname);
%     im = imread(fullfile(img_path, img_dir(ii).name));
%     nper_img(ii) = prod(size(im));%计�??��???????�?��
end