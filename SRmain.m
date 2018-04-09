img_path = 'testTD';
type = '*.tiff';
time=zeros(1,3);
img_dir = dir(fullfile(img_path, type));

img_num = length(img_dir);
% nper_img = zeros(1, img_num);
TD = 3;
for ii = 1:length(img_dir),
    imname = img_dir(ii).name;
    my = Demo_mySR(time,img_path,TD,imname);
%     im = imread(fullfile(img_path, img_dir(ii).name));
%     nper_img(ii) = prod(size(im));
end
TD = 2;
for ii = 1:length(img_dir),
    imname = img_dir(ii).name;
    my = Demo_mySR(time,img_path,TD,imname);
%     im = imread(fullfile(img_path, img_dir(ii).name));
%     nper_img(ii) = prod(size(im));
end