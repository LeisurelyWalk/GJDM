img_path = 'testTD';
type = '*.tiff';
img_dir = dir(fullfile(img_path, type));%?·å???»¶??bmp

img_num = length(img_dir);%?¾ç?ä¸??
% nper_img = zeros(1, img_num);
TD = 1;
for ii = 1:length(img_dir),
    imname = img_dir(ii).name;
    my = Demo_mySR(TD,imname);
%     im = imread(fullfile(img_path, img_dir(ii).name));
%     nper_img(ii) = prod(size(im));%è®¡ç??°ç???????ä¹?§¯
end
TD = 2;
for ii = 1:length(img_dir),
    imname = img_dir(ii).name;
    my = Demo_mySR(TD,imname);
%     im = imread(fullfile(img_path, img_dir(ii).name));
%     nper_img(ii) = prod(size(im));%è®¡ç??°ç???????ä¹?§¯
end