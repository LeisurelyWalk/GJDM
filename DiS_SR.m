img_path = 'testTD';
type = '*.tiff';
img_dir = dir(fullfile(img_path, type));%?��???��??bmp

img_num = length(img_dir);%?��?�??
% nper_img = zeros(1, img_num);
TD = 1;%1��256��2��512��3��1024�� 4��2048
res = zeros(img_num,4);%psnr Bic,Jmd/ScSR,GSCSR,GJMD
for ii = 1:length(img_dir),
    imname = img_dir(ii).name;
    res1 = zeros(1,4);
    res1 = Demo_mySR(TD,imname);
    res(ii,:) = res1(:);
%     im = imread(fullfile(img_path, img_dir(ii).name));
%     nper_img(ii) = prod(size(im));%计�??��???????�?��
end
dict_path = ['newfig/res_' num2str(TD) '.mat' ];
save(dict_path,'res');
TD = 2;
res = zeros(img_num,4);%psnr Bic,Jmd/ScSR,GSCSR,GJMD
for ii = 1:length(img_dir),
    imname = img_dir(ii).name;
    res1 = zeros(1,4);
    res1 = Demo_mySR(TD,imname);
    res(ii,:) = res1;
%     im = imread(fullfile(img_path, img_dir(ii).name));
%     nper_img(ii) = prod(size(im));%计�??��???????�?��
end
dict_path = ['newfig/res_' num2str(TD) '.mat' ];
save(dict_path,'res');
TD = 3;
res = zeros(img_num,4);%psnr Bic,Jmd/ScSR,GSCSR,GJMD
for ii = 1:length(img_dir),
    imname = img_dir(ii).name;
    res1 = zeros(1,4);
    res1 = Demo_mySR(TD,imname);
    res(ii,:) = res1;
%     im = imread(fullfile(img_path, img_dir(ii).name));
%     nper_img(ii) = prod(size(im));%计�??��???????�?��
end
dict_path = ['newfig/res_' num2str(TD) '.mat' ];
save(dict_path,'res');
TD = 4;
res = zeros(img_num,4);%psnr Bic,Jmd/ScSR,GSCSR,GJMD
for ii = 1:length(img_dir),
    imname = img_dir(ii).name;
    res1 = zeros(1,4);
    res1 = Demo_mySR(TD,imname);
    res(ii,:) = res1;
%     im = imread(fullfile(img_path, img_dir(ii).name));
%     nper_img(ii) = prod(size(im));%计�??��???????�?��
end
dict_path = ['newfig/res_' num2str(TD) '.mat' ];
save(dict_path,'res');

TD = 5;
res = zeros(img_num,4);%psnr Bic,Jmd/ScSR,GSCSR,GJMD
for ii = 1:length(img_dir),
    imname = img_dir(ii).name;
    res1 = zeros(1,4);
    res1 = Demo_mySR(TD,imname);
    res(ii,:) = res1;
%     im = imread(fullfile(img_path, img_dir(ii).name));
%     nper_img(ii) = prod(size(im));%计�??��???????�?��
end
dict_path = ['newfig/TDres_' num2str(TD) '.mat' ];
save(dict_path,'res');

res = zeros(img_num,1);
for ii = 1:length(img_dir),
    imname = img_dir(ii).name;
    res1 = Jmd(imname);
    res(ii,1) = res1;
%     im = imread(fullfile(img_path, img_dir(ii).name));
%     nper_img(ii) = prod(size(im));%计�??��???????�?��
end
dict_path = ['newfig/res_jmd' '.mat' ];
save(dict_path,'res');