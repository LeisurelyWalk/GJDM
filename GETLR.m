img_path = 'newtest';
type = '*.jpg';
%type = '*.tiff';
img_dir = dir(fullfile(img_path, type));%?��???��??bmp

img_num = length(img_dir);%?��?�??

res = zeros(img_num,4);%psnr Bic,Jmd/ScSR,GSCSR,GJMD
for ii = 7:length(img_dir),
    imname = img_dir(ii).name; 
    Demo_myLR(img_path,imname);
end