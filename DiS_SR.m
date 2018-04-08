function DiS_SR(img_path,type)
%img_path = 'testTD';
%img_path = 'newtest';
%type = '*.jpg';
%type = '*.tiff';
img_dir = dir(fullfile(img_path, type));%?��???��??bmp
%dict_path = ['newres/img_dir' '.txt' ];
%save(dict_path,img_dir);
img_num = length(img_dir);%?��?�??
% nper_img = zeros(1, img_num);
TD = 1;%1��256��2��512��3��1024�� 4��2048
res = zeros(img_num,4);%psnr Bic,Jmd/ScSR,GSCSR,GJMD
times=zeros(1,3);
% for ii = 7:length(img_dir),
%     imname = img_dir(ii).name;
%     res1= zeros(1,4);
%     [res1 ,times]= Demo_mySR(times,img_path,TD,imname);
%     res(ii,:) = res1(:);
% %     im = imread(fullfile(img_path, img_dir(ii).name));
% %     nper_img(ii) = prod(size(im));%计�??��???????�?��
% end
% dict_path = ['newres/dz_s2_res_256_' num2str(TD) '.mat' ];
% save(dict_path,'res');
% TD = 2;
% res = zeros(img_num,4);%psnr Bic,Jmd/ScSR,GSCSR,GJMD
% for ii = 7:length(img_dir),
%     imname = img_dir(ii).name;
%     res1 = zeros(1,4);
%     [res1 ,times]= Demo_mySR(times,img_path,TD,imname);
%     res(ii,:) = res1;
% %     im = imread(fullfile(img_path, img_dir(ii).name));
% %     nper_img(ii) = prod(size(im));%计�??��???????�?��
% end
% dict_path = ['newres/dz_s2_res_512_' num2str(TD) '.mat' ];
% save(dict_path,'res');
TD = 2;
res = zeros(img_num,4);%psnr Bic,Jmd/ScSR,GSCSR,GJMD
times= zeros(img_num,3);%psnr Bic,Jmd/ScSR,GJMD
for ii = 1:length(img_dir),
    
    imname = img_dir(ii).name;
    if imname(1)=='.'
        continue
    end
    res1 = zeros(1,4);
    [res1 ,time1]= Demo_mySR(times,img_path,TD,imname);
    res(ii,:) = res1;
    times(ii,:)=time1;
%     im = imread(fullfile(img_path, img_dir(ii).name));
%     nper_img(ii) = prod(size(im));%计�??��???????�?��
end
dict_path = ['newres/dz_s2_res_1024_' num2str(TD) '.mat' ];
save(dict_path,'res');
dict_path = ['newres/timesSc,SJ,BC,1024_' num2str(TD) type '.mat' ];
save(dict_path,'times');
% TD = 4;
% res = zeros(img_num,4);%psnr Bic,Jmd/ScSR,GSCSR,GJMD
% for ii = 7:length(img_dir),
%     imname = img_dir(ii).name;
%     res1 = zeros(1,4);
%     [res1 ,times]= Demo_mySR(times,img_path,TD,imname);
%     res(ii,:) = res1;
% %     im = imread(fullfile(img_path, img_dir(ii).name));
% %     nper_img(ii) = prod(size(im));%计�??��???????�?��
% end
% dict_path = ['newres/dz_s2_res_2048_' num2str(TD) '.mat' ];
% save(dict_path,'res');

