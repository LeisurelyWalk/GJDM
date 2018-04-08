function  Demo_myLR(path,imname)
    im = imread(fullfile(path,imname));
    im_l = imresize(im,1/2);
    bname = ['lr_' '_' imname];
    newpath = 'newres';
    imwrite(im_l,fullfile(newpath,bname));