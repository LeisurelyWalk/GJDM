function im_h= MyBpC(imh,iml)
addpath('Data/Testing');
addpath('testTD');
addpath('Results');
addpath('newfig');
newpath = 'newfig';
alfa = 1.0;
lambda = 0.01;
im_h = GualAsent(imh,alfa,iml,lambda);
