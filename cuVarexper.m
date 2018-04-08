function res=cuVarexper(imh,iml)
addpath('Data/Testing');
addpath('testTD');
addpath('Results');
addpath('newfig');
newpath = 'newfig';
for alfa=1:20
    for lambda=1:20
        al=alfa/10;
        la=lambda/100;
        im_h = GualAsent(imh,al,iml,la);
        
        [nrow, ncol] = size(im_h);
        
        
    end
end