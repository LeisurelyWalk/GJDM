function [hf] = presetimage(f)
addpath('Training');
addpath('RegularizedSC');
load par;
p = par;
p.LR = f;
p.scale = 1;
p.B = Set_blur_matrix(p);
hf = Superresolution(f,p);