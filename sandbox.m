clear all
close all

load('X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\405B_2.mat');
load('X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\405B_2mnr.mat');


[HE_reg, nnr, F_rgb ] = select_common_points( im1a, imHE, mprob );


M = mean( [ range(yl)/range(yr) range(xl)/range(xr) ] );
HE_sm = imresize(imHE.I, M);

imIn = rgb2gray(HE_sm);
Ps = [xl' yl'];
Pd = [xr' yr'].*M;

imOut = warpItFun(rgb2gray(HE_sm), [xl' yl'], [xr' yr'].*M);

HE_reg = cat(3,imOut,imOut,imOut);

addpath('X:\#6 - Code\Matlab Code\TOOLBOXES\CPD2\')




% 
% L = [linspace(20,100,128) linspace(100,20,128)]';
% 
% R = 1 - (1 ./ (1 + exp(-1 .* (linspace(-10,5,256)))));
% G = (1 ./ (1 + exp(-1 .* (linspace(-5,10,256)))));
% B =  zeros(256,1);
% plot(1:256,R,1:256,G,1:256,B); xlim([0 256])
% 
% Lab = applycform(([R' G' (R'+G')-1]),makecform('srgb2lab'));
% Lab(:,1) = (R+G).*54.1;
% 
% rgb = applycform(Lab,makecform('lab2srgb'));




bot = imread('X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 5 (F501-F512) DEX\F505\Pearl WB\WP_20151113_13_52_38_Pro (2).jpg');
top_t = bfopen('X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 5 (F501-F512) DEX\F505\Pearl WB\F505 whole brain\0065622_01\0065622_01_800.TIF');
top = double(top_t{1}{1});
bot_rz = imresize(bot,size(top));
imagesc(top)

I1 = imcrop(top, [617 347 500 400]);
I2 = imcrop(bot_rz, [600 300 500 400]);

I1_n = nm(I1); I1_n = (I1_n.* 0.9).^1.1; I1_n(1,1) = 1;
% imagesc(5.*cat(3,zeros(size(I1)),I1,zeros(size(I1))) + im2double(I2)); axis image

options = fuseset( 'model', 'logistic', 'cmap_type', 'sRGB', 'im2_type', 'sRGB', ...
                                'L', 0.8, 'x0', 0.4, 'k', 10);
F = fuse_images( nm(I1_n), I2, goodmap('cyans'), options );

figure; imagesc(F); axis image; axis off; colorbar; colormap(goodmap('cube1'))

