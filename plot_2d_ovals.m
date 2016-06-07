close all
clear all


cd('X:\#6 - Code\Matlab Code\# # Coregistration\');
addpath('X:\#6 - Code\Matlab Code\TOOLBOXES\nonrigid_version23\');

% Bottom image
pathname = 'X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\';

% load first image - res of odyssey usually 42.
im1 = load_images(pathname);

% load second image
im2 = load_images(im1.pathname, 3);

close all;
T = [-5 0];
[ im1a, im2a ] = imalign( im1, im2, 0, T );

im_blend = cat( 3, ( nm(im2a.I_a)) .* 0.3,zeros(size(im1a.I_a)), 1.9 .* nm(im1a.I_a) );
imagesc(im_blend);


he_path = 'X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\F98 EGFR_pos B04\hico\'
he_file = '407A_2'



figure('color','white');
subplot(1,2,1)
imagesc(imread([he_path he_file '.TIF'])); axis image;

subplot(1,2,2)
imagesc(nm( im1a.I_a )); caxis([0 0.6]); axis image;
title('healthy')
% circle out the ROI
hFH = imfreehand();
pos = hFH.getPosition();
brain_mask = hFH.createMask();    

subplot(1,2,2)
imagesc(nm( im1a.I_a )); caxis([0 0.6]); axis image;

% circle out the ROI
hFH = imfreehand();
pos = hFH.getPosition();
tumor_mask = hFH.createMask();


mf = [2 2]; 

% unnormalized predictors?
P1 = medfilt2(im1a.I_a, mf);
P2 = medfilt2(im2a.I_a, mf);

% normalize predictors first?
P1 = medfilt2(nm(im1a.I_a), mf);
P2 = medfilt2(nm(im2a.I_a), mf);

% build data matrix
X1 = [P1(tumor_mask) P2(tumor_mask)];
Y1 = ones(sum(tumor_mask(:)),1);

X = [X1; [P1(brain_mask) P2(brain_mask)]];
Y = [Y1; zeros(sum(brain_mask(:)),1)];

svmStruct = svmtrain(X,Y,'ShowPlot',true)

% fit a logit function to the date.

[mp, prob1, r1, t01] = get_logistic( P1(tumor_mask), P1(brain_mask), 100 );
[mp, prob2, r2, t02] = get_logistic( P2(tumor_mask), P2(brain_mask), 100 );


v = linspace(0.01,1,101);
lb = v(1:100)';
ub = v(2:101)';
mp = median([lb ub]')';

for i = 1:100
    p1t = P1(tumor_mask);
    p1n = P1(brain_mask);
    
    prob1(i) = sum( (p1t >= lb(i)) & (p1t < ub(i)) ) / ...
        (sum( (p1t >= lb(i)) & (p1t < ub(i)) ) + sum( (p1n >= lb(i)) & (p1n < ub(i)) ));

end

paramguess=[100,0.15];
params=fminsearch(@penalty_logistic,paramguess,[],mp(~isnan(prob1)),prob1(~isnan(prob1)));


[B, dev, stats] = mnrfit(X,Y+1);





