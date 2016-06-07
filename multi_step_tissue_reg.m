function [imHE, smpl_pts] = multi_step_tissue_reg(im1a, im2a, im3a, imHE)
%MULTI_STEP_TISSUE_REG Summary of this function goes here
%   Detailed explanation goes here

addpath(['.' filesep 'icp_finite']);

% 1. Need to resize the imHE image.
figure('Name','Draw line spanning midline');  imagesc(im1a.I_a);
l1 = imline; p = wait(l1);
len1 = sqrt((p(1,1) - p(1,2))^2 + ( p(2,1) - p(2,2) )^2);
close;

figure('Name','Draw line spanning midline');  imagesc(imHE.I);
l2 = imline; p = wait(l2);
len2 = sqrt((p(1,1) - p(1,2))^2 + ( p(2,1) - p(2,2) )^2);
close;

M = len1 ./ len2;
I_s = imresize(imHE.I,M);
close;

% circle out the ROI
figure('Name','Draw boundary around slice'); 
h2 = imagesc(I_s); axis image; axis off; drawnow;
hFH = imfreehand();
pos = hFH.getPosition();
bin_I = hFH.createMask();

% attempt to resize based on im1a.
[imHE.I_c, pos_c] = cropfromfh( I_s, pos );
imHE.BW = logical( cropfromfh( double(bin_I), pos ) );
imHE.BN = pos_c;

% tighten the boundary on imHE
imHE = tighten_boundary( imHE );
close;

% Trace tumors to figure out center of mass.
figure('color','white','Name','Trace around tumor (for center of mass)');
imagesc( imHE.I_c ); axis image;
% imagesc( im1a.I_a ); axis image;
title('Outline Tumor')

% Circle the tumor
hFH = imfreehand();
pos = hFH.getPosition();
mask.tumor = hFH.createMask();    
close;

t_com = mean(pos);
t_pos = pos;

% 2. Corregister boundaries.

[HE_reg, nnr, F ] = select_common_points( im1a, imHE, '', 8 );
clear HE_reg;

for i = 1:length(nnr.xl)
    d = sqrt( (imHE.BN(:,1)-nnr.xr(i)).^2 + (imHE.BN(:,2)-nnr.yr(i)).^2 );  
    min_d(i) = min(d);
    clear d;
end

ext_cps = logical( min_d < 15 );

tform = fitgeotrans([nnr.xr(ext_cps)' nnr.yr(ext_cps)'],[nnr.xl(ext_cps)' nnr.yl(ext_cps)'], 'affine');
B_t = imwarp(imHE.I_c,tform,'FillValues',squeeze(imHE.I_c(1,1,:)),'OutputView',imref2d(size(im1a.I_a)));
BW_t = logical( imwarp(double(imHE.BW),tform,'FillValues',0,'OutputView',imref2d(size(im1a.I_a))) );
tum_t = logical( imwarp(double(mask.tumor),tform,'FillValues',0,'OutputView',imref2d(size(im1a.I_a))) );

[t_X,t_Y] = tformfwd(maketform('affine',tform.T), nnr.xr, nnr.yr);

[t_BN(:,1) t_BN(:,2)] = tformfwd(maketform('affine',tform.T), imHE.BN(:,1), imHE.BN(:,2));
[t_pos(:,1) t_pos(:,2)] = tformfwd(maketform('affine',tform.T), pos(:,1), pos(:,2));
t_com = mean(t_pos);
%transform tumor region and then use for c of mass.


%t_BN = [imHE.BN + repmat([-0.2 2].*tform.T(3,1:2),length(imHE.BN),1); ] * tform.T(1:2,1:2)  

[r0,c0,~] = size(imHE.I_c);
[r1,c1,~] = size(B_t);

% B = B_t(floor(abs(c0-c1)/2)+1:end-floor(abs(c0-c1)/2),ceil(abs(r0-r1)/2):end-floor(abs(r0-r1)/2),:);
B=B_t;

figure;
subplot(2,2,1);
imagesc(imHE.I_c); axis image; axis off; hold on; plot(nnr.xr(ext_cps),nnr.yr(ext_cps)); hold off
subplot(2,2,2);
imagesc(im1a.I_a); axis image; axis off; hold on; plot(nnr.xl(ext_cps),nnr.yl(ext_cps),'r','LineWidth',2); hold off
subplot(2,2,3);
imagesc(B); axis image; axis off; hold on; plot(nnr.xr(ext_cps),nnr.yr(ext_cps),'.--b','LineWidth',2); plot(nnr.xl(ext_cps),nnr.yl(ext_cps),'.--r','LineWidth',2); hold off
subplot(2,2,4);
imagesc(B); axis image; axis off; hold on; plot(t_X,t_Y,'.b','LineWidth',2);

% use ray projection to determine the boundary CPs
no_rays = 18;
th = linspace(0,pi,no_rays);
k = linspace(-1000,1000,1000);
cp_xr = []; cp_xl = []; cp_yr = []; cp_yl = [];

for i = 1:no_rays
    A = [t_com(1) + cos( th(i) ) t_com(2) + sin( th(i) )]
    ry(:,1) = t_com(1) + k .* cos( th(i) );
    ry(:,2) = t_com(2) + k .* sin( th(i) );
%      plot(ry(:,1),ry(:,2),'-'); hold on;
    
    % find two closest points to the line.
    for j = 1:length(k);
        d1(j,:) = sqrt( (im1a.BN(:,1)-ry(j,1)).^2 + (im1a.BN(:,2)-ry(j,2)).^2 );
        d2(j,:) = sqrt( (t_BN(:,1)-ry(j,1)).^2 + (t_BN(:,2)-ry(j,2)).^2 );
    end
    [ ~, I_d1 ] = min_two( d1 );
    [ ~, I_d2 ] = min_two( d2 );
    
    cp_xr = [cp_xr; ry(I_d2(1,1),1); ry(I_d2(2,1),1)];
    cp_yr = [cp_yr; ry(I_d2(1,1),2); ry(I_d2(2,1),2)];
    
    cp_xl = [cp_xl; ry(I_d1(1,1),1); ry(I_d1(2,1),1)];
    cp_yl = [cp_yl; ry(I_d1(1,1),2); ry(I_d1(2,1),2)];
    
%     plot(A(1),A(2),'o'); hold on;
end

% reject any control points that aren't close to each other
d_cp = sqrt( (cp_xr-cp_xl).^2 + (cp_yr-cp_yl).^2 );  
il_i = reshape([2:2:length(d_cp); 1:2:length(d_cp)],1,[]); 
d_cp_f = sqrt( (cp_xr-cp_xl(il_i)).^2 + (cp_yr-cp_yl(il_i)).^2 );     

if sum(d_cp_f)*1.25 < sum(d_cp) 
    cp_xl = cp_xl(il_i);
    cp_yl = cp_yl(il_i);
end
% done!

d_cp = sqrt( (cp_xr-cp_xl).^2 + (cp_yr-cp_yl).^2 );  
% select out only below thresh
thresh = 60;
cp_xr = cp_xr(d_cp < thresh);
cp_xl = cp_xl(d_cp < thresh);
cp_yr = cp_yr(d_cp < thresh);
cp_yl = cp_yl(d_cp < thresh);

figure('Name','Refining the boundary control points');
imagesc(B); axis image; axis off; 
hold on;
plot(cp_xr,cp_yr,'ro',cp_xl,cp_yl,'bo')
plot(im1a.BN(:,1),im1a.BN(:,2),'b-','LineWidth',1)
plot(t_BN(:,1),t_BN(:,2),'r-','LineWidth',1)

% 3. Do the whole pressure thing!
% Still plan to do this!

% tform = fitgeotrans([cp_xr cp_yr; [t_X' t_Y']],[cp_xl cp_yl; [nnr.xr' nnr.yr']], 'affine');
tform = fitgeotrans([cp_xr cp_yr],[cp_xl cp_yl], 'affine');
B_tt = imwarp(B,tform,'FillValues',squeeze(B(1,1,:)),'OutputView',imref2d(size(B)));
BW_tt = logical( imwarp(double(BW_t),tform,'FillValues',0,'OutputView',imref2d(size(BW_t))) );
tum_tt = logical( imwarp(double(tum_t),tform,'FillValues',0,'OutputView',imref2d(size(tum_t))) );
[cp_xr_t,cp_yr_t] = tformfwd(maketform('affine',tform.T), cp_xr, cp_yr);
[tt_BN(:,1),tt_BN(:,2)] = tformfwd(maketform('affine',tform.T), t_BN(:,1),t_BN(:,2));

figure;
imagesc(B_tt); axis image; axis off; 
hold on;
plot(cp_xr_t,cp_yr_t,'ro',cp_xl,cp_yl,'bo')
plot(im1a.BN(:,1),im1a.BN(:,2),'b--','LineWidth',1)
plot(tt_BN(:,1),tt_BN(:,2),'r-','LineWidth',1)

% perform non-rigid registration based on control points/landmarks.
options.Verbose = true;
% [O_trans,Spacing]=point_registration(size(B_tt),[cp_xr cp_yr; [t_X' t_Y']],[cp_xl cp_yl; [nnr.xr' nnr.yr']],options);
 [O_trans,Spacing]=point_registration(size(B),[cp_xr_t cp_yr_t],[cp_xl cp_yl],options);
% [O_trans,Spacing]=point_registration(size(B),[cp_xl cp_yl; [t_X' t_Y']],[cp_xr cp_yr; [nnr.xr' nnr.yr']],options);
% transform H&E image.

HE_reg = bspline_transform(O_trans,B_tt,Spacing,0);
BW_reg_t = bspline_transform(O_trans,BW_tt,Spacing,3);
BW_reg = bwareaopen(logical(BW_reg_t >0.5), 100);
tum_reg_t = bspline_transform(O_trans,tum_tt,Spacing,3);
tum_reg = bwareaopen(logical(tum_reg_t >0.5), 100);

% create fusion image of fluorescence map (im1a.I_a)
options = fuseset( 'model', 'logistic', 'cmap_type', 'sRGB', 'im2_type', 'sRGB', ...
                                'L', 0.6, 'x0', 0.5, 'k', 10);
F = fuse_images( nm(im2a.I_a), HE_reg(1:size(im2a.I_a,1),1:size(im2a.I_a,2),:), goodmap('cube1'), options );
imagesc(F);
hold on;

% create a margin region.
nx = 10;
mar_1 = imerode( tum_reg, strel('disk', 10, 8));

class_reg = (BW_reg + tum_reg + mar_1) .* imerode(BW_reg, strel('disk',10,8));
% randomly sample spots inside tissue
smpl_pts = [];
rX = randi([1 size(BW_reg,1)],1000,1);
rY = randi([1 size(BW_reg,2)],1000,1);
counter = [0 0 0];
j = 1;
for i = 1:1000
    switch (class_reg(rX(i),rY(i)))
        
        case 1
        if counter(1) <= nx
            smpl_pts(j,:) = [rY(i) rX(i) (class_reg(rX(i),rY(i)))];
            j = j + 1;
            counter(1) = counter(1) + 1;
        end
        
        case 2
        if counter(2) <= nx
            smpl_pts(j,:) = [rY(i) rX(i) (class_reg(rX(i),rY(i)))];
            j = j + 1;
            counter(2) = counter(2) + 1;
        end
        
        case 3
        if counter(3) <= nx
            smpl_pts(j,:) = [rY(i) rX(i) (class_reg(rX(i),rY(i)))];
            j = j + 1;
            counter(3) = counter(3) + 1;
        end
    end
end

% save('temp.mat');


