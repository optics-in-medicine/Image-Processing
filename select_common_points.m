function [HE_reg, nnr, F ] = select_common_points( im1a, imHE, xmap, no_pts )
%SELECT_COMMON_POINTS Summary of this function goes here
%   Detailed explanation goes here

% load dependencies.
addpath('X:\#6 - Code\Matlab Code\# # DATAVIS\Adaptive Colormap Evaluation');
rmpath('X:\#6 - Code\Matlab Code\TOOLBOXES\nonrigid_version23\')
addpath(genpath('X:\#6 - Code\Matlab Code\# # Coregistration\Modified Nonrigid Registration\'));

% Option of overlaying a different map than just intensity of im1a.
if nargin == 2
    xmap = im1a.I_a;
    no_pts = 8;
elseif nargin == 3
    no_pts = 8;
    if isempty(xmap)
        xmap = im1a.I_a;
    end
elseif nargin == 4
    if isempty(xmap)
        xmap = im1a.I_a;
    end
end

if ~isfield(imHE,'I_c')
    imHE.I_c = imHE.I;
end

% images to select control points.
figure('color','white'); 
h1 = subplot(1,2,1)
imagesc(nm(im1a.I_a)); caxis([0 0.5]); axis image; axis off
h2 = subplot(1,2,2)
imagesc(imHE.I_c); axis image; axis off
hold on;

j=1;
for i = 1:(no_pts*2)
    
    if mod(i,2) % odd  
        
        axes(h1)
        set(h2,'Visible','off')
        
        title(['Select point no. ' int2str(j) ' from left'])
        [xl(j),yl(j)]=ginput(1)
        
        % convert to coordinates.
        hold on;
        plot(xl(1:j),yl(1:j),'k+');
        plot(xl(j),yl(j),'r+');
        str2 = [int2str(j)];
        text(xl(j),yl(j),str2,'BackgroundColor',[1 1 1])

        title('')
        set(h2,'Visible','on')
        
    % set to correct axes
    else %
        axes(h2)
        title(['Select point no. ' int2str(j) ' from right'])
        [xr(j),yr(j)]=ginput(1);
        
        % convert to coordinates.
        hold on;
        plot(xr(1:j),yr(1:j),'k+');
        plot(xr(j),yr(j),'r+');
        str2 = [int2str(j)];
        text(xr(j),yr(j),str2,'BackgroundColor',[1 1 1])
        title('')
                j = j + 1;
    end
end
hold off
close(gcf)

BW = poly2mask(im1a.BN(:,1),im1a.BN(:,2),size(im1a.I_a,1),size(im1a.I_a,2));
[X, Y] = size(im1a.I_a);

% calculate density of control points.
for i = 1:X
    for j = 1:Y
        cpd(i,j) = sum( sqrt( (yl-i).^2 + (xl-j).^2 ) );
    end
end

cpd(~BW) = nan;

% plot density map of control points (landparks).
subplot(1,3,1);
imagesc(nm(cpd)); colormap(goodmap('kryptonite'));
hold on;
plot(xl,yl,'r+');
plot(im1a.BN(:,1),im1a.BN(:,2),'b','LineWidth',2)        
hold off;
axis image
axis off
title('Density of control point information');

% convert the CPs to the proper image units.

% first roughly resize the images according to the variance of the points.
M = mean( [ range(yl)/range(yr) range(xl)/range(xr) ] );
HE_sm = imresize(imHE.I_c, M);
xr_adj = xr.*M; yr_adj = yr.*M;
HE_sm(size(im1a.I_a,1),size(im1a.I_a,2),3) = 0;

% then roughly transform the exterior of the 

% perform non-rigid registration based on control points/landmarks.
options.Verbose = true;
[O_trans,Spacing]=point_registration(size(HE_sm),[xl(:) yl(:)],[xr_adj(:) yr_adj(:)],options);

% transform H&E image.
HE_reg = bspline_transform(O_trans,HE_sm,Spacing,3);

% create fusion image of fluorescence map (im1a.I_a)
options = fuseset( 'model', 'logistic', 'cmap_type', 'sRGB', 'im2_type', 'sRGB', 'L', 0.8, 'x0', 0.1, 'k', 20);
F = fuse_images( nm(xmap .* BW), HE_reg(1:size(im1a.I_a,1),1:size(im1a.I_a,2),:), goodmap('cube1'), options );

blk_m = logical( F(:,:,1) <= 0.2 & F(:,:,2) <= 0.2 & F(:,:,3) <= 0.2 );
blk_md = imdilate(blk_m,strel('square',5));

F_mask = F;
F_mask(cat(3,blk_md,blk_md,blk_md)) = 1;

    subplot(1,3,2);
    imagesc(F_mask); axis image; axis off;


% create image showing grids.
Igrid = make_grid_image(Spacing, size(HE_sm));
[X,Y] = meshgrid(1:size(HE_sm,1),1:size(HE_sm,2));

Ireg = bspline_transform( O_trans, Igrid, Spacing, 3 );
X_dots = X(Ireg>0.8);
Y_dots = Y(Ireg>0.8);


HE_r = HE_reg(:,:,1); HE_r(Ireg>0.7) = 0; HE_r(blk_md) = 1;
HE_g = HE_reg(:,:,2); HE_g(Ireg>0.7) = 0.4; HE_g(blk_md) = 1;
HE_b = HE_reg(:,:,3); HE_b(Ireg>0.7) = 0; HE_b(blk_md) = 1;
HE_grid = cat(3,HE_r, HE_g, HE_b);

       subplot(1,3,3);
        imagesc(HE_grid);
        axis image; axis off
% save all registration info into structure, nrr
nnr.xl = xl;
nnr.xr = xr;
nnr.yl = yl;
nnr.yr = yr;
nnr.HE_grid = HE_grid;
nnr.O_trans = O_trans;
nnr.Spacing = Spacing;

% hold on;for i=1:10, plot([y1(i) y2(i)],[x1(i) x2(i)],'b'); end

%restore references to the global TOOLBOX.
addpath('X:\#6 - Code\Matlab Code\TOOLBOXES\nonrigid_version23\')
rmpath('X:\#6 - Code\Matlab Code\# # Coregistration\Modified Nonrigid Registration\');

