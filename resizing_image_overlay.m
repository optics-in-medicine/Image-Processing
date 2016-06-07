 %% Load Histopathology Data

 im1a = imread('beige.tif');
 imHE = imread('purple.tif');
 
    figure('Name','Draw line spanning midline'); imagesc(im1a);
    l1 = imline;
    p = wait(l1);
    len1 = sqrt((p(1,1) - p(1,2))^2 + ( p(2,1) - p(2,2) )^2);
    close;

    figure('Name','Draw line spanning midline'); imagesc(imHE);
    l2 = imline;
    p = wait(l2);
    len2 = sqrt((p(1,1) - p(1,2))^2 + ( p(2,1) - p(2,2) )^2);
    close;

    M = len1 ./ len2;

    imHE_s = imresize(imHE,M);

    figure('Name','Draw a circle around the slice');
    h2 = imagesc(imHE_s); axis image; axis off; drawnow;

    % circle out the ROI
    hFH = imfreehand();
    pos = hFH.getPosition();

    binaryImage = hFH.createMask();

    % attempt to resize based on im1a.
    imHE_c = cropfromfh( imHE_s, pos );

    % attempt to corregister H&E sections
    [HE_reg, nnr, F ] = select_common_points( im1a, imHE, '', 8 );

    close all;
    % save('401B_test_case.mat');
    %% recplace below with

    %%%
    [imHE, smpl_pts] = multi_step_tissue_reg(im1a,im2a,im3a, imHE);
    %%%

    stats = sample_from_points( smpl_pts, im1a, im2a, im3a );

    save([textpath imHE.fname '.mat'])

%end


% 
% 
% 
% 
% 
% run optimization to minimize energy, difference in length distribution,
% and difference in location of points.
options.Verbose = true;

[O_trans,Spacing]=point_registration(size(imHE_c),[nnr.xl(:) nnr.yl(:)],[nnr.xr(:) nnr.yr(:)],options);
% transform H&E image.
HE_reg = bspline_transform(O_trans,imHE_c,Spacing,3);
imagesc(rgb2gray(HE_reg(1:225,1:352,:))+0.005.*im1a(1:225,1:352))


imHE_a = HE_reg;
figure; imagesc(imHE_a);

%% Draw user masks based on H&E
figure('color','white');
% imagesc( imHE_a ); axis image;
imagesc( im1a ); axis image;
title('Outline Tumor')

% Circle the tumor
hFH = imfreehand();
pos = hFH.getPosition();
mask.tumor = hFH.createMask();    

% Circle the section outline
% imagesc( imHE_a); axis image
imagesc( im1a); axis image
title('Outline Section')
hFH = imfreehand();
pos = hFH.getPosition();
mask.section = hFH.createMask();    
close(gcf);

mask.brain = logical( mask.section .* imcomplement(mask.tumor) );


% Load IHC Data
ihc1
ihc2
ihc3



save([pathname 'processed\' he_file 'wrp.mat'],'im1a','im2a','im3a','imHE','mask','im_blend');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 