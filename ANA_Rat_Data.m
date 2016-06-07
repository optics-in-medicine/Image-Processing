%Adjusted version of LOOP_CORREGISTRATION.m
%Uses one fluorescent image (800 nm) and one H&E (or similar) image

%BEFORE RUNNING CODE: Make sure to edit the code based on the orientation
%of the two images

%Upload text file specifying location of the two images
[filename, textpath] = uigetfile([cd filesep '*.txt']);

% load the filepaths for all the available files.
% textpath = 'X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\Corregistered Data Files\File Path Info\';
% filename = '404A1.txt';

fileID = fopen([textpath filename],'r');

% determine how many different paths are saved in the text file.
i=1;
tline = fgetl(fileID);
while ischar(tline)
pth{i} = tline;
tline = fgetl(fileID);
i = i + 1;
end
fclose(fileID);

if (length(pth) >= 2) % has to be at least 2 images to corregister - H&E and fluorescent
he_path = 'X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\F98 EGFR_pos B04\hico\'
imHE.fname = filename(1:end-4);


%%%If the H&E image is rotated compared to the fluor. image, 
%%%uncomment the first line in the following chunk of text (and comment 2nd). 
%%%Otherwise, uncomment the second line and comment the first. 
%%%When using Ana's data, the images should already be aligned (uncomment 2nd line)

imHE.I = fliplr( flipud(imread([pth{2} '.TIF'])) );
%imHE.I = (imread([pth{2} '.TIF']));
h2 = imagesc(imHE.I); axis image; axis off; drawnow;

im1 = load_images(pth{1});

close all;
 
[ im1a,im2a ] = imalign( im1, im1, 0 );
im_blend = cat( 3,4.* ( nm(im2a.I_a)) .* 0.5,zeros(size(im1a.I_a)), 5.* nm(im1a.I_a) );
imagesc(im_blend); axis image; axis off;

[ ~, im3a ] = imalign( im1, im1, 0 );

% create 3-channel fusion.

im_blend = cat(3, 2.2.*nm(im3a.I_a), 1.8.* nm(im1a.I_a), 0.7 .* nm(im2a.I_a).^1.5);
imagesc(im_blend); axis image; axis off
 end

%% Load Histopathology Data
figure('Name','Draw line spanning midline'); imagesc(im1a.I_a);
l1 = imline;
p = wait(l1);
len1 = sqrt((p(1,1) - p(1,2))^2 + ( p(2,1) - p(2,2) )^2);
close;

figure('Name','Draw line spanning midline'); imagesc(imHE.I);
l2 = imline;
p = wait(l2);
len2 = sqrt((p(1,1) - p(1,2))^2 + ( p(2,1) - p(2,2) )^2);
close;

M = len1 ./ len2;

imHE.I_s = imresize(imHE.I,M);

figure('Name','Draw a circle around the slice');
h2 = imagesc(imHE.I_s); axis image; axis off; drawnow;

% circle out the ROI
hFH = imfreehand();
pos = hFH.getPosition();

binaryImage = hFH.createMask();

% attempt to resize based on im1a.
imHE.I_c = cropfromfh( imHE.I_s, pos );

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

[O_trans,Spacing]=point_registration(size(imHE.I_c),[nnr.xl(:) nnr.yl(:)],[nnr.xr(:) nnr.yr(:)],options);
% transform H&E image.
HE_reg = bspline_transform(O_trans,imHE.I_c,Spacing,3);
imagesc(rgb2gray(HE_reg(1:225,1:352,:))+0.005.*im1a.I_a(1:225,1:352))


imHE.I_a = HE_reg;
figure; imagesc(imHE.I_a);

%% Draw user masks based on H&E
figure('color','white');
% imagesc( imHE.I_a ); axis image;
imagesc( im1a.I_a ); axis image;
title('Outline Tumor')

% Circle the tumor
hFH = imfreehand();
pos = hFH.getPosition();
mask.tumor = hFH.createMask();    

% Circle the section outline
% imagesc( imHE.I_a); axis image
imagesc( im1a.I_a); axis image
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
