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
T = [0 0];
[ im1a, im2a ] = imalign( im1, im2, 0, T );

im_blend = cat( 3, ( nm(im2a.I_a)) .* 1.1,zeros(size(im1a.I_a)), 1.6*nm(im1a.I_a) );
imagesc(im_blend);

im3 = load_images(im1.pathname);

T = [-10 6];
[ ~, im3a ] = imalign( im1, im3, 0, T );

% blend

im_blend = cat(3, 0.9.*nm(im2a.I_a), 2 .* nm(im1a.I_a).^0.9, 0.9 .* nm(im3a.I_a));
imagesc(im_blend);

% MRI
im4 = load_images(im1.pathname);




% OUTLINE TUMOR BY USING H&E AS GUIDE (EVENTUALLY WILL CORREGISTER LIKE
% OTHERS)

imagesc(nm( im1a.I_a )); caxis([0 0.2]); axis image;

% circle out the ROI
hFH = imfreehand();
pos = hFH.getPosition();
tumor_mask = hFH.createMask();


imagesc(nm( im1a.I_a )); caxis([0 0.2]); axis image;

% circle out the ROI
hFH = imfreehand();
pos = hFH.getPosition();
brain_mask = hFH.createMask();    

mf = [2 2]; 
P1 = medfilt2(im1a.I_a, mf);
P2 = medfilt2(im2a.I_a, mf);
P3 = medfilt2(im3a.I_a, mf);

figure('color','white')';
subplot(2,2,1)

% build data matrix
X1 = [P1(tumor_mask) P2(tumor_mask) P3(tumor_mask)];
Y1 = ones(sum(tumor_mask(:)),1);

X = [X1; [P1(brain_mask) P2(brain_mask) P3(brain_mask)]];
Y = [Y1; zeros(sum(brain_mask(:)),1)];

u_Y = unique(Y);    
m_vect = zeros( length(u_Y), 3 );
for i = 1:length(u_Y)
    m_vect(i,:) = [median(P1(Y==u_Y(i))) median(P2(Y==u_Y(i))) median(P3(Y==u_Y(i)))];
end

% plot scatter

plot3(X(Y==0,1),X(Y==0,2),X(Y==0,3),'g.', 'MarkerSize',2)
hold on
plot3(X(Y==1,1),X(Y==1,2),X(Y==1,3),'r.', 'MarkerSize',2)

plot3(m_vect(1,1),m_vect(1,2),m_vect(1,3),'ko', 'MarkerSize',10)
plot3(m_vect(2,1),m_vect(2,2),m_vect(2,3),'ko', 'MarkerSize',10)

hold off




% TRY AN UNSUPERVISED CLASSIFICATION!

[X,Y] = size(im1a.I_a);
data = [nm( im1a.I_a(:) ) nm( im2a.I_a(:) ) nm( im3a.I_a(:) )];
mask = imresize(im1.BW,0.84);

[cluster_idx cluster_center] = kmeans(data(mask==1), 3);

idx_f = zeros(X*Y,1);
idx_f(mask==1) = cluster_idx;

I_idx = reshape(idx_f,[X Y 1]);
figure; imagesc(I_idx)


ha = tight_subplot(2,3,[.01 .01],[.1 .1],[.1 .1]);
axes(ha(1)); imagesc(nm(im2a.I_a)); axis image; colormap(goodmap('cube1')); title('ABY029'); axis off; caxis ([0 0.2])
axes(ha(2)); imagesc(nm(im1a.I_a)); axis image; colormap(goodmap('cube1')); title('IRDye680'); axis off; caxis ([0 0.8])
axes(ha(3));  imagesc(nm(im3a.I_a)); axis image; colormap(goodmap('cube1')); title('PpIX'); axis off
axes(ha(4));  imagesc(im_blend); axis image; colormap(goodmap('cube1')); axis off
axes(ha(5));  imagesc(I_idx); title('kmeans'); axis image; axis off
axes(ha(6));  axis off; 
plot3(nm(im1a.I_a(mask==2)),nm(im2a.I_a(mask==2)),nm(im3a.I_a(mask==2)),'r.', ...
      nm(im1a.I_a(mask==1)),nm(im2a.I_a(mask==1)),nm(im3a.I_a(mask==1)),'b.'); 



% TRY A SUPERVISED CLASSIFICATION!






