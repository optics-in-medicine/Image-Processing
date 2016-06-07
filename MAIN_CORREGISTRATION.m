%MAIN_LOGISTIC_ANALYSIS Load in multiparametric data and H&E sections,
%create masks and perform logistic analysis to generate tumor probability
%map.


close all
clear all

cd('X:\#6 - Code\Matlab Code\# # Coregistration\');
addpath('X:\#6 - Code\Matlab Code\TOOLBOXES\nonrigid_version23\');

% Bottom image
pathname = 'X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\';

%% Image load and corregister

% load first image (ABY029) - res of odyssey usually 42.
im1 = load_images(pathname);

% load second image (PpIX)
im2 = load_images(im1.pathname, 0);

close all;
T = [0 -5];
[ im1a, im2a ] = imalign( im1, im2, 0, T );

im_blend = cat( 3, ( nm(im2a.I_a)) .* 0.93,zeros(size(im1a.I_a)), 1.9 .* nm(im1a.I_a) );
imagesc(im_blend);

% load third image (IRD680)
im3 = load_images(im1.pathname);

T = [30 -5];
[ ~, im3a ] = imalign( im1, im3, 0, T );

% create 3-channel fusion.

im_blend = cat(3, 2.0.*nm(im3a.I_a), 1.7 .* nm(im1a.I_a).^0.8, 0.7 .* nm(im2a.I_a));
imagesc(im_blend); axis image; axis off

%% Load Histopathology Data

% YOU LEFT OFF HERE, CORREGISTER WITH HE 410A_2

he_path = 'X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\F98 EGFR_pos B04\hico\'
he_file = '406A_2';
imHE.fname = he_file
imHE.I = fliplr( flipud(imread([he_path he_file '.TIF'])) );
% imHE.I = ( flipud(imread([he_path he_file '.TIF'])) );
% imHE.I = (imread([he_path he_file '.TIF']));

%% Load IHC Data
ihc1
ihc2
ihc3

% attempt to corregister H&E sections
[HE_reg, nnr, F ] = select_common_points( im1a, imHE, '', 12 );

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


save([pathname 'processed\' he_file 'wrp.mat'],'im1a','im2a','im3a','imHE','mask','im_blend');



%% Get the logistic analsyis part
[stats, mprob] = MAIN_LOGISTIC_ANALYSIS( im1a, im2a, im3a, imHE, mask, 'panel' );
set(gcf,'Position',[488   113   890   649]);
saveas(gcf,[pathname 'processed\' he_file '.fig'])
saveas(gcf,[pathname 'processed\' he_file '.png'])

% warped image
options = fuseset( 'model', 'logistic', 'cmap_type', 'sRGB', 'im2_type', 'sRGB', ...
                                'L', 0.6, 'x0', 0.3, 'k', 20);
F = fuse_images( nm(mprob), imHE.I_a(1:size(im1a.I_a,1),1:size(im1a.I_a,2),:), goodmap('comet'), options );

figure; imagesc(F); axis image; axis off; colorbar; colormap(goodmap('comet'))


fnames = {'405B_2','405C_2','407A_2','407B_1','410B_2', '408B_1','406A_2'};

for i = 1:length(fnames)
% pick control points, number them after each click.
load(['X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\',...
    fnames{i} 'wrp.mat'])
plot3( median( im1a.I_a(mask.tumor) ), ...
       median( im2a.I_a(mask.tumor) ), ...
       median( im3a.I_a(mask.tumor) ), ...
      'ro','MarkerSize',10);
hold on;
plot3( median( im1a.I_a(mask.brain) ), ...
       median( im2a.I_a(mask.brain) ), ...
       median( im3a.I_a(mask.brain) ), ...
      'go','MarkerSize',10);

% plot3( im1a.I_a(mask.tumor), ...
%        im2a.I_a(mask.tumor), ...
%        im3a.I_a(mask.tumor), ...
%       'r.','MarkerSize',1);
% plot3( im1a.I_a(mask.brain), ...
%        im2a.I_a(mask.brain), ...
%        im3a.I_a(mask.brain), ...
%       'g.','MarkerSize',1);
end

fnames = {'405B_2','405C_2','407B_1','410B_2', '408B_1','406A_2'};

h_s = figure;

for i = 1:length(fnames);
     load(['X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\',...
     fnames{i} 'wrp.mat']);

     [c, c_nm, bnds] = create_profile( im1a, im2a, im3a, mask )
     
     subplot(3,2,i);
    plot(1:400,c_nm(:,1),'g',1:400,c_nm(:,2),'b',1:400,c_nm(:,3),'r',bnds(1)*[1 1],[0 1],'k--',bnds(2)*[1 1],[0 1],'k--')     
    
end

for i = 1:length(fnames)
    
end

% 
% h1 = figure;
% h2 = figure;
% hold on;
% for i = 1:length(fnames)
%      load(['X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\',...
%      fnames{i} 'wrp.mat']);
% %      plot( median( im1a.I_a(mask.tumor) ), median( im2a.I_a(mask.tumor) ), 'ro','MarkerSize',10 );
% %      hold on
% %      plot( median( im1a.I_a(mask.brain) ), median( im2a.I_a(mask.brain) ), 'go','MarkerSize',10 );
%      figure(h1)
%      hold on
%      errorbarxy( median( im1a.I_a(mask.tumor) ), median( im2a.I_a(mask.tumor)), ...
%          abs(median( im1a.I_a(mask.tumor) ) - prctile( im1a.I_a(mask.tumor) ,25 )), ...
%          abs(median( im2a.I_a(mask.tumor) ) - prctile( im2a.I_a(mask.tumor) ,25 )),{'go--', 'g', 'g'});
%      hold on;
%      errorbarxy( median( im1a.I_a(mask.brain) ), median( im2a.I_a(mask.brain)), ...
%          abs(median( im1a.I_a(mask.brain) ) - prctile( im1a.I_a(mask.brain) ,25 )), ...
%          abs(median( im2a.I_a(mask.brain) ) - prctile( im2a.I_a(mask.brain) ,25 )),{'ro--', 'r', 'r'});
%      hold on;
%      
%      figure(h2)
%      hold on
%      errorbarxy( median( im1a.I_a(mask.tumor) ), median( im3a.I_a(mask.tumor)), ...
%          abs(median( im1a.I_a(mask.tumor) ) - prctile( im1a.I_a(mask.tumor) ,25 )), ...
%          abs(median( im3a.I_a(mask.tumor) ) - prctile( im3a.I_a(mask.tumor) ,25 )),{'go--', 'g', 'g'});
%      hold on;
%      errorbarxy( median( im1a.I_a(mask.brain) ), median( im3a.I_a(mask.brain)), ...
%          abs(median( im1a.I_a(mask.brain) ) - prctile( im1a.I_a(mask.brain) ,25 )), ...
%          abs(median( im3a.I_a(mask.brain) ) - prctile( im3a.I_a(mask.brain) ,25 )),{'ro--', 'r', 'r'});
%      hold on;
%      
% end

h1 = figure;
hold on;
for i = 1:length(fnames)
     load(['X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\',...
     fnames{i} 'wrp.mat']);
     
 figure(h1)
 covariance = cov([im1a.I_a(mask.tumor) im2a.I_a(mask.tumor)]);
 [V, eigenval ] = eig(covariance);
 theta = atan2(V(2,1),V(1,1))
 
 plot( median( im1a.I_a(mask.tumor) ), median( im2a.I_a(mask.tumor) ), 'k^','MarkerFaceColor','black','MarkerSize',5 );
 ellipse(abs(median( im2a.I_a(mask.tumor) ) - prctile( im2a.I_a(mask.tumor) ,32 )), ...
         abs(median( im1a.I_a(mask.tumor) ) - prctile( im1a.I_a(mask.tumor) ,32 )), 0,...%theta, ...
         median( im1a.I_a(mask.tumor) ), median( im2a.I_a(mask.tumor) ),'r');
 
     
 hold on
 covariance = cov([im1a.I_a(mask.brain) im2a.I_a(mask.brain)]);
 [V, eigenval ] = eig(covariance);
 theta = atan2(V(2,1),V(1,1));
 
 plot( median( im1a.I_a(mask.brain) ), median( im2a.I_a(mask.brain) ), 'ko','MarkerFaceColor','black','MarkerSize',5 );
 ellipse(abs(median( im2a.I_a(mask.brain) ) - prctile( im2a.I_a(mask.brain) ,32 )), ...
         abs(median( im1a.I_a(mask.brain) ) - prctile( im1a.I_a(mask.brain) ,32 )), 0,...%theta, ...
         median( im1a.I_a(mask.brain) ), median( im2a.I_a(mask.brain) ),'g');
     
end

h2 = figure;
hold on;
for i = 1:length(fnames)
     load(['X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\',...
     fnames{i} 'wrp.mat']);
     
 figure(h2)
 
 hold on
 plot(im1a.I_a(mask.tumor),im3a.I_a(mask.tumor),'g.','MarkerSize',1);

 
 hold on
 covariance = cov([im1a.I_a(mask.tumor) im3a.I_a(mask.tumor)]);
 [V, eigenval ] = eig(covariance);
 theta = atan2(V(2,1),V(1,1));
 
 plot( median( im1a.I_a(mask.tumor) ), median( im3a.I_a(mask.tumor) ), 'k^','MarkerFaceColor','black','MarkerSize',5 );
 ellipse( abs(median( im3a.I_a(mask.tumor) ) - prctile( im3a.I_a(mask.tumor) ,5 )), ...
          abs(median( im1a.I_a(mask.tumor) ) - prctile( im1a.I_a(mask.tumor) ,5 )), theta,...
          median( im1a.I_a(mask.tumor) ), median( im3a.I_a(mask.tumor) ),'k');
 
hold on
plot(im1a.I_a(mask.brain),im3a.I_a(mask.brain),'g.','MarkerSize',1);


 hold on
 covariance = cov([im1a.I_a(mask.brain) im3a.I_a(mask.brain)]);
 [V, eigenval ] = eig(covariance);
 theta = atan2(V(2,1),V(1,1));
 
 plot( median( im1a.I_a(mask.brain) ), median( im3a.I_a(mask.brain) ), 'ko','MarkerFaceColor','black','MarkerSize',5 );
 ellipse(abs(median( im3a.I_a(mask.brain) ) - prctile( im3a.I_a(mask.brain) ,5 )), ...
         abs(median( im1a.I_a(mask.brain) ) - prctile( im1a.I_a(mask.brain) ,5 )), theta, ...
         median( im1a.I_a(mask.brain) ), median( im3a.I_a(mask.brain) ),'k');

end

%   figure(h1)
%      hold on
%      errorbarxy( median( im1a.I_a(mask.tumor) ), median( im2a.I_a(mask.tumor)), ...
%          abs(median( im1a.I_a(mask.tumor) ) - prctile( im1a.I_a(mask.tumor) ,25 )), ...
%          abs(median( im2a.I_a(mask.tumor) ) - prctile( im2a.I_a(mask.tumor) ,25 )),{'go--', 'g', 'g'});
%      hold on;
%      errorbarxy( median( im1a.I_a(mask.brain) ), median( im2a.I_a(mask.brain)), ...
%          abs(median( im1a.I_a(mask.brain) ) - prctile( im1a.I_a(mask.brain) ,25 )), ...
%          abs(median( im2a.I_a(mask.brain) ) - prctile( im2a.I_a(mask.brain) ,25 )),{'ro--', 'r', 'r'});
%      hold on;
%      
%      figure(h2)
%      hold on
%      errorbarxy( median( im1a.I_a(mask.tumor) ), median( im3a.I_a(mask.tumor)), ...
%          abs(median( im1a.I_a(mask.tumor) ) - prctile( im1a.I_a(mask.tumor) ,25 )), ...
%          abs(median( im3a.I_a(mask.tumor) ) - prctile( im3a.I_a(mask.tumor) ,25 )),{'go--', 'g', 'g'});
%      hold on;
%      errorbarxy( median( im1a.I_a(mask.brain) ), median( im3a.I_a(mask.brain)), ...
%          abs(median( im1a.I_a(mask.brain) ) - prctile( im1a.I_a(mask.brain) ,25 )), ...
%          abs(median( im3a.I_a(mask.brain) ) - prctile( im3a.I_a(mask.brain) ,25 )),{'ro--', 'r', 'r'});
%      hold on;
%      
end
     

figure(h1)
xlabel('ABY029 (RFU)'); ylabel('PpIX (RFU)');

figure(h2)
xlabel('ABY029 (RFU)'); ylabel('IRDye680 (RFU)');


xlabel('ABY029 (RFU)');
ylabel('PpIX (RFU)');
zlabel('IRD680 (RFU)');
% do some ROC analysis,
[X,Y,T,AUC] = perfcurve(double(2.*mask.brain(:) + 1.*mask.tumor(:)),mprob(:),1);



for i = 1:10
    













