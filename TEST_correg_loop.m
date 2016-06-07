    %Load previous set
    [filename, textpath] = uigetfile([cd filesep '*.txt']);
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

    if (length(pth) >= 2)     % has to be at least 4 images to corregister.

        he_path = 'X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\F98 EGFR_pos B04\hico\'
        imHE.fname = filename(1:end-4);
    
        
        % CHANGE IF H&E IMAGE IS ROTATED COMPARED WITH FLUORESCENCE IMAGE.
        imHE.I = fliplr( flipud(imread([pth{2} '.TIF'])) );
        %   imHE.I = flipud(imread([pth{4} '.TIF']));
%            imHE.I = ((imread([pth{4} '.TIF'])) );
%         imHE.I = rot90(fliplr(imread([pth{4} '.TIF'])));
        h2 = imagesc(imHE.I); axis image; axis off; drawnow;

        % load first image (ABY029) - res of odyssey usually 42.
        im1 = load_images(pth{2});

        % load second image (PpIX) - res of typhoon usually 50.
        im2 = load_images(pth{3}, 0);

        close all;
        
        % should automate this by putting in a loop and having user input
        % 'T'.
        
        T = [-4 -10]; %<- translation, since it is usually still off a little (no idea why!).
        [ im1a,im2a ] = imalign( im1, im2, 0, T );

        % fusion of the first two fluorescence images.
        im_blend = cat( 3,4.* ( nm(im2a.I_a)) .* 0.5,zeros(size(im1a.I_a)), 5.* nm(im1a.I_a) );
        imagesc(im_blend); axis image; axis off;

        
        % load third image (IRD680)
        im3 = load_images(pth{1});

        
        % fusion of the first three fluorescence images.
        T = [0 -6]; %<- translation, since it is usually still off a little (no idea why!).
        [ ~, im3a ] = imalign( im1, im3, 0, T );

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
% 
% 
% 



% 
% 
% %% Get the logistic analsyis part
% [stats, mprob] = MAIN_LOGISTIC_ANALYSIS( im1a, im2a, im3a, imHE, mask, 'panel' );
% set(gcf,'Position',[488   113   890   649]);
% saveas(gcf,[pathname 'processed\' he_file '.fig'])
% saveas(gcf,[pathname 'processed\' he_file '.png'])
% 
% % warped image
% options = fuseset( 'model', 'logistic', 'cmap_type', 'sRGB', 'im2_type', 'sRGB', ...
%                                 'L', 0.6, 'x0', 0.3, 'k', 20);
% F = fuse_images( nm(mprob), imHE.I_a(1:size(im1a.I_a,1),1:size(im1a.I_a,2),:), goodmap('comet'), options );
% 
% figure; imagesc(F); axis image; axis off; colorbar; colormap(goodmap('comet'))
% 
% 
% fnames = {'405B_2','405C_2','407A_2','407B_1','410B_2', '408B_1','406A_2'};
% 
% for i = 1:length(fnames)
% % pick control points, number them after each click.
% load(['X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\',...
%     fnames{i} 'wrp.mat'])
% plot3( median( im1a.I_a(mask.tumor) ), ...
%        median( im2a.I_a(mask.tumor) ), ...
%        median( im3a.I_a(mask.tumor) ), ...
%       'ro','MarkerSize',10);
% hold on;
% plot3( median( im1a.I_a(mask.brain) ), ...
%        median( im2a.I_a(mask.brain) ), ...
%        median( im3a.I_a(mask.brain) ), ...
%       'go','MarkerSize',10);
% 
% % plot3( im1a.I_a(mask.tumor), ...
% %        im2a.I_a(mask.tumor), ...
% %        im3a.I_a(mask.tumor), ...
% %       'r.','MarkerSize',1);
% % plot3( im1a.I_a(mask.brain), ...
% %        im2a.I_a(mask.brain), ...
% %        im3a.I_a(mask.brain), ...
% %       'g.','MarkerSize',1);
% end
% 
% fnames = {'405B_2','405C_2','407B_1','410B_2', '408B_1','406A_2'};
% 
% h_s = figure;
% 
% for i = 1:length(fnames);
%      load(['X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\',...
%      fnames{i} 'wrp.mat']);
% 
%      [c, c_nm, bnds] = create_profile( im1a, im2a, im3a, mask )
%      
%      subplot(3,2,i);
%     plot(1:400,c_nm(:,1),'g',1:400,c_nm(:,2),'b',1:400,c_nm(:,3),'r',bnds(1)*[1 1],[0 1],'k--',bnds(2)*[1 1],[0 1],'k--')     
%     
% end
% 
% for i = 1:length(fnames)
%     
% end
% 
% % 
% % h1 = figure;
% % h2 = figure;
% % hold on;
% % for i = 1:length(fnames)
% %      load(['X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\',...
% %      fnames{i} 'wrp.mat']);
% % %      plot( median( im1a.I_a(mask.tumor) ), median( im2a.I_a(mask.tumor) ), 'ro','MarkerSize',10 );
% % %      hold on
% % %      plot( median( im1a.I_a(mask.brain) ), median( im2a.I_a(mask.brain) ), 'go','MarkerSize',10 );
% %      figure(h1)
% %      hold on
% %      errorbarxy( median( im1a.I_a(mask.tumor) ), median( im2a.I_a(mask.tumor)), ...
% %          abs(median( im1a.I_a(mask.tumor) ) - prctile( im1a.I_a(mask.tumor) ,25 )), ...
% %          abs(median( im2a.I_a(mask.tumor) ) - prctile( im2a.I_a(mask.tumor) ,25 )),{'go--', 'g', 'g'});
% %      hold on;
% %      errorbarxy( median( im1a.I_a(mask.brain) ), median( im2a.I_a(mask.brain)), ...
% %          abs(median( im1a.I_a(mask.brain) ) - prctile( im1a.I_a(mask.brain) ,25 )), ...
% %          abs(median( im2a.I_a(mask.brain) ) - prctile( im2a.I_a(mask.brain) ,25 )),{'ro--', 'r', 'r'});
% %      hold on;
% %      
% %      figure(h2)
% %      hold on
% %      errorbarxy( median( im1a.I_a(mask.tumor) ), median( im3a.I_a(mask.tumor)), ...
% %          abs(median( im1a.I_a(mask.tumor) ) - prctile( im1a.I_a(mask.tumor) ,25 )), ...
% %          abs(median( im3a.I_a(mask.tumor) ) - prctile( im3a.I_a(mask.tumor) ,25 )),{'go--', 'g', 'g'});
% %      hold on;
% %      errorbarxy( median( im1a.I_a(mask.brain) ), median( im3a.I_a(mask.brain)), ...
% %          abs(median( im1a.I_a(mask.brain) ) - prctile( im1a.I_a(mask.brain) ,25 )), ...
% %          abs(median( im3a.I_a(mask.brain) ) - prctile( im3a.I_a(mask.brain) ,25 )),{'ro--', 'r', 'r'});
% %      hold on;
% %      
% % end
% 
% h1 = figure;
% hold on;
% for i = 1:length(fnames)
%      load(['X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\',...
%      fnames{i} 'wrp.mat']);
%      
%  figure(h1)
%  covariance = cov([im1a.I_a(mask.tumor) im2a.I_a(mask.tumor)]);
%  [V, eigenval ] = eig(covariance);
%  theta = atan2(V(2,1),V(1,1))
%  
%  plot( median( im1a.I_a(mask.tumor) ), median( im2a.I_a(mask.tumor) ), 'k^','MarkerFaceColor','black','MarkerSize',5 );
%  ellipse(abs(median( im2a.I_a(mask.tumor) ) - prctile( im2a.I_a(mask.tumor) ,32 )), ...
%          abs(median( im1a.I_a(mask.tumor) ) - prctile( im1a.I_a(mask.tumor) ,32 )), 0,...%theta, ...
%          median( im1a.I_a(mask.tumor) ), median( im2a.I_a(mask.tumor) ),'r');
%  
%      
%  hold on
%  covariance = cov([im1a.I_a(mask.brain) im2a.I_a(mask.brain)]);
%  [V, eigenval ] = eig(covariance);
%  theta = atan2(V(2,1),V(1,1));
%  
%  plot( median( im1a.I_a(mask.brain) ), median( im2a.I_a(mask.brain) ), 'ko','MarkerFaceColor','black','MarkerSize',5 );
%  ellipse(abs(median( im2a.I_a(mask.brain) ) - prctile( im2a.I_a(mask.brain) ,32 )), ...
%          abs(median( im1a.I_a(mask.brain) ) - prctile( im1a.I_a(mask.brain) ,32 )), 0,...%theta, ...
%          median( im1a.I_a(mask.brain) ), median( im2a.I_a(mask.brain) ),'g');
%      
% end
% 
% h2 = figure;
% hold on;
% for i = 1:length(fnames)
%      load(['X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\',...
%      fnames{i} 'wrp.mat']);
%      
%  figure(h2)
%  
%  hold on
%  plot(im1a.I_a(mask.tumor),im3a.I_a(mask.tumor),'g.','MarkerSize',1);
% 
%  
%  hold on
%  covariance = cov([im1a.I_a(mask.tumor) im3a.I_a(mask.tumor)]);
%  [V, eigenval ] = eig(covariance);
%  theta = atan2(V(2,1),V(1,1));
%  
%  plot( median( im1a.I_a(mask.tumor) ), median( im3a.I_a(mask.tumor) ), 'k^','MarkerFaceColor','black','MarkerSize',5 );
%  ellipse( abs(median( im3a.I_a(mask.tumor) ) - prctile( im3a.I_a(mask.tumor) ,5 )), ...
%           abs(median( im1a.I_a(mask.tumor) ) - prctile( im1a.I_a(mask.tumor) ,5 )), theta,...
%           median( im1a.I_a(mask.tumor) ), median( im3a.I_a(mask.tumor) ),'k');
%  
% hold on
% plot(im1a.I_a(mask.brain),im3a.I_a(mask.brain),'g.','MarkerSize',1);
% 
% 
%  hold on
%  covariance = cov([im1a.I_a(mask.brain) im3a.I_a(mask.brain)]);
%  [V, eigenval ] = eig(covariance);
%  theta = atan2(V(2,1),V(1,1));
%  
%  plot( median( im1a.I_a(mask.brain) ), median( im3a.I_a(mask.brain) ), 'ko','MarkerFaceColor','black','MarkerSize',5 );
%  ellipse(abs(median( im3a.I_a(mask.brain) ) - prctile( im3a.I_a(mask.brain) ,5 )), ...
%          abs(median( im1a.I_a(mask.brain) ) - prctile( im1a.I_a(mask.brain) ,5 )), theta, ...
%          median( im1a.I_a(mask.brain) ), median( im3a.I_a(mask.brain) ),'k');
% 
% end
% 
% %   figure(h1)
% %      hold on
% %      errorbarxy( median( im1a.I_a(mask.tumor) ), median( im2a.I_a(mask.tumor)), ...
% %          abs(median( im1a.I_a(mask.tumor) ) - prctile( im1a.I_a(mask.tumor) ,25 )), ...
% %          abs(median( im2a.I_a(mask.tumor) ) - prctile( im2a.I_a(mask.tumor) ,25 )),{'go--', 'g', 'g'});
% %      hold on;
% %      errorbarxy( median( im1a.I_a(mask.brain) ), median( im2a.I_a(mask.brain)), ...
% %          abs(median( im1a.I_a(mask.brain) ) - prctile( im1a.I_a(mask.brain) ,25 )), ...
% %          abs(median( im2a.I_a(mask.brain) ) - prctile( im2a.I_a(mask.brain) ,25 )),{'ro--', 'r', 'r'});
% %      hold on;
% %      
% %      figure(h2)
% %      hold on
% %      errorbarxy( median( im1a.I_a(mask.tumor) ), median( im3a.I_a(mask.tumor)), ...
% %          abs(median( im1a.I_a(mask.tumor) ) - prctile( im1a.I_a(mask.tumor) ,25 )), ...
% %          abs(median( im3a.I_a(mask.tumor) ) - prctile( im3a.I_a(mask.tumor) ,25 )),{'go--', 'g', 'g'});
% %      hold on;
% %      errorbarxy( median( im1a.I_a(mask.brain) ), median( im3a.I_a(mask.brain)), ...
% %          abs(median( im1a.I_a(mask.brain) ) - prctile( im1a.I_a(mask.brain) ,25 )), ...
% %          abs(median( im3a.I_a(mask.brain) ) - prctile( im3a.I_a(mask.brain) ,25 )),{'ro--', 'r', 'r'});
% %      hold on;
% %      
%   end
%      
% 
% figure(h1)
% xlabel('ABY029 (RFU)'); ylabel('PpIX (RFU)');
% 
% figure(h2)
% xlabel('ABY029 (RFU)'); ylabel('IRDye680 (RFU)');
% 
% 
% xlabel('ABY029 (RFU)');
% ylabel('PpIX (RFU)');
% zlabel('IRD680 (RFU)');
% % do some ROC analysis,
% [X,Y,T,AUC] = perfcurve(double(2.*mask.brain(:) + 1.*mask.tumor(:)),mprob(:),1);
% 
% 
% 
% for i = 1:10
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




