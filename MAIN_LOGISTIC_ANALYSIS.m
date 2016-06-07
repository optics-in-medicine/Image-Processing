function [stats, mprob] = MAIN_LOGISTIC_ANALYSIS( im1a, im2a, im3a, imHE, mask, wplot )
%MAIN_LOGISTIC_ANALYSIS Load in corregistered data and perform logistic 
% analysis to generate tumor probability map.

% clear all
% close all
% 
% % First, run MAIN_CORREGISTRATION program and then load the saved data.
% 
% filename = ['X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\' ...
%             '405B_2.mat'];
% 
% load( filename );

% Statistical analysis on classified data

mf = [2 2]; 

% unnormalized predictors?
% P1 = medfilt2(im1a.I_a, mf);
% P2 = medfilt2(im2a.I_a, mf);
% P3 = medfilt2(im3a.I_a, mf);

s.min1 = min(im1a.I_a(:));
s.min2 = min(im2a.I_a(:));
s.min3 = min(im3a.I_a(:));

s.max1 = max(im1a.I_a(:));
s.max2 = max(im2a.I_a(:));
s.max3 = max(im3a.I_a(:));

% % normalize predictors first?
P1 = medfilt2(nm(im1a.I_a), mf);
P2 = medfilt2(nm(im2a.I_a), mf);
P3 = medfilt2(nm(im3a.I_a), mf);

% build data matrix
X1 = [P1(mask.tumor) P2(mask.tumor) P3(mask.tumor)];
Y1 = ones(sum(mask.tumor(:)),1);

X = [X1; [P1(mask.brain) P2(mask.brain)  P3(mask.brain)]];
Y = [Y1; zeros(sum(mask.brain(:)),1)];

% fit a logit function to the data.
[mp, prob1, r1, t01] = get_logistic( P1(mask.tumor), P1(mask.brain), 100 );
[mp, prob2, r2, t02] = get_logistic( P2(mask.tumor), P2(mask.brain), 100 );
[mp, prob3, r3, t03] = get_logistic( P3(mask.tumor), P3(mask.brain), 100 );

beta0 = -(r1*t01 + r2*t02 + r3*t03);

 mprob = 1 ./ (1 + ( exp(-( (beta0 + r1.*nm(im1a.I_a) + r2.*nm(im2a.I_a) + r3.*nm(im3a.I_a))))));
% 
% mprob = (logistic_func(nm(im1a.I_a),1,r1,t01) + ...
%          logistic_func(nm(im2a.I_a),1,r2,t02)) ./ 2;
     
 
% mprob = 1 ./ (1 + ( exp(-( (-14.75 + 75.9.*nm(im1a.I_a) + 11.5.*nm(im2a.I_a) + 21.6.*nm(im3a.I_a))))));    

stats.prob1 = prob1;
stats.prob2 = prob2;
stats.prob3 = prob3;
stats.beta0 = beta0;
stats.mprob = mprob;
stats.r1 = r1;
stats.r2 = r2; 
stats.r3 = r3;



% 
% 
% 
% 
% save([filename(1:end-4) 'mnr.mat'], 'prob1','prob2','prob3','mp','r1','r2','r3','t01','t02','t03','X','Y','mprob');

if strcmpi(wplot,'panel')

    figure('color','white');
    subplot('Position',[0.05 0.75 0.3 0.2])
    imagesc(nm(im1a.I_a)); colormap(goodmap('cube1')); caxis([0 0.6]); axis image; axis off; title('ABY029'); freezeColors
    subplot('Position',[0.35 0.75 0.3 0.2])
    imagesc(nm(im2a.I_a)); colormap(goodmap('cube1')); caxis([0 1.0]); axis image; axis off; title('PpIX'); freezeColors
    subplot('Position',[0.65 0.75 0.3 0.2])
    imagesc(nm(im3a.I_a)); colormap(goodmap('cube1')); caxis([0 0.4]); axis image; axis off; title('IRD680'); freezeColors

    subplot('Position',[0.08 0.38 0.4 0.35])
    im_blend = cat(3, 1.0.*nm(im3a.I_a), 1.7 .* nm(im1a.I_a).^0.8, 0.7 .* nm(im2a.I_a));
    imagesc(im_blend); axis image; axis off; 
    text(10,10,'ABY029','color','green','FontWeight','bold')
    text(10,30,'PpIX','color',[0.5 0.4 1],'FontWeight','bold')
    text(10,50,'IRD680','color','red','FontWeight','bold')

    subplot('Position',[0.58 0.45 0.32 0.28])
%     plot(mp,logistic_func(mp,1,r1,t01),'g',mp,logistic_func(mp,1,r2,t02),'b',mp,logistic_func(mp,1,r3,t03),'r'); 
    plot(mp,logistic_func([1 r1 t01],mp),'g',mp,logistic_func([1 r2 t02],mp),'b',mp,logistic_func([1 r3 t03],mp),'r'); 
    
    legend('ABY029','PpIX','IRD680');

    hold on
    plot(mp, prob1, '.g',...
         mp, prob2, '.b',...
         mp, prob3, '.r', 'MarkerSize',3)
    hold off

    xlabel('Relative Intensity'); ylabel('Tumor Probability')

    subplot('Position',[0.08 0.0 0.4 0.35])
    imagesc(imHE.I); axis image; axis off; 

    subplot('Position',[0.55 0.05 0.4 0.3])
    imagesc(mprob); colormap(goodmap('newdawn')); freezeColors; colorbar( 'FontWeight','bold'); axis image; axis off; 
    text(size(mprob,2)./3.2, 7,'Tumor Probability Map', 'FontWeight','bold')
    

    hold on
    plot(im1a.BN(:,1),im1a.BN(:,2),'LineWidth',3, 'color',[0.5 0.1 1])
    hold off

else%if strcmpi(wplot,'fusion')
    imagesc(mprob);    
    colormap(goodmap('vascular'))

end
% 
% % 
% [rows,cols] = size(P1);
% 
% [B, dev, stats] = mnrfit(X,Y+1);
% 
% act_el = logical(mask.tumor) + logical(mask.brain);
% 
% TPM = zeros(rows,cols);
% TPM(logical(act_el)) = pihat(:,2);
% 
% 
% 
% B0 = -44.56;
% B1 = 169.81;
% B2 = 19.55;
% B3 = 9.8;
% 
% mprob_avg = 1 ./ (1 + ( exp(-( (B0 + B1.*nm(im1a.I_a) + B2.*nm(im2a.I_a) + B3.*nm(im3a.I_a))))));


