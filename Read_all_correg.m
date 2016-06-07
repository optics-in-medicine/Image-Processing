clear all
close all


files = {'405C','405B','404A','505B','508A','407A','409A'};
figure('color','white')
for kk = 1:length(files)
    load(strcat('X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\Corregistered Data Files\Corregistered Data\',char(files(kk))),'im*','B');
    %subplot(6,4,4*kk-3)
    subaxis(length(files),4,4*kk-3, 'Spacing', 0.0, 'Padding', 0, 'Margin', 0);
    pcolor(double(real((flipud(nm(im1a.I_a)))))); shading flat; axis image; axis off; colormap(goodmap('greens')); 
    if kk <= 5
        caxis([0 0.6])
    else
        caxis([0.1 1.1]);
    end
    freezeColors;
    subaxis(length(files),4,4*kk-2, 'Spacing', 0.0, 'Padding', 0, 'Margin', 0);
%     subaxis(5,5,i, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
    pcolor(double(real((flipud(nm(im2a.I_a)))))); shading flat; axis image; axis off; colormap(goodmap('blues')); caxis([0 1.5]);
    freezeColors;
    subaxis(length(files),4,4*kk-1, 'Spacing', 0.0, 'Padding', 0, 'Margin', 0);
%     subaxis(5,5,i, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
    pcolor(double(real((flipud(nm(im3a.I_a)))))); shading flat; axis image; axis off; colormap(goodmap('reds')); caxis([0.25 1]);
  
    if kk == 7
        caxis([0.25 0.5])
    else
        caxis([0.25 1]);
    end
    
    subaxis(length(files),4,4*kk, 'Spacing', 0.0, 'Padding', 0, 'Margin', 0);
%     subaxis(5,5,i, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
    imagesc(B); axis image; axis off;
end


clear all
close all

files = {'405C','405B','407A','406A'};%,'501B','508A','406A','509B','505B'};

files = {'501B','508A','509B','505B'};


lr_vect = [];
for kk = 1:length(files)
   load(strcat('X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\Corregistered Data Files\Corregistered Data\',char(files(kk))),'stats','smpl_pts');
   lr_vect = [lr_vect; [stats.mean smpl_pts(:,3)~=1]];
end

xout = [];
for i = 1:3
    
%     X(:,i) = linspace(prctile(lr_vect(:,i),5),prctile(lr_vect(:,i),95),20)';
%     N = hist(lr_vect(lr_vect(:,4)==0,i),X(:,i))';
%     P = hist(lr_vect(lr_vect(:,4)==1,i),X(:,i))';
%     
    [xdata,I] = sort( lr_vect(:,i), 1, 'ascend' );
    ydata = lr_vect(I,4);
    x = lsqcurvefit(@logistic_func, [0.5 0.5 1], nm(xdata), ydata);
    hold on;
    plot(nm(xdata),logistic_func(x,nm(xdata)),'-');%, nm(xdata), ydata,'o');
    xout = [xout [nm(xdata) logistic_func(x,nm(xdata))]];
    
    
end




    
clear all
close all


files = {'410B','405C','405B','407A','406A','501B','508A','406A','509B','505B'};

% files = {'501B','508A','509B','505B'};

for kk = 1:length(files)
    load(strcat('X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\Corregistered Data Files\Corregistered Data\',char(files(kk))),'im*','tum_reg','BW_reg');
    mask.tumor = logical(tum_reg);
    mask.brain = logical(abs(BW_reg - tum_reg));
    [log_stats{kk}, mprob{kk}] = MAIN_LOGISTIC_ANALYSIS( im1a, im2a, im3a, imHE, mask, 'panel' );
    
    [X,Y,T,AUC(kk)] = perfcurve( double(mask.tumor(:)), mprob{kk}(:) ,1 )

end

r1_c = []; r2_c = []; r3_c = []; beta_c = [];
for kk = 1:length(files)
r1_c = [r1_c; log_stats{kk}.r1];
r2_c = [r2_c; log_stats{kk}.r2];
r3_c = [r3_c; log_stats{kk}.r3];
beta_c = [beta_c; log_stats{kk}.beta0];
end


clear all
close all
avgs = [];
allpts = [];
%files = {'405C','405B','407A','410B','501B','508A','406A','509B','505B'};
% files = {'404A1','409A','507B','509B','505B'};
files = {'405C','406A','405B','410B','404A','404A1','501B','508A','505B','506B','507B','502A','509B','F8','407A','409A'};
for kk = 1:length(files)
    load(strcat('X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\Corregistered Data Files\Corregistered Data\',char(files(kk))),'stats','smpl_pts');
    im1a_avg(kk) = mean(stats.mean(smpl_pts(:,3) == 1,1));
    im2a_avg(kk) = mean(stats.mean(smpl_pts(:,3) == 1,2));
    im3a_avg(kk) = mean(stats.mean(smpl_pts(:,3) == 1,3));
    
    im1a_tum(kk) = mean(stats.mean(smpl_pts(:,3) ~= 1,1));
    im2a_tum(kk) = mean(stats.mean(smpl_pts(:,3) ~= 1,2));
    im3a_tum(kk) = mean(stats.mean(smpl_pts(:,3) ~= 1,3));
    
    
%     allpts = [allpts; [stats.mean smpl_pts(:,3)]];
    avgs = [avgs; [stats.mean(smpl_pts(:,3) == 1,1) stats.mean(smpl_pts(:,3) == 1,2) stats.mean(smpl_pts(:,3) == 1,3)]];
    % ABY029 marg TBR 
    grp_mean(kk,:) = [mean( stats.mean(smpl_pts(:,3) == 2,1)./im1a_avg(kk) ); mean( stats.mean(smpl_pts(:,3) == 3,1)./im1a_avg(kk) );
                mean( stats.mean(smpl_pts(:,3) == 2,2)./im2a_avg(kk) ); mean( stats.mean(smpl_pts(:,3) == 3,2)./im2a_avg(kk) );
                mean( stats.mean(smpl_pts(:,3) == 2,3)./im3a_avg(kk) ); mean( stats.mean(smpl_pts(:,3) == 3,3)./im3a_avg(kk) )]';
            
    grp_SEM(kk,:) =  [std( stats.mean(smpl_pts(:,3) == 2,1)./im1a_avg(kk) ); std( stats.mean(smpl_pts(:,3) == 3,1)./im1a_avg(kk) );
                std( stats.mean(smpl_pts(:,3) == 2,2)./im2a_avg(kk) ); std( stats.mean(smpl_pts(:,3) == 3,2)./im2a_avg(kk) );
                std( stats.mean(smpl_pts(:,3) == 2,3)./im3a_avg(kk) ); std( stats.mean(smpl_pts(:,3) == 3,3)./im3a_avg(kk) )]' ./sqrt(10);    
end

roc_vect = [];
for kk = 1:length(files)
    load(strcat('X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\Corregistered Data Files\Corregistered Data\',char(files(kk))));

    roc_vect = [roc_vect; [stats.mean smpl_pts(:,3)]];
    % ABY029 marg TBR 
   
end

[X,Y,T,AUC] = perfcurve( double(roc_vect(1:end,4)>1), roc_vect(1:end,1) ,1 )

[mp, prob1, r1, t01] = get_logistic( roc_vect(smpl_pts(:,3)>1,1:3), roc_vect(smpl_pts(:,3)==1,1:3), 100 );


[stats, mprob] = MAIN_LOGISTIC_ANALYSIS( im1a, im2a, im3a, imHE, mask, wplot )




% load('X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\Corregistered Data Files\Corregistered Data\407A')


    figure('color','white');
    im_blend = cat(3, 1.5.*nm(im3a.I_a), 3.5 .* nm(im1a.I_a) .^ 1.1, 0.5 .* nm(im2a.I_a));
    imagesc(im_blend); axis image; axis off
    
    hold on;
    xx = 1;
    for i = 1:30
        switch smpl_pts(i,3)
            case 1
                plot(smpl_pts(i,1),smpl_pts(i,2),'wo');
            case 2
                plot(smpl_pts(i,1),smpl_pts(i,2),'ro');
            case 3
                plot(smpl_pts(i,1),smpl_pts(i,2),'go');
        end
        text(smpl_pts(i,1)+5,smpl_pts(i,2)-5,num2str(i),'BackgroundColor',[1 1 1])
    end
    
%     plot(smpl_pts(smpl_pts(:,3) == 2,1),smpl_pts(smpl_pts(:,3) == 2,2),'ro')
%     plot(smpl_pts(smpl_pts(:,3) == 3,1),smpl_pts(smpl_pts(:,3) == 3,2),'go')
%     hold off;
    
    im1a_avg = mean(stats.mean(smpl_pts(:,3) == 1,1))
    im2a_avg = mean(stats.mean(smpl_pts(:,3) == 1,2))
    im3a_avg = mean(stats.mean(smpl_pts(:,3) == 1,3))
    
    figure('color','white');
    subplot(1,3,1)
    boxplot(stats.mean(:,1)./im1a_avg,smpl_pts(:,3),'labels',{'normal','margin','core'});
    title('ABY029');
    subplot(1,3,2)
    boxplot(stats.mean(:,2)./im2a_avg,smpl_pts(:,3),'labels',{'normal','margin','core'});
    title('PpIX');
    subplot(1,3,3)
    boxplot(stats.mean(:,3)./im3a_avg,smpl_pts(:,3),'labels',{'normal','margin','core'});
    title('IRDye680RD');
    
    
    
    load('X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\Corregistered Data Files\Corregistered Data\501B')
    
    %correlation
    plot(stats.mean(smpl_pts(:,3)==1,2),stats.mean(smpl_pts(:,3)==1,1),'ko')
    hold on
    plot(stats.mean(smpl_pts(:,3)==2,2),stats.mean(smpl_pts(:,3)==2,1),'go')
    plot(stats.mean(smpl_pts(:,3)==3,2),stats.mean(smpl_pts(:,3)==3,1),'ro')
    
    
    load('X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\Corregistered Data Files\Corregistered Data\407A')
    
    %correlation
    plot(stats.mean(smpl_pts(:,3)==1,2),stats.mean(smpl_pts(:,3)==1,1),'ko')
    hold on
    plot(stats.mean(smpl_pts(:,3)==2,2),stats.mean(smpl_pts(:,3)==2,1),'go')
    plot(stats.mean(smpl_pts(:,3)==3,2),stats.mean(smpl_pts(:,3)==3,1),'ro')
    
    load('X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\Corregistered Data Files\Corregistered Data\405C')
    
    %correlation
    plot(stats.mean(smpl_pts(:,3)==1,2),stats.mean(smpl_pts(:,3)==1,1),'ko')
    hold on
    plot(stats.mean(smpl_pts(:,3)==2,2),stats.mean(smpl_pts(:,3)==2,1),'go')
    plot(stats.mean(smpl_pts(:,3)==3,2),stats.mean(smpl_pts(:,3)==3,1),'ro')
    
    
    for i = 1:3
        
        
        
        
        
        
        
    end
    
    
    
    
%     KS test
k = 1;

        [h,p,ks2stat] = kstest2( interp1(KS(:,k),KS(:,k+1),linspace(KS(1,k),KS(end,k),100),'pchip'), ...
                                 interp1(KS(:,k+2),KS(:,k+3),linspace(KS(1,k+2),KS(end,k+2),100),'pchip'));
                             
                             
      plot( interp1(KS(:,k),KS(:,k+1),linspace(KS(1,k),KS(end,k),100),'pchip'), ...
                                 interp1(KS(:,k+2),KS(:,k+3),linspace(KS(1,k+2),KS(end,k+2),100),'pchip'))

plot(interp1(KS(:,k),KS(:,k+1),linspace(KS(1,k),KS(end,k),100),'pchip'),'.'); hold on;

plot(interp1(KS(:,k+2),KS(:,k+3),linspace(KS(1,k+2),KS(end,k+2),100),'pchip'));

