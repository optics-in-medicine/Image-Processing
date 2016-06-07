 clear all
 close all

[filename,filepath] = uigetfile({'*.mat;*.dat;*.txt','All Data Files';...
          '*.*','All Files' },'mytitle',...
          'X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\');
      
load([filepath filename]);      
      
im1a.pathname

imHE.fname
TBR
TBR.mean(1) = mean(im1a.I_a(mask.tumor))./mean(im1a.I_a(mask.brain));
TBR.mean(2) = mean(im2a.I_a(mask.tumor))./mean(im2a.I_a(mask.brain));
TBR.mean(3) = mean(im3a.I_a(mask.tumor))./mean(im3a.I_a(mask.brain));


% 
% [X,Y,T,AUC(1)] = perfcurve(double(2.*mask.brain(:) + 1.*mask.tumor(:)),im1a.I_a(:),1);
% [X,Y,T,AUC(2)] = perfcurve(double(2.*mask.brain(:) + 1.*mask.tumor(:)),im2a.I_a(:),1);
% [X,Y,T,AUC(3)] = perfcurve(double(2.*mask.brain(:) + 1.*mask.tumor(:)),im3a.I_a(:),1);
% [X,Y,T,AUC(4)] = perfcurve(double(2.*mask.brain(:) + 1.*mask.tumor(:)),mprob(:),1);





