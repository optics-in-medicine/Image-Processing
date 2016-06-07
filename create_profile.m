function [c, c_nm, bnds] = create_profile( im1a, im2a, im3a, mask )
%CREATE_PROFILE Summary of this function goes here
%   Detailed explanation goes here

h2 = figure;
im_blend = cat(3, 2.0.*nm(im3a.I_a), 1.7 .* nm(im1a.I_a).^0.8, 0.7 .* nm(im2a.I_a));
imagesc(im_blend); axis image; axis off

hl = imline;
position = wait(hl);
close(h2);
M = (position(2,2)-position(2,1)) / (position(1,2)-position(1,1))

if abs(M) > 2
    x1 = 1; x2 = -1; y1 = 0; y2 = 0;        
elseif abs(M) <= 0.5
    x1 = 0; x2 = 0; y1 = 1; y2 = -1;    
else
    x1 = 1; x2 = 0; y1 = 0; y2 = +1;    
end

position = position';
% close all
n = 400;

c1 = (1/3) .* ( improfile(im1a.I_a,position(1,:),position(2,:),n,'bicubic') + ...
     improfile(im1a.I_a,position(1,:)+x1,position(2,:)+y1,n,'bicubic') + ...
     improfile(im1a.I_a,position(1,:)+x2,position(2,:)+y2,n,'bicubic'));
     
c2 = (1/3) .* ( improfile(im2a.I_a,position(1,:),position(2,:),n,'bicubic') + ...
     improfile(im2a.I_a,position(1,:)+x1,position(2,:)+y1,n,'bicubic') + ...
     improfile(im2a.I_a,position(1,:)+x2,position(2,:)+y2,n,'bicubic'));
     
c3 = (1/3) .* ( improfile(im3a.I_a,position(1,:),position(2,:),n,'bicubic') + ...
     improfile(im3a.I_a,position(1,:)+x1,position(2,:)+y1,n,'bicubic') + ...
     improfile(im3a.I_a,position(1,:)+x2,position(2,:)+y2,n,'bicubic'));
 
 ct =  improfile(double(mask.tumor),position(1,:),position(2,:),n,'bicubic') + ...
     improfile(double(mask.tumor),position(1,:)+x1,position(2,:)+y1,n,'bicubic') + ...
     improfile(double(mask.tumor),position(1,:)+x2,position(2,:)+y2,n,'bicubic');
     
c = [c1 c2 c3];
c_nm = [nm(c1) nm(c2) nm(c3)];
% c_sm = [smooth(nm(c1),10,'sgolay',6) smooth(nm(c2),10,'sgolay',6) smooth(nm(c3),10,'sgolay',6)];
% c_nm = [( c1 - im1a.LL ) ./ (im1a.UL - im1a.LL) ...
%         ( c2 - im2a.LL ) ./ (im2a.UL - im2a.LL) ...
%         ( c3 - im3a.LL ) ./ (im3a.UL - im3a.LL)];

% finds the bounds of the tumor (assumes only two interfaces).
f1 = fdo(length(ct),length(ct))*ct;
[~,minpt] = min(f1);
[~,maxpt] = max(f1);
bnds = [min([minpt maxpt]) max([minpt maxpt])];

% plot(1:400,nm(c1),'g',1:400,nm(c2),'b',1:400,nm(c3),'r',bnds(1)*[1 1],[0 1],'k--',bnds(2)*[1 1],[0 1],'k--')
end

