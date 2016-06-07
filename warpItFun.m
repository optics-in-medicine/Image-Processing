function imOut = warpItFun(imIn,Ps,Pd)
% The function a very simple Thin-Plane spline warping of the two images. 
% Image imIn is warped to imOut.  The warping is defined by a 
% set of reference points Ps = [Xp Yp] and a set of destination points 
% Pd = [Xd Yd].
% The warping will move points in Ps to Pd exactly using the thin plate
% splines algorithm.
% 
% Reference: F.L. Bookstein, "Principal Warps: Thin-Plate splines and the 
% decomposition % of deformations", IEEE Tr. on Pattern Analysis and 
% Machine Intel, vol. 11, No. 6, June 1989
% Originally created by Alex Gloria Ossadtchi at Halfacre, aos@usc.edu 
% 
% Modified by Yuqing Wan 10/2010
% 
% Additional modifications by Alex Hartov for ENGG 111 03/04/12

NPs = size(Ps,1);
Xs = Ps(:,1)';  % Source points
Ys = Ps(:,2)'; 
Xd = Pd(:,1)';  % Desitnation points
Yd = Pd(:,2)';
rXs = repmat(Xs(:),1,NPs);  % Built only for computing r in K matrix
rYs = repmat(Ys(:),1,NPs);  % Same

% Compute the distance
wR = sqrt((rXs-rXs').^2 + (rYs-rYs').^2);
% matrix K.  Note the use of 1e-20 to avoice log(0)
wK = 2*(wR.^2).*log(wR.^2+1e-20);
% matrix P
wP = [ones(NPs,1) Xs(:) Ys(:)];
% matrix L
wL = [wK wP;wP' zeros(3,3)];
% matrix Y
wY = [Xd(:) Yd(:); zeros(3,2)];
% matrix W
wW = inv(wL)*wY;
img1_tmp = imIn;

% point coordinates in the old image
[X Y] = meshgrid(1:size(img1_tmp,1),1:size(img1_tmp,2)); 
X = X(:)'; 
Y = Y(:)'; 
NWs = length(X);

% calculate the new coordinate for points in the old image
rX = repmat(X,NPs,1); 
rY = repmat(Y,NPs,1);
rXs = repmat(Xs(:),1,NWs); 
rYs = repmat(Ys(:),1,NWs);
wR = sqrt((rXs-rX).^2 + (rYs-rY).^2); 
wK = 2*(wR.^2).*log(wR.^2+1e-20); 
wP = [ones(NWs,1) X(:) Y(:)]'; 
wL = [wK;wP]';

% new coordinates
Xw = wL*wW(:,1); 
Yw = wL*wW(:,2);
Xw = round(max(min(Xw,size(imIn,1)),1)); 
Yw = round(max(min(Yw,size(imIn,2)),1));
iw = sub2ind(size(imIn),Xw,Yw); 
ip = sub2ind(size(imIn),X',Y');

% new image after warping
imOut = uint8(zeros(size(imIn)));

% put in the warped pixels
%img1w(round(iw)) = img1_tmp(round(ip)); 

% Use reversed indices to avoid artifacts resulting from gaps in pixel
% assigments resulting from expanding.  In order to use this approach, the
% source and destination points must also be reversed.  Works in most
% situations.  Not tested for the general case.  (AH mod)
imOut(round(ip)) = img1_tmp(round(iw));


