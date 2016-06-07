function [ stats ] = sample_from_points( varargin )
%SAMPLE_FROM_POINTS Summary of this function goes here
%   varargin, first should be smpl_pts

smpl_pts = varargin{1};

for i = 1 : nargin-1
    im{i} = varargin{i+1};
end
 

    
r = 5;
[Xgr,Ygr] = meshgrid(1:size(im{1}.I_a,2),1:size(im{1}.I_a,1));

for i = 1:length(smpl_pts)
    
    % make the logical operator to cut out around the roi point.
    roi = logical( ((Xgr - smpl_pts(i,1)).^2 +  (Ygr - smpl_pts(i,2)).^2 ) <= r.^2 );
    for j = 1:length(im)
        cutout = im{j}.I_a(roi);
        stats.mean(i,j) = mean( cutout );
        stats.sd(i,j) = std( cutout );
        stats.med(i,j) = median( cutout );
    end
end
