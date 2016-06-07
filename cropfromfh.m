function [ I_c, pos_c ] = cropfromfh( I, pos )
%CROPFROMFH Summary of this function goes here
%   Detailed explanation goes here

% crop image to +- 10% of maximum and minimum values.
min_px = min(pos);
max_px = max(pos);
mean_px = mean(pos);

% find boundary

% determine the distance between x's and y's
dist_px = ceil( abs( min_px - max_px) );

% define the [xmin ymin ...] by adding 10% to each dim.
xmin = fix(min_px(1) - 0.1 .* dist_px(1));
ymin = fix(min_px(2) - 0.1 .* dist_px(2));
width = fix(dist_px(1) + dist_px(1) .* 0.2);
height = fix(dist_px(2) + dist_px(2) .* 0.2);

base_dim = max([width height]);

I_c = imcrop( I, [xmin ymin width height] );

[rows cols chan] = size(I_c);

% if the image is too close to the edge, it will want to add extra space so
% that it is centered in the frame and has +- 10% on each side.

c_val = I_c(1,1,:);

% check to see if cropped image corner is negative x value.
if xmin < 0

    for i = 1:chan
        I_c_t(:,:,i) = [c_val(i) .* ones([rows abs(xmin)],'like',I_c) I_c(:,:,i)];
    end
    I_c = I_c_t;
end

% check to see if cropped image corner is negative y value.
if ymin < 0

    for i = 1:chan
        I_c_t(:,:,i) = [c_val(i) .* ones([abs(ymin) cols],'like',I_c); I_c(:,:,i)];
    end
    I_c = I_c_t;
end

% update the points specifying the free hand contour drawn in previous
% steps.
pos_c(:,1) = pos(:,1) + abs(xmin);
pos_c(:,2) = pos(:,2) - abs(ymin);

end

