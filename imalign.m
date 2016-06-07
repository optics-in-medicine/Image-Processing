function [ im1a, im2a ] = imalign( im1, im2, var_size, T )
%IMALIGHN Summary of this function goes here
%   Detailed explanation goes here
% apply the point cloud cooregistration

if nargin == 2
    % default is that size can be variable, and is the smallest dimensions
    % across both images.
    var_size = 1;
    T = [0 0];
elseif nargin == 3
    T = [0 0];
end

im1a = im1;
im2a = im2;

MF1 = im1.res ./ 50;
MF2 = im2.res ./ 50;

im1a.BN = im1.BN .* MF1;
im2a.BN = im2.BN .* MF2;

I1_t = imresize(im1.I_c, MF1);
I2_t = imresize(im2.I_c, MF2);

[X1, Y1] = size( I1_t )
[X2, Y2] = size( I2_t )

xcrop = abs( X1 - X2 );
ycrop = abs( Y1 - Y2 );

% crop x
if xcrop >= 2
    if X1 < X2
        I1_tt = I1_t(1:X1,:);
        I2_tt = I2_t( fix(xcrop/2) : end-fix(xcrop/2) - 1, : );
        im2a.BN(:,2) = im2.BN(:,2) - fix(xcrop/2);
    elseif X2 < X1
        if var_size
            % force size to be the smallest image.
            I2_tt = I2_t(1:X2,:);
            I1_tt = I1_t( fix(xcrop/2) : end-fix(xcrop/2) - 1, : );
            im1a.BN(:,2) = im1.BN(:,2) - fix(xcrop/2);
        else
            %force size to be equal to im1
            I1_tt = I1_t(1:X1,:);
            I2_tt = I2_t;
            I2_tt(X1,:) = 0;

        end
    end
else
    I2_tt = I2_t;
    I1_tt = I1_t;
end

% crop y
if ycrop >= 2
    if Y1 < Y2
        im1.I_cr = I1_tt(:,1:Y1);
        im2.I_cr = I2_tt(:, fix(ycrop/2) : end-fix(ycrop/2) - 1 );
        im2a.BN(:,1) = im2.BN(:,1) - fix(ycrop/2);
    elseif Y2 < Y1
        
        if var_size        
            im2.I_cr = I2_tt(:,1:Y2);
            im1.I_cr = I1_tt(:, fix(ycrop/2) : end-fix(ycrop/2) - 1 );
            im1a.BN(:,1) = im1.BN(:,1) - fix(ycrop/2);
        else
            %force size to be equal to im1
            im1.I_cr = I1_tt(:,1:Y1);
            im2.I_cr = I2_tt;
            im2.I_cr(:,Y1) = 0;
        end
    end
else
    im2.I_cr = I2_tt;
    im1.I_cr = I1_tt;
end

[X1, Y1] = size( im1.I_cr );
[X2, Y2] = size( im2.I_cr );

im2a.I_cr = im2.I_cr(1:min([X1 X2]),1:min([Y1 Y2]));
im1a.I_cr = im1.I_cr(1:min([X1 X2]),1:min([Y1 Y2]));

[TR, TT, ER, t] = icp( [im1a.BN zeros(length(im1a.BN),1)]', [im2a.BN zeros(length(im2a.BN),1)]');

M_rot = eye(3,3); M_rot(1:2,1:2) = TR(1:2,1:2);
M_tra = eye(3,3); M_tra(3,1) = TT(1); M_tra(3,2) = TT(2);

BN_t = im2a.BN;
size(BN_t)
BN_tt = BN_t * TR(1:2,1:2);
size(BN_tt)
im2a.BN = [BN_tt(:,1) + TT(1)   BN_tt(:,2) + TT(2) ];

tform1 = affine2d(M_rot);
tform2 = affine2d(M_tra);

I_w1 = imwarp(im2a.I_cr,tform1);
% I_w2 = imwarp(I_w1,tform2);
I_w2 = imtranslate(I_w1,[TT(1)+ T(1), TT(2)+ T(2)]);

[X_t, Y_t] = size(I_w2);
[X, Y] = size(im1a.I_cr);

xcrop = abs( X - X_t );
ycrop = abs( Y - Y_t );

if (xcrop > 1 | ycrop > 1)
    im2a.I_t = I_w2( fix(xcrop/2) : end-fix(xcrop/2) - 1, fix(ycrop/2) : end-fix(ycrop/2) - 1 );
else
    im2a.I_t = I_w2;
end

I_temp = im2a.I_t;
clear im2.I_t


[X, Y] = size(im1.I_cr);
im1a.I_a(1:X,1:Y) = im1.I_cr(1:X,1:Y);
im2a.I_a(1:X,1:Y) = I_temp(1:X,1:Y);

im1a.I_a = im1.I_cr(1:X,1:Y);
im2a.I_a = I_temp(1:X,1:Y);
end

