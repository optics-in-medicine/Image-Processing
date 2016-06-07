% This is a script to demonstrate warping an image using the Thin Plate 
% Spline method presented by Bookstein (see ref).
%
% Alex Hartov, ENGG 111, 03/04/12

clear



%% Example with the outline of a guy in a binary image
imIn=imread('Outline.tiff');    % Read a test image
R=double(imIn(:,:,1));
R=R/max(R(:));                  % Double valued red plane in range 0 to 1
G=zeros(size(R));
B=zeros(size(R));


% Modify the picture by adding a grid so we can see distortion
step=20;
[nr nc]=size(R);
for r=1:step:nr
    imIn(r,:)=~imIn(r,:);
end
for c=1:step:nc
    imIn(:,c)=~imIn(:,c);
end

% Show what we have
imIn=im2bw(imIn(:,:,1));        % Single plane B & W image.
[nr nc]=size(imIn);


%% Make a skinny guy
figure(1);
subplot(1,3,1);
imshow(imIn);
drawnow;

% Source reference points for outline image - Note Row Col coordinates
Ps=[224 170     % Upper left
    332 170     % Next down
    420 170     % Next down
    494 170     % Next down
    224 217     % Upper right
    332 217     % Next down
    420 217     % Next down
    494 217];   % Next down

% Destination reference points X +/- 5 going towards the medial axis for
% outline image
Pd1=[224 170    % No change
     332 175    % To the right towards the median
     420 175    % Same
     494 170    % No change
     224 217    % No change
     332 212    % To the left towards the median
     420 212    % Same
     494 217];  % No change

hold on
plot(Ps(:,2),Ps(:,1),'rx');
plot(Pd1(:,2),Pd1(:,1),'bx');
hold off

%imOut=warpItFun(imIn,Ps,Pd1);
imOut=warpItFun(imIn,Pd1,Ps);  % Use with mod warpItFun to remove artifacts
subplot(1,3,2);
imshow(imOut,[]);
G(imOut>0)=1;
B(R~=G)=1;
subplot(1,3,3);
imRGB(:,:,1)=R;
imRGB(:,:,2)=G;
imRGB(:,:,3)=B;
imshow(imRGB);


%% Make a fat guy
figure(2);
subplot(1,3,1);
imshow(imIn);

% Destination reference points X +/- 5 going away from the medial axis for
% the outline image
Pd2=[224 170    % No change from starting position
     332 165    % To the left away from median
     420 165    % Same
     494 170    % No change
     224 217    % No Change
     332 222    % To the right away from median
     420 222    % Same
     494 217];  % No change

hold on
plot(Ps(:,2),Ps(:,1),'rx');
plot(Pd2(:,2),Pd2(:,1),'bx');
hold off

%imOut=warpItFun(imIn,Ps,Pd2);
imOut=warpItFun(imIn,Pd2,Ps);  % Use with mod warpItFun to remove artifacts
subplot(1,3,2);
imshow(imOut,[]);
G(imOut>0)=1;
B(R~=G)=1;
subplot(1,3,3);
imRGB(:,:,1)=R;
imRGB(:,:,2)=G;
imRGB(:,:,3)=B;
imshow(imRGB);

