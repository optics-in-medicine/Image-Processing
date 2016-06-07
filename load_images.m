% clear all
% close all
% 
% cd('X:\#6 - Code\Matlab Code\# # Coregistration\');
% 
% % Bottom image
% PathName = 'X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\F405\';

function [im] = load_images( pathname, rot, im )

if nargin == 1
    [ im.I, imstat, im.pathname, im.filename ] = pick_image( pathname );
    if strcmpi(im.filename(end-2:end),'gel')
        im.I = flipud(im.I);
    end
    rot = 0;
elseif nargin == 2
    [ im.I, imstat, im.pathname, im.filename ] = pick_image( pathname );

    if strcmpi(im.filename(end-2:end),'gel')
        im.I = flipud(im.I);
    end
    
    if rot ~= 0
        im.I = rot90(im.I,rot);
    end
else
    imstat = 1;
end

if imstat
    im.res = str2num( char( inputdlg('What is the image resolution in um?')) );
    
    hp = figure;
    h = imagesc(im.I); axis image; colormap(goodmap('cube1'));
    hfig = imcontrast(h);
    waitfor(hfig);

    I1_adj = getimage(h);
    im.LL = mean(im.I(I1_adj == 0));
    im.UL = mean(im.I(I1_adj == 1));

    figure(hp);
    imagesc(im.I); caxis([im.LL im.UL]); axis image;

    % circle out the ROI
    hFH = imfreehand();
    pos = hFH.getPosition();
    binaryImage = hFH.createMask();
    
    im.I_c = cropfromfh( im.I, pos );
   
    
    hp = figure('color','white');
    subplot(1,2,1);
    imagesc(im.I_c); caxis([im.LL im.UL]); colormap(goodmap('cube1')); axis image;
    h2 = subplot(1,2,2);
    
    im.I_n = (im.I_c - im.LL) ./ (im.UL - im.LL);
    
    lvl = linspace(0,1,100);
    % Clever little code 
    for i = 1:length(lvl)
%         I_test = im2bw(im.I_c ./ prctile(im.I_c(:),i));
        I_test = im2bw(im.I_n, lvl(i));
        obj_fct(i) = abs( sum( binaryImage(:) > 0 ) - sum( I_test(:) > 0 ) );
        
    end
    
    
    [min_of, index] = min(obj_fct);
    BW = imfill( im2bw(im.I_n, lvl(index)), 'holes') ; 
    im.BW = BW;
    
    quit_boundary = 0;
    
    while( ~quit_boundary )
        imagesc(BW);axis image;

        hFH = imfreehand(h2,'Closed',0);
        pos = hFH.getPosition();
        wait(hFH);

        % basically, need to reinterpolate this so that 
    %      edge_i = [linspace(pos(1,2),pos(end,2),100)' linspace(pos(1,1),pos(end,1),100)'];
         edge_i = [interp1(linspace(0,1,length(pos)),pos(:,1),linspace(0,1,100),'pchip')' ...
                   interp1(linspace(0,1,length(pos)),pos(:,2),linspace(0,1,100),'pchip')']

        j = 1;
        boundary = [];
        while (isempty(boundary) & j <= length(edge_i))
            boundary = bwtraceboundary(BW,[fix(edge_i(j,2)), fix(edge_i(j,1))],'S');
            j = j + 1;

        end

        if isempty(boundary)
            if ( strcmpi(questdlg('Boundary not found! Make sure that a line is drawn in the south-north direction across the boundary. Try drawing again?', 'Repeat trace?','Yes','No','Yes'), 'no') );
                warndlg('No boundary coordinates have been saved.');
                quit_boundary = 1;
                im.BN = -1;
            end
        else
            quit_boundary = 1;
            im.BN = [boundary(:,2) boundary(:,1)];
            figure(hp);
            subplot(1,2,[1 2]);
            imagesc(im.I_c); hold on; plot(im.BN(:,1),im.BN(:,2)); hold off; axis image;
        end
    end
else
    im = [];
end




