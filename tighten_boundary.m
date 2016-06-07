function [ im ] = tighten_boundary( im )
%TIGHTEN_BOUNDARY Summary of this function goes here
%   Detailed explanation goes here
    
    im_gs = nm( double( rgb2gray( im.I_c) ) );  
    figure;
    h = imagesc(im_gs); axis image; colormap(goodmap('hot'));
   
    hfig = imcontrast(h);
    waitfor(hfig);

    I1_adj = getimage(h);
    I_n = I1_adj;
    
    hp = figure;
    subplot(1,2,1);
    imagesc(im.I_c); axis image; colormap(goodmap('santa'));
    
    h2 = subplot(1,2,2);
    
    BW = imfill( im2bw(abs(1-I_n)), 'holes') ; 
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

end

