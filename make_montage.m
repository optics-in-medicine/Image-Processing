


clear all
close all

name_vec = {'A_1','A_2','B_1','B_2','C_1','C_2','D_1','D_2'};
file_dir = 'X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 3 (F301-F310) WT\Histopathology B02-03\Histopathology\hico\'


I = [];
for i = 1:4:42
    R = [];
    for j = 1:4
        try
        P1 = imread( [file_dir 'slide' int2str(i+j-1) '-1.tif'] );
        catch
            P1 = ones(1200,1600,3);
        end
        
        try
            P2 = imread( [file_dir 'slide' int2str(i+j-1) '-2.tif'] );
        catch
            P2 = ones(1200,1600,3);
        end
%         P = imread( [file_dir '4' num2str(i,'%02.0f') name_vec{j} '.tif'] );
        R = [R; [P1; P2]]; 
    end
    I = [I R];
end

imwrite(I,'montage_b02_03.png')



%% Fluorescence image montage
clear all
close all

file_dir = 'X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 5 (F501-F512) DEX\';

cx = 350 .* [2 2 2 2 2 2 2 2 2 2];
cy = 125 .* [2 2 2 2 2 2 2 2 2 2]; 

I7 = []; I8 = [];
for i = [1 2 4 5 6 7 8 9 10 3]
    R7 = [];
    R8 = [];
    for j = 1:2
        subdir = ['F5' num2str(i,'%02.0f') '\Odyssey\'] 
        files = dir([file_dir subdir])
        for k = 3:length(files)
            filename = files(k).name;
            if strcmpi( filename(findstr(filename,'.')-3:findstr(filename,'.')-1), '700' )
                P7 = double( adapthisteq( imread( [file_dir subdir filename] ) ) );
                P8 = double( adapthisteq( imread( [file_dir subdir strrep(filename, '_700','_800')] ) ) );
  
                [X, Y] = size(P7); 
                a = prctile(P7(:),20);
                b = prctile(P8(:),20);
             
                PP7 = [a.* ones(200,Y+400); [a.* ones(X,200) P7 a.* ones(X,200)]; a.* ones(200,Y+400)];
                PP8 = [b.* ones(200,Y+400); [b.* ones(X,200) P8 b.* ones(X,200)]; b.* ones(200,Y+400)];
                
                figure; imagesc(PP7); h = impoint(gca,[]); pos = wait(h);
                PPP7 =  PP7( fix(pos(2)) - cx(i) : fix(pos(2)) + cx(i), fix(pos(1)) - cy(i): fix(pos(1)) + cy(i) );
                PPP8 =  PP8( fix(pos(2)) - cx(i) : fix(pos(2)) + cx(i), fix(pos(1)) - cy(i): fix(pos(1)) + cy(i) ); 
               
                R7 = [R7; flipud(imresize(PPP7,[cx(1)*2 cy(1)*2]))];
                R8 = [R8; flipud(imresize(PPP8,[cx(1)*2 cy(1)*2]))];
                
            end
        end
    end
    I7 = [I7 R7];
    I8 = [I8 R8];
end

% set lim
I7adj = (I7 - 2.25e4) ./ (3.25e4 - 2.25e4);
I7adj(I7adj>1) = 1; I7adj(I7adj<0) = 0;
imagesc(I7adj); colormap(goodmap('cube'))

[ mont7 ] = gen_cmap_image( I7adj, linspace(0,1,256)'*[1 1 1] );

I8adj = (I8 - 2.25e4) ./ (3.75e4 - 2.25e4);
I8adj(I8adj>1) = 1; I8adj(I8adj<0) = 0;
imagesc(I8adj); colormap(goodmap('cube'))

[ mont8 ] = gen_cmap_image( I8adj,  linspace(0,1,256)'*[1 1 1] );

imwrite(mont7,'montage_700_dex.tiff')
imwrite(mont8,'montage_800_dex.tiff')




%% Fluorescence image montage
clear all
close all

file_dir = 'X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 5 (F501-F512) DEX\';

cx = 350 .* [2 2 2 2 2 2 2 2 2 2];
cy = 125 .* [2 2 2 2 2 2 2 2 2 2]; 

I7 = []; 
for i = [1 2 4 5 6 7 8 9 10 3]
    R7 = [];
    for j = 1:2
        subdir = ['F5' num2str(i,'%02.0f') '\Typhoon\']
        files = dir([file_dir subdir])
        for k = 3:length(files)
            filename = files(k).name;
            if strcmpi( filename(findstr(filename,'.')+1:findstr(filename,'.')+3), 'gel' )
                P7 = double( adapthisteq( imread( [file_dir subdir filename] ) ) );
              
                [X, Y] = size(P7); 
                if X > Y
                    P7 = rot90(P7);
                    [X, Y] = size(P7);
                end
                
                a = prctile(P7(:),20);
             
                PP7 = [a.* ones(200,Y+400); [a.* ones(X,200) P7 a.* ones(X,200)]; a.* ones(200,Y+400)];
                
                figure; imagesc(PP7); h = impoint(gca,[]); pos = wait(h);
                PPP7 =  PP7( fix(pos(2)) - cx(i) : fix(pos(2)) + cx(i), fix(pos(1)) - cy(i): fix(pos(1)) + cy(i) );
               
                R7 = [R7; imresize(PPP7,[cx(1)*2 cy(1)*2])];
            end
        end
    end
    I7 = [I7 R7];
end


% I7f = [flipud(I7(1:2800,:)); flipud(I7(2801:end,:))];
% 

% set lim
I7adj = (I7 - 1.0e4) ./ (5.5e4 - 1e4);
I7adj(I7adj>1) = 1; I7adj(I7adj<0) = 0;
imagesc(I7adj)

[ mont9 ] = gen_cmap_image( I7adj, goodmap('ppix') );

imwrite(mont9,'montage_ppix_dex.tiff')


