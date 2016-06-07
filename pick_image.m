function [ I, imstat, PathName, FileName ] = pick_image( PathName )
%PICK_IMAGE Summary of this function goes here
%   Detailed explanation goes here

image_types = {'Matlab Datafile (.mat)', 'Pearl 22-bit (.tiff)', 'Odyssey 16-bit (.tiff)', 'Typhoon (.gel)', 'Zeiss HIS (.tiff)', 'Zeiss JPEG (.jpg)', 'MRI DICOM', 'Other Image Format', 'DHMC FGR structure'};
image_ext = {'.mat','.TIF','.TIF','.gel','.tiff','*.jpg', '*.*', '*.*', '*.mat'};
[s,v] = listdlg('PromptString','Choose Image Type:', 'SelectionMode','single','ListString', image_types);

if v
    [FileName, PathName, filterIndex] = uigetfile(image_ext{s}, 'Open ...', PathName);
    
    if filterIndex ~= 0
        
        switch s
            case {5, 6, 8}
                I = imread([PathName FileName]);    
            case {2, 3, 4}
                I_bf = bfopen([PathName FileName]);
                I = double( I_bf{1}{1} );
            case 7
                %dicom reader
            case 1
                data_in = load([PathName FileName]);
                var_names = fieldnames( data_in );
                [s1, v1] = listdlg('PromptString','Select variable:', 'SelectionMode','single','ListString', var_names);
                I = data_in.( var_names{s1} );
    %             var_choice = strcat('data_in.', var_names(s1));
    %             eval( strcat('I = data_in.', var_names(s1), ';') );
            case 9
                data_in = load([PathName FileName]);
                var_names = fieldnames( data_in );
                [s1, v1] = listdlg('PromptString','Select variable:', 'SelectionMode','single','ListString', var_names);
                I = data_in.( var_names{s1} );
                [X, Y, Z] = size(I);

                if Z ~= 3
                    I_t = I;
                    clear I;
                    imn = 1;
                    if Z ~= 1
                        imn = inputdlg(['Which image in this series of ' int2str(Z) ' would you like?']);
                    end
                    I = squeeze( I_t(:,:,imn) );
                end
        end
        imstat = 1;
    else
        imstat = 0; I = []; PathName = [];
    end
else
    imstat = 0; I = []; PathName = [];
end

end

