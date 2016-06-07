function [no_images] = improve_contrast( file_dir, FL )
%RESAVE_USCOPE_TIFFS Loads the images from the Borwell microscope, applies
%a field correction (specified, stored from before as .mat, or ignore), and
%then stretches the lightness channel to utilize the full dynamic range
%thereby enhancing the contrast of the image.
%
% Author: Prof. Jonathan T. Elliott
%         Thayer School of Engineering at Dartmouth
%         <jte@dartmouth.edu>
%
% Date:   2015-Oct-22
%

% if no input arguments are put in, then define file_dir as empty string.
if nargin == 0
    file_dir = [];
end

% if you put an empty string '[]' for file_dir, a dialog will pop up.
if isempty(file_dir)
    file_dir = uigetdir('C:\Users\Clare\Documents\MATLAB\# # Coregistration (2)');
else
    % the code assumes that there is no terminal "filesep" (i.e. \ or /).
    % This removes if there is one.
    if strcmpi(file_dir(end), filesep)
        file_dir = file_dir(1:end-1);
    end
end

%select directory
files = dir(file_dir);

% load FL image
if nargin == 1
    % if no field correction is available, ignore the nonuniformity.
    LPF = zeros(size(imread( [file_dir filesep files(3).name] )));
else
    if ischar(FL)
        if strcmpi( FL(end-2:end), 'mat' )
            % if the data is already stored in a .mat file, load it.
            data_in = load(FL);
            LPF = data_in.FL;
            clear data_in
        else strcmpi( FL(end-2:end), 'tif' ) | strcmpi( FL(end-2:end), 'iff')
            % assumes the same 12-bit tiff images that the microscope
            % takes.
            LPF = uint8( imread( FL ) .* (2^12 / 2^16) );
        end
    else
        % if the variable 'FL' specifies an actual image element.
        [X, Y, Z] = size( FL );
        if Z == 1   % make sure its an RGB image, if not, assume monochrom.
            LPF = cat(3, FL, FL, FL);
        elseif Z == 3   
            LPF = FL;
        else       % Can't catch them all!
            LPF = cat(3, FL(:,:,1),FL(:,:,1),FL(:,:,1));
        end
    
    end
end
% convert the light-field image to CIELAB space and separate channels into
% 3 variables.
LPF_lab = squeeze( applycform( im2double( LPF ), makecform('srgb2lab') ));
LPF_L = squeeze( LPF_lab(:,:,1) );
LPF_a = squeeze( LPF_lab(:,:,2) );
LPF_b = squeeze( LPF_lab(:,:,3) );

% make the directory "hico" and return status to proceed.
mk_dir_stat = mkdir([file_dir filesep 'hico'])
no_images = 0; % initiate image counter.

if mk_dir_stat
    for i = 3:length(files)
        % convert 16-bit to 8-bit for color processing
        filename = files(i).name;
        if strcmpi( filename(end-2:end), 'tif' ) | strcmpi( filename(end-2:end), 'iff')
            
            % keep track of the number of images processed.
            no_images = no_images + 1;

            % load the i'th image.
            F = uint8( imread( [file_dir filesep filename] ) .* (2^12 / 2^16) );
            
            % convert to labspace.
            F_lab = applycform( im2double(F ), makecform('srgb2lab'));

            % separate CIELAB channels into separate variables.
            L = squeeze(F_lab(:,:,1)); 
            a = squeeze(F_lab(:,:,2));
            b = squeeze(F_lab(:,:,3));

            % normalize illumination field for L channel, and then do some
            % autocontrast enhancement.
            L_temp = L - LPF_L + mean( LPF_L(:) );
            L_ad = adapthisteq((L_temp ./ 100)); %assumes between 0 and 1.
            L_c = L_ad .* 100;  % multiply back to 0 - 100.

            % concatenate the enhanced CIELAB channels.
            Lab_adj = cat(3, L_c, ...
                             a - LPF_a + mean(LPF_a(:)), ...
                             b - LPF_b + mean(LPF_b(:)));

            % convert back to RGB space
            rgb_adj = applycform( Lab_adj, makecform('lab2srgb'));
    
            % save image as 8-bit tiff to simplify viewing on PC/Mac.
            imwrite(rgb_adj, [file_dir filesep 'hico' filesep filename]);
        end
    end
    % end with a msgbox stating how many images were processed.
    msgbox(['Contrast enhancement complete. Processed and saved ' int2str(no_images) ' images.'], 'Processing Complete');
else
    % warn that the processing failed because new directory was not
    % created.
    warndlg('Could not create new directory. Images not processed');
end

