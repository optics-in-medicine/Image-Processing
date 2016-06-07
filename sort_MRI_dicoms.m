function [stat, pat_name,prot_name] = sort_MRI_dicoms( folder_path )
%SORT_MRI_DICOMS Will take a folder of DICOMS and sort them into subfolders according to Patient ID and Protocol Name.
%Uses PatientID and ProtocolName dicom header fields. These must be 
%specified during the scan set-up.
%
%   Author: Jonathan T. Elliott <jte@dartmouth.edu>
%   Date:   Sept. 30, 2015


% list directory contents in structure, 'files'.
files = dir(folder_path);

try
    for i = 3:length(files)
        
        % load dicom header to pull patient_ID and protocol_name
        info = dicominfo([folder_path filesep files(i).name]);
        
        % store the labels for output
        prot_name{i-2} = info.ProtocolName;
        pat_name{i-2} = info.PatientID;
        
        % copy the dicom file to the new folder.
        copyfile([folder_path filesep files(i).name '*'], [folder_path filesep pat_name{i-2} filesep prot_name{i-2}]);
    end
    
    stat = 1;           % made it through without error.
catch
    stat = 0;           % flag error so that user knows.
end




