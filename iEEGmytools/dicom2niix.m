ft_defaults


% here is where i stored the folders with the dicom data
path_dir = '/Users/alex/Dropbox/Projects/Predcision/Raw/sub_01/anat/dicom/' 

% lets list the folders within the dicom folder to find the dicoms inside
fnames = get_subfolders(path_dir)

for ixfname = 1 : length(fnames) % explore first level of the folder jerarchy
    
    sub_path_dir =  fullfile(path_dir, fnames(ixfname).name)    
    ffnames = get_subfolders(sub_path_dir)
    for ixxfname = 1 : length(ffnames) % explore second level of the folder jerarchy
        
        scan_path_dir = fullfile(sub_path_dir, ffnames(ixxfname).name);
        dicom_files = dir(scan_path_dir);
        dicom_path = fullfile(scan_path_dir, dicom_files(3).name);
        
        %% Reading dicom and writting nifti files in folders
        
        % some files might fail (probably not relevant scans like the pre-imp MRI and the CT). Check the final result
        try % Skipping on errors 
        [mri] = ft_read_mri(dicom_path, 'dataformat', 'dicom')       
        scan_name = mri.hdr(1).SeriesDescription(~isspace(mri.hdr(1).SeriesDescription)) 
        scan_name = strrep(scan_name,'.','') % removing points to avoid conflicts
        cfg = []
        cfg.parameter     =  'anatomy'
        cfg.filename      =  fullfile(sub_path_dir, scan_name)
        cfg.filetype      = 'nifti'
        ft_volumewrite(cfg, mri)
        end
    end
end
