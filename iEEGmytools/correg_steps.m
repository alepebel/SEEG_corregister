ft_defaults

subjID = 'sub_01'

path_anats = ['/Users/alex/Dropbox/Projects/Predcision/Raw/',subjID,'/anat/']; 

acpc_mri_path = fullfile(path_anats,  'PreMR_acpc.nii')

%% Preprocessing MRI
if ~exist(acpc_mri_path)
    
    p_mri = ft_read_mri(fullfile(path_anats, 'PreMRI.nii'))      
    %Determine what axis is the horizontal
    %Additionally write down if the values on the left?right axis increase
    %to the right (indicated by a ?+? sign), then the scan has a left-to-right
    %orientation. If the values on the left?right axis increase to the left,
    %then the scan has a right-to-left orientation.
    
    % visualizing
    ft_determine_coordsys(p_mri)  % x goes from left to right (positive values at left side of the head)
    
    % Selecting, what is the right hemisphere (R), anterior (A), posterior comisure
    % (P) and top part of the brain (Z). Quit (Q)
    cfg = [];
    cfg.method = 'interactive';
    cfg.coordsys = 'acpc';
    p_mri_acpc = ft_volumerealign(cfg, p_mri);
    
    cfg = [];
    cfg.filename = acpc_mri_path;
    cfg.filetype = 'nifti';
    cfg.parameter = 'anatomy';
    ft_volumewrite(cfg, p_mri_acpc);
else
    p_mri_acpc = ft_read_mri(acpc_mri_path)      
end

%% Running freesurfer parcellation
fshome = '/Applications/freesurfer'; 
subdir = path_anats; 
mrfile = acpc_mri_path; 

system(['export FREESURFER_HOME=' fshome '; ' ...
'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
'mri_convert -c -oc 0 0 0 ' mrfile ' ' [subdir '/tmp.nii'] '; ' ...
'recon-all -i ' [subdir '/tmp.nii'] ' -s ' 'freesurfer' ' -sd ' subdir ' -all'])


ft_sourceplot(p_mri_acpc)

%% Preprocessing CT acan

acpc_ct_path = fullfile(path_anats ,'ct_acpc.nii')

if ~exist(acpc_ct_path)
    
    ct = ft_read_mri(fullfile(path_anats, 'CT.nii'));
    
    cfg = [];
    cfg.method = 'interactive';
    cfg.coordsys = 'ctf';
    ct_ctf = ft_volumerealign(cfg, ct);
    
    ct_acpc = ft_convert_coordsys(ct_ctf, 'acpc',  [], 0); % if you dont have spm installed it seems that you need to use method 0
    
    ft_volumerealign(cfg, ct_acpc)
    
    cfg = [];
    cfg.filename = acpc_ct_path;
    cfg.filetype = 'nifti';
    cfg.parameter = 'anatomy';
    ft_volumewrite(cfg, ct_acpc);
    
else
   ct_acpc = ft_read_mri(acpc_ct_path)    
end

% fusion
cfg = []; 
cfg.method = 'spm'; 
cfg.spmversion = 'spm12'; 
cfg.coordsys = 'acpc'; 
cfg.viewresult = 'yes';
ct_acpc_f = ft_volumerealign(cfg, ct_acpc, p_mri_acpc); 