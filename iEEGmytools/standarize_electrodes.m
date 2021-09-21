% Alexis Pérez Bellido (2021)
% To use this script you only need to loop over the participants folders loading the file with the localised electrodes and the preimplantation T1 anat scans.


ft_defaults


subjs = {'S2', 'S3', 'S5'};

path_anats = ['~/Insync/alexisperez@ub.edu/OneDrive Biz/Projects/corregisLudo/'];


for isubj = 1 : length(subjs)
    
    subj_path = fullfile(path_anats, char(subjs(isubj)));
    
    
    % loading acpc alligned preimplantation T1
    p_mri_acpc = ft_read_mri(fullfile(subj_path, 'MR_acpc.nii'));
    
    % loading acpc alligned localized electrodes
    
    electrodes_path = fullfile(subj_path, 'elec_acpc_f.mat');
    load(electrodes_path)
    
    
    % Loading labels info in order to
    iEEGcorr_path = fullfile(subj_path, 'iEEG_correg.mat');
    iEEGcorr(isubj) = load(iEEGcorr_path);
    
    % Here you can check whether the electrodes are correctly alligned to
    % the native space T1 image
    
    %figure;
    %ft_plot_ortho(p_mri_acpc.anatomy, 'transform', p_mri_acpc.transform, 'style', 'intersect');
    %ft_plot_sens(elec_acpc_f, 'label', 'on', 'fontcolor', 'w', 'fontsize', 6, 'style','r', 'elecsize', 10);
    
    % Standarizing electrodes
    
    p_mri_acpc.coordsys = 'acpc';
    
    cfg            = [];
    cfg.elec       = elec_acpc_f;
    cfg.method     = 'mni';
    cfg.mri        = p_mri_acpc;
    cfg.spmversion = 'spm12';
    cfg.spmmethod  = 'new';
    cfg.nonlinear  = 'yes';
    elec_mni_fstd(isubj) = ft_electroderealign(cfg);
    elec_data =  elec_mni_fstd(isubj);
    elect_savepath = fullfile( subj_path, 'elec_mni_fstd.mat');
    
    save(elect_savepath, 'elec_data')
    
    
end



%% visualize the results
[ftver, ftpath] = ft_version;

% in surface pial

load([ftpath filesep 'template/anatomy/surface_pial_left.mat']);
template_lh = mesh; %clear mesh;

load([ftpath filesep 'template/anatomy/surface_pial_right.mat']);
template_rh = mesh; %clear mesh;

% in mni volume

mni_path = fullfile(ftpath , 'template/anatomy/single_subj_T1_1mm.nii');
mni = ft_read_mri(mni_path);


% Define colors for each participant electrodes...
cmap = [0, 0.5, 0; 0.5, 0.25, 0 ; 0, 0.5, 1; 0.8, 0.5, 0];


%% Visualize in 3D
figure;
ft_plot_mesh(template_lh, 'facealpha', 0.4);
ft_plot_mesh(template_rh, 'facealpha', 0.4);

% When you load the electrodes, you can use the information from iEEGcorr
% in order to filter only a specific set of a electrode channels.

for isubj = 1 : length(subjs)
    ft_plot_sens(elec_mni_fstd(isubj), 'fontcolor' , 'white', 'fontsize', 10, 'style',cmap(isubj,:), 'elecsize',10);
end

view([-90 20]);
material dull;
lighting gouraud;
camlight;


%% Visualize in volume space
figure;
ft_plot_ortho(mni.anatomy, 'transform', mni.transform, 'style', 'intersect');
for isubj = 1 : length(subjs)
    ft_plot_sens(elec_mni_fstd(isubj), 'label', 'on','fontcolor' , 'white', 'fontsize', 5, 'style',cmap(isubj,:), 'elecsize',15);
end
  

%ft_plot_sens(elec_mni_frv, 'label', 'on', 'fontcolor', 'w','style', 'b');
