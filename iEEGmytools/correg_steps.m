ft_defaults

subjID = 'sub_01'

path_anats = ['/Users/alex/Dropbox/Projects/Predcision/Raw/',subjID,'/anat/']; 

%path_anats = ['~/Insync/alexisperez@ub.edu/OneDrive Biz/Projects/corregisLudo/S2']; 

path_anats = ['/Users/alex/Library/CloudStorage/OneDrive-UniversitatdeBarcelona/Projects/IDIBAPS/Exp_1_WM_tms/Raw/S01/anat']; 

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
fshome = '/Applications/freesurfer/7.2.0'; 
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
    
    ct_acpc = ft_convert_coordsys(ct_ctf, 'acpc',  0); % if you dont have spm installed it seems that you need to use method 0
    
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




% Loading header
%subj_dir = fullfile('/users' ,'alex', 'Dropbox', 'Projects', 'Predcision','Raw', subjID, 'func')%
subj_dir = fullfile('/Users/alex/Library/CloudStorage/OneDrive-UniversitatdeBarcelona/Projects/IDIBAPS/Exp_1_WM_tms/Raw/S01/func')%

fileList = dir([subj_dir, '/*.edf'])

hdr = ft_read_header(fullfile(subj_dir,fileList(1).name))


% Electrodes placement

cfg         = [];
cfg.channel = hdr.label;
elec_acpc_f = ft_electrodeplacement(cfg, ct_acpc_f, p_mri_acpc);
load('elecS01.mat')
save('elecS01', 'elec_acpc_f')

% Checking correct corregister
ft_plot_ortho(p_mri_acpc.anatomy, 'transform',p_mri_acpc.transform, 'style', 'intersect');
ft_plot_sens(elec_acpc_f, 'label', 'on', 'fontcolor', 'w');





%% Creating table with electrodes locations

% loading atlases and localizing electrodes
aparcaseg = ft_read_atlas('/Users/alex/Library/CloudStorage/OneDrive-UniversitatdeBarcelona/Projects/IDIBAPS/Exp_1_WM_tms/Raw/S01/anat/freesurfer/mri/aparc+aseg.mgz');
aparcaseg.coordsys = 'acpc';



cfg               = [];
cfg.roi           = elec_acpc_f.chanpos; % from elec_nat
cfg.inputcoord    = 'acpc';
cfg.atlas         = aparcaseg;
cfg.output        = 'multiple'; % since v2
labels_aparcaseg = ft_volumelookup(cfg, aparcaseg);


[~, indx] = max(labels_aparcaseg.count);
labels_aparcaseg.name(indx)


% loading atlas aparc 2009 and localizing electrodes
aparc2009 = ft_read_atlas('/Users/alex/Library/CloudStorage/OneDrive-UniversitatdeBarcelona/Projects/IDIBAPS/Exp_1_WM_tms/Raw/S01/anat/freesurfer/mri/aparc.a2009s+aseg.mgz');
aparc2009.coordsys = 'acpc';



cfg               = [];
cfg.roi           = elec_acpc_f.chanpos; % from elec_nat
cfg.inputcoord    = 'acpc';
cfg.atlas         = aparc2009;
cfg.output        = 'multiple'; % since v2
labels_aparc2009 = ft_volumelookup(cfg, aparc2009);


[~, indx] = max(labels_aparc2009.count);
labels_aparc2009.name(indx)


% loading atlas DSK and localizing electrodes
DKT = ft_read_atlas('/Users/alex/Library/CloudStorage/OneDrive-UniversitatdeBarcelona/Projects/IDIBAPS/Exp_1_WM_tms/Raw/S01/anat/freesurfer/mri/aparc.DKTatlas+aseg.mgz');
DKT.coordsys = 'acpc';



cfg               = [];
cfg.roi           = elec_acpc_f.chanpos; % from elec_nat
cfg.inputcoord    = 'acpc';
cfg.atlas         = DKT;
cfg.output        = 'multiple'; % since v2
labels_DKT = ft_volumelookup(cfg, DKT);


tabledata = {}

for e = 1 : length(labels_DKT)
    ch_name = elec_acpc_f.label{e}
    chanpos = elec_acpc_f.chanpos(e,:)

    [cnt, idx] = max(labels_DKT(e).count);
    DKTlab  = char(labels_DKT(e).name(idx)); % anatomical label

    [cnt, idx] = max(labels_aparcaseg(e).count);
    aparclab  = char(labels_aparcaseg(e).name(idx)); % anatomical label

    [cnt, idx] = max(labels_aparc2009(e).count);
    aparc2009lab  = char(labels_aparc2009(e).name(idx)); % anatomical label

    aux = vertcat([{ch_name},DKTlab, aparclab , aparc2009lab, chanpos(1), chanpos(2), chanpos(3)]);
    tabledata = vertcat(tabledata,aux);
    %disp(lab)
end

T = cell2table(tabledata, "VariableNames",["Chan" "DKT" "Aparc" "Aparc2009" "posx" "posy" "posz"])

writetable(T,'s01iEEG_correg.xlsx')


%% Plotting electrodes and ROI



pial_lh = ft_read_headshape('freesurfer/surf/lh.pial');
pial_lh.coordsys = 'acpc';
ft_plot_mesh(pial_lh);
lighting gouraud;
camlight;


pial_rh = ft_read_headshape('freesurfer/surf/lh.pial');
pial_rh.coordsys = 'acpc';
ft_plot_mesh([pial_rh,pial_lh]);
lighting gouraud;
camlight;


atlas = ft_read_atlas('freesurfer/mri/aparc+aseg.mgz');
atlas.coordsys = 'acpc';
cfg            = [];
cfg.inputcoord = 'acpc';
cfg.atlas      = atlas;
cfg.roi        = {'Left-Hippocampus', 'Right-Hippocampus'};
mask_rha = ft_volumelookup(cfg, atlas);


seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain = mask_rha;
cfg             = [];
cfg.method      = 'iso2mesh';
cfg.radbound    = 2;
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 1000;
cfg.smooth      = 3;
cfg.spmversion  = 'spm12';
mesh_rha = ft_prepare_mesh(cfg, seg);



cfg              = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim  = [-.5 .5];
cfg.method       = 'cloud';
cfg.slice        = '3d';
cfg.nslices      = 2;
cfg.facealpha    = .25;
ft_sourceplot(cfg, [], mesh_rha);
view([120 40]);
lighting gouraud;
camlight;


% Creating hippocampus masks
atlas = ft_read_atlas('freesurfer/mri/aparc.a2009s+aseg.mgz');
atlas.coordsys = 'acpc';
cfg            = [];
cfg.inputcoord = 'acpc';
cfg.atlas      = atlas;
cfg.roi        = { 'Right-Hippocampus'};
mask_rha = ft_volumelookup(cfg, atlas);


seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain = mask_rha;
cfg             = [];
cfg.method      = 'iso2mesh';
cfg.radbound    = 2;
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 500;
cfg.smooth      = 3;
cfg.spmversion  = 'spm12';
mesh_rha = ft_prepare_mesh(cfg, seg);


atlas.coordsys = 'acpc';
cfg            = [];
cfg.inputcoord = 'acpc';
cfg.atlas      = atlas;
cfg.roi        = { 'Left-Hippocampus'};
mask_lha = ft_volumelookup(cfg, atlas);


seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain = mask_lha;
cfg             = [];
cfg.method      = 'iso2mesh';
cfg.radbound    = 2;
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 500;
cfg.smooth      = 3;
cfg.spmversion  = 'spm12';
mesh_lha = ft_prepare_mesh(cfg, seg);




figure('color','w')
fspial_lh = ft_read_headshape([fshome '/subjects/fsaverage/surf/lh.pial']);
fspial_lh.coordsys = 'acpc';
fspial_rh = ft_read_headshape([fshome '/subjects/fsaverage/surf/rh.pial']);
fspial_rh.coordsys = 'acpc';
ft_plot_mesh(fspial_lh, 'facealpha', 0.6, 'facecolor', 'skin');
ft_plot_mesh(fspial_rh, 'facealpha', 0.6, 'facecolor', 'skin');
ft_plot_mesh(mesh_rha,'facecolor', "#80B3FF", 'edgecolor', 'none', 'facealpha', 0.75);
ft_plot_mesh(mesh_lha,'facecolor', "#80B3FF", 'edgecolor', 'none', 'facealpha', 0.75);
ft_plot_sens(elec_acpc_a );
view([-90 -20]);
material dull;
lighting gouraud;
camlight;

elec_mni_fstd

elec_acpc_a = elec_acpc_f
elec_acpc_a.label = elec_mni_fstd(1).label
elec_acpc_a.elecpos = elec_mni_fstd(1).elecpos
elec_acpc_a.chanpos = elec_mni_fstd(1).chanpos
elec_acpc_a.chantype = elec_mni_fstd(1).chantype
elec_acpc_a.tra = elec_acpc_f.tra(1:2,1:2)
ft_plot_sens(elec_acpc_f);
 
elec_acpc_f




%% visualize the results like Ludo
rmpath(genpath('/media/ludovico/DATA/iEEG_Ludo/spm12')) %remove spm 12 to avoid conflicts with fieldtrip
[ftver, ftpath] = ft_version;

% in surface pial

load([ftpath filesep 'template/anatomy/surface_pial_left.mat']);
template_lh = mesh; %clear mesh;

load([ftpath filesep 'template/anatomy/surface_pial_right.mat']);
template_rh = mesh; %clear mesh;

% in mni volume

mni_path = fullfile(ftpath , 'template/anatomy/single_subj_T1_1mm.nii');
mni = ft_read_mri(mni_path);

%% create subcortical space
%right hippomcapus
atlas = ft_read_atlas(fullfile(ftpath, 'template/atlas/aal/ROI_MNI_V4.nii'));
atlas = ft_convert_units(atlas, 'mm');
atlas.coordsys = 'mni';
cfg            = [];
cfg.inputcoord = 'mni';
cfg.atlas      = atlas;
cfg.roi        = {'Hippocampus_R'}; %change name here if other area is interesting
mask_rha = ft_volumelookup(cfg, atlas);


rmpath(genpath('/media/ludovico/DATA/iEEG_Ludo/spm12'))
seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain = mask_rha;
cfg             = [];
cfg.method      = 'iso2mesh';
cfg.radbound    = 2;
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 1000;
cfg.smooth      = 3;
cfg.spmversion  = 'spm12';
mesh_rha1 = ft_prepare_mesh(cfg, seg);

%left hippomcapus
atlas = ft_read_atlas(fullfile(ftpath, 'template/atlas/aal/ROI_MNI_V4.nii'));
atlas = ft_convert_units(atlas, 'mm');
atlas.coordsys = 'mni';
cfg            = [];
cfg.inputcoord = 'mni';
cfg.atlas      = atlas;
cfg.roi        = {'Hippocampus_L'};%change name here if other area is interesting
mask_rha = ft_volumelookup(cfg, atlas);


rmpath(genpath('/media/ludovico/DATA/iEEG_Ludo/spm12'))
seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain = mask_rha;
cfg             = [];
cfg.method      = 'iso2mesh';
cfg.radbound    = 2;
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 1000;
cfg.smooth      = 3;
cfg.spmversion  = 'spm12';
mesh_rha2 = ft_prepare_mesh(cfg, seg);

% Define colors for each participant electrodes...

cmap = brewermap(10, 'Paired');

%% Visualize in 3D whole brain
figure;
ft_plot_mesh(template_lh, 'facealpha', 0.1);
ft_plot_mesh(template_rh, 'facealpha', 0.1);
ft_plot_mesh(mesh_rha1, 'facecolor', [0 0.7 1], 'facealpha', 0.2, 'edgealpha', 0)
ft_plot_mesh(mesh_rha2, 'facecolor', [0 0.7 1], 'facealpha', 0.2, 'edgealpha', 0)


% When you load the electrodes, you can use the information from iEEGcorr
% in order to filter only a specific set of a electrode channels.
for subj = 1 : length(subjc)
    ft_plot_sens(elec_mni_fstd(subj), 'fontcolor' , 'white', 'fontsize', 10, 'style',cmap(subj,:), 'elecsize',30);
end

view([110 20]); %change orientation here. [-110 20] to see from the other side
material dull;
lighting gouraud;
camlight;
set(gcf, 'color', 'white');
set(gcf, 'renderer', 'Painters')
set(gcf, 'color', 'white');



%% Visualize in 3D only hipp
figure;
% ft_plot_mesh(template_lh, 'facealpha', 0.1);
% ft_plot_mesh(template_rh, 'facealpha', 0.1);
ft_plot_mesh(mesh_rha1, 'facecolor', [0 0.7 1], 'facealpha', 0.2, 'edgealpha', 0)
ft_plot_mesh(mesh_rha2, 'facecolor', [0 0.7 1], 'facealpha', 0.2, 'edgealpha', 0)

% When you load the electrodes, you can use the information from iEEGcorr
% in order to filter only a specific set of a electrode channels.
for subj = 1 : length(subjc)
    ft_plot_sens(elec_mni_fstd(subj), 'fontcolor' , 'white', 'fontsize', 10, 'style',cmap(subj,:), 'elecsize',30);
end

view([-120 40]);
material dull;
lighting gouraud;
camlight;
set(gcf, 'renderer', 'Painters')
set(gcf, 'color', 'white');
%
% %% Visualize in volume space
% figure;
% ft_plot_ortho(mni.anatomy, 'transform', mni.transform, 'style', 'intersect');
% for subj = 1 : length(subjc)
%     ft_plot_sens(elec_mni_fstd(subj), 'label', 'on','fontcolor' , 'white', 'fontsize', 10, 'style',cmap(subj,:), 'elecsize',50);
% end
%
%
% %ft_plot_sens(elec_mni_frv, 'label', 'on', 'fontcolor', 'w','style', 'b');