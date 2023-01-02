clear
cd ('/media/ludovico/DATA/iEEG_Ludo')
addpath ('/media/ludovico/DATA/iEEG_Ludo/toolbox2.0')
addpath('/media/ludovico/DATA/iEEG_Ludo/fieldtrip-20200607')
addpath('/media/ludovico/DATA/iEEG_Ludo/ElecPlacementToolbox')
addpath ('/media/ludovico/DATA/iEEG_Ludo/spm12')
%If you do not have Paris aptients comment this
ft_defaults;
%% important variables
%For each subject change those. 

Subject                     = 19;
config.mriname              = 'T1.nii'
config.CTname               = 'CT.nii'

%% main script

sub = {'07', '08', '09', '10', '12', '13', '14', '15' '16', '17' '18' '19'};
sublinic = [1 2 3 4 5 6 7 8 9 10 11 13; 7 8 9 10 12 13 14 15 16 17 18 19]'; 

%Root file of the raw data of clinic and Paris patients, ignore paris
%patients if you do not have them
config.Home             = '/media/ludovico/DATA/iEEG_Ludo/';
config.HomeClinic       = '/media/ludovico/DATA/iEEG_Ludo/Globus/HospitalClinic';
config.SubjnameBIDS     = sprintf('sub-%s', sub{find(sublinic(:, 2) == Subject)})
config.EEG              = sprintf('/media/ludovico/DATA/iEEG_Ludo/BIDS/%s/ses-Day_1/eeg/%s_ses-Day_1_task-MemSeqIntra_eeg.EDF',...
    config.SubjnameBIDS, config.SubjnameBIDS)
%Where to save your file
config.Homesave         = '/media/ludovico/DATA/iEEG_Ludo/Results';
config.Structural       = 'Structural'
config.SubjID           = sprintf ('Subject_%d', Subject);
config.subjname         = config.SubjID;
config.repreproc        = 0;
config.ResDir           = 'Results';
%get the names of the CT and MRI



config.HomeSubj             = fullfile(config.Home, config.ResDir, config.SubjID, config.Structural);
config.fullEEG              = fullfile(config.EEG);
config.mrifile              = fullfile(config.Home, config.ResDir, config.SubjID,  config.Structural, config.mriname);

savename                        = 'MR_acpc.nii';
Savedir                         = fullfile(config.Homesave, config.SubjID, config.Structural, savename);
%if mri has already been processed load it otherwise process it or
%reprocess it if config.repreproc = 1
if exist(Savedir) > 0 || config.repreproc == 1
    fsmri_acpc                      = ft_read_mri(Savedir);
    fsmri_acpc.coordsys             = 'acpc';
else
    fsmri_acpc                      = MRIMakeSave(config)
end

config.fsmri_acpc                   = fsmri_acpc;

config.freesurferdir = 'freesurfer'
Savedir = fullfile(config.Homesave, config.subjname, config.Structural, config.freesurferdir);

if exist(Savedir)<1
    
    fshome          = '/usr/local/freesurfer/';
    subdir          = fullfile(config.Homesave, config.subjname, config.Structural);
    mrifilename     = 'MR_acpc.nii'
    mrfile          = fullfile(subdir, mrifilename);
    system(['export FREESURFER_HOME=' fshome '; ' ...
        'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
        'mri_convert -c -oc 0 0 0 ' mrfile ' ' [subdir '/tmp.nii'] '; ' ...
        'recon-all -i ' [subdir '/tmp.nii'] ' -s ' 'freesurfer' ' -sd ' subdir ' -all'])
end

%% control freesurfer
%troubleshooting fo freeesurfer volumes.3d reconstruction
Pathlhpial      = 'freesurfer/surf/lh.pial';
Pathrhpial      = 'freesurfer/surf/rh.pial';

path2freesurferpial = fullfile(config.Homesave, config.subjname, config.Structural, Pathlhpial)
pial_lh = ft_read_headshape(path2freesurferpial);
pial_lh.coordsys = 'acpc';
ft_plot_mesh(pial_lh);
lighting gouraud;
camlight;

path2freesurferpial = fullfile(config.Homesave, config.subjname, config.Structural, Pathrhpial)
pial_rh = ft_read_headshape(path2freesurferpial);
pial_rh.coordsys = 'acpc';
ft_plot_mesh(pial_rh);
lighting gouraud;
camlight;

%% Ct scan

config.CTfile             = fullfile(config.Homesave, config.subjname, config.Structural, config.CTname);


savename        = 'MRI_CT_acpc.nii';
Savedir         = fullfile(config.Homesave, config.subjname, config.Structural, savename);

if exist(Savedir)>0 || config.repreproc == 1
    ct_acpc_f               = ft_read_mri(Savedir);
    fsmri_acpc.coordsys     = 'acpc';
else
    ct_acpc_f               = CTMakeSave(config)
end

config.ct_acpc_f            = ct_acpc_f;


%Import header from iEEG file
hdr = ft_read_header(config.fullEEG)

%Localize electrodes
savename        = 'elec_acpc_f.mat';
Savedir         = fullfile(config.Homesave, config.subjname, config.Structural, savename); 
if exist(Savedir)>0
    load(Savedir)
    
    ft_plot_ortho(config.fsmri_acpc.anatomy, 'transform', config.fsmri_acpc.transform, 'style', 'intersect');
    ft_plot_sens(elec_acpc_f, 'label', 'on', 'fontcolor', 'w');
else
    cfg         = [];
    cfg.channel = hdr.label;
    elec_acpc_f = ft_electrodeplacement(cfg, ct_acpc_f, config.fsmri_acpc);

    %step 19
    ft_plot_ortho(config.fsmri_acpc.anatomy, 'transform', config.fsmri_acpc.transform, 'style', 'intersect');
    ft_plot_sens(elec_acpc_f, 'label', 'on', 'fontcolor', 'w');

    %step 20

    save(Savedir, 'elec_acpc_f');
end



%% atlases
% 
% atlas = ft_read_atlas('/media/ludovico/DATA/iEEG_Ludo/Globus/TeamValeroCabre/Subject_5/Structural/freesurfer/mri/aparc.DKTatlas+aseg.mgz');
% atlas.coordsys = 'acpc'
% for i = 1:length(elec_acpc_fr.label)
%     cfg            = [];
%     cfg.roi        = elec_acpc_fr.chanpos(match_str(elec_acpc_fr.label, elec_acpc_fr.label(i)),:);
%     cfg.atlas      = atlas;
%     cfg.inputcoord = 'acpc';
%     cfg.output     = 'label';
%     labels = ft_volumelookup(cfg, atlas);
%     
%     [~, indx] = max(labels.count);
%     labels.name(indx)
%     labels(i, 1) = elec_acpc_fr.label(i);
%     labels(i, 2) = labels.name
% end

% Importing the different atlases in participants native space
elec_acpc_fr = elec_acpc_f;
DKTname                     = 'aparc.DKTatlas+aseg.mgz';
aparc_asegname              = 'aparc+aseg.mgz';
aparc2009name               = 'aparc.a2009s+aseg.mgz';
% Hippname                    = 'lh.hippoAmygLabels-T1-T2.v21.mgz'
fsfolder                    = 'freesurfer/mri/';
locationatlas               = fullfile(config.Homesave, config.subjname, config.Structural, fsfolder, DKTname);
atlas_nat(1).atlas          = ft_read_atlas(locationatlas);
atlas_nat(1).atlas.coordsys = 'acpc';
locationatlas               = fullfile(config.Homesave, config.subjname, config.Structural, fsfolder, aparc_asegname);
atlas_nat(2).atlas          = ft_read_atlas(locationatlas);
atlas_nat(2).atlas.coordsys = 'acpc';
locationatlas               = fullfile(config.Homesave, config.subjname, config.Structural, fsfolder, aparc2009name);
atlas_nat(3).atlas          = ft_read_atlas(locationatlas);
atlas_nat(3).atlas.coordsys = 'acpc';
% locationatlas               = fullfile(config.Homesave, config.subjname, config.Structural, fsfolder, Hippname);
% atlas_nat(4).atlas          = ft_read_atlas(locationatlas);
% atlas_nat(4).atlas.coordsys = 'acpc';
string_labels               = elec_acpc_fr.label;


for iatlas = 1 : length(atlas_nat)
  % here use the elect location and atlases in native space
  cfg      = [];
  cfg.roi    = elec_acpc_f.chanpos(match_str(elec_acpc_f.label,string_labels),:); %[17.7 ,0.0  ,-17.0]
  cfg.atlas   = atlas_nat(iatlas).atlas;
  cfg.inputcoord = 'acpc'; %'mni'; or 'tal' % tailarach or mni for volume based. In freesurfer atlases use acpc
  cfg.output   = 'multiple' %'single';
  labels = ft_volumelookup(cfg, atlas_nat(iatlas).atlas);
  % assigning anatomical labels to each electrode.
  anat_label = {};
  for ix = 1 : length(labels)
    [~, indx] = max(labels(ix).count);
    elec_acpc_f.atlas(iatlas).anat_label{ix} = char(labels(ix).name(indx));
  end     
end
% Use multiple atlases and save the data
atlases_tab = cell2table(horzcat(elec_acpc_f.label , ...
elec_acpc_f.atlas(1).anat_label',elec_acpc_f.atlas(2).anat_label',...
elec_acpc_f.atlas(3).anat_label'),'VariableNames',{'Chan' ,'DKT', 'Aparc', 'Aparc2009'})
% atlases_tab = cell2table(horzcat(elec_acpc_f.label,...
%     elec_acpc_f.atlas(1), anat_label',anat_label, anat_label'),'VariableNames',{'Chan' ,'DKT', 'Aparc', 'Aparc2009' })
elec_pos = array2table(elec_acpc_f.elecpos,'VariableNames',{'posx' 'posy' 'posz' } );
iEEG_correg = [atlases_tab, elec_pos];
writetable(iEEG_correg,'SEEG_correg_elec_native.xlsx','Sheet',1,'Range','A1')

savename        = 'iEEG_correg';
Savedir         = fullfile(config.Homesave, config.subjname, config.Structural, savename); 
save(Savedir, 'iEEG_correg')