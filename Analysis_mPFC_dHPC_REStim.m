%% Analysis

% Script to analyze all animals
% - The code relies on the following functions : -> Open_Files.m
%                                                -> Extracting_LFPs_and_events_from_all.m
%                                                -> ... other desired functions

% by Flavio Mourao.
% email: mourao.fg@illinois.edu
% Maren Lab -  Beckman Institute for Advanced Science and Technology
% University of Illinois Urbana-Champaign

% Started in:  12/2023
% Last update: 07/2025

%%

tic
set(groot,'DefaultFigureColormap',jet)

[files] = Open_files();

for ms = 1:length(files.FilesLoaded{1, 1}) % Loop over files path

    % Settings
    id          = files.id(ms).name;
    FilesLoaded = {files.FilesLoaded{1,1}(ms).name};
    Path        = files.FilesLoaded{1,1}(ms).folder;

    %% Open and extract files

    %  Load function - raw files
    % [data,parameters] = Extracting_LFPs_and_events(id,FilesLoaded,Path);
    
    % Load .*mat files
    name = strcat(Path,'/',id);
    fprintf('\n loading data... \n');
    %load(name,'data','parameters'); % pre processed files
    load(name);                      % analised data

    %% Optogenetic Stim Conterbalanced Order

    % Habituation -> to set plots, pre-processing and future stats

    % First Cohort
    % OptoSTIM{1} = [8 4]; % 8Hz -> 4Hz - R02 - chR2
    % OptoSTIM{2} = [8 4]; % 8Hz -> 4Hz - R07 - chR2

    % Second Cohort
    % mCherries
    % OptoSTIM{1} = [4 8]; % 4Hz -> 8Hz - R01 - mCherry
    % OptoSTIM{2} = [4 8]; % 4Hz -> 8Hz - R07 - mCherry
    % OptoSTIM{3} = [8 4]; % 8Hz -> 4Hz - R08 - mCherry
    % OptoSTIM{4} = [8 4]; % 8Hz -> 4Hz - R09 - mCherry

    % ChR2s
    % OptoSTIM{1} = [8 4]; % 8Hz -> 4Hz - R10 - chR2
    % OptoSTIM{2} = [4 8]; % 4Hz -> 8Hz - R11 - chR2


    %% Scripts pre-processing

    %    Pre_processing_mPCF_dHPC_Habituation_noCS_4Hz_8Hz_Stim_v2;
    %    Pre_processing_mPCF_dHPC_Extinction_CS_8Hz_stim
    %    Pre_processing_mPCF_dHPC_Retrieval_CS_NoStim
    %    Pre_processing_mPCF_dHPC_Renewal_noCS_8Hz_Stim

    %    CorCov_mPFC_dHPC       % Correlation and Covariance Matrices betwwen channels
    %    Main_Plots_mPCF_dHPC_  % plot to check raw data

    %% Habituation
    %  4Hz / 8Hz opto stimulation

    % Load. Choose data from data.lfp according pre_processing.m and main plots script
    % name = strcat(Path,'/',id);
    % fprintf('\n loading data... \n');
    % 
    % load(name);

    % 1) Baseline timeepochs without noise

    % B_clean{1} = dsearchn(data.timev_decimated',[10 170]');   % R02 - chR2
    % B_clean{2} = dsearchn(data.timev_decimated',[10 170]');   % R07 - chR2
    % B_clean{3} = dsearchn(data.timev_decimated',[10 170]');   % R10 - chR2
    % B_clean{4} = dsearchn(data.timev_decimated',[10 170]');   % R11 - chR2

    % B_clean{1} = dsearchn(data.timev_decimated',[10 170]');   % R01 - mCherry
    % B_clean{2} = dsearchn(data.timev_decimated',[10 170]');   % R07 - mCherry
    % B_clean{3} = dsearchn(data.timev_decimated',[10 170]');   % R08 - mCherry
    % B_clean{4} = dsearchn(data.timev_decimated',[10 170]');   % R09 - mCherry

    % 2) Selecting  trials

    % Stim ON
    % CSIT{1} = [1:5];   % R02 - R01
    % CSIT{2} = [1:5];   % R07 - R07
    % CSIT{3} = [1:5];   % R10 - R08
    % CSIT{4} = [1:5];   % R11 - R09
    
    % Stim OFF
    % CSIT_1{1} = [1:5];   % R02 - R01
    % CSIT_1{2} = [1:5];   % R07 - R07
    % CSIT_1{3} = [1:5];   % R10 - R08
    % CSIT_1{4} = [1:5];   % R11 - R09


    % 3) Scripts

    % p_welch
    % p_welch_final_plots
    % wavelets_spec
    % PLV_Hilbert_
    % PLV_Hilbert_final_plots
    % PLV_spec_timecourse

    %% Extinction
    %  8Hz opto stimulation + CS-Tones

    % Load. Choose data from data.lfp according pre_processing.m and main plots script
    % name = strcat(Path,'/',id);
    % fprintf('\n loading data... \n');
    %
    % load(name);

    % 1) Baseline timeepochs without noise

    % B_clean{1} = dsearchn(data.timev_decimated',[10 170]');   % R02 - chR2
    % B_clean{2} = dsearchn(data.timev_decimated',[10 170]');   % R07 - chR2
    % B_clean{3} = dsearchn(data.timev_decimated',[10 170]');   % R10 - chR2
    % B_clean{4} = dsearchn(data.timev_decimated',[10 170]');   % R11 - chR2

    % B_clean{1} = dsearchn(data.timev_decimated',[10 170]');   % R01 - mCherry
    % B_clean{2} = dsearchn(data.timev_decimated',[10 170]');   % R07 - mCherry
    % B_clean{3} = dsearchn(data.timev_decimated',[10 170]');   % R08 - mCherry
    % B_clean{4} = dsearchn(data.timev_decimated',[10 170]');   % R09 - mCherry
    

    % 2) Selecting  trials

    % Stim ON
    % CSIT{1} = [1:45];   % R02 - R01
    % CSIT{2} = [1:45];   % R07 - R07
    % CSIT{3} = [1:45];   % R10 - R08
    % CSIT{4} = [1:45];   % R11 - R09
    
    % Stim OFF
    % CSIT_1{1} = [1:45];   % R02 - R01
    % CSIT_1{2} = [1:45];   % R07 - R07
    % CSIT_1{3} = [1:45];   % R10 - R08
    % CSIT_1{4} = [1:45];   % R11 - R09


    % 3) Scripts

    % p_welch
    % p_welch_final_plots
    % wavelets_spec
    % PLV_Hilbert_
    % PLV_Hilbert_final_plots
    % PLV_spec_timecourse

    %% Clear
    if ms < length(files.FilesLoaded{1, 1})
        %  clear('FilesLoaded','Path','data','parameters','newStr1','path', 'name' )
        clear('FilesLoaded','Path','newStr1','path', 'name' )

    else
        clear('id','FilesLoaded','Path','ms')
    end

end
toc

fprintf('\n Done. \n');

%% last update 07/17/2025
%  listening: This will destroy you - quiet