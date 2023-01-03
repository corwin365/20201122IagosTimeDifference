clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Impact of climate-dynamical processes on flight times
%
%script to run analyses
%
%Corwin Wright, c.wright@bath.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%what settings are we using?
SettingsID = 'test567';

%load them
Settings = load(['data/',SettingsID,'.mat']);
clear SettingsID  %it's contained in the file we loaded

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% import and generate raw data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % %generate airport geolocation
% % % generate_geolocation(Settings);
% % % 
% % % %generate IAGOS data
% % % if sum(ismember(Settings.Choices.DataSets,1)); import_iagos_data(Settings); end
% % % 
% % % %generate data from Ed?
% % % if sum(ismember(Settings.Choices.DataSets,2)); import_ed_data(   Settings); end
% % % 
% % % %generate merged dataset
% % % generate_merged_data(Settings);

%generate climate indices for the data
generate_climate_indices_v2(Settings);
