function generate_merged_data(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%merge all datasets used
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/01/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Merging data')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load IAGOS?
if sum(ismember(Settings.Choices.DataSets,1)) == 1
  IAGOS = load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_iagos.mat'],'FlightData');
  Data.set1 = IAGOS.FlightData;
end

%Ed's data?
if sum(ismember(Settings.Choices.DataSets,2)) == 1
  EdG = load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_edG.mat'],'FlightData');
  Data.set2 = EdG.FlightData;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% merge data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDataSet=1:1:numel(Settings.Choices.DataSets)
  ThisDataSet = Data.(['set',num2str(Settings.Choices.DataSets(iDataSet))]);

  if ~exist('FlightData'); FlightData = ThisDataSet;
  else                     FlightData = [FlightData;ThisDataSet];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


save([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_merged.mat'],'FlightData')

disp('--------------------------')
disp('Datasets merged')
disp('--------------------------')

end

