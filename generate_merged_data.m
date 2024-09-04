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
if sum(ismember(Settings.Choices.DataSets,1)) == 1 | sum(ismember(Settings.Choices.DataSets,3)) == 1
  IAGOS = load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_iagos.mat'],'FlightData');
  Data.set1 = IAGOS.FlightData;
  clear IAGOS
end

%Ed's data?
if sum(ismember(Settings.Choices.DataSets,2)) == 1
  EdG = load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_edG.mat'],'FlightData');
  Data.set2 = EdG.FlightData;
  clear EdG
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if we're using the 'eddays' option,
% subset the IAGOS data to just these dates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sum(ismember(Settings.Choices.DataSets,3)) == 1

  %this is a clunky approach that is computationally slow, but the dataset is small
  DaysToChoose = load('data/eddays.mat','eddays');
  DaysToChoose = DaysToChoose.eddays;
  
  DaysUsed = floor(Data.set1.Date);
  idx = [];
  for iDay=1:1:numel(DaysUsed)
    jdx = sum(ismember(DaysToChoose,DaysUsed(iDay)));
    if jdx > 0; idx = [idx,iDay]; end
  end

  Data.set1 = Data.set1(idx,:);

  %override dataset list
  Settings.Choices.DataSets(Settings.Choices.DataSets == 3) = 1;

  clear DaysToChoose DaysUsed iDay jdx idx
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


%add an index column
FlightData(:,size(FlightData,2)+1) = array2table(transpose(1:1:size(FlightData,1)));
FlightData = renamevars(FlightData,'Var9','FlightIndex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% deduplicate data 
% (this happens in several situations, eg if IAGOS and subsetted-IAGOS 
% are combined as options or if we're using both old and new format data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,idx] = unique(FlightData.Date); %assumes no two flights start at *exactly* the same time
FlightData = FlightData(idx,:);
clear idx


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tidy data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%a small number of flights have the same DEP and ARR airport. Remove these.
Good = find(FlightData.Dep ~= FlightData.Arr);
FlightData = FlightData(Good,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




save([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_merged.mat'],'FlightData')

disp('--------------------------')
disp('Datasets merged')
disp('--------------------------')

end

