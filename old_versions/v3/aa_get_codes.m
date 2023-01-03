function [] = aa_get_codes(Paths)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get a list of aircraft codes present in the IAGOS dataset
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/11/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get all flight departure and arrival airports
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%find files
Flights = wildcardsearch(Paths.AeolusData ,'*.nc');

%find starts and ends
Codes  = {};
Coords = {};

textprogressbar('Identifying unique airport codes in IAGOS dataset ')
for iFile=1:1:numel(Flights)
  try
  %get data
  Data = rCDF(Flights{iFile});
  
  %get metadata
  MetaData = Data.MetaData.Attributes.Global;

  %hence, identify start and end location
  NFields = numel(MetaData);
  for iField=1:1:NFields
    
    if strcmp(MetaData(iField).Name,'departure_airport')
      Dep = MetaData(iField).Value;
    elseif strcmp(MetaData(iField).Name,'arrival_airport')
      Arr = MetaData(iField).Value;
% %     elseif strcmp(MetaData(iField).Name,'departure_coord')
% %       Coords{end+1} = MetaData(iField).Value;      
% %     elseif strcmp(MetaData(iField).Name,'arrival_coord')
% %       Coords{end+1} = MetaData(iField).Value;  
    end

  end

  
  %trim down
  Dep = Dep(1:3);
  Arr = Arr(1:3);
  
  %and store
  Codes{end+1} = Dep;
  Codes{end+1} = Arr;

  
  if mod(iFile,100); textprogressbar(iFile./numel(Flights).*100); end
  catch;end
end
textprogressbar(100);textprogressbar('!')


%reduce to unique set
[Codes,ia] = unique(Codes);
% Coords = Coords(ia);

clearvars -except Codes Coords Paths
disp('Unique airport codes in IAGOS data identified')


%save
save([Paths.StoreDir,'/airport_codes_iagos.mat'],'Codes')

%and return
return
