function [] = bb_merge_data(Paths,WhichData,Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fmerge original flavour IAGOS data with new Ed-supplied data (or not)
%
%Corwin Wright, c.wright@bath.ac.uk, 2022/07/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the datasets
if WhichData ~=2; 
  IAGOS = load([Paths.StoreDir,'/iagos_flight_data_',Paths.SourceIdentifier,'.mat']); 
end
if WhichData ~= 1
  Extra = load([Paths.StoreDir,'/extra_flight_data_',Paths.SourceIdentifier,'.mat']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create empty output fields. we'll drop the first element later
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fields   = {'Dep','Arr','PlaneID','InstID','t','Date'};
Flights = struct();
for iF=1:1:numel(Fields);
  Flights.(Fields{iF}) = cell(1);
end; clear  iF
Flights.t = Flights.t{1};
Flights.Date = Flights.Date{1};

PathFields = {'Lat','Lon','U','V','T'};
Flights.Paths = struct();
for iF=1:1:numel(PathFields);
  Flights.Paths.(PathFields{iF}) = cell(1);
end; clear  iF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add IAGOS data?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if WhichData == 0 | WhichData == 1; 

  %copy main data over
  for iF=1:1:numel(Fields)
    Flights.(Fields{iF}) = cat(1,Flights.(Fields{iF}),IAGOS.Flights.(Fields{iF}));
  end

  %retain individual flight traces
  for iF=1:1:numel(PathFields)
    Flights.Paths.(PathFields{iF}) = cat(2,Flights.Paths.(PathFields{iF}),IAGOS.Flights.Paths.(PathFields{iF}));
  end

  Settings.IAGOS = IAGOS.Settings;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add extra data?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if WhichData == 0 | WhichData == 2; 

  %transpose, to match other dataset
  %do it
  for iF=1:1:numel(Fields)
    Flights.(Fields{iF}) = cat(1,Flights.(Fields{iF}),Extra.Flights.(Fields{iF}));
  end

  %retain individual flight traces. That are actually blank...
  for iF=1:1:numel(PathFields)
    Flights.Paths.(PathFields{iF}) = cat(2,Flights.Paths.(PathFields{iF}),Extra.Flights.Paths.(PathFields{iF}));
  end


  Settings.Extra = Extra.Settings;  
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%drop empty first field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iF=1:1:numel(Fields)
  if strcmp(Fields{iF},'t') | strcmp(Fields{iF},'Date') ; continue; end
  a = Flights.(Fields{iF});
  a = a(2:end);
  Flights.(Fields{iF}) = a;
end
for iF=1:1:numel(PathFields)
  a = Flights.Paths.(PathFields{iF});
  a = a(2:end);
  Flights.Paths.(PathFields{iF}) = a;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([Paths.StoreDir,'/flight_data_',Paths.SourceIdentifier,'.mat'],'Flights','Settings')

%and return
return
