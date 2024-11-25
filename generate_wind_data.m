function [] = generate_wind_data(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate ERA5 wind data used in subsequent analyses
%
%wind data will be produced in two forms:
%  1. pointwise U and V along each flight track
%  2. maps of U and V every three hours
%
%tropopause height and stratopause pressure will also be extracted, as this 
%is computationally efficient to get here and might be needed
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/09/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Generating wind data')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get flight data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_merged.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get wind data based on this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%generate a list of dates that will include every flight in the dataset
FullDateRange = (floor(min(FlightData.Date))-1):1:(ceil(max(FlightData.Date))+1);

%create a store for the flight traces
FlightStore = {};

%loop over dates
textprogressbar('Importing contextual ERA5 wind ')
for iDay=1:1:numel(FullDateRange) %note that by definition this starts one day before the first flight


  %% generate list of points needed for MAP
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %first, generate a list of points to get a map for this day
  %get the data on a six-hourly mesh
  t = FullDateRange(iDay)+(0:6:18)./24;
  [xi,yi,ti] = meshgrid(Settings.Choices.WindMap.LonScale,Settings.Choices.WindMap.LatScale,t);

  %put these into a line, but retain the size so we can stitch it back together later
  sz = size(xi);
  List.Lon  = xi(:);
  List.Lat  = yi(:);
  List.Time = ti(:);
  List.Fid  = zeros(size(List.Time)); %'file ID' - 0 for map points
  List.P    = ones( size(List.Time)).*250; %pressure level for the maps
  clear t xi yi ti

  %% generate list of points needed for FLIGHTS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %identify all flights starting on this day
  FlightsThisDay = find(floor(FlightData.Date) == FullDateRange(iDay));

  if numel(FlightsThisDay) > 0;

    %add all flights starting on this day into memory
    for iFile=1:1:numel(FlightsThisDay)

      if FlightData.DataSource(FlightsThisDay(iFile)) ~= "IAGOS"; continue; end

      %I generated the flight list on my laptop but want to run on 0184,
      %so this line just switches filepaths if needed
      FilePath = FlightData.FilePath(FlightsThisDay(iFile));
      FilePath = strrep(FilePath,'D:\Data\',LocalDataDir);
      FilePath = strrep(FilePath,'\','/');

      %load data and store
      Data =  rCDF(FilePath);
      List.Lat  = [List.Lat; Data.lat];
      List.Lon  = [List.Lon; Data.lon];
      List.Time = [List.Time;FullDateRange(iDay)+Data.UTC_time/86400];
      List.Fid  = [List.Fid;ones(size(Data.lon)).*FlightData.FlightIndex(FlightsThisDay(iFile))];
      List.P    = [List.P;Data.air_press_AC./100]; % Use values from the aircraft barometer for pressure level

    end
  else
    continue; 
  end
  clear FlightsThisDay iFile FilePath Data

  %% get winds at every point. 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
  %drop bad points for the wind interpolant, but put back later
  Bad = find(List.Lon < -180 | List.Lon > 180  ...
           | List.Lat <  -90 | List.Lat >  90  ...
           | List.Time < FullDateRange(iDay)-1 ...
           | List.Time > FullDateRange(iDay)+2 ...
           | List.P    > 1200                  ...
           | isnan(List.Lon) | isnan(List.Lat) | isnan(List.Time) | isnan(List.P));
  Good = 1:1:numel(List.Lon); Good(Bad) = [];
  List = reduce_struct(List,Good,{},0);

  %get the winds
  Wind = get_context(List.Lon,List.Lat,'TimePoints',List.Time,'Pressure',List.P,'Wind',true,'PressureAsPoints',true,'Pauses',true);

  %if we removed any points, put them back as gaps in the data
  if numel(Bad) > 0
    NewList = spawn_uniform_struct({'Lon','Lat','Time','Fid'},[numel(Good)+numel(Bad),1]);
    NewWind = spawn_uniform_struct({'U','V','TP','SP'},       [numel(Good)+numel(Bad),1]);
    NewList.Lon( Good) = List.Lon;
    NewList.Lat( Good) = List.Lat;
    NewList.Time(Good) = List.Time;
    NewList.Fid        = List.Fid; %this is used for indexing below so must stay the same as the original
    NewWind.U(   Good) = Wind.U;
    NewWind.V(   Good) = Wind.V;
    NewWind.TP(  Good) = Wind.Tropopause;
    NewWind.SP(  Good) = Wind.Stratopause;    

    List = NewList; clear NewList
    Wind = NewWind; clear NewWind
  end

  % store the map of wind
  if ~exist('WindStore')
    WindStore.U   = NaN([sz,numel(FullDateRange)]);
    WindStore.V   = NaN([sz,numel(FullDateRange)]);
    WindStore.t   = FullDateRange;
    WindStore.h   = 0:6:18;
    WindStore.Lat = Settings.Choices.WindMap.LatScale;
    WindStore.Lon = Settings.Choices.WindMap.LonScale;
  end
  WindStore.U( :,:,:,iDay) = reshape(Wind.U( List.Fid == 0),sz);
  WindStore.V( :,:,:,iDay) = reshape(Wind.V( List.Fid == 0),sz);
  WindStore.TP(:,:,:,iDay) = reshape(Wind.Tropopause( List.Fid == 0),sz);
  WindStore.Sp(:,:,:,iDay) = reshape(Wind.Stratopause(List.Fid == 0),sz);  

  %pull out flight traces and store
  if any(List.Fid ~= 0);
    Fids = unique(List.Fid(List.Fid ~= 0));
    for iFlight=1:1:numel(Fids)
      idx = find(List.Fid == Fids(iFlight));
      UStore{ Fids(iFlight)} = Wind.U( idx);
      VStore{ Fids(iFlight)} = Wind.V( idx);
      TPStore{Fids(iFlight)} = Wind.Tropopause(idx);
      SPStore{Fids(iFlight)} = Wind.Stratopause(idx);
    end; clear iFlight idx
  end

  clear Wind List Good Bad Fids sz
catch; end





  textprogressbar(iDay./numel(FullDateRange).*100)
end; clear iDay FullDateRange
textprogressbar('!')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([Settings.Paths.DataDir,'/',Settings.ID,'_flightwinds.mat'],'UStore','VStore','TPStore','SPStore')
save([Settings.Paths.DataDir,'/',Settings.ID,'_windmaps.mat'],'WindStore')




disp('--------------------------')
disp('Wind data generated')
disp('--------------------------')
