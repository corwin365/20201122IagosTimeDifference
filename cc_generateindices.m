function [] = cc_generateindices(Paths,Indices,Settings,IndexRange)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load files containing climate indices and compute them for each flight
%Corwin Wright, c.wright@bath.ac.uk, 2021/12/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load flight data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load flights
FlightData = load([Paths.StoreDir,'/flight_data_',Paths.SourceIdentifier,'.mat'],'Flights');
FlightData = FlightData.Flights;

%store (1) flight times and (2) create a linear daily scale from the first to the last date. 
TimeScale.Flights = FlightData.Date';
TimeScale.Daily = floor(min(TimeScale.Flights)):1:ceil(max(TimeScale.Flights));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for each requested index:
%1. interpolate it to the time of the flight
%2. interpolate it to the daily time scale
%2. store it in an array named liked the index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IndexStore = struct();
Root = Paths.Indices;


for iIndex=1:1:numel(Indices)
  
  switch Indices{iIndex}
    case 'QBO'
      QBO = load([Root,'/QBO.mat']);
      a = interp1(QBO.Time,QBO.QBO,TimeScale.Flights);
      b = interp1(QBO.Time,QBO.QBO,TimeScale.Daily); 
      clear QBO
    case 'ENSO'
      ENSO = load([Root,'/nino34.mat']);
      a = interp1(ENSO.Time,ENSO.Nino34,TimeScale.Flights);
      b = interp1(ENSO.Time,ENSO.Nino34,TimeScale.Daily);      
      clear ENSO
    case 'HadCRUT'
      HadCRUT = rCDF([Root,'/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc']);
      HadCRUT.MatlabTime = datenum(1850,1,HadCRUT.time);
      HadCRUT.NH = squeeze(nanmean(HadCRUT.tas_mean(:,HadCRUT.latitude > 0,:),[1,2]));
      a = interp1(HadCRUT.MatlabTime,HadCRUT.NH,TimeScale.Flights);
      b = interp1(HadCRUT.MatlabTime,HadCRUT.NH,TimeScale.Daily); 
      clear HadCRUT
    case 'NAM'
      NAM = load([Root,'/daily_nam.mat']);
      a = interp1(NAM.Time,NAM.NAM,TimeScale.Flights);
      b = interp1(NAM.Time,NAM.NAM,TimeScale.Daily); 
      clear NAM
    case 'NAO'
      NAO = load([Root,'/nao.mat']);
      a = interp1(NAO.Date,NAO.NAO,TimeScale.Flights);
      b = interp1(NAO.Date,NAO.NAO,TimeScale.Daily); 
      clear NAO      
    case 'TSI'
      TSI = load([Root,'/tsi.mat']);
      TSI.TSI(TSI.TSI < 1358) = NaN; %very noisy - remove an extreme outlier
      TSI.TSI = smoothn(TSI.TSI,[31,1]);%very noisy - smooth to make it fair
      a = interp1(TSI.Time,TSI.TSI,TimeScale.Flights);
      b = interp1(TSI.Time,TSI.TSI,TimeScale.Daily); 
      clear TSI
    case 'Time'
      a = TimeScale.Flights;
      b = TimeScale.Daily; 
    case 'SeaIce'
      SeaIce = readmatrix([Root,'/N_seaice_extent_daily_v3.0.csv']);
      t = datenum(SeaIce(:,1),SeaIce(:,2),SeaIce(:,3));
      a = interp1(t,SeaIce(:,4),TimeScale.Flights);
      b = interp1(t,SeaIce(:,4),TimeScale.Daily); 
      clear SeaIce t
    otherwise
      disp(['Index ',Settings.Indices{iIndex},' not specified; stopping'])
      stop
  end
 
  %store
  IndexStore.Flights.(Indices{iIndex}) = a;
  IndexStore.Daily.(  Indices{iIndex}) = b;  

  disp(['Loaded ',Indices{iIndex}])

  clear a b
end; clear iIndex Root 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% produce a deseasonalised form of each index
%this will not be valid if we're using less than a few years of data
%but that's a choice that can be made later
%
%also, several will not be physically meaningful. Again, deal with this downstream.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find day-of-year for each point
dd   = floor(date2doy(TimeScale.Daily));
dd_f = date2doy(TimeScale.Flights);


for iIndex=1:1:numel(Indices);


  %get index
  a = IndexStore.Daily.(Indices{iIndex});

  %find day-of-yearly averages
  DV = NaN(366,1);
  for iDay=1:1:366;
    DV(iDay) = nanmean(a(dd == iDay));
  end; clear iDay

  %smooth a bit in time
  DV = smoothn([DV;DV;DV],[Settings.DSSmooth,1]); 


  %...daily...
  for iDay=1:1:366;
    a(dd == iDay) = a(dd == iDay)-DV(iDay+366); %the +366 is so we use the 'middle' (i.e. correctly smoothed at ends) time series
  end; clear iDay 

  %and flight-based data
  b = IndexStore.Flights.(Indices{iIndex}) - interp1(1:1:numel(DV),DV,dd_f+366,'linear','extrap');

  %and store
  IndexStore.DS.Daily.(  Indices{iIndex})= a;
  IndexStore.DS.Flights.(Indices{iIndex})= b;


end; clear iIndex a b DV dd dd_f



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tidy variable space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%flight-based
FlightIndices.Time      = TimeScale.Flights;
FlightIndices.Raw       = IndexStore.Flights;
FlightIndices.DS        = IndexStore.DS.Flights;

%date-based
DateIndices.Time       = TimeScale.Daily;
DateIndices.Raw        = IndexStore.Daily;
DateIndices.DS         = IndexStore.DS.Daily;

clearvars -except Paths IndexRange DateIndices FlightIndices Indices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% normalise the data (based on the daily data as these are less biased)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




for iForm=1:2; %1 is raw, 2 is deseasonalised

  %get the data

  if iForm==1;
    Daily = DateIndices.Raw;
    Flights = FlightIndices.Raw;
  elseif iForm ==2;
    Daily = DateIndices.DS;
    Flights = FlightIndices.DS;
  end
  RangeStore = struct();

  for iIndex =1:1:numel(Indices)

    %compute the necessary stats from the daily data
    RangeStore.(Indices{iIndex}) = prctile(Daily.(Indices{iIndex}),IndexRange);

    %apply the stats to normalise the data
    a = Daily.(Indices{iIndex});
    a = a./range(RangeStore.(Indices{iIndex}));
    a = a-min(a);
    a= (a.*2)-1;
    Daily.(Indices{iIndex}) = a;

    a = Flights.(Indices{iIndex});
    a = a./range(RangeStore.(Indices{iIndex}));
    a = a-min(a);
    a= (a.*2)-1;
    Flights.(Indices{iIndex}) = a;


  end; clear iIndex

 
  %and store
  if iForm==1;
    DateIndices.Raw   = Daily;    DateIndices.Raw.Ranges   = RangeStore;
    FlightIndices.Raw = Flights;  FlightIndices.Raw.Ranges = RangeStore;
  elseif iForm ==2;
    DateIndices.DS   = Daily;    DateIndices.DS.Ranges   = RangeStore;
    FlightIndices.DS = Flights;  FlightIndices.DS.Ranges = RangeStore;
  end  
  
end; clear iForm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



save([Paths.StoreDir,'/indices_',Paths.SourceIdentifier,'.mat'],'FlightIndices','DateIndices')


%and return
return
