function [] = ff_generateindices(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load files containing climate indices and compute them for each flight
%also deseasonalise, if requested
%also delinearise, if requested <-- TO DO
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/12/27 - 2021/12/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load flight data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load flights
FlightData = load('data/flight_data.mat');

%get flight dates
TimeScale = FlightData.Results.Date;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for each requested index:
%1. interpolate it to the date of the flight
%2. store it in an array named liked the index
%
%all indices will be normalised to a range of 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Indices = struct();


for iIndex=1:1:numel(Settings.Indices)
  
  switch Settings.Indices{iIndex}
    case 'QBO'
      QBO = load('indices/QBO.mat');
      a = interp1(QBO.Time,QBO.QBO,TimeScale);
      clear QBO
    case 'ENSO'
      ENSO = load('indices/nino34.mat');
      a = interp1(ENSO.Time,ENSO.Nino34,TimeScale);
      clear ENSO
    case 'HadCRUT'
      HadCRUT = rCDF('indices/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc');
      HadCRUT.MatlabTime = datenum(1850,1,HadCRUT.time);
      HadCRUT.NH = squeeze(nanmean(HadCRUT.tas_mean(:,HadCRUT.latitude > 0,:),[1,2]));
% %       HadCRUT = rCDF('indices/HadCRUT.4.6.0.0.median.nc');
% %       HadCRUT.MatlabTime = datenum(1850,1,HadCRUT.time);
% %       HadCRUT.NH = squeeze(nanmean(HadCRUT.temperature_anomaly(:,HadCRUT.latitude > 0,:),[1,2]));
      a = interp1(HadCRUT.MatlabTime,HadCRUT.NH,TimeScale);
      clear HadCRUT
    case 'NAM'
      NAM = load('indices/daily_nam.mat');
      a = interp1(NAM.Time,NAM.NAM,TimeScale);
      clear NAM
    case 'NAO'
      NAO = load('indices/nao.mat');
      a = interp1(NAO.Date,NAO.NAO,TimeScale);
      clear NAO      
    case 'TSI'
      TSI = load('indices/tsi.mat');
      TSI.TSI(TSI.TSI < 1358) = NaN; %very noisy - remove an extreme outlier
      TSI.TSI = smoothn(TSI.TSI,[31,1]);%very noisy - smooth to make it fair
      a = interp1(TSI.Time,TSI.TSI,TimeScale);
      clear TSI
    case 'Time'
      a = TimeScale;
    case 'SeaIce'
      SeaIce = readmatrix('indices/N_seaice_extent_daily_v3.0.csv');
      t = datenum(SeaIce(:,1),SeaIce(:,2),SeaIce(:,3));
      a = interp1(t,SeaIce(:,4),TimeScale);
      clear SeaIce t
    otherwise
      disp(['Index ',Settings.Indices{iIndex},' not specified; stopping'])
      stop
  end


  %store the original range, in case we need to refer to it.
  Indices.OriginalRange.(Settings.Indices{iIndex}) = minmax(a(:));
 
  %normalise
  a = a./range(a);
  a = a-min(a);
  a= (a.*2)-1;

  
  %store
  Indices.(Settings.Indices{iIndex}) = a;
  
  clear a
end; clear iIndex


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% we may want to regress data with lags. 
%if so, pregenerate the lagged data here
%this piece of code is very inefficient with file loading,
%but the files are small and it saves a full rewrite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.Reg.Lag == 1;
  
  Indices.Lagged = struct();
  LagScale = Settings.Reg.Steps;
  Indices.Lagged.LaggedScale = LagScale;
  
  for iIndex=1:1:numel(Settings.Indices)   
    Store = NaN(numel(LagScale),numel(TimeScale));
    
    for iLag=1:1:numel(LagScale)
      
      switch Settings.Indices{iIndex}
        case 'QBO'
          QBO = load('indices/QBO.mat');
          QBO.Time = floor(QBO.Time); %shift from noon to midnight to make the logic easier - on a 91-day smoothing this is very minor...
          a = interp1(QBO.Time+LagScale(iLag),QBO.QBO,TimeScale);
          clear QBO
        case 'ENSO'
          ENSO = load('indices/nino34.mat');
          a = interp1(ENSO.Time+LagScale(iLag),ENSO.Nino34,TimeScale);
          clear ENSO
        case 'HadCRUT'
          HadCRUT = rCDF('indices/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc');
          HadCRUT.MatlabTime = datenum(1850,1,HadCRUT.time);
          HadCRUT.NH = squeeze(nanmean(HadCRUT.tas_mean(:,HadCRUT.latitude > 0,:),[1,2]));
          % %       HadCRUT = rCDF('indices/HadCRUT.4.6.0.0.median.nc');
          % %       HadCRUT.MatlabTime = datenum(1850,1,HadCRUT.time);
          % %       HadCRUT.NH = squeeze(nanmean(HadCRUT.temperature_anomaly(:,HadCRUT.latitude > 0,:),[1,2]));
          a = interp1(HadCRUT.MatlabTime,HadCRUT.NH,TimeScale);
          clear HadCRUT
        case 'NAM'
          NAM = load('indices/daily_nam.mat');
          a = interp1(NAM.Time+LagScale(iLag),NAM.NAM,TimeScale);
          clear NAM
        case 'NAO'
          NAO = load('indices/nao.mat');
          a = interp1(NAO.Date+LagScale(iLag),NAO.NAO,TimeScale);
          clear NAO
        case 'TSI'
          TSI = load('indices/tsi.mat');
          TSI.TSI(TSI.TSI < 1358) = NaN; %very noisy - remove an extreme outlier
          TSI.TSI = smoothn(TSI.TSI,[15,1]);%very noisy - smooth to make it fair
          a = interp1(TSI.Time,TSI.TSI,TimeScale);
          clear TSI
        case 'Time'
          a = TimeScale+LagScale(iLag);
        otherwise
          disp(['Index ',Settings.Indices{iIndex},' not specified; stopping'])
          stop
      end
      
      %normalise
      a = a./range(a);
      a = a-min(a);
      a= (a.*2)-1;
  
      %store
      Store(iLag,:) = a;
      
    end

    Indices.Lagged.(Settings.Indices{iIndex}) = Store;
    clear Store
  end; clear iIndex
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% deseasonalise data and indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.DS == 1;

  %find day-of-year for each point
  dd = floor(date2doy(TimeScale));

  for iIndex=0:1:numel(Settings.Indices); %the 0th entry is to do the raw data with the same logic

    %skip non-annual-cycle indices
    if iIndex ~= 0;
      if strcmp(Settings.Indices{iIndex},   'ENSO'); continue; end
      if strcmp(Settings.Indices{iIndex},    'QBO'); continue; end
      if strcmp(Settings.Indices{iIndex},    'TSI'); continue; end
      if strcmp(Settings.Indices{iIndex},'HadCRUT'); continue; end
      if strcmp(Settings.Indices{iIndex},   'Time'); continue; end
    end
  
    %get index
    if iIndex == 0; a = FlightData.Results.tRel;
    else            a = Indices.(Settings.Indices{iIndex});
    end

    %find day-of-yearly averages
    DV = NaN(366,1);
    for iDay=1:1:365;
      DV(iDay) = nanmean(a(dd == iDay));
    end; clear iDay

    %smooth by a week
    DV = smoothn([DV;DV;DV],[7,1]); DV = DV(367:366*2);

    %and remove from the raw data
    for iDay=1:1:365;
      a(dd == iDay) = a(dd == iDay)-DV(iDay);
    end; clear iDay DV

    %finally, renormalise

    if iIndex == 0;
    else; 
      b = Indices.OriginalRange.(Settings.Indices{iIndex})
      b = b ./(2./range(a))
    end

    a = a./range(a);
    a = a-min(a);
    a= (a.*2)-1;



    %store index
    if iIndex == 0; FlightData.Results.tRel           = a;
    else            Indices.(Settings.Indices{iIndex})= a;
    end    


  end; clear iIndex


  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('data/indices.mat','Indices')

Airports = FlightData.Airports; Results = FlightData.Results; Settings = FlightData.Settings;
save('data/flight_data.mat','Airports','Results','Settings')