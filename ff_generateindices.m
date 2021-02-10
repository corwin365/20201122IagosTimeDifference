function [] = ff_generateindices(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load files containing climate indices and compute them for each flight
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/12/27
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
      QBO.Time = floor(QBO.Time); %shift from noon to midnight to make the logic easier - on a 91-day smoothing this is very minor...
      a = interp1(QBO.Time,QBO.QBO,TimeScale);
      clear QBO
    case 'ENSO'
      ENSO = load('indices/nino34.mat');
      a = interp1(ENSO.Time,ENSO.Nino34,TimeScale);
      clear ENSO
    case 'HadCRUT'
      HadCRUT = rCDF('indices/HadCRUT.4.6.0.0.median.nc');
      HadCRUT.MatlabTime = datenum(1850,1,HadCRUT.time);
      HadCRUT.NH = squeeze(nanmean(HadCRUT.temperature_anomaly(:,HadCRUT.latitude > 0,:),[1,2]));
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
      TSI.TSI = smoothn(TSI.TSI,[15,1]);%very noisy - smooth to make it fair
      a = interp1(TSI.Time,TSI.TSI,TimeScale);
      clear TSI
    case 'Time'
      a = TimeScale;
    otherwise
      disp(['Index ',Settings.Indices{iIndex},' not specified; stopping'])
      stop
  end
  
 
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
          HadCRUT = rCDF('indices/HadCRUT.4.6.0.0.median.nc');
          HadCRUT.MatlabTime = datenum(1850,1,HadCRUT.time);
          HadCRUT.NH = squeeze(nanmean(HadCRUT.temperature_anomaly(:,HadCRUT.latitude > 0,:),[1,2]));
          a = interp1(HadCRUT.MatlabTime+LagScale(iLag),HadCRUT.NH,TimeScale);
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
%% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('data/indices.mat','Indices')