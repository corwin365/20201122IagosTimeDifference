function process_climate_indices_v2(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load climate indices and interpolate them to flight times
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/01/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Generating climate indices')
disp('+++++++++++++++++++++++++++')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load flight data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat'])

%generate a continuous timescale at daily cadence
TimeScale  = min(Settings.Choices.TimeRange):1:max(Settings.Choices.TimeRange);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for each requested index:
%1. interpolate it to the daily time scale
%2. store it in an array named liked the index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IndexStore = struct();
Root = Settings.Paths.Indices;


for iIndex=1:1:numel(Settings.Indices.List)
  
  switch Settings.Indices.List{iIndex}
    case 'QBO'
      QBO = load([Root,'/QBO.mat']);
      b = interp1(QBO.Time,QBO.QBO,TimeScale); 
      clear QBO
    case 'ENSO'
      ENSO = load([Root,'/nino34.mat']);
      b = interp1(ENSO.Time,ENSO.Nino34,TimeScale);      
      clear ENSO
    case 'Fuel'
      Fuel = load([Root,'/jet_fuel_price.mat']);
      b = interp1(Fuel.Time,Fuel.Price,TimeScale);      
      clear Fuel
    case 'HadCRUT'
      HadCRUT = rCDF([Root,'/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc']);
      HadCRUT.MatlabTime = datenum(1850,1,HadCRUT.time);
      HadCRUT.NH = squeeze(mean(HadCRUT.tas_mean(:,HadCRUT.latitude > 0,:),[1,2],'omitnan'));
      b = interp1(HadCRUT.MatlabTime,HadCRUT.NH,TimeScale); 
      clear HadCRUT
    case 'NAM'
      NAM = load([Root,'/daily_nam.mat']);
      b = interp1(NAM.Time,NAM.NAM,TimeScale); 
      clear NAM
    case 'NAO'
      NAO = load([Root,'/nao.mat']);
      b = interp1(NAO.Time,NAO.NAO,TimeScale); 
      clear NAO      
    case 'SSTs'
      SSTs = load([Root,'/ssts.mat']);
      b = interp1(SSTs.Time,SSTs.SSTs,TimeScale); 
      clear SSTs      
    case 'TSI'
      TSI = load([Root,'/tsi.mat']);
      TSI.TSI(TSI.TSI < 1358) = NaN; %very noisy - remove an extreme outlier
      TSI.TSI = inpaint_nans(TSI.TSI); %fill the new gap
      TSI.TSI = smoothn(TSI.TSI,[31,1]);%very noisy - smooth to make it fair
      b = interp1(TSI.Time,TSI.TSI,TimeScale); 
      clear TSI
    case 'Time'
      b = TimeScale; 
    case 'SeaIce'
      SeaIce = readmatrix([Root,'/N_seaice_extent_daily_v3.0.csv']);
      t = datenum(SeaIce(:,1),SeaIce(:,2),SeaIce(:,3));
      b = interp1(t,SeaIce(:,4),TimeScale); 
      clear SeaIce t
    case 'AMO'
      AMO = load([Root,'/AMO.mat']);
      b = interp1(AMO.Time,AMO.AMO,TimeScale); 
      clear AMO
    case 'Annual'
      b = round(sin(pi/2 + date2doy(TimeScale)./366.*2.*pi),3); % 1 is Jan 1st, -1 is Jul 1
    otherwise
      disp(['Index ',Settings.Indices{iIndex},' not specified; stopping'])
      stop
  end
 
  %store
  IndexStore.(  Settings.Indices.List{iIndex}) = b;  

  disp(['Loaded ',Settings.Indices.List{iIndex}])

  clear b
end; clear iIndex Root 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% produce a deseasonalised form of each index
%this will not be valid if we're using less than a few years of data
%but that's a choice that can be made later
%
%also, several will not be physically meaningful. Again, deal with this downstream.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find day-of-year for each point
dd   = floor(date2doy(TimeScale));


for iIndex=1:1:numel(Settings.Indices.List);


  %get index
  a = IndexStore.(Settings.Indices.List{iIndex});

  %find day-of-yearly averages
  DV = NaN(366,1);
  for iDay=1:1:366;
    DV(iDay) = mean(a(dd == iDay),'omitnan');
  end; clear iDay


  %if we don't have a value for day 366 due to no leap year data, interpolate one so the smoothing is safe
  if sum(isnan(DV)) == 1 && isnan(DV(366)); DV(366) = mean(DV([1,365])); end
  
  %smooth a bit in time
  DV = smoothn([DV;DV;DV],[Settings.Indices.DSSmooth,1]); 


  %extract daily values
  for iDay=1:1:366;
    a(dd == iDay) = a(dd == iDay)-DV(iDay+366); %the +366 is so we use the 'middle' (i.e. correctly smoothed at ends) time series
  end; clear iDay 

  %and store
  IndexStore.Raw.(Settings.Indices.List{iIndex})= IndexStore.(Settings.Indices.List{iIndex});
  IndexStore.DS.( Settings.Indices.List{iIndex})= a;
  IndexStore = rmfield(IndexStore,Settings.Indices.List{iIndex});


end; clear iIndex a b DV dd

IndexStore.TimeScale = TimeScale; clear TimeScale


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% normalise the data (based on the daily data as these are less biased than the flight data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iForm=1:2; %1 is raw, 2 is deseasonalised

  %get the data

  if iForm==1;      Daily = IndexStore.Raw;
  elseif iForm ==2; Daily = IndexStore.DS;
  end
  RangeStore = struct();

  for iIndex =1:1:numel(Settings.Indices.List)

    %compute the necessary stats from the daily data
    RangeStore.(Settings.Indices.List{iIndex}) = prctile(IndexStore.Raw.(Settings.Indices.List{iIndex}),Settings.Indices.IndexRange);
    
    %apply the stats to normalise the data
    a = Daily.(Settings.Indices.List{iIndex});
    a = a./range(RangeStore.(Settings.Indices.List{iIndex}));
    a = a-min(a);
    a= (a.*2)-1;
    Daily.(Settings.Indices.List{iIndex}) = a;

  end; clear iIndex

 
  %and store
  if iForm==1;
    DateIndices.Raw   = Daily;    DateIndices.Raw.Ranges   = RangeStore;
  elseif iForm ==2;
    DateIndices.DS   = Daily;    DateIndices.DS.Ranges   = RangeStore;
  end  

  clear Daily Flights a
end; clear iForm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reformat the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%create a data table to keep the values in
VarNames = ["Date",Settings.Indices.List];
VarTypes = repmat("double",numel(VarNames),1);

%create table
Out.Raw = table('Size',[numel(IndexStore.TimeScale),numel(VarNames)], ...
                'VariableTypes',VarTypes,          ...
                'VariableNames',VarNames           );
Out.DS = Out.Raw;

%fill table
Out.Raw{:,1} = IndexStore.TimeScale';
Out.DS{ :,1} = IndexStore.TimeScale';
for iIndex=2:1:numel(VarNames); 
  Out.Raw{:,iIndex} = DateIndices.Raw.(VarNames{iIndex})'; 
  Out.DS{ :,iIndex} = DateIndices.DS.( VarNames{iIndex})'; 
end

IndicesDS = Out.DS;
IndicesRaw  = Out.Raw;

clear Out iIndex IndexStore DateIndices


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now we're going to alter some of the indices based on the requests
% in the settings file. 
%
%To do this, we'll take the 'raw' data as the baseline, then progressively
%replace it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%copy over raw to 'output', for clarity
Output = IndicesRaw; clear IndicesRaw

%% which datasets do we want to be deseasonalised?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);

%and replace those we want to be deseasonalised
for iIndex=1:1:numel(Settings.Indices.DS)
  if ~isTableCol(Output,Settings.Indices.DS{iIndex}); continue; end %occurs if we're asking to deseasonalise an index we never created.
  Output.(Settings.Indices.DS{iIndex}) = IndicesDS.(Settings.Indices.DS{iIndex});
end; clear iIndex

clear IndicesDS

%% which datasets do we want to be delinearised?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iIndex=1:1:numel(Settings.Indices.DL)
  if ~isTableCol(Output,Settings.Indices.DL{iIndex}); continue; end  %occurs if we're asking to delinearise an index we never created.
  p = polyfit(Output.Date,Output.(Settings.Indices.DL{iIndex}),1);
  tD = polyval(p,Output.Date);
  Output.(Settings.Indices.DL{iIndex}) = Output.(Settings.Indices.DL{iIndex}) - tD;
  
end; clear iIndex tD p

%% smooth indices?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%smooth?
if Settings.Indices.SmoothLength ~= 0;

  for iIndex=1:1:numel(Settings.Indices.List)
    if contains(Settings.Indices.List{iIndex},Settings.Indices.DoNotSmooth) == 1; continue; end

    %daily data
    a = smoothn(Output.(Settings.Indices.List{iIndex}),[Settings.Indices.SmoothLength,1]);
    Output.(Settings.Indices.List{iIndex}) = a;    

  end

end
clear a iIndex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finally, interpolate to the individual flight times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = Output.Date;
xq = FlightData.Date;

%create table
FlightIndices = table('Size',[numel(xq),numel(VarNames)+1],     ...
                      'VariableTypes',["double";VarTypes],    ...
                      'VariableNames',["FlightIndex",VarNames]);

FlightIndices.Date = xq;
for iIndex=2:1:numel(VarNames); 
  FlightIndices.(VarNames{iIndex}) = interp1(x,Output.(VarNames{iIndex}),xq);
end

FlightIndices.FlightIndex = FlightData.FlightIndex; %we'll use this for table-joins later

clear x xq VarNames VarTypes

DateIndices   = Output; clear Output




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



save([Settings.Paths.DataDir,'/',Settings.ID,'_indices.mat'],'FlightIndices','DateIndices','RangeStore')

disp('--------------------------')
disp('Climate indices generated')
disp('--------------------------')

end

