function process_climate_indices_v3(Settings,LagMode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load climate indices and interpolate them to flight times
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/01/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%lagmode lets us use this same function if we want the data lagged
%it's off by default, and only the lag calculator function sets it
if ~exist('LagMode','var'); LagMode = 0; end

if LagMode == 0;

  disp('+++++++++++++++++++++++++++')
  disp('Generating climate indices')
  disp('+++++++++++++++++++++++++++')

else

  disp('+++++++++++++++++++++++++++')
  disp('Lagging climate indices')
  disp('+++++++++++++++++++++++++++')

  load([Settings.Paths.DataDir,'/',Settings.ID,'_optimallags.mat'])

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load flight data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat'])

%generate a continuous timescale at daily cadence. Add a year on each end if possible so we can circshift them for lagging
TimeScale  = min(Settings.Choices.TimeRange)-365:1:max(Settings.Choices.TimeRange)+365;

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
    case 'u1060'
      U_1060 = load([Root,'/u_1060.mat']);
      b = interp1(U_1060.Time,U_1060.U,TimeScale); 
      clear U_1060     
    case 'TSI'
      TSI = load([Root,'/tsi.mat']);
      TSI.TSI(TSI.TSI < 1358) = NaN; %very noisy - remove an extreme outlier
      TSI.TSI = inpaint_nans(TSI.TSI); %fill the new gap
      TSI.TSI = smoothn(TSI.TSI,[31,1]);%very noisy - smooth to make it fair
      b = interp1(TSI.Time,TSI.TSI,TimeScale); 
      clear TSI
    case 'Time'
      %this looks weird, as we needs a small correction to put the correct part into -1 to 1, since
      %we allowed some extra for circshifting on the master timescale defined above
      b = interp1(minmax(TimeScale),Settings.Choices.TimeRange,TimeScale); 
    case 'SeaIce'
      SeaIce = readmatrix([Root,'/N_seaice_extent_daily_v3.0.csv']);
      t = datenum(SeaIce(:,1),SeaIce(:,2),SeaIce(:,3));
      b = interp1(t,SeaIce(:,4),TimeScale); 
      clear SeaIce t
    case 'AMO'
      AMO = load([Root,'/amo_noaa.mat']);      
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
%% if we want to apply lags, do it now
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if LagMode == 1;

  %shift the time series by the lag for this index
  for iIndex=1:1:numel(Settings.Indices.List)
    for iDir=1:1:numel(Settings.Choices.Directions)
      dir = Settings.Choices.Directions{iDir};
      index = Settings.Indices.List{iIndex};
      IndexStore.(dir).(index) = circshift(IndexStore.(index),-Lags.(dir).(index));
    end
    IndexStore = rmfield(IndexStore,index);
  end
  clear dir index iIndex iDir

else

  %no lag used 
  %make copies of the index for each direction, so we can use the same code 
  %as for the lagged version
  for iIndex=1:1:numel(Settings.Indices.List)
    for iDir=1:1:numel(Settings.Choices.Directions)
      dir = Settings.Choices.Directions{iDir};
      index = Settings.Indices.List{iIndex};
      IndexStore.(dir).(index) = IndexStore.(index);
    end
    IndexStore = rmfield(IndexStore,index);
  end
  clear dir index iIndex iDir
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% produce a deseasonalised form of each index
%this will not be valid if we're using less than a few years of data
%but that's a choice that can be made later
%
%also, several will not be physically meaningful. Again, deal with this downstream.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find day-of-year for each point
dd   = floor(date2doy(TimeScale));

for iDir=1:1:numel(Settings.Choices.Directions)
  for iIndex=1:1:numel(Settings.Indices.List);

    %get index
    a = IndexStore.(Settings.Choices.Directions{iDir}).(Settings.Indices.List{iIndex});

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
    IndexStore.Raw.(Settings.Choices.Directions{iDir}).(Settings.Indices.List{iIndex})= IndexStore.(Settings.Choices.Directions{iDir}).(Settings.Indices.List{iIndex});
    IndexStore.DS.( Settings.Choices.Directions{iDir}).(Settings.Indices.List{iIndex})= a;


  end
  IndexStore = rmfield(IndexStore,Settings.Choices.Directions{iDir});
end; clear iIndex a b DV dd iDir

IndexStore.TimeScale = TimeScale; clear TimeScale


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% normalise the data (based on the daily data as these are less biased than the flight data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%store original data ranges
RangeStore = struct();

for iDir=1:1:numel(Settings.Choices.Directions)
  for iForm=1:2; %1 is raw, 2 is deseasonalised

    %get the data
    if iForm==1;      Daily = IndexStore.Raw.(Settings.Choices.Directions{iDir});
    elseif iForm ==2; Daily = IndexStore.DS.(Settings.Choices.Directions{iDir});
    end
    RangeStore.(Settings.Choices.Directions{iDir}) = struct();

    for iIndex =1:1:numel(Settings.Indices.List)

      %compute the necessary stats from the daily data
      RangeStore.(Settings.Choices.Directions{iDir}).(Settings.Indices.List{iIndex}) = prctile(IndexStore.Raw.(Settings.Choices.Directions{iDir}).(Settings.Indices.List{iIndex}),Settings.Indices.IndexRange);

      %choose normalisation type
      if strcmp(Settings.Indices.NormType,'r') %range normalisation

        %apply the stats to normalise the data
        a = Daily.(Settings.Indices.List{iIndex});
        a = a./range(RangeStore.(Settings.Choices.Directions{iDir}).(Settings.Indices.List{iIndex}));
        a = a-min(a);
        a= (a.*2)-1;
        Daily.(Settings.Choices.Directions{iDir}).(Settings.Indices.List{iIndex}) = a;

      elseif strcmp(Settings.Indices.NormType,'z'); %z-score normalise

        a = Daily.(Settings.Indices.List{iIndex});
        m = mean(a,'omitnan');
        s = std( a,'omitnan');
        Daily.(Settings.Choices.Directions{iDir}).(Settings.Indices.List{iIndex}) = ((a - m)./s);

      else

        disp('Error: invalid normalisation option set for indices')
        stop
      end
    end; clear iIndex


    %and store
    if iForm==1;
      DateIndices.Raw.(Settings.Choices.Directions{iDir})   = Daily.(Settings.Choices.Directions{iDir});   
      DateIndices.Raw.Ranges.(Settings.Choices.Directions{iDir})   = RangeStore;
    elseif iForm ==2;
      DateIndices.DS.(Settings.Choices.Directions{iDir})   = Daily.(Settings.Choices.Directions{iDir});    
      DateIndices.DS.Ranges.(Settings.Choices.Directions{iDir})   = RangeStore;
    end

    clear Daily Flights a
  end; clear iForm

end; clear iDir


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reformat the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a data table to keep the values in
VarNames = ["Date",Settings.Indices.List];
VarTypes = repmat("double",numel(VarNames),1);

%create table
for iDir=1:1:numel(Settings.Choices.Directions)

  Out.(Settings.Choices.Directions{iDir}).Raw = table('Size',[numel(IndexStore.TimeScale),numel(VarNames)], ...
                                                      'VariableTypes',VarTypes,          ...
                                                      'VariableNames',VarNames           );
  Out.(Settings.Choices.Directions{iDir}).DS = Out.(Settings.Choices.Directions{iDir}).Raw;

  %fill table
  Out.(Settings.Choices.Directions{iDir}).Raw{:,1} = IndexStore.TimeScale';
  Out.(Settings.Choices.Directions{iDir}).DS{ :,1} = IndexStore.TimeScale';
  for iIndex=2:1:numel(VarNames);
    Out.(Settings.Choices.Directions{iDir}).Raw{:,iIndex} = DateIndices.Raw.(Settings.Choices.Directions{iDir}).(VarNames{iIndex})';
    Out.(Settings.Choices.Directions{iDir}).DS{ :,iIndex} = DateIndices.DS.(Settings.Choices.Directions{iDir}).( VarNames{iIndex})';
  end
  IndicesDS.(Settings.Choices.Directions{iDir}) = Out.(Settings.Choices.Directions{iDir}).DS;
  IndicesRaw.(Settings.Choices.Directions{iDir})  = Out.(Settings.Choices.Directions{iDir}).Raw;

end
clear Out iIndex IndexStore DateIndices iDir


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
  for iDir=1:1:numel(Settings.Choices.Directions)
    if ~isTableCol(Output.(Settings.Choices.Directions{iDir}),Settings.Indices.DS{iIndex}); continue; end %occurs if we're asking to deseasonalise an index we never created.

    Output.(Settings.Choices.Directions{iDir}).(Settings.Indices.DS{iIndex}) = IndicesDS.(Settings.Choices.Directions{iDir}).(Settings.Indices.DS{iIndex});
  end
end; clear iIndex iDir

clear IndicesDS


%% which datasets do we want to be delinearised?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iIndex=1:1:numel(Settings.Indices.DL)
  for iDir=1:1:numel(Settings.Choices.Directions)
    if ~isTableCol(Output.(Settings.Choices.Directions{iDir}),Settings.Indices.DL{iIndex}); continue; end  %occurs if we're asking to delinearise an index we never created.
    p = polyfit(Output.(Settings.Choices.Directions{iDir}).Date,Output.(Settings.Choices.Directions{iDir}).(Settings.Indices.DL{iIndex}),1);
    tD = polyval(p,Output.(Settings.Choices.Directions{iDir}).Date);
    Output.(Settings.Choices.Directions{iDir}).(Settings.Indices.DL{iIndex}) = Output.(Settings.Choices.Directions{iDir}).(Settings.Indices.DL{iIndex}) - tD;
  end; clear iDir
end; clear iIndex tD p


%% smooth indices?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%smooth?
if Settings.Indices.SmoothLength ~= 0;

  for iIndex=1:1:numel(Settings.Indices.List)
    if contains(Settings.Indices.List{iIndex},Settings.Indices.DoNotSmooth) == 1; continue; end

    %daily data
    for iDir=1:1:numel(Settings.Choices.Directions)
      a = smoothn(Output.(Settings.Choices.Directions{iDir}).(Settings.Indices.List{iIndex}),[Settings.Indices.SmoothLength,1]);
      Output.(Settings.Choices.Directions{iDir}).(Settings.Indices.List{iIndex}) = a;
    end; clear iDir
  end

end
clear a iIndex


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finally, interpolate to the individual flight times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDir=1:1:numel(Settings.Choices.Directions)

  x = Output.(Settings.Choices.Directions{iDir}).Date;
  xq = FlightData.Date;

  %create table
  FlightIndices.(Settings.Choices.Directions{iDir}) = table('Size',[numel(xq),numel(VarNames)+1],     ...
                                                            'VariableTypes',["double";VarTypes],    ...
                                                            'VariableNames',["FlightIndex",VarNames]);

  FlightIndices.(Settings.Choices.Directions{iDir}).Date = xq;
  for iIndex=2:1:numel(VarNames);
    FlightIndices.(Settings.Choices.Directions{iDir}).(VarNames{iIndex}) = interp1(x,Output.(Settings.Choices.Directions{iDir}).(VarNames{iIndex}),xq);
  end

  FlightIndices.(Settings.Choices.Directions{iDir}).FlightIndex = FlightData.FlightIndex; %we'll use this for table-joins later

  
end; clear iDir
clear x xq VarNames VarTypes
DateIndices = Output; clear Output




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if LagMode == 0; save([Settings.Paths.DataDir,'/',Settings.ID,'_indices.mat'],      'FlightIndices','DateIndices','RangeStore')
else             save([Settings.Paths.DataDir,'/',Settings.ID,'_laggedindices.mat'],'FlightIndices','DateIndices','RangeStore')
end

disp('--------------------------')
disp('Climate indices generated')
disp('--------------------------')

end

