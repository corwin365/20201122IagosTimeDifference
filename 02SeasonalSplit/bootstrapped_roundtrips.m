clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyse IAGOS data for effect of ENSO, NAO, QBO and HadCRUT on
%trans-atlantic flight times, winter only
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/11/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file containing airports and flight times
Settings.DataFile = '../01InitialTesting/flightpairs.mat';

%indices to use
Settings.Indices = {'QBO','ENSO','HadCRUT','NAM','TSI','NAO'};

%minimum points for comparison
Settings.MinPoints = 10;

%seasons - days of year in each
Settings.Seasons.DJF = date2doy(datenum(2000,12,1):datenum(2001, 3,1)-1);
Settings.Seasons.MAM = date2doy(datenum(2000, 3,1):datenum(2000, 6,1)-1);
Settings.Seasons.JJA = date2doy(datenum(2000, 6,1):datenum(2000, 9,1)-1);
Settings.Seasons.SON = date2doy(datenum(2000, 9,1):datenum(2000,12,1)-1);

Settings.Seasons.All = 1:1:366;

%how many bands for each index?
Settings.NBands = 4;

%bootstrap properties
Settings.BS.Straps  = 10000;
Settings.BS.Samples = 2000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data = load(Settings.DataFile);

%add day-of-year to the data
Data.Results.DoY = floor(date2doy(Data.Results.Date));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% results arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%master data array
Used = sum(~isnan(Data.Results.t));
Results = NaN(Used,7+numel(Settings.Indices));
clear Used


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find each set of pairwise airports
% also normalise flight time for each pair of airports to their median
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 0;
for iDep = 1:1:numel(Data.Airports)
  for iArr=1:1:numel(Data.Airports)
   
    ThisPair = find(Data.Results.Dep == iDep ...
                  & Data.Results.Arr == iArr);

    if numel(ThisPair) < Settings.MinPoints;
      Data.Results.t(ThisPair) = NaN;
      continue;
    end
                
    MedianTime = nanmedian(Data.Results.t(ThisPair));
    Data.Results.Duration(ThisPair) =   Data.Results.t(ThisPair)./60./60; %flight time in hours
    Data.Results.t(ThisPair) = Data.Results.t(ThisPair)./MedianTime;
              
    %is this flight eastbound or westbound?
    Dep = Data.Airports(iDep);
    Arr = Data.Airports(iArr);
    
    if     ismember(Dep,Data.Settings.NA)  && ismember(Arr,Data.Settings.Eur);
      EW = 1;
    elseif ismember(Dep,Data.Settings.Eur) && ismember(Arr,Data.Settings.NA);
      EW = 2;
    else disp('Error'); stop; end
    
    for iFlight =1:1:numel(ThisPair)
    
      k = k+1;
      %add results to big list
      Line = [iDep,iArr, ...
              Data.Results.Date(ThisPair(iFlight)),    ...
              Data.Results.t(   ThisPair(iFlight)),    ...
              EW,                                      ...
              Data.Results.DoY(ThisPair(iFlight)),     ...
              Data.Results.Duration(ThisPair(iFlight))];
      Results(k,1:7) = Line;
    end
  end
end

%there is a single monster outlier. remove it
Bad = find(Results(:,4) > 3);
Results(Bad,:) = NaN;

clearvars -except TimeScale Results Settings Data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load indices and interpolate to data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the indices, and interpolate to each data point
TimeScale = Results(:,3);

for iIndex=1:1:numel(Settings.Indices)
  
  switch Settings.Indices{iIndex}
    case 'QBO'
      QBO = load([LocalDataDir,'/Miscellany/QBO.mat']);
      QBO.Time = floor(QBO.Time); %shift from noon to midnight to make the logic easier - on a 91-day smoothing this is very minor...
      a = interp1(QBO.Time,QBO.QBO,TimeScale);
      clear QBO
    case 'ENSO'
      ENSO = load([LocalDataDir,'/Miscellany/nino34.mat']);
      a = interp1(ENSO.Time,ENSO.Nino34,TimeScale);
      clear ENSO
    case 'HadCRUT'
      HadCRUT = rCDF([LocalDataDir,'/Miscellany/HadCRUT.4.6.0.0.median.nc']);
      HadCRUT.MatlabTime = datenum(1850,1,HadCRUT.time);
      HadCRUT.NH = squeeze(nanmean(HadCRUT.temperature_anomaly(:,HadCRUT.latitude > 0,:),[1,2]));
      a = interp1(HadCRUT.MatlabTime,HadCRUT.NH,TimeScale);
      clear HadCRUT
    case 'NAM'
      NAM = load([LocalDataDir,'/Miscellany/daily_nam.mat']);
      a = interp1(NAM.Time,NAM.NAM,TimeScale);
      clear NAM
    case 'NAO'
      NAO = load([LocalDataDir,'/Miscellany/nao.mat']);
      a = interp1(NAO.Date,NAO.NAO,TimeScale);
      clear NAO      
    case 'TSI'
      TSI = load([LocalDataDir,'/Miscellany/tsi.mat']);
      a = interp1(TSI.Time,TSI.TSI,TimeScale);
      clear TSI
  end
  Results(:,7+iIndex) = a;
  clear a
end; clear iIndex TimeScale




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bootstrap the dependence of relative flight time on each index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Seasons = fieldnames(Settings.Seasons);
BSOut    = NaN(numel(Settings.Indices),numel(Seasons),Settings.NBands-1,5); %5 is the summary stats
BandsOut = NaN(numel(Settings.Indices),Settings.NBands);


for iIndex=1:1:numel(Settings.Indices);
  
  %pull out values
  Index = Results(:,7+iIndex);
  
  %split up into bands, then loop over them
  Bands = linspace(prctile(Index,5),prctile(Index,95),Settings.NBands);
  BandsOut(iIndex,:) = Bands;
  
  %loop over season
  for iSeason=1:1:numel(Seasons)
    
    %find flights in season
    ThisSeason = find(ismember(Results(:,6),Settings.Seasons.(Seasons{iSeason})));
    
    %loop over bands
    for iBand=1:1:numel(Bands)-1
      
      %find ponts in this band
      InBand = find(Index >= Bands(iBand) & Index <= Bands(iBand+1));
      
      %and in season
      InBand = intersect(InBand,ThisSeason);
      if numel(InBand) < Settings.MinPoints; clear InBand; continue; end

      %convert to relative flight time
      FlightTimes = Results(InBand,4);
      
      %bootstrap!
      %%%%%%%%%%%%
      
      %create index array
      Indices = randi(numel(InBand),Settings.BS.Straps,Settings.BS.Samples);
      
      %create value array
      Values = FlightTimes(Indices); clear Indices
      
      %take statistics of samples
      Means = nanmean(Values,2);
      
      %and the characterise this distribution
      Distrib = prctile(Means,[2.5,18,50,82,97.5]);
      
      %and store!
      BSOut(iIndex,iSeason,iBand,:) = Distrib;
      
      clear InBand FlightTimes Indices Values Means Distrib
      
    end; clear iBand
    
    
  end; clear iSeason ThisSeason  Bands
end; clear iIndex Index Seasons


%tidy up
clear Results
Results = BSOut;
clear BSOut

%shift bin edges to bin centres
BandsOut = BandsOut(:,1:end-1)+ mean(diff(BandsOut,1,2),2)./2;

%scale times to %
Results = (Results-1).*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
set(gcf,'color','w')

k = 0;

Seasons = fieldnames(Settings.Seasons);
 
for iSeason=1:1:numel(Seasons)
  
  for iIndex=1:1:numel(Settings.Indices)
        
    %create panel
    k = k+1;
    subplot(numel(fieldnames(Settings.Seasons)),numel(Settings.Indices),k)
      
    %x-scale of bands
    xs = [BandsOut(iIndex,:),BandsOut(iIndex,end:-1:1)];
    
    %2 stdev band
    ys = [squeeze(Results(iIndex,iSeason,:,5));squeeze(Results(iIndex,iSeason,end:-1:1,1))]';
    patch(xs,ys,[1,1,1].*0.9,'edgecolor','none'); hold on
    
    %1 stdev band
    ys = [squeeze(Results(iIndex,iSeason,:,4));squeeze(Results(iIndex,iSeason,end:-1:1,2))]';
    patch(xs,ys,[1,1,1].*0.6,'edgecolor','none')    
    
    plot(BandsOut(iIndex,:),squeeze(Results(iIndex,iSeason,:,3)),'k-');
    hold on
    
    if iSeason == 1; title(Settings.Indices{iIndex}); end
    if iIndex  == 1; ylabel(Seasons{iSeason}); end
    
  end; clear iSeason
  
end; clear iIndex