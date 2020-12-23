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

%minimum points for comparison. must have this many flights EACH WAY.
%(there are two checks - the latter one is both way, the first one being both isn't a bug!)
Settings.MinPoints = 2;

%seasons - days of year in each
Settings.Seasons.DJF = date2doy(datenum(2000,12,1):datenum(2001, 3,1)-1);
Settings.Seasons.MAM = date2doy(datenum(2000, 3,1):datenum(2000, 6,1)-1);
Settings.Seasons.JJA = date2doy(datenum(2000, 6,1):datenum(2000, 9,1)-1);
Settings.Seasons.SON = date2doy(datenum(2000, 9,1):datenum(2000,12,1)-1);

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
end; clear iIndex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate histograms
%combine east and west, reweighting the distribution to make them equal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ListOfSeasons = fieldnames(Settings.Seasons);
xscale = linspace(0.85,1.15,33);
store  = NaN(numel(ListOfSeasons),numel(xscale),numel(Settings.Indices),2,2);

for iSeason=1:1:numel(ListOfSeasons)
  
  Season = Settings.Seasons.(ListOfSeasons{iSeason});
  ThisSeason = find(ismember(Results(:,6),Season));
  clear Season
  
  
  for EW=1:2;
    
    %find this direction and season
    ThisDir = intersect(find(Results(:,5) == EW),ThisSeason);
    
    for iIndex=1:1:numel(Settings.Indices)
      
      %extract data
      x = Results(ThisDir,4);
      y = Results(ThisDir,7+iIndex);
      
      %split into bottom third and top third, then generate and store scaled histograms
      for Range=[1,2];
        
        switch Range
          case 1; idx = find(y < prctile(y,25));
          case 2; idx = find(y > prctile(y,75));
        end
        
        [ys,xs] = hist(x(idx),xscale);
        
        ys(1) = 0;
        ys(end) = 0;
        
        store(iSeason,:,iIndex,Range,EW) = ys;
      end 
    end
  end
end

%normalise east/west trips to be the same number total
sz = size(store);
for iX=1:1:sz(1)
  for iY=1:1:sz(2)
    for iZ = 1:1:sz(3);
      for iA = 1:1:sz(4)
        
        
        %how many flights each way?
        NE = store(iX,iY,iZ,iA,1);
        NW = store(iX,iY,iZ,iA,2);
        
        
        %is the minimum numnber of flights each way between each airport met?
        if NE < Settings.MinPoints | NW < Settings.MinPoints; 
          store(iX,iY,iZ,iA,:) = 0;  
          continue;
        end
        %scale the data
        Scale = NE ./ NW;
        store(iX,iY,iZ,iA,2) = store(iX,iY,iZ,iA,2).*Scale;
      end
    end
  end
end

store  = sum(store,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot histograms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
k = 0;
for iSeason=1:1:numel(ListOfSeasons)

  for iIndex=1:1:numel(Settings.Indices)
    k = k+1;
    subplot(numel(ListOfSeasons),numel(Settings.Indices),k)
    for Range=1:1:2;
      switch Range
        case 1; Colour = 'r'; case 2; Colour = 'b';
      end
    
      ys = squeeze(store(iSeason,:,iIndex,Range));
      ys = smoothn(ys,5);
      
      
      patch([0.79,xs,1.21],[0,ys,0],Colour,'facealpha',0.3,'edgecolor',Colour,'linewi',1)
      hold on
      xlim([0.85 1.15])
      
      if iIndex == 1;
        ylabel(ListOfSeasons{iSeason},'fontsize',24)
      end
      if iSeason == 1;
        xlabel(Settings.Indices{iIndex},'fontsize',24)
        set(gca,'xaxislocation','top')
      end
    end
  end
end