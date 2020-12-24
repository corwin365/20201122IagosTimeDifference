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
Settings.DataFile = '../03CleanerFlights/flightpairs.mat';

%indices to use
Settings.Indices = {'CCFast','CCSlow'};

%minimum points for comparison
Settings.MinPoints = 20;

%outlier definition - flights this far off the median will be excluded
Settings.Outlier = [0.9,1.1]; %proportion of median time

%seasons - days of year in each
Settings.Seasons.DJF = date2doy(datenum(2000,12,1):datenum(2001, 3,1)-1);
Settings.Seasons.MAM = date2doy(datenum(2000, 3,1):datenum(2000, 6,1)-1);
Settings.Seasons.JJA = date2doy(datenum(2000, 6,1):datenum(2000, 9,1)-1);
Settings.Seasons.SON = date2doy(datenum(2000, 9,1):datenum(2000,12,1)-1);

Settings.Seasons.All = 1:1:366;


% % % %bootstrap properties
% % % Settings.BS.Straps  = 2000;
% % % Settings.BS.Samples = 2000;

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
    Data.Results.t(ThisPair) = Data.Results.t(ThisPair)./MedianTime; %normalised flight time
              
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

% % % %there is a single monster outlier. remove it
% % % Bad = find(Results(:,4) > 3);
% % % Results(Bad,:) = NaN;

%outlier removal
Bad = find(Results(:,4) < Settings.Outlier(1).*nanmedian(Results(:,4)) ...
         | Results(:,4) > Settings.Outlier(2).*nanmedian(Results(:,4)));
Results(Bad,:) = NaN;
disp([num2str(numel(Bad)),' outlier time series removed'])
clear Bad

clearvars -except TimeScale Results Settings Data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load indices and interpolate to data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the indices, and interpolate to each data point
TimeScale = Results(:,3);

for iIndex=1:1:numel(Settings.Indices)
  
  switch Settings.Indices{iIndex}
    case 'CCFast'
      HadCRUT = rCDF([LocalDataDir,'/Miscellany/HadCRUT.4.6.0.0.median.nc']);
      HadCRUT.MatlabTime = datenum(1850,1,HadCRUT.time);
      HadCRUT.NH = smoothn(squeeze(nanmean(HadCRUT.temperature_anomaly(:,HadCRUT.latitude > 0,:),[1,2])),13*3);
      a = interp1(HadCRUT.MatlabTime,HadCRUT.NH,TimeScale);
      clear HadCRUT
    case 'CCSlow'
      HadCRUT = rCDF([LocalDataDir,'/Miscellany/HadCRUT.4.6.0.0.median.nc']);
      HadCRUT.MatlabTime = datenum(1850,1,HadCRUT.time);
      HadCRUT.NH = squeeze(nanmean(HadCRUT.temperature_anomaly(:,HadCRUT.latitude > 0,:),[1,2]));
      HadCRUT.NH = HadCRUT.NH - smoothn(HadCRUT.NH,13*3);
      a = interp1(HadCRUT.MatlabTime,HadCRUT.NH,TimeScale);
      clear HadCRUT
  end
  Results(:,7+iIndex) = a;
  clear a
end; clear iIndex TimeScale

%retain the average time of a flight, in minutes
AverageFlightTime = nanmean(Results(:,7)).*60;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do linear regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Seasons = fieldnames(Settings.Seasons);

for iSeason=1:1:numel(Seasons)
  for EW = 1:2
  
    %which points are in this season and direction?
    ThisSeason = find(ismember(Results(:,6),Settings.Seasons.(Seasons{iSeason})));
    Indices = intersect(ThisSeason,find(Results(:,5) == EW));
    
    %pull out relative time taken
    TimeTaken = Results(Indices,4);
    
    %and the regression time series
    RegSeries = Results(Indices,8:end);
    
    %do the regression
    mdl = fitlm(RegSeries,TimeTaken);
    Coefs = table2array(mdl.Coefficients);
    
% % %     %print out results
% % %     disp('===========================')
% % %     if     EW == 1; disp([Seasons{iSeason},' Eward'])
% % %     elseif EW == 2; disp([Seasons{iSeason},' Wward']);
% % %     end
% % %     mdl.Coefficients 
    
    %and store the outputs
    Reg.Est(iSeason,EW,:) = Coefs(2:end,1);
    Reg.SE( iSeason,EW,:) = Coefs(2:end,2);
    Reg.T(  iSeason,EW,:) = Coefs(2:end,3);
    Reg.P(  iSeason,EW,:) = Coefs(2:end,4);
    

% Indices = Results(:,8:end);
  end
end
clear mdl TimeTaken RegSeries ThisSeason Indices EW Coefs
clear Results iSeason EW


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot regression results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
set(gcf,'color','w')

Symbols = 'sd<o^>';
Colours = 'rgbycm';

for iSeason=1:1:numel(Seasons)
  for EW = 1:2
    
    %generate panel
    subplot(numel(Seasons),2,(iSeason-1)*2+EW)
    
    %plot each regression coefficient
    for iCoef=1:1:size(Reg.Est,3)
      
      %get value
      Value = Reg.Est(iSeason,EW,iCoef);
      
      %convert value from fraction of flight to minutes for mean flight
      Value = Value.*AverageFlightTime;
      
      
      %get significane
      Sig   = Reg.P(iSeason,EW,iCoef);
      if Sig < 0.05; Alpha = 1; else Alpha = 0; end
      
      
      h = plot(Value,iCoef,'marker',Symbols(iCoef), ...
               'color','k','markerfacecolor',Colours(iCoef), ...
               'markersize',20);
      setMarkerColor(h,Colours(iCoef),Alpha);
      hold on
      
      %indicate coefficient standard error
      StErr = Value + [-1,1] .* AverageFlightTime .* Reg.SE(iSeason,EW,iCoef);
      plot(StErr,[1,1].*iCoef,'k-','linewi',3)
      
      
    end
    
    
    %tidy
    set(gca,'ytick',[])
    ylabel(Seasons{iSeason},'fontsize',24);
    if EW == 2; set(gca,'yaxislocation','right'); end
    plot([0,0],[0,size(Reg.Est,3)+1],'k-','linewi',3)
    xlim([-1,1].*1.05.*max(Reg.Est(:)).*AverageFlightTime)
    ylim([0,size(Reg.Est,3)+1])
    
    if iSeason == 1 | iSeason == numel(Seasons);
      if EW == 1; title('Eastward'); else; title('Westward'); end
      if iSeason == 1; set(gca,'xaxislocation','top'); end
      end
    
  end
end

%% key - relative to last plot
for iCoef=1:1:size(Reg.Est,3)
  plot(-20+3.*iCoef,-1.5,'clipping','off','markersize',20, ...
       'marker',Symbols(iCoef),'color','k','markerfacecolor',Colours(iCoef))
  text(-20+3.*iCoef,-2,Settings.Indices{iCoef},'fontsize',12,'horizontalalignment','center')  
  
end