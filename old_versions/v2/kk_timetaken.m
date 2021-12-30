function [] = kk_timetaken(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot time series of time taken over the instrument record
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/02/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load route and flight data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FlightData = load('data/flight_data.mat');
RouteData  = load('data/routes.mat');
Indices    = load('data/indices.mat');

%list of seasons
Seasons = fieldnames(Settings.Seasons);

%flight pairs and median time. This is a bare load rather than a load to
%a struct for backwards compatability below to an older version that computed 
%the data locally
load('data/relative_times.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find monthly-mean notional delay for each plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NotionalDelay = (FlightData.Results.tRel-1) .* MedianFlightTime;

%time scale
[yy1,~,~]  = datevec(min(FlightData.Results.Date));
[yy2,~,~]  = datevec(max(FlightData.Results.Date));
Scale = datenum(yy1,1:2:(yy2-yy1+1).*12,1);
clear yy1 yy2

%indices for planes
[PlaneNames,~,ic]= unique(FlightData.Results.PlaneID);

%storage array
Store = NaN(numel(PlaneNames),numel(Scale),2);

for iEW=1:2
  for iPlane=1:1:numel(PlaneNames)
    ThisGroup = find(ic == iPlane & FlightData.Results.Eward == iEW);
    
    Store(iPlane,:,iEW) = bin2matN(1,FlightData.Results.Date(ThisGroup), ...
                                     NotionalDelay(ThisGroup),           ...
                                     Scale,'@nanmean');
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
clf
set(gcf,'color','w')

subplot = @(m,n,p) subtightplot (m, n, p, 0.04, 0.1, [0.07 0.2]);
Colours = flipud(cbrewer('qual','Set2',numel(PlaneNames)));



for iEW=1:2;
  
  
  %make panel
  subplot(2,1,iEW)
  
  %get data
  Data = Store(:,:,iEW);
 
  
  %plot plane time series
  for iPlane=1:1:numel(PlaneNames)
    plot(Scale,Data(iPlane,:),'o-','color',Colours(iPlane,:),'linewi',2);%,'markerfacecolor',Colours(iPlane,:))
    hold on
    if iEW == 1;
      text(datenum(2021,1,1),20-15.*(iPlane-1),PlaneNames(iPlane),'color',Colours(iPlane,:),'fontweight','bold')
    end
  end
  hold on
  
  %linear trend across all data
  a = nanmean(Data,1);
  Good = find(~isnan(a));
  p = polyfit(Scale(Good),a(Good),1);
  plot(Scale,polyval(p,Scale),'k--')
  
  
  %tidy up
  datetick
  axis([datenum([1994,2020],1,1),-30,50])
  if iEW == 1; set(gca,'xaxislocation','top'); end
  set(gca,'tickdir','out')
  set(gca,'xtick',datenum(1995:3:2020,1,1),'xticklabel',datestr(datenum(1995:3:2020,1,1),'yyyy'))
  ylabel('Delay [min]')
  
  
end
