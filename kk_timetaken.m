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

%number of airports
NPorts = numel(RouteData.Airports);

%list of seasons
Seasons = fieldnames(Settings.Seasons);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute relative flight times, and split into Eward and Wward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FlightData.Results.tRel  = NaN.*FlightData.Results.t;
FlightData.Results.Eward = NaN.*FlightData.Results.t;

for iDep=1:1:NPorts
  for iArr=1:1:NPorts
   for iSeason=1:1:numel(Seasons)
    
     %find all flights meeting these criteria
     Flights = squeeze(RouteData.Flights(iSeason,iDep,iArr,:));
     Flights = Flights(~isnan(Flights));
     if numel(Flights) == 0; continue; end

     %what is the relative time taken by each flights?
     FlightData.Results.tRel(Flights) = FlightData.Results.t(Flights)./nanmedian(FlightData.Results.t(Flights));
     
     %are we heading east or west?
     if ismember(RouteData.Airports{iDep},Settings.NA); E = 1; else E = 2; end  %1 eastward, 2 westward
     FlightData.Results.Eward(Flights) = E;
     
   end
  end
end
clear iDep iArr iSeason Flights NPorts E

%also compute the mean overall flight time, so we can convert our
%coefficients from relative time to minutes
MedianFlightTime = nanmedian(FlightData.Results.t(:))./60;


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
