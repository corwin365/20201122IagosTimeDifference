function [] = hh_histograms(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate and plot histograms of relative flight time under different
%climate indices
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/12/27
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
%% find the top and bottom X% against each index, and split 
%the data into those groups for each direction/season
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hists   = NaN(numel(Seasons),2,numel(Settings.Indices),2,numel(Settings.HistBins)-1);
CutOffs = NaN(numel(Seasons),2,numel(Settings.Indices),2);


for iSeason=1:1:numel(Seasons)
  for iEW=1:2; %1 eastward, 2 westward
    
    %find all flights in this season
    Flights = flatten(RouteData.Flights(iSeason,:,:,:));
    Flights = Flights(~isnan(Flights));
    
    %and direction
    Flights = Flights(FlightData.Results.Eward(Flights) == iEW);
    
    %now, loop over indices
    for iIndex=1:1:numel(Settings.Indices)
      
      %pull out all the values for this set of flights
      Index = Indices.Indices.(Settings.Indices{iIndex});
      Index = Index(Flights);
      
      %work out the cutoff values for the two histograms
      CutOff = prctile(Index,[Settings.Frac,1-Settings.Frac].*100);
      CutOffs(iSeason,iEW,iIndex,:) = CutOff;
      
      %hence, generate the two histograms
      tLow  = FlightData.Results.tRel(Flights(find(Index < CutOff(1))));
      tHigh = FlightData.Results.tRel(Flights(find(Index > CutOff(2)))); 
      yLow  = histcounts( tLow,Settings.HistBins);
      yHigh = histcounts(tHigh,Settings.HistBins);

      Hists(iSeason,iEW,iIndex,1,:) = yLow;
      Hists(iSeason,iEW,iIndex,2,:) = yHigh;
    end
    
    
    
  end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%convert bin *edges* in #frac flight time# to bin *centres* in #minutes#
BinCent = Settings.HistBins(1:end-1) + mean(diff(Settings.HistBins))./2;
BinCent = (Settings.HistBins.*MedianFlightTime)-MedianFlightTime;

%compute start and end bin locations (added for patch object)
a = (Settings.HistBins(  1).*MedianFlightTime)-MedianFlightTime;
b = (Settings.HistBins(end).*MedianFlightTime)-MedianFlightTime;

%smooth the data
Hists = smoothn(Hists,[1,1,1,1,Settings.HistSmooth]);

for iEW=1:2
  
  %create figure
  figure
  set(gcf,'color','w')
  subplot = @(m,n,p) subtightplot (m, n, p, 0.03, [0.05,0.15], 0.05);

  k = 0;
  for iSeason=1:1:numel(Seasons)
    for iIndex=1:1:numel(Settings.Indices)
      
      %prepare panel
      k = k+1;
      subplot(numel(Seasons),numel(Settings.Indices),k)
      
      %plot the histograms
      xs = [a,BinCent,b];
      
      %small values
      ysA = [0,squeeze(Hists(iSeason,iEW,iIndex,1,:))',0,0];
      patch(xs,ysA,'b','facealpha',0.7,'edgecolor','b')
      hold on
      
      %large values
      ysB = [0,squeeze(Hists(iSeason,iEW,iIndex,2,:))',0,0];
      patch(xs,ysB,'r','facealpha',0.7,'edgecolor','r')
      
      %overlap
      ysC = ysB;
      ysC(:) = 0;
      for iY=1:1:numel(ysC)
        ysC(iY) = min([ysB(iY),ysA(iY)]);
      end
      patch(xs,ysC,[178,102,204]./255,'facealpha',0.7,'edgecolor','none')
      
      
      %tidy up
      if iSeason == 1; title(Settings.Indices{iIndex},'fontsize',25); end
      if iIndex  == 1; ylabel(Seasons{iSeason},'fontsize',25); end
      
      plot([0,0],[0,max([ysA,ysB,ysC])],'k--')
      axis tight
      
      drawnow
      
    end
  end
  
  if iEW==1; sgtitle('Eastwards','fontsize',25)
  else       sgtitle('Westwards','fontsize',25)
  end
  
end
