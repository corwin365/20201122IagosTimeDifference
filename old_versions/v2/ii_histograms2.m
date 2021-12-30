function [] = ii_histograms2(Settings)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate and plot histograms of relative flight time under different
%climate indices - composite bootstrapped round-trips
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/12/28
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
%% produce data on the pairs, and generate histogram data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hists   = NaN(numel(Seasons),numel(Settings.Indices),2,numel(Settings.HistBins)-1);
CutOffs = NaN(numel(Seasons),numel(Settings.Indices),2);

for iSeason=1:1:numel(Seasons)
  
  %get list of pairs
  ThisSeason = squeeze(FlightPairs(iSeason,:,:));
  ThisSeason = ThisSeason(find(~isnan(ThisSeason(:,1))),:);
  
  %get relative times and sum for the round trip
  tRel = sum(FlightData.Results.tRel(ThisSeason),2);
  
  %get a mean value for each index over both legs
  %should be very similar, this is just me being fussy!
  IndexVals = NaN(numel(tRel),numel(Settings.Indices));
  for iIndex=1:1:numel(Settings.Indices)
     ThisIndex = Indices.Indices.(Settings.Indices{iIndex});     
     IndexVals(:,iIndex) = nanmean(ThisIndex(ThisSeason),2);
  end
  clear iIndex ThisIndex
  
  %generate histograms
  for iIndex=1:1:numel(Settings.Indices)

    %get index
    Index = IndexVals(:,iIndex);
    
    %work out the cutoff values for the two histograms
    CutOff = prctile(Index,[Settings.Frac,1-Settings.Frac].*100);
    CutOffs(iSeason,iIndex,:) = CutOff;
    
    %hence, generate the two histograms
    tLow  = tRel(find(Index < CutOff(1)));
    tHigh = tRel(find(Index > CutOff(2)));
  
    yLow  = histcounts( tLow,Settings.HistBins); 
    yHigh = histcounts(tHigh,Settings.HistBins);
    
    Hists(iSeason,iIndex,1,:) = yLow;
    Hists(iSeason,iIndex,2,:) = yHigh;
  end; clear iIndex CutOff tLow tHigh yLow yHigh
  

end; clear iSeason

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the histograms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%convert bin *edges* in #frac flight time# to bin *centres* in #minutes#
BinCent = Settings.HistBins(1:end-1) + mean(diff(Settings.HistBins))./2;
BinCent = (BinCent - 2).*MedianFlightTime;

%compute start and end bin locations (added for patch object)
a = BinCent(  1)-mean(diff(BinCent))./2;
b = BinCent(end)+mean(diff(BinCent))./2;

%smooth the data
Hists = smoothn(Hists,[1,1,1,Settings.HistSmooth]);


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
    ysA = [0,squeeze(Hists(iSeason,iIndex,1,:))',0];
    patch(xs,ysA,'b','facealpha',0.7,'edgecolor','b')
    hold on
    
    %large values
    ysB = [0,squeeze(Hists(iSeason,iIndex,2,:))',0];
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
sgtitle('Round-Trips','fontsize',25)