function [] = ij_boxplots2(Settings)


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
%% plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%smooth the data
Hists = smoothn(Hists,[1,1,1,Settings.HistSmooth]);

%convert bin *edges* in #frac flight time# to bin *centres* in #minutes#
BinCent = Settings.HistBins(1:end-1) + mean(diff(Settings.HistBins))./2;
BinCent = (Settings.HistBins.*MedianFlightTime)-2.*MedianFlightTime;

%compute start and end bin locations (added for patch object)
a = (Settings.HistBins(  1).*MedianFlightTime)-MedianFlightTime;
b = (Settings.HistBins(end).*MedianFlightTime)-MedianFlightTime;

Seasons = fieldnames(Settings.Seasons);

k = 0;
figure
set(gcf,'color','w')
clf
subplot = @(m,n,p) subtightplot (m, n, p, 0.02, 0.1, 0.2);

for iSeason=1:1:numel(Seasons);
  
  %create panel
  k = k+1; subplot(2,numel(Seasons),k)
  hold on
  
  %create axes
  axis([ -50 50 0.5 numel(Settings.Indices)+0.5])
  set(gca,'tickdir','out','ticklen',[1,1].*0.02)
  
  for iIndex = 1:1:numel(Settings.Indices)
    
    for iSegment = [2,1];
      for iPass=1:2 %first pass fill shades, second pass plot lines
        
        switch iSegment;
          case 2; Colour = [255,128,000]./255;
          case 1; Colour = [000,255,255]./255;
        end
        
        %get data
        Data = [0;squeeze(Hists(iSeason,iIndex,iSegment,:))];%the 0 lines it up with BinCent
        
        %plot box
        x = [BinCent,BinCent(end:-1:1)];
        y = [Data;-Data(end:-1:1)]'; y = y./max(y)./2.25; y = y+iIndex;
        if iPass == 1;
          patch(x,y,Colour,'edgecolor','none','facealpha',0.5)
        else
          plot(x,y,'color','k','linewi',0.25)
        end
        
      end
    end
    %overplot distribution maxima
    for iSegment=1:2
      switch iSegment;
        case 2; Colour = [255,128,000]./255;
        case 1; Colour = [000,000,200]./255;
      end
      
      %get data
      Data = [0;squeeze(Hists(iSeason,iIndex,iSegment,:))];%the 0 lines it up with BinCent
      x = [BinCent,BinCent(end:-1:1)];
      y = [Data;-Data(end:-1:1)]'; y = y./nansum(abs(y)).*10; y = y+iIndex;
      [~,idx] = max(y(1:floor(numel(y)./2)));
      plot([1,1].*x(idx),[-1,1].*0.33+mean(y),'color',Colour,'linewi',2)
      
    end
    
    
    
  end
  
  %tidy panel
  if     iSeason ==              1; set(gca,'ytick',1:1:numel(Settings.Indices),'yticklabel',Settings.Indices,'yaxislocation', 'left','fontsize',15);
  elseif iSeason == numel(Seasons); set(gca,'ytick',1:1:numel(Settings.Indices),'yticklabel',Settings.Indices,'yaxislocation','right','fontsize',15);
  else                              set(gca,'yticklabel',[]);
  end
%   if iEW == 1; set(gca,'xaxislocation','top'); else set(gca,'xaxislocation','bottom'); end
  set(gca,'xtick',-100:25:100)
  set(gca,'ydir','reverse','fontsize',12)
  title(Seasons{iSeason},'fontsize',28)
  xlabel('Delay [minutes]','fontsize',12)
  plot([0,0],[0,numel(Settings.Indices)+0.5],'k-')
  
  box on
  grid on
  
end
