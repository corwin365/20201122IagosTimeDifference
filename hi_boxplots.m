function [] = hi_boxplots(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate and plot box plots of relative flight time under different
%climate indices
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
MedianFlightTime = nanmedian(FlightData.Results.t(:))./60; %MINUTES


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

%smooth the data
Hists = smoothn(Hists,[1,1,1,1,Settings.HistSmooth]);

%convert bin *edges* in #frac flight time# to bin *centres* in #minutes#
BinCent = Settings.HistBins(1:end-1) + mean(diff(Settings.HistBins))./2;
BinCent = (Settings.HistBins.*MedianFlightTime)-MedianFlightTime;

%compute start and end bin locations (added for patch object)
a = (Settings.HistBins(  1).*MedianFlightTime)-MedianFlightTime;
b = (Settings.HistBins(end).*MedianFlightTime)-MedianFlightTime;

Seasons = fieldnames(Settings.Seasons);

k = 0;
figure
set(gcf,'color','w')
clf
subplot = @(m,n,p) subtightplot (m, n, p, 0.02, 0.1, 0.2);

for iEW=1:2 %loop over east west
  
  for iSeason=1:1:numel(Seasons);
    
    %create panel
    k = k+1; subplot(2,numel(Seasons),k)
    hold on
    
    %create axes
    axis([ -50 50 0.5 numel(Settings.Indices)+0.5])
    set(gca,'tickdir','out','ticklen',[1,1].*0.02)
    
    for iIndex = 1:1:numel(Settings.Indices)
      
      for iSegment = [1,2];
        for iPass=1:2 %first pass fill shades, second pass plot lines
        
          switch iSegment; 
            case 2; Colour = [255,128,000]./255;
            case 1; Colour = [000,255,255]./255;
          end
        
        %get data
        Data = [0;squeeze(Hists(iSeason,iEW,iIndex,iSegment,:))];%the 0 lines it up with BinCent
        
        %plot box
        x = [BinCent,BinCent(end:-1:1)];
        y = [Data;-Data(end:-1:1)]'; y = y./nansum(abs(y)).*8; y = y+iIndex;
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
        Data = [0;squeeze(Hists(iSeason,iEW,iIndex,iSegment,:))];%the 0 lines it up with BinCent
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
    if iEW == 1; set(gca,'xaxislocation','top'); else set(gca,'xaxislocation','bottom'); end
    set(gca,'xtick',-100:25:100)
    set(gca,'ydir','reverse','fontsize',12)
    if iEW == 1;
      xlabel(Seasons{iSeason},'fontsize',28)
    else
      xlabel('Delay [minutes]','fontsize',12)
    end
    plot([0,0],[0,numel(Settings.Indices)+0.5],'k-')
    
    box on
    grid on
    drawnow
    
  end
end
