function [] = gg_times(Paths,Settings)

Indices  = Settings.Indices;
Seasons  = Settings.Seasons;
RTRange  = Settings.RelativeTime;
Colours  = Settings.Colours;
xPercent = Settings.CutOff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot time series of time taken over the instrument record
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/02/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load route and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%all flight data
Flights = load([Paths.StoreDir,'/flight_data_',Paths.SourceIdentifier,'.mat']);
Flights = Flights.Flights;

%which flights are actually used?
load([Paths.StoreDir,'/routes_',Paths.SourceIdentifier,'_',Paths.PPIdentifier,'.mat']);

%load indices
Index = load([Paths.StoreDir,'/indices_',Paths.SourceIdentifier,'_',Paths.PPIdentifier,'.mat']);
FlightIndices = Index.Flight;
clear Index

%get season names, and add 'all' 
SeasonNames = fieldnames(Seasons);
SeasonNames{end+1} = 'All';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% split out the top and bottom x % by each index, then plot KDFs of the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %prepare figure
%   figure(iDir)
  clf
  set(gcf,'color','w')
  k= 0;

for iDir=1:1:3 %1 is eastward, 2 is westward, 3 is round-trip



 


  
  for iSeason=1:1:numel(SeasonNames)

    k = k+1;
    
    subplot(3,numel(SeasonNames),k)
    set(gca,'color',[1,1,1].*0.8)
    axis([-1,1,-1,1].*22); axis square manual
    hold on
    if iDir == 1; title(SeasonNames{iSeason}); end
    if iSeason == 1; ylabel(iDir);end
    box on
    plot([-1,1].*999,[0,0],'k-')
    plot([0,0],[-1,1].*999,'k-')
%     plot([-1,1].*999,[-1,1].*999,'k--')
    plot([-1,1].*999,[1,-1].*999,'k--')

    %select flights in this season and direction
    %also load their flight times


    if iDir == 1;
      ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1 ...
                     & Working.Eastward                        == true);
      TotalFlightTime = Working.tRel.(SeasonNames{iSeason}); 
      TotalFlightTime = TotalFlightTime(ThisGroup);
    elseif iDir == 2;
      ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1 ...
                     & Working.Eastward                        == false);
      TotalFlightTime = Working.tRel.(SeasonNames{iSeason}); 
      TotalFlightTime = TotalFlightTime(ThisGroup);
    elseif iDir == 3;
      ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1 ...
                     & ~isnan(Working.Pair));
      TotalFlightTime = Working.tRel.(SeasonNames{iSeason}); 
      TotalFlightTime = TotalFlightTime(ThisGroup)+TotalFlightTime(Working.Pair(ThisGroup));
    end


    %convert flight time to delay in minutes
    OverallMedian = 500.*60
    TotalFlightTime = TotalFlightTime.*OverallMedian;
    if iDir == 3; TotalFlightTime = TotalFlightTime-OverallMedian.*2;
    else          TotalFlightTime = TotalFlightTime-OverallMedian;
    end
    TotalFlightTime = TotalFlightTime./60;


    for iIndex=1:1:numel(Indices)

      %select index
      Index = FlightIndices.(Indices{iIndex});

      %subset to just these flights
      Index = Index(ThisGroup);

      %pull out the relative flight times for the top and bottom x percent of each index
      Bottom = TotalFlightTime(find(Index <= prctile(Index,    xPercent)));
      Top    = TotalFlightTime(find(Index >= prctile(Index,100-xPercent)));

      %generate KDFs of the data
      x = -80:1:80;%
      a = pdf(fitdist(   Top,'kernel','Kernel','normal'),x); a = a./nansum(a(:));
      b = pdf(fitdist(Bottom,'kernel','Kernel','normal'),x); b = b./nansum(b(:));


      %test similarity of the two datasets
      [h,p,kstat]= kstest2(Top,Bottom,'Alpha',0.05);      


      %is the data statistically significant?
      if h == 1; 
        Alpha = 1; LineStyle = '-';
      else 
        Alpha = 0.5;  LineStyle = '--';
      end

      %plot mean top and bottom value for index
      TopVal = mean(Top,'omitnan');
      BottomVal = mean(Bottom,'omitnan');

%       %line
%       plot([0,mean(Top,'omitnan')],[0,mean(Bottom,'omitnan')],'-','color',Colours.(Indices{iIndex}),'linestyle',LineStyle)

      %value
      scatter(BottomVal,TopVal,60,'k', ...
              'markerfacecolor',Colours.(Indices{iIndex}), ...
              'marker','o','markerfacealpha',Alpha)
      drawnow

      xlabel('Index Min Delay'); ylabel('Index Max Delay')


    end; clear iIndex   TotalFlightTime


  end; clear iSeason
end; clear iDir