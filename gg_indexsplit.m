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



for iDir=1:1:3 %1 is eastward, 2 is westward, 3 is round-trip


  %prepare figure
  figure(iDir)
  clf
  set(gcf,'color','w')
  k = 0;
 


  
  for iSeason=1:1:numel(SeasonNames)

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
      Top    = TotalFlightTime(find(Index <= prctile(Index,    xPercent)));
      Bottom = TotalFlightTime(find(Index >= prctile(Index,100-xPercent)));


      %generate KDFs of the data
      x = -80:1:80;%
% % 
% %       if    iDir == 3; x = 2.*min(RTRange):0.005:max(RTRange).*2;
% %       else             x =    min(RTRange):0.005:max(RTRange);
% %       end
      a = pdf(fitdist(   Top,'kernel','Kernel','normal'),x); a = a./nansum(a(:));
      b = pdf(fitdist(Bottom,'kernel','Kernel','normal'),x); b = b./nansum(b(:));


% %       %convert scale from relative times to minutes of delay
% %       if iDir == 3;  x = x.*(2.*OverallMedian) - 2.*OverallMedian;
% %       else           x = x.*OverallMedian - OverallMedian;
% %       end
% %       x = x./60;

      %test similarity of the two datasets
      [h,p,kstat]= kstest2(Top,Bottom,'Alpha',0.05);      

      %find overlap region
      c = min([a;b],[],1);

      %plot
      k = k+1;
      subplot(numel(SeasonNames),numel(Indices),k)
      cla
      hold on
      for y=0:0.01:0.04; plot(minmax(x),[1,1].*y,'-','linewi',0.5,'color',[1,1,1].*0.7); end

      %histograms
      patch([x,max(x),min(x)],[a,0,0],         'r','facealpha',0.5,'edgecolor','none')
      patch([x,max(x),min(x)],[b,0,0],         'b','facealpha',0.5,'edgecolor','none')
      
      %shade the common region dark grey if significant, and light grey if not
      if h == 1; patch([x,max(x),min(x)],[c,0,0],[1,1,1].*0.7,'facealpha',1,'edgecolor','none')
      else       patch([x,max(x),min(x)],[c,0,0],[1,1,1].*0.9,'facealpha',1,'edgecolor','none')
      end

      %mean values
      scatter(mean(Top),0,60,[1,1,1].*0.6,'markerfacecolor','r','markerfacealpha',0.6,'marker','^')
      scatter(mean(Bottom),0,60,[1,1,1].*0.6,'markerfacecolor','b','markerfacealpha',0.6,'marker','v')


      %print test results
      xpos = 55;
      if h == 1; Label = ['\it\bf{{p=',num2str(round(p,2,'significant')),'}}'];
      else       Label = [    '\it{p=',num2str(round(p,2,'significant')),'}' ];%\it{n/s}';
      end
      text(xpos,0.037,Label,                             'horizontalalignment','right','fontsize',10);
      text(xpos,0.032,['\it{n=',num2str(numel(Top)),'}'],'horizontalalignment','right','fontsize',10);
      

      %tidy up panel
      axis([minmax(x) 0 0.04])
      grid off
      axis manual
      plot([1,1].*min(x),[-1,1].*999,'w-','linewi',3)
      set(gca,'ytick',0:0.01:0.04,'yticklabel',{' ','1%','2%','3%','4%'})
      plot([0,0],[0,4].*1e-2,'-','color',[1,1,1].*0)

      %add a swatch of the index's signature colour
      a = -55; b = -42; c= 0.032; d =0.038;
      patch([a,b,b,a,a],[c,c,d,d,c],Colours.(Indices{iIndex}),'edgecolor','k','clipping','off')

      if  iIndex == 1; ylabel(SeasonNames{iSeason}); end
      if iSeason == 1; title(Indices{iIndex},'fontsize',17); end


      drawnow
      clear Top Bottom x a b h c y idx Index Label kstat p


    end; clear iIndex   TotalFlightTime


  end; clear iSeason
end; clear iDir