function [] = clim_bb(Paths,Settings)

Seasons  = Settings.Seasons;
xPercent = Settings.CutOff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%split histograms of different index segments, for climate subpaper
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
figure()
clf
set(gcf,'color','w','position',[97 226 1493 646])
k = 0;



for iDir=[2,3,1] %1 is eastward, 2 is westward, 3 is round-trip




  
  for iSeason=1:1:numel(SeasonNames)

    %select flights in this season and direction
    %also load their flight times


    if iDir == 1;
      ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1 ...
                     & Working.Eastward                        == true);
      TotalFlightTime = Working.tRel.(SeasonNames{iSeason}); 
      TotalFlightTime = TotalFlightTime(ThisGroup);
      OverallMedian = OverallMedians.(SeasonNames{iSeason}).East;
    elseif iDir == 2;
      ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1 ...
                     & Working.Eastward                        == false);
      TotalFlightTime = Working.tRel.(SeasonNames{iSeason}); 
      TotalFlightTime = TotalFlightTime(ThisGroup);
      OverallMedian = OverallMedians.(SeasonNames{iSeason}).West;
    elseif iDir == 3;
      ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1 ...
                     & ~isnan(Working.Pair));
      TotalFlightTime = Working.tRel.(SeasonNames{iSeason}); 
      TotalFlightTime = TotalFlightTime(ThisGroup)+TotalFlightTime(Working.Pair(ThisGroup));
      OverallMedian = (OverallMedians.(SeasonNames{iSeason}).East+OverallMedians.(SeasonNames{iSeason}).West)/2;
    end


    %convert flight time to delay in minutes
    TotalFlightTime = TotalFlightTime.*OverallMedian;
    if iDir == 3; TotalFlightTime = TotalFlightTime-OverallMedian.*2;
    else          TotalFlightTime = TotalFlightTime-OverallMedian;
    end
    TotalFlightTime = TotalFlightTime./60;



      %select index
      Index = FlightIndices.HadCRUT;

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

      %find overlap region
      c = min([a;b],[],1);

      %plot
      k = k+1;
      subplot(3,numel(SeasonNames),k)
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

      %print test results
      xpos = 77;
      if h == 1; Label = ['\it\bf{{p=',num2str(round(p,2,'significant')),'}}'];
      else       Label = [    '\it{p=',num2str(round(p,2,'significant')),'}' ];%\it{n/s}';
      end
      text(xpos,0.037,Label,                             'horizontalalignment','right','fontsize',14);
      text(xpos,0.032,['\it{n=',num2str(numel(Top)),'}'],'horizontalalignment','right','fontsize',14);
      
      %overplot the distribution edges
      plot(x,a,'r','linewi',0.5)
      plot(x,b,'b','linewi',0.5)

      %tidy up panel
      axis([minmax(x) 0 0.04])
      grid off
      axis manual
      plot([1,1].*min(x),[-1,1].*999,'w-','linewi',3)

      if iSeason == 1 | iSeason == numel(SeasonNames);
        set(gca,'ytick',0:0.01:0.04,'yticklabel',{' ','1%','2%','3%','4%'})
      else
        set(gca,'ytick',0:0.01:0.04,'yticklabel',{' '})
      end
      if iSeason == numel(SeasonNames); 
        set(gca,'yaxislocation','right'); 
        plot([1,1].*max(x),[-1,1].*999,'w-','linewi',3);
      end
      plot([0,0],[0,4].*1e-2,'-','color',[1,1,1].*0)

      if iDir == 2; title(SeasonNames{iSeason},'fontsize',20); end
      if iDir == 1; xlabel('Delay [minutes]'); end
% 

      drawnow
      clear Top Bottom x a b h c y idx Index Label kstat p


clear   TotalFlightTime


  end; clear iSeason
 
end; clear iDir


