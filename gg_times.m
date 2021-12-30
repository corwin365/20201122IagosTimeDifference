function [] = gg_times(Paths,Indices,Seasons,RTRange)

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

%identify direction of each individual flight
Working.Eward = NaN.*Working.Used;
for iRoute=1:1:size(RouteInfo,1)
  OnThisRoute = find(Working.Route == RouteInfo{iRoute,1});
  Working.Eward(OnThisRoute) = RouteInfo{iRoute,5};
end; clear iRoute OnThisRoute

%load indices
Index = load([Paths.StoreDir,'/indices_',Paths.SourceIdentifier,'_',Paths.PPIdentifier,'.mat']);
FlightIndices = Index.Flight;
clear Index


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% split out the top and bottom x % by each index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xPercent = 20;



for iDir=3%1:1:2 %1 is eastward, 2 is westward, 3 is round-trip


  %prepare figure
  figure(iDir)
  set(gcf,'color','w')
  k = 0;
 

  SeasonNames = fieldnames(Seasons);
  
  for iSeason=1:1:numel(SeasonNames)
    for iIndex=1:1:numel(Indices)

      %prepare data
      Index = FlightIndices.(Indices{iIndex});
      Total = Working.tRel; %this will be overwritten for the paired trips
      if iDir == 1; %eastward
        Index(find(Working.(SeasonNames{iSeason}).Used ~= 1 | Working.Eward ~= 2)) = NaN;
      elseif iDir == 2; %westward
        Index(find(Working.(SeasonNames{iSeason}).Used ~= 1 | Working.Eward ~= 1)) = NaN;
      else

        %round trip. More complicated!
        %first, find all westbound flights in the season which have an eastbound pair
        Paired = find(Working.(SeasonNames{iSeason}).Used == 1 & ~isnan(Working.Pair));
        
        %now, add the relative flight times of the outbounds to the inbounds
        Total = NaN(numel(Working.Used),1);
        for iOut=1:1:numel(Paired);
          Total(Paired(iOut)) = Working.tRel(Paired(iOut)) + Working.tRel(Working.Pair(Paired(iOut)));
        end


        
        %now, pull out the index array, then we're done
        i2 = Index.*NaN;
        i2(Paired) = Index(Paired);
        Index = i2;
        clear Paired iOut i2

      end

      %pull out the relative flight times for the top and bottom x percent of each index
      Top    = Total(find(Index <= prctile(Index,    xPercent)));
      Bottom = Total(find(Index >= prctile(Index,100-xPercent)));

      %generate KDFs of the data

      if    iDir == 3; x = 2.*min(RTRange):0.005:max(RTRange).*2;
      else             x =    min(RTRange):0.005:max(RTRange);
      end
      a = pdf(fitdist(   Top,'kernel','Kernel','normal'),x); a = a./nansum(a(:));
      b = pdf(fitdist(Bottom,'kernel','Kernel','normal'),x); b = b./nansum(b(:));

% %       %convert scale from relative times to minutes of delay
% %       if iDir == 3;  x = x.*(2.*OverallMedian) - 2.*OverallMedian;
% %       else           x = x.*OverallMedian - OverallMedian;
% %       end
% %       x = x./60;

      %test similarity of the two datasets
      h= kstest2(Top,Bottom,'Alpha',0.05);      

      %find overlap region
      c = min([a;b],[],1);

      %plot
      k = k+1;
      subplot(numel(SeasonNames),numel(Indices),k)
      cla
      hold on
      for y=0:0.02:0.08; plot(minmax(x),[1,1].*y,'-','linewi',0.5,'color',[1,1,1].*0.7); end

      %histograms
      patch([x,max(x),min(x)],[a,0,0],         'r','facealpha',0.5,'edgecolor','none')
      patch([x,max(x),min(x)],[b,0,0],         'b','facealpha',0.5,'edgecolor','none')
      
      %shade the common region dark grey if significant, and light grey if not
      if h == 1; patch([x,max(x),min(x)],[c,0,0],[1,1,1].*0.7,'facealpha',1,'edgecolor','none')
      else       patch([x,max(x),min(x)],[c,0,0],[1,1,1].*0.9,'facealpha',1,'edgecolor','none')
      end

      %central values
      [~,idx] = max(a); scatter(x(idx),0,40,[1,1,1].*0.6,'markerfacecolor','r','markerfacealpha',0.6)
      [~,idx] = max(b); scatter(x(idx),0,40,[1,1,1].*0.6,'markerfacecolor','b','markerfacealpha',0.6)


      %print test results
      if h == 1; text(1,0.01,'\Delta','horizontalalignment','center');
      else       text(1,0.01,'n/s','horizontalalignment','center');
      end

      axis([minmax(x) 0 0.08])
      grid off
      axis manual
      plot([1,1].*min(x),[-1,1].*999,'w-','linewi',3)
      set(gca,'ytick',0:0.02:0.08,'yticklabel',{' ','2%','4%','6%','8%'})


      if  iIndex == 1; ylabel(SeasonNames{iSeason}); end
      if iSeason == 1; title(Indices{iIndex}); end


      drawnow
      clear Top Bottom x a b h c y idx Index Total


    end; clear iIndex


  end; clear iSeason
end; clear iDir