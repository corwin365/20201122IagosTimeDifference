function [] = dd_routeplot(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot metadata about the airports used and routes between them
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/12/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create new figure for each season
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get list of seasons
Seasons = fieldnames(Settings.Seasons);

for iSeason=1:1:numel(Seasons)
  
  
  figure
  clf
  set(gcf,'color','w')
  subplot = @(m,n,p) subtightplot (m, n, p, 0.08, 0.05, 0.05);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% load route and flight data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  FlightData = load('data/flight_data.mat');
  RouteData  = load('data/routes.mat');
  
  NPorts = numel(RouteData.Airports);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% generate a table showing how many flights we have
  % between each airport
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %prepare axes
  subplot(1,2,2)
  axis([0 NPorts 0 NPorts])
  axis square
  hold on
  set(gca,'xtick',[],'ytick',[])
  
  
  %plot lines and label airports
  for iX=1:1:NPorts
    plot([1,1].*iX,[0 numel(RouteData.Airports)],'k-')
    if ismember(RouteData.Airports{NPorts-iX+1},Settings.NA); Colour = 'b';
    else                                                      Colour = 'r';
    end
    text(-0.25,iX-0.5,['--> ',RouteData.Airports{NPorts-iX+1}],'horizontalalignment','right','color',Colour)
  end
  for iY=1:1:NPorts
    plot([0 numel(RouteData.Airports)],[1,1].*iY,'k-')
    if ismember(RouteData.Airports{iY},Settings.NA); Colour = 'b';
    else                                                      Colour = 'r';
    end
    text(iY-0.5,NPorts+0.2,RouteData.Airports{iY},'horizontalalignment','center','color',Colour)
  end
  
  %print number of flights in pair
  for iY=1:1:NPorts
    for iX=1:1:NPorts
      
      %how many flights?
      N = RouteData.NFlights(iSeason,iX,NPorts-iY+1);
      
      %are we travelling east or west?
      Dep = RouteData.Airports{iX};
      if ismember(Dep,Settings.NA); Colour = 'b';
      else                          Colour = 'r';
      end
      
      
      if iX == NPorts-iY+1; Value = '---'; Colour = 'k';
      elseif N == 0;        Value = ' ';
      else                  Value = num2str(N);
      end
      
      text(iX-0.5,iY-0.5,Value,'horizontalalignment','center','color',Colour,'fontweight','bold','fontsize',11)
      
    end
  end
  
  %key
  text(       0.5,-0.25,'from North America (Eastward)','color','b','horizontalalignment','left')
  text(NPorts-0.5,-0.25,'from Europe (Westward)','color','r','horizontalalignment','right')
  text(NPorts./2,NPorts+.6,Seasons{iSeason},'fontsize',24,'fontweight','bold','horizontalalignment','center')

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% generate maps showing the airports
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  for iSide=1:2;
    
    if iSide == 1;
      subplot(2,2,1)
      m_proj('robinson','lat',[30,70],'lon',[-20,25]);
      %     set(gca,'yaxislocation','right')
      title('Europe')
    else
      subplot(2,2,3)
      m_proj('robinson','lat',[25,65],'lon',[-100,-55]);
      title('North America')
    end
    
    m_coast('patch',[1,1,1].*0.9,'edgecolor','none');
    hold on
    
    for iAirport=1:1:NPorts
      
      Port = RouteData.Airports{iAirport};
      
      %colour?
      if ismember(Port,Settings.NA); Colour = 'b';
      else                           Colour = 'r';
      end
      
      if iSide == 1 & strcmp(Colour,'b'); continue; end
      if iSide == 2 & strcmp(Colour,'r'); continue; end
      
      %location?
      ID = find(contains(FlightData.Airports.IDs,Port));
      
      m_text(FlightData.Airports.Lons(ID), ...
        FlightData.Airports.Lats(ID), ...
        Port, ...
        'color',Colour,'horizontalalignment','center','fontsize',8,'fontweight','bold')
      
      
    end
    m_grid('xtick',[-185:10:180],'ytick',[-90:10:90])
    
  end
  drawnow
end
  
