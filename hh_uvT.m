function [] = hh_uvT(Paths,Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot time series of time taken over the instrument record
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/02/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
SeasonNames = fieldnames(Settings.Seasons);
SeasonNames{end+1} = 'All';


disp('*****Properties being manually set in child routine - move to parent at end*****')
Properties = {'U','V','T'};
NBins = 50;
x.U = [-50,150];
x.V= [-80,80];
x.T = [190,260];
% Percentiles = [0,2.5,18,50,82,97.5,100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% produce a histogram of each property for each flight
%also store some key percentiles for each flight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create storage
HistStore = NaN(numel(Properties),numel(Flights.Date),NBins);
% PercStore = NaN(numel(Properties),numel(Flights.Date),numel(Percentiles)); 


%loop over properties
for iProp=1:1:numel(Properties)
  ThisProp = Flights.Paths.(Properties{iProp});

  %loop over individual flights
  for iFlight=1:1:numel(Flights.Date)

    %generate and bin the histogram
    PP = ThisProp{iFlight};
    if sum(PP,'all','omitnan') == 0; continue; end %some flights don't contain any of these data

    xx = linspace(min(x.(Properties{iProp})), ...
                  max(x.(Properties{iProp})), ...
                  NBins);

    HistStore(iProp,iFlight,:) = hist(PP,xx);
%     PercStore(iProp,iFlight,:) = prctile(PP,Percentiles);

  end; clear iFlight
end; clear iProp
clear PP xx

%end bins can be a dumping ground - remove
HistStore(:,:,[1,end]) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% divide data by index and show distribution for top and bottom fractions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%   figure
figure
clf
set(gcf,'color','w','position',[3.8 234.2 1531.6 508])
k=0;

  %loop over properties
for iProp=1:1:numel(Properties)
  %loop over directions

  for iDir=3%1:1:3; %1 is eastward, 2 is westward, 3 is round-trip
    %loop over seasons
    for iSeason=numel(SeasonNames);

      %find all flights in this season and direction
      if iDir == 1;
        ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1 ...
                       & Working.Eastward                        == true);
      elseif iDir == 2;
        ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1 ...
                       & Working.Eastward                        == false);
      elseif iDir == 3;
        ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1 ...
                       & ~isnan(Working.Pair));
      end

      %loop over indices
      for iIndex=1:1:numel(Settings.Indices)

        %select index
        Index = FlightIndices.(Settings.Indices{iIndex});

        %subset to just these flights
        Index = Index(ThisGroup);

        %pull out the relative flight times for the top and bottom x percent of each index
        Bottom = find(Index <= prctile(Index,    Settings.CutOff));
        Top    = find(Index >= prctile(Index,100-Settings.CutOff));


        %compress the distribution for these sets of flights
        Distrib.Bottom = squeeze(nansum(HistStore(iProp,Bottom,:),2));
        Distrib.Bottom = Distrib.Bottom./sum(Distrib.Bottom,'omitnan');
        
        Distrib.Top    = squeeze(nansum(HistStore(iProp,Top   ,:),2));
        Distrib.Top    = Distrib.Top   ./sum(Distrib.Top   ,'omitnan');


        Distrib.Top = smoothn(Distrib.Top,5);
        Distrib.Bottom = smoothn(Distrib.Bottom,5);


        %get values of bins
        xx = linspace(min(x.(Properties{iProp})), ...
                      max(x.(Properties{iProp})), ...
                      NBins);

        k = k+1;
        subplot(numel(Properties),numel(Settings.Indices),k)
        cla
        hold on




        %plot
        Delta = Distrib.Top-Distrib.Bottom;
        Delta = Delta.*100; %percent
        plot(xx,Delta,'-','color',[1,1,1].*0.3)
        plot(xx,zeros(size(xx)),'k-')
        hold on

        A = Delta; A(find(Delta > 0)) = 0; A = A';
        B = Delta; B(find(Delta < 0)) = 0; B = B';


        patch([xx,xx(end),min(xx)],[A,0,0],'b','facealpha',0.5,'edgecolor','none')
        patch([xx,xx(end),min(xx)],[B,0,0],'r','facealpha',0.5,'edgecolor','none')       
%         patch([xx,0,0],[A,0,0],'r')


% % plot(xx,Distrib.Top./Distrib.Bottom,'r');

% %         plot(xx,Distrib.Top,'r');
% %         hold on
% %         plot(xx,Distrib.Bottom,'b')

       xlim(minmax(xx))
%        ylim([-1,1].*1)


        title(Settings.Indices{iIndex})
        ylabel(Properties{iProp})

        drawnow
      end; clear iIndex

    end; clear iSeason


  end; clear iIndex
end; clear iProp

