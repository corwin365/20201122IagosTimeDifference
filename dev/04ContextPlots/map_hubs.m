clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot maps of routes taken between each airport pair in each season
%and direction
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/12/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%file containing airports and flight times
Settings.DataFile = '../03CleanerFlights/flightpairs_identified_withpaths.mat';

%selection criteria - should be the same as those used in the main analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%minimum points for comparison
Settings.MinPoints = 0;

%outlier definition - flights this far off the median will be excluded
Settings.Outlier = [0.9,1.1]; %proportion of median time

%seasons - days of year in each
Settings.Seasons.DJF = date2doy(datenum(2000,12,1):datenum(2001, 3,1)-1);
Settings.Seasons.MAM = date2doy(datenum(2000, 3,1):datenum(2000, 6,1)-1);
Settings.Seasons.JJA = date2doy(datenum(2000, 6,1):datenum(2000, 9,1)-1);
Settings.Seasons.SON = date2doy(datenum(2000, 9,1):datenum(2000,12,1)-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data = load(Settings.DataFile);

%add day-of-year to the data
Data.Results.DoY = floor(date2doy(Data.Results.Date));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find each set of pairwise airports
% also normalise flight time for each pair of airports to their median
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 0;
for iDep = 1:1:numel(Data.Airports)
  for iArr=1:1:numel(Data.Airports)
   
    ThisPair = find(Data.Results.Dep == iDep ...
                  & Data.Results.Arr == iArr);

    if numel(ThisPair) < Settings.MinPoints;
      Data.Results.t(ThisPair) = NaN;
      continue;
    end
                
    MedianTime = nanmedian(Data.Results.t(ThisPair));
    Data.Results.Duration(ThisPair) =   Data.Results.t(ThisPair)./60./60; %flight time in hours
    Data.Results.t(ThisPair) = Data.Results.t(ThisPair)./MedianTime; %normalised flight time
              
    %is this flight eastbound or westbound?
    Dep = Data.Airports(iDep);
    Arr = Data.Airports(iArr);
    
    if     ismember(Dep,Data.Settings.NA)  && ismember(Arr,Data.Settings.Eur);
      EW = 1;
    elseif ismember(Dep,Data.Settings.Eur) && ismember(Arr,Data.Settings.NA);
      EW = 2;
    else disp('Error'); stop; end
    
    for iFlight =1:1:numel(ThisPair)
    
      k = k+1;
      %add results to big list
      Line = [iDep,iArr, ...
              Data.Results.Date(ThisPair(iFlight)),     ...
              Data.Results.t(   ThisPair(iFlight)),     ...
              EW,                                       ...
              Data.Results.DoY(ThisPair(iFlight)),      ...
              Data.Results.Duration(ThisPair(iFlight)), ...
              ThisPair(iFlight)];
      Results(k,1:8) = Line;
    end
  end
end

%outlier removal
Bad = find(Results(:,4) < Settings.Outlier(1).*nanmedian(Results(:,4)) ...
         | Results(:,4) > Settings.Outlier(2).*nanmedian(Results(:,4)));
Results(Bad,:) = NaN;
disp([num2str(numel(Bad)),' outlier time series removed'])
clear Bad

clearvars -except TimeScale Results Settings Data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot map of hubs we have used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%prepare plot
clf
% subplot(1,2,1)

m_proj('lambert','lat',[30,60],'lon',[-90,25]);
m_coast('patch','g');
hold on


%definition of trans-Atlantic flight:
%between one of these 
Settings.NA = {'ATL','BOS','BWI','CDW','CLE','CLT','CVG','DRM','DTW','EWR','FOK','IAD','JFK','LUK','MKE','MRB','ORD','PHL','PNE','YMX','YQB','YUL','YYR','YYZ','YZD'};
%and one of these
Settings.Eur = {'AGA','AGP','AHO','AMM','AMS','ATH','AYT','BCN','BEY','BOD','BRE','BRU','BTS','BUD','CAI','CDG','CGN','CIA','CRL','DBV','DLM','DME','DRS','DUS','ESB','FCO','FKB','FRA','GHF','GRO','HAJ','HAM','HEL','HER','HSK','IST','LCA','LEI','LEJ','LGW','LHR','LIS','LNZ','LYS','MAD','MAN','MLA','MRS','MUC','MXP','NCE','NUE','ORY','OST','OTP','PMI','PRG','PSA','PUY','RHO','RIX','RLG','SDV','SKG','SNN','SPM','STN','SXB','SZG','SZW','TLS','TLV','TOJ','TXL','UTC','VIE','ZNV','ZQL','ZRH'};
%in either direction

Codes = load('../01InitialTesting/airport_codes.mat');

%get list of valid airports, and plot them all as codes
Good = find(~isnan(Results(:,1)+Results(:,2)));
Good = unique([Results(Good,1),Results(Good,2)]);
for iAirport=1:1:numel(Good)
  
  
  %get the airport ID
  Airport = Data.Airports{Good(iAirport)};
  
  %is it in North America or Europe?
  if ismember(Airport,Settings.NA); Colour = 'r';
  else                              Colour = 'b';
  end

  %where is it?
  ID = find(contains(Codes.Codes,Airport));
  Coords = strsplit(Codes.Coords{ID},' ');
  
  %plot it!
  m_text(str2num(Coords{1}),str2num(Coords{2}),Airport, ...
         'horizontalalign','left','color',Colour,'fontsize',6)

  
end
         
m_grid
          
          

% % % % % 
% % % % % for iA = 1:1:numel(Data.Airports)
% % % % %   for iB=1:1:numel(Data.Airports)
% % % % %    
% % % % %     %get the list of flights
% % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %     
% % % % %     ThisPairDir1 = find(Results(:,1) == iA & Results(:,2) == iB);
% % % % %     ThisPairDir2 = find(Results(:,1) == iB & Results(:,2) == iA);
% % % % % 
% % % % %     if numel(ThisPairDir1) < Settings.MinPoints;
% % % % %       clear ThisPairDir1 ThisPairDir 2
% % % % %       continue;
% % % % %     end
% % % % %     if numel(ThisPairDir2) < Settings.MinPoints;
% % % % %       clear ThisPairDir1 ThisPairDir 2
% % % % %       continue;
% % % % %     end    
% % % % %     
% % % % %     %plot the routes for each season
% % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %     
% % % % %     figure
% % % % %     clf
% % % % %     set(gcf,'color','w')
% % % % %     sgtitle([Data.Airports{iA},'  --  ',Data.Airports{iB}])
% % % % %     
% % % % % 
% % % % %     for iDir=1:2;
% % % % %       
% % % % %       %which direction?
% % % % %       if iDir == Results(ThisPairDir1(1),5); ThisDir = ThisPairDir1;
% % % % %       else                                   ThisDir = ThisPairDir2;
% % % % %       end
% % % % %       
% % % % %       if iDir == 1; Dir = 'E'; else Dir = 'W'; end
% % % % %       
% % % % %       %pull out info on the flights
% % % % %       Flights = Results(ThisDir,:);
% % % % %       
% % % % %       %loop over seasons
% % % % %       Seasons = fieldnames(Settings.Seasons);
% % % % %       for iSeason=1:1:numel(Seasons)
% % % % % 
% % % % %         ThisSeason = Flights(find(ismember(Flights(:,6),Settings.Seasons.(Seasons{iSeason}))),:);
% % % % %         
% % % % %         %plot the flights!
% % % % %         subplot(2,numel(Seasons),iSeason+(numel(Seasons).*(iDir-1)))
% % % % %         m_proj('lambert','lat',[24,80],'lon',[-100,25]);
% % % % %         m_coast('patch','g');    
% % % % %           hold on        
% % % % % 
% % % % %         for iFlight=1:1:size(ThisSeason,1)
% % % % %           
% % % % %           ID = ThisSeason(iFlight,8);
% % % % %           
% % % % %           m_plot(Data.Results.Path.Lon{ID}, ...
% % % % %                  Data.Results.Path.Lat{ID}, ...
% % % % %                  '.-','color',[1,1,1].*0.7,'linewi',0.25);
% % % % % 
% % % % %           
% % % % %         end; clear iFlight
% % % % %        
% % % % %         m_grid('xtick',[],'ytick',[]);
% % % % % 
% % % % %         title([Seasons{iSeason},'    ',Dir,'ward'])
% % % % %         drawnow
% % % % %         
% % % % %       end; clear iSeason Seasons
% % % % %       
% % % % %       
% % % % %       
% % % % %     end; clear iDir
% % % % %     
% % % % %     
% % % % %     
% % % % %     
% % % % %     
% % % % %     %done!
% % % % %     
% % % % %     clear ThisPairDir1 ThisPairDir 2
% % % % %   end; clear iB
% % % % % end; clear iA