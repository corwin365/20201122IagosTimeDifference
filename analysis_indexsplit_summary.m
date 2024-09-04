function analysis_indexsplit_summary(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%KDF analysis of merged aircraft data - summary results
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/01/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Index-split KDF analysis')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat'])

%load indices
load([Settings.Paths.DataDir,'/',Settings.ID,'_indices.mat'])

%join indices onto flight data, and drop/re
Data = innerjoin(FlightData,FlightIndices);

clear DateIndices FlightData FlightIndices RangeStore RouteData

Seasons= fieldnames(Settings.Seasons);

Percentiles =  Settings.Choices.KDFPc.Percentiles; %move this to settings file when script complete

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KDFStore = NaN(numel(Settings.Seasons.List),          ...
               numel(Settings.Choices.Directions),    ...
               numel(Settings.Choices.KDFSplit.Bins), ...
               numel(Settings.Indices.List),          ...
               2); %upper and lower region

PCStore = NaN(numel(Settings.Seasons.List),          ...
               numel(Settings.Choices.Directions),   ...
               numel(Percentiles),                   ...
               numel(Settings.Indices.List),         ...
               2); %upper and lower region

NStore = KDFStore(:,:,1,1); %number of flights contributing

PStoreK = NStore; %statistical test result
PStoreT = NStore; %statistical test result

for iSeason=1:1:numel(Settings.Seasons.List)
  for iDirection=1:1:numel(Settings.Choices.Directions)

    %find flights in this season and direction
    InThisSeason    = find(table2array(Data.InSeasons(:,iSeason)) == 1);
    InThisDirection = find(Data.Direction == Settings.Choices.Directions{iDirection});
    InThisSet = intersect(InThisDirection,InThisSeason);

    %hence, cut down the master table to just these entries
    FlightsInUse = Data(InThisSet,:);


    %get flight delays and convert to minutes (used for this plot)
    t = table2array(FlightsInUse.Delay)./60;


    for iIndex=1:1:numel(Settings.Indices.List)
   
      %get index
      Index = FlightsInUse.(Settings.Indices.List{iIndex});
      
      %find top and bottom X% by this index
      CutOff = prctile(Index,[Settings.Choices.KDFSplit.CutOff,100-Settings.Choices.KDFSplit.CutOff]);

      %generate KDFs - these are only used for the statistical tests
      stop
      a = pdf(fitdist(t(Index < CutOff(1)),'kernel','Kernel','normal'),Settings.Choices.KDFSplit.Bins); a = a./nansum(a(:));
      b = pdf(fitdist(t(Index > CutOff(2)),'kernel','Kernel','normal'),Settings.Choices.KDFSplit.Bins); b = b./nansum(b(:));
      
      %do statistical tests on the distributions
      [~,pk,~] = kstest2(t(Index < CutOff(1)),t(Index > CutOff(2)),'Alpha',Settings.Choices.KDFSplit.Alpha);
      [~,pt,~] = ttest2( t(Index < CutOff(1)),t(Index > CutOff(2)),'Alpha',Settings.Choices.KDFSplit.Alpha);

      %find interesting percentiles
      PCb = prctile(t(Index < CutOff(1)),Percentiles);
      PCa = prctile(t(Index > CutOff(2)),Percentiles);

      %and store
      KDFStore(iSeason,iDirection,:,iIndex,1) = a;
      KDFStore(iSeason,iDirection,:,iIndex,2) = b;
      PCStore( iSeason,iDirection,:,iIndex,1) = PCa;
      PCStore( iSeason,iDirection,:,iIndex,2) = PCb;
      NStore(  iSeason,iDirection,  iIndex)   = numel(Index);
      PStoreK( iSeason,iDirection,  iIndex)   = pk;
      PStoreT( iSeason,iDirection,  iIndex)   = pt;

    end
  end
end

clear iSeason iDirection InThisSeason InThisDirection InThisSet FlightsInUse t iIndex Index CutOff a b pt pk


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
clf
for iDirection=1:1:numel(Settings.Choices.Directions)
  
  %create panel
  subplot(1,numel(Settings.Choices.Directions), iDirection)
  xlim([minmax(Settings.Choices.KDFSplit.Bins)])
  title(Settings.Choices.Directions{iDirection})
  hold on

  for iIndex=1:1:numel(Settings.Indices.List)


    %get plotting variables
    x = Settings.Choices.KDFSplit.Bins;
    Bottom = squeeze(PCStore(:,iDirection,:,iIndex,1));
    Top    = squeeze(PCStore(:,iDirection,:,iIndex,2));
    N      = squeeze(NStore(  :,iDirection,  iIndex));
    PK     = squeeze(PStoreK( :,iDirection,  iIndex));
    PT     = squeeze(PStoreT( :,iDirection,  iIndex));

    %patch in the colour of this dataset, mixed with white
    shift  =  size(Bottom,1).*iIndex;    
    Colour = colorGradient(Settings.Indices.Colours.(Settings.Indices.List{iIndex}),[1,1,1],5);
    N = size(Bottom,1);
    patch([-1,1,1,-1,-1].*max(Settings.Choices.KDFSplit.Bins),shift+[0,0,N,N,0]+0.5, ...
          Colour(4,:),'edgecolor','k')
 
    patch([-1,1,1,-1,-1].*max(Settings.Choices.KDFSplit.Bins),shift+[0,0,1,1,0]+0.5, ...
          Colour(3,:),'edgecolor','none')
   
    for iSeason=1:1:size(Bottom,1)
      width = ceil(numel(Percentiles)./2)+4;
      
      for iPair=1:1:ceil(numel(Percentiles)./2)

        if iPair == 1; idx = ceil(numel(Percentiles)./2);                   MarkerSize = 20;   
        else           idx = ceil(numel(Percentiles)./2)+[-1,1].*(iPair-1); MarkerSize = 0.001;
        end

        plot(Bottom(iSeason,idx),ones(numel(idx),1).*iSeason+0.1+shift,'.-','color','b','linewi',(width-2*iPair),'markersize',MarkerSize)
        plot(Top(   iSeason,idx),ones(numel(idx),1).*iSeason-0.1+shift,'.-','color','r','linewi',(width-2*iPair),'markersize',MarkerSize)
        plot(minmax(Settings.Choices.KDFSplit.Bins),[1,1].*iSeason+0.5+shift,'k:')

        %add a star to the line if the difference is statistically significant
        if PK(iSeason) < Settings.Choices.KDFSplit.Alpha %...
        % && PT(iSeason) < Settings.Choices.KDFSplit.Alpha;
          text(max(Settings.Choices.KDFSplit.Bins)-1,iSeason+shift+0.15,'*','verticalalignment','middle', ...
               'horizontalalignment','right','fontweight','bold','fontsize',20)
        end

      end; clear iPair
      text(min(Settings.Choices.KDFSplit.Bins)+1,iSeason+shift,Seasons{iSeason}, ...
          'verticalalignment','middle','horizontalalignment','left')
    end; clear iSeason


  end; clear iIndex
  ylim([size(Bottom,1)+0.5 shift+size(Bottom,1)+0.5])
  set(gca,'ytick',[])
  box on
  set(gca,'ydir','reverse','layer','top')
  xlabel('Delay [minutes]')

  for iX=-1000:20:1000;plot([1,1].*iX,[-1,1].*999,':','color',[1,1,1].*0.5); end
  plot([0,0],[-1,1].*999,'-','color',[1,1,1].*0.5)
  drawnow

end; clear iDirection


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% output results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------------------------')
disp('Index-split KDF summary complete')
disp('--------------------------')

end

