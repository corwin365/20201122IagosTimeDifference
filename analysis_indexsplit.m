function analysis_indexsplit(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%KDF analysis of merged aircraft data
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
if Settings.Choices.ApplyLags == 0; load([Settings.Paths.DataDir,'/',Settings.ID,'_indices.mat'])
else;                               load([Settings.Paths.DataDir,'/',Settings.ID,'_laggedindices.mat'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KDFStore = NaN(numel(Settings.Seasons.List),          ...
               numel(Settings.Choices.Directions),    ...
               numel(Settings.Choices.KDFSplit.Bins), ...
               numel(Settings.Indices.List),          ...
               2); %upper and lower region

NStore = KDFStore(:,:,1,1); %number of flights contributing

PStoreK = NStore; %statistical test result
PStoreT = NStore; %statistical test result

for iSeason=1:1:numel(Settings.Seasons.List)
  for iDirection=1:1:numel(Settings.Choices.Directions)

    %join indices onto flight data
    Data = innerjoin(FlightData,FlightIndices.(Settings.Choices.Directions{iDirection}));

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

      %generate KDFs
      % a = hist(t(Index < CutOff(1)),Settings.Choices.KDFSplit.Bins); a = a./nansum(a(:));
      % b = hist(t(Index > CutOff(1)),Settings.Choices.KDFSplit.Bins); b = b./nansum(b(:));
      a = pdf(fitdist(t(Index < CutOff(1)),'kernel','Kernel','normal'),Settings.Choices.KDFSplit.Bins); a = a./nansum(a(:));
      b = pdf(fitdist(t(Index > CutOff(2)),'kernel','Kernel','normal'),Settings.Choices.KDFSplit.Bins); b = b./nansum(b(:));
      
      %do statistical tests on the distributions
      [~,pk,~] = kstest2(t(Index < CutOff(1)),t(Index > CutOff(2)),'Alpha',Settings.Choices.KDFSplit.Alpha);
      [~,pt,~] = ttest2( t(Index < CutOff(1)),t(Index > CutOff(2)),'Alpha',Settings.Choices.KDFSplit.Alpha);

      %and store
      KDFStore(iSeason,iDirection,:,iIndex,1) = a;
      KDFStore(iSeason,iDirection,:,iIndex,2) = b;
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


Letters = 'abcdefghijklmnopqrstuvwxyz';

for iSeason = 1%1:1:numel(Settings.Seasons.List)
  % figure; 
  clf; set(gcf,'color','w'); 
  % sgtitle(Settings.Seasons.List{iSeason})
  % subplot = @(m,n,p) subtightplot (m, n, p, [0.1,0.01], 0.08, 0.2);
  subplot = @(m,n,p) subtightplot (m, n, p, [0.1,0.02], 0.08, 0.15);
  k = 0; %panel count
  for iDirection=1:1:numel(Settings.Choices.Directions)
    for iIndex=1:1:numel(Settings.Indices.List)

      %get plotting variables
      x = Settings.Choices.KDFSplit.Bins;
      Bottom = squeeze(KDFStore(iSeason,iDirection,:,iIndex,1));
      Top    = squeeze(KDFStore(iSeason,iDirection,:,iIndex,2));
      N      = squeeze(NStore(  iSeason,iDirection,  iIndex));
      PK     = squeeze(PStoreK( iSeason,iDirection,  iIndex));
      PT     = squeeze(PStoreT( iSeason,iDirection,  iIndex));

      %test similarity of the two datasets
      if PK < Settings.Choices.KDFSplit.Alpha ...
      && PT < Settings.Choices.KDFSplit.Alpha; h = 1; else h = 0; end       

      %find overlap region
      c = min([Top,Bottom],[],2);      

      %generate plot panel
      k = k+1;
      subplot(numel(Settings.Choices.Directions),numel(Settings.Indices.List),k)
      hold on
      % set(gca,'color',Settings.Indices.Colours.(Settings.Indices.List{iIndex}))
      for y=0:0.01:0.05; plot(minmax(x),[1,1].*y,'-','linewi',0.5,'color',[1,1,1].*0.9); end      

      %histograms
      patch([x,max(x),min(x)],[Top',0,0],    'r','facealpha',0.5,'edgecolor','none')
      patch([x,max(x),min(x)],[Bottom',0,0], 'b','facealpha',0.5,'edgecolor','none')

      %shade the common region dark grey if significant, and light grey if not
      if h == 1; patch([x,max(x),min(x)],[c',0,0],[1,1,1].*0.7,'facealpha',1,'edgecolor','none')
      else       patch([x,max(x),min(x)],[c',0,0],[1,1,1].*0.9,'facealpha',1,'edgecolor','none')
      end

      %overplot the distribution edges
      plot(x,Top,'r','linewi',0.5)
      plot(x,Bottom,'b','linewi',0.5)

      %print KS test results
      if PK < 1e-3; KString = '<0.001';
      else          KString  = ['=',num2str(round(PK,2,'significant'))];
      end
      if h == 1; Label = ['\it\bf{{PK',KString,'}}'];
      else       Label = [    '\it{PK',KString,'}' ];
      end
      text(max(x)-2,0.05,Label,'horizontalalignment','right','verticalalignment','top','fontsize',10);

      %print T-test test results   
      if PT < 1e-3; TString = '<0.001';
      else          TString  = ['=',num2str(round(PT,2,'significant'))];
      end      
      if h == 1; Label = ['\it\bf{{PT',TString,'}}'];
      else       Label = [    '\it{PT',TString,'}' ];
      end
      text(max(x)-2,0.044,Label,'horizontalalignment','right','fontsize',10);

      % %print number of samples
      % text(xpos,0.042,['\it{n=',num2str(N),'}'],'horizontalalignment','right','fontsize',12);      

      %tidy up panel
      axis([minmax(x) 0 0.05])
      grid off
      axis manual
      plot([1,1].*min(x),[-1,1].*999,'w-','linewi',3);plot([1,1].*max(x),[-1,1].*999,'w-','linewi',3)
      if     iIndex ==                           1;  set(gca,'ytick',0:0.01:0.05,'yticklabel',{' '}) 
      elseif iIndex == numel(Settings.Indices.List); set(gca,'ytick',0:0.01:0.05,'yticklabel',{'0','1%','2%','3%','4%','5%'},'yaxislocation','right')
      else                                           set(gca,'ytick',0:0.01:0.05,'yticklabel',{' '})
      end
      plot([0,0],[0,5].*1e-2,'-','color',[1,1,1].*0)

      %panel lettering
      text(min(x),0.05,['(',Letters(k),')'],'verticalalignment','top')




      %labelling
      if  iIndex == 1; 
        if     Settings.Choices.Directions{iDirection} == 'E'; DirName = 'Eastwards';
        elseif Settings.Choices.Directions{iDirection} == 'W'; DirName = 'Westwards';
        elseif Settings.Choices.Directions{iDirection} == 'R'; DirName = 'Round Trip';
        end
        ylabel(DirName,'fontsize',20);
      end

      if iDirection == 1;
        title(Settings.Indices.List{iIndex},'fontsize',30,'color',Settings.Indices.Colours.(Settings.Indices.List{iIndex}));
      end
      if iDirection == numel(Settings.Choices.Directions); xlabel('Delay [minutes]'); end      

      % % %indicate interquartile spread of values
      % % stop
      % % a = cumsum(Top); pc = prctile(a,[45,55]); Left = x(closest(a,pc(1))); Right = x(closest(a,pc(2))); 
      % % plot([Left,Right],[1,1].*-0.0005,'r','linewi',2,'clipping','off')
      % % a = cumsum(Bottom); pc = prctile(a,[45,55]); Left = x(closest(a,pc(1))); Right = x(closest(a,pc(2))); 
      % % plot([Left,Right],[1,1].*-0.0010,'b','linewi',2,'clipping','off')

      %indicate peak values
      [~,MaxTop] = max(Top);     MaxTop = x(MaxTop);
      [~,MaxBot] = max(Bottom);  MaxBot = x(MaxBot);
      plot(MaxTop,-0.0005,'^','color','k','markerfacecolor','r','markersize',10,'clipping','off')
      plot(MaxBot,-0.0005,'v','color','k','markerfacecolor','b','markersize',10,'clipping','off')



      drawnow
    end
  end
end

clear iSeason iDirection iIndex


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% output results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------------------------')
disp('Index-split KDF analysis complete')
disp('--------------------------')

end

