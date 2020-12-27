clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyse IAGOS data for effect of ENSO, NAO, QBO and HadCRUT on
%trans-atlantic flight times
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/11/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file containing airports and flight times
Settings.DataFile = 'flightpairs.mat';

%indices to use
Settings.Indices = {'QBO','ENSO','HadCRUT','NAM','TSI','NAO'};

%minimum points for comparison
Settings.MinPoints = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data = load(Settings.DataFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% results array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Used = sum(~isnan(Data.Results.t));
Results = NaN(Used,5+numel(Settings.Indices));
clear Used



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
    Data.Results.t(ThisPair) = Data.Results.t(ThisPair)./MedianTime;
              
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
              Data.Results.Date(ThisPair(iFlight)), ...
              Data.Results.t(   ThisPair(iFlight)), ...
              EW];
      Results(k,1:5) = Line;
    end
  end
end

%there is a single monster outlier. remove it
Bad = find(Results(:,4) > 3);
Results(Bad,:) = NaN;

clearvars -except TimeScale Results Settings Data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load indices and interpolate to data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the indices, and interpolate to each data point
TimeScale = Results(:,3);

for iIndex=1:1:numel(Settings.Indices)
  
  switch Settings.Indices{iIndex}
    case 'QBO'
      QBO = load([LocalDataDir,'/Miscellany/QBO.mat']);
      QBO.Time = floor(QBO.Time); %shift from noon to midnight to make the logic easier - on a 91-day smoothing this is very minor...
      a = interp1(QBO.Time,QBO.QBO,TimeScale);
      clear QBO
    case 'ENSO'
      ENSO = load([LocalDataDir,'/Miscellany/nino34.mat']);
      a = interp1(ENSO.Time,ENSO.Nino34,TimeScale);
      clear ENSO
    case 'HadCRUT'
      HadCRUT = rCDF([LocalDataDir,'/Miscellany/HadCRUT.4.6.0.0.median.nc']);
      HadCRUT.MatlabTime = datenum(1850,1,HadCRUT.time);
      HadCRUT.NH = squeeze(nanmean(HadCRUT.temperature_anomaly(:,HadCRUT.latitude > 0,:),[1,2]));
      a = interp1(HadCRUT.MatlabTime,HadCRUT.NH,TimeScale);
      clear HadCRUT
    case 'NAM'
      NAM = load([LocalDataDir,'/Miscellany/daily_nam.mat']);
      a = interp1(NAM.Time,NAM.NAM,TimeScale);
      clear NAM
    case 'NAO'
      NAO = load([LocalDataDir,'/Miscellany/nao.mat']);
      a = interp1(NAO.Date,NAO.NAO,TimeScale);
      clear NAO      
    case 'TSI'
      TSI = load([LocalDataDir,'/Miscellany/tsi.mat']);
      a = interp1(TSI.Time,TSI.TSI,TimeScale);
      clear TSI
  end
  Results(:,5+iIndex) = a;
  clear a
end; clear iIndex

% % % % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % %% finally, plot results
% % % % % % % % % % % % %  do one plot combining all (normalised) routes reaching threshold Npoints
% % % % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % 
% % % % % % % % % % % % clf
% % % % % % % % % % % % set(gcf,'color','w')
% % % % % % % % % % % % subplot = @(m,n,p) subtightplot (m, n, p, 0.05, 0.05, 0.05);
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % k = 0;
% % % % % % % % % % % % for EW=1:2;
% % % % % % % % % % % %   
% % % % % % % % % % % %   %find this direction
% % % % % % % % % % % %   ThisDir = find(Results(:,5) == EW);
% % % % % % % % % % % %   if EW == 1; Dir = 'E'; else; Dir = 'W'; end
% % % % % % % % % % % %   
% % % % % % % % % % % %   for iIndex=1:1:numel(Settings.Indices)
% % % % % % % % % % % %     k = k+1;
% % % % % % % % % % % %     subplot(2,numel(Settings.Indices),k)
% % % % % % % % % % % %   
% % % % % % % % % % % %     %plot data 
% % % % % % % % % % % %     x = Results(ThisDir,4);
% % % % % % % % % % % %     y = Results(ThisDir,5+iIndex);
% % % % % % % % % % % %     
% % % % % % % % % % % %     plot(x,y,'ko')
% % % % % % % % % % % %     hold on
% % % % % % % % % % % %     
% % % % % % % % % % % %     %linear fit
% % % % % % % % % % % %     Good = find(~isnan(x+y));
% % % % % % % % % % % %     if numel(Good) > 10;
% % % % % % % % % % % %       p = polyfit(x(Good),y(Good),1);
% % % % % % % % % % % %       y2 = polyval(p,x);
% % % % % % % % % % % %       plot(x,y2,'r-')
% % % % % % % % % % % %     end
% % % % % % % % % % % %     
% % % % % % % % % % % %     %correlation
% % % % % % % % % % % %     r = corrcoef(x(Good),y(Good));
% % % % % % % % % % % %     r = round(r(2).*100)./100;
% % % % % % % % % % % %     
% % % % % % % % % % % %     title([Settings.Indices{iIndex},', ',Dir,', r=',num2str(r)])
% % % % % % % % % % % %     
% % % % % % % % % % % %     axis square
% % % % % % % % % % % %     
% % % % % % % % % % % %   end
% % % % % % % % % % % % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finally, plot histograms
%  do one plot combining all (normalised) routes reaching threshold Npoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, 0.02, 0.05, 0.05);



k = 0;
for EW=1:2;
  
  %find this direction
  ThisDir = find(Results(:,5) == EW);
  if EW == 1; Dir = 'flights W->E'; else; Dir = 'flight E->W'; end
  
  for iIndex=1:1:numel(Settings.Indices)
    k = k+1;
    subplot(2,numel(Settings.Indices),k)
  
    %extract data 
    x = Results(ThisDir,4);
    y = Results(ThisDir,5+iIndex);
    
    %split into bottom third and top third, then generate and plot histograms
    for Range=[1,2];
      
      switch Range
        case 1; idx = find(y < prctile(y,25)); Colour = 'r';
        case 2; idx = find(y > prctile(y,75)); Colour = 'b';
      end
      
      [ys,xs] = hist(x(idx),linspace(0.9,1.1,100));
      
      ys(1) = 0;
      ys(end) = 0;
      ys = smoothn(ys,5);
      
      
      patch([0.79,xs,1.21],[0,ys,0],Colour,'facealpha',0.3,'edgecolor',Colour)
      hold on

      title([Settings.Indices{iIndex},', ',Dir])
      xlim([0.9 1.1])
      xlabel('Relative travel time')
      
    end
    axis square
    
  end
end



% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %% finally, plot results
% % % % % %  do a separate plot for each index and airport pair
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % 
% % % % % Pairs = zeros(numel(Data.Airports),numel(Data.Airports));
% % % % % 
% % % % % for iA=1:1:numel(Data.Airports);
% % % % %   for jA=1:1:numel(Data.Airports);
% % % % %    
% % % % %     ThisPair = find(Data.Results.Dep == iA ...
% % % % %                   & Data.Results.Arr == jA);
% % % % %     if numel(ThisPair) > Settings.MinPoints
% % % % %       Pairs(iA,jA) = numel(ThisPair);
% % % % %     end
% % % % %     
% % % % %   end
% % % % % end
% % % % % UniquePairs = find(Pairs > 0);
% % % % % [iA,iB] = ind2sub(size(Pairs),UniquePairs);
% % % % % 
% % % % % 
% % % % % PlotsPerPage = 5;
% % % % % Pages = ceil(numel(UniquePairs)./PlotsPerPage);
% % % % % 
% % % % % for iPlot=1:1:numel(UniquePairs)
% % % % %   
% % % % %   %new page?
% % % % %   if mod(iPlot,PlotsPerPage) == 1;
% % % % %     %new figure
% % % % %     figure(iPlot)
% % % % %     clf
% % % % %     set(gcf,'color','w')
% % % % %     k = 0;
% % % % %   end
% % % % %   
% % % % %   %extract data
% % % % %   ThisPair.A = iA(iPlot);
% % % % %   ThisPair.B = iB(iPlot);
% % % % %   
% % % % %   Dep = Data.Airports(ThisPair.A); Dep = Dep{1};
% % % % %   Arr = Data.Airports(ThisPair.B); Arr = Arr{1};
% % % % %   
% % % % %   ThisPair.Rows = find(Results(:,1) == ThisPair.A ...
% % % % %                      & Results(:,2) == ThisPair.B);
% % % % %   
% % % % %   
% % % % %   %loop over indices
% % % % %   for iIndex=1:1:numel(Settings.Indices)
% % % % %     
% % % % %     %prepare panel
% % % % %     k = k+1;
% % % % %     subplot(PlotsPerPage,numel(Settings.Indices),k)
% % % % %     
% % % % %     %plot travel time vs this index
% % % % %     plot(Results(ThisPair.Rows,4),Results(ThisPair.Rows,iIndex+4),'kx')
% % % % %     hold on
% % % % %     
% % % % %     %linear fit
% % % % %     x = Results(ThisPair.Rows,4);
% % % % %     y = Results(ThisPair.Rows,iIndex+4);
% % % % %     Good = find(~isnan(x+y));
% % % % %     if numel(Good) > 10;
% % % % %       p = polyfit(x,y,1);
% % % % %       y2 = polyval(p,x);
% % % % %       plot(x,y2,'r-')
% % % % %     end
% % % % %       
% % % % %     %label
% % % % % %     xlabel('Relative Time')
% % % % %     ylabel(Settings.Indices{iIndex})
% % % % %     title([Dep,' -> ',Arr])
% % % % %     
% % % % %     
% % % % %     drawnow
% % % % %   end
% % % % % end
% % % % % 
% % % % % % % % for iPage=1;
% % % % % % % %   
% % % % % % % %   %create figure
% % % % % % % %   figure(iPage)
% % % % % % % %   clf
% % % % % % % %   set(gcf,'color','w')
% % % % % % % %   
% % % % % % % %   %loop over pairs
% % % % % % % %   for iPair = (iPage-1).*PlotsPerPage:1:
% % % % % % % %   
% % % % % % % %   
% % % % % % % %   
% % % % % % % % end
% % % % % % % % 
% % % % % % % % 
% % % % % % % % 
% % % % % % % % 
