function analysis_cost(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cost analysis of delays implied by regression analysis
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/11/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++++++')
disp('Regression imputed cost analysis')
disp('+++++++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load flight data
load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat'])
clear RouteData

%load indices
if Settings.Choices.ApplyLags == 0;
      load([Settings.Paths.DataDir,'/',Settings.ID,'_indices.mat'])
else; load([Settings.Paths.DataDir,'/',Settings.ID,'_laggedindices.mat'])
end
clear RangeStore OverallMedianTimes

%load regression results and scale from seconds to minutes
%resulting dimension order is [Directions,Seasons,Indices.List] 
Reg = load([Settings.Paths.DataDir,'/',Settings.ID,'_regressioncoefficients.mat']);
RegSE = Reg.Reg.SE ./60;
Reg   = Reg.Reg.Est./60;

%load cost results
MasterCost = load([Settings.Paths.DataDir,'/',Settings.ID,'_fuel_and_emissions.mat']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% specify costs for each flight, in units of [per minute of delay scaled to all flights that day]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create table cols needed
FlightData.Cost_Co2        = FlightData.FlightIndex .* NaN;
FlightData.Cost_Co2_Fixed  = FlightData.FlightIndex .* NaN;
FlightData.Cost_Fuel       = FlightData.FlightIndex .* NaN;
FlightData.Cost_Fuel_Fixed = FlightData.FlightIndex .* NaN;

%interpolate total number of flights to date of each flight
NFlights.W = interp1(MasterCost.NFlights.Time,MasterCost.NFlights.Val(:,1),FlightData.Date);
NFlights.E = interp1(MasterCost.NFlights.Time,MasterCost.NFlights.Val(:,2),FlightData.Date);

%interpolate fuel price to each flight
FuelPrice = interp1(MasterCost.Fuel.Time,MasterCost.Fuel.Price,FlightData.Date);

%find cost of fuel and number of flights at a fixed date
FixedDate      = Settings.Choices.PriceFixDate;
FixedFuelPrice = interp1(MasterCost.Fuel.Time,MasterCost.Fuel.Price,FixedDate);
FixedN         = mean(interp1(MasterCost.NFlights.Time,MasterCost.NFlights.Val,FixedDate));

for iF=1:1:size(FlightData,1)

  if FlightData.Direction(iF) == 'W' | FlightData.Direction(iF) == 'E'

    %one-way flight
    n = NFlights.(string(FlightData.Direction(iF))); n = n(iF);
    FlightData.Cost_Co2(       iF) = n                        .* MasterCost.Co2PerMinute(    MasterCost.PlaneType(iF));
    FlightData.Cost_Fuel(      iF) = n      .* FuelPrice(iF)  .* MasterCost.GallonsPerMinute(MasterCost.PlaneType(iF));
    FlightData.Cost_Co2_Fixed( iF) = FixedN                   .* MasterCost.Co2PerMinute(    MasterCost.PlaneType(iF));
    FlightData.Cost_Fuel_Fixed(iF) = FixedN .* FixedFuelPrice .* MasterCost.GallonsPerMinute(MasterCost.PlaneType(iF));   

  else
    %round-trip flight - take mean of each leg
    n = (NFlights.E(iF) + NFlights.W(iF))./2;
    FlightData.Cost_Co2(        iF) = n                        .* mean(MasterCost.Co2PerMinute(    MasterCost.PlaneType([FlightData.OriginalW(iF),FlightData.OriginalE(iF)])));
    FlightData.Cost_Fuel(       iF) = n      .* FuelPrice(iF)  .* mean(MasterCost.GallonsPerMinute(MasterCost.PlaneType([FlightData.OriginalW(iF),FlightData.OriginalE(iF)])));
    FlightData.Cost_Co2_Fixed(  iF) = FixedN                   .* mean(MasterCost.Co2PerMinute(    MasterCost.PlaneType([FlightData.OriginalW(iF),FlightData.OriginalE(iF)])));
    FlightData.Cost_Fuel_Fixed( iF) = FixedN .* FixedFuelPrice .* mean(MasterCost.GallonsPerMinute(MasterCost.PlaneType([FlightData.OriginalW(iF),FlightData.OriginalE(iF)])));
  end
  clear n

end; clear iF

clear NFlights FixedDate FixedN FuelPrice


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for each each flight compute the price due to the indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%specify cost types
CostTypes = {'Cost_Co2','Cost_Fuel','Cost_Co2_Fixed','Cost_Fuel_Fixed'};
CostUnits = {'kT CO2','Million USD','kT CO2','Million USD'};
[yy,~,~] = datevec(Settings.Choices.PriceFixDate);
CostTypeNames = {'CO2','Actual USD',[num2str(yy),'-level CO2'],[num2str(yy),'-level USD']};
clear yy

%create storage arrays
Cost   = NaN([size(FlightData,1),size(Reg),numel(CostTypes)]);
CostSE = Cost;

%compute values for each flight in each [dir,season] for each index with a simple scaling
for iF=1:1:size(FlightData,1)
  for iCostType=1:1:numel(CostTypes)
    c  = FlightData.(CostTypes{iCostType});
    Cost(  iF,:,:,:,iCostType) = c(iF) .* Reg;
    CostSE(iF,:,:,:,iCostType) = c(iF) .* RegSE;
    clear c
  end; clear iCostType
end; clear iF
clear Reg RegSE 

%now remove combinations that didn't actually happen for that specific flight:

%----> plane didn't fly in this direction
for iDir=1:1:numel(Settings.Choices.Directions)
  DidntHappen = find(FlightData.Direction ~= Settings.Choices.Directions(iDir));
  Cost(  DidntHappen,iDir,:,:,:) = NaN;
  CostSE(DidntHappen,iDir,:,:,:) = NaN; 
  clear DidntHappen
end; clear iDir

%----> plane didn't fly in this season
SeasFlags = NaN.*Cost;
sz = size(Cost);
for iF = 1:1:size(FlightData,1)

  %create an array which is 1 if the flight was in this season and NaN otherwise
  ThisFlightSeasFlags = single(table2array(FlightData.InSeasons(iF,:)));
  ThisFlightSeasFlags(ThisFlightSeasFlags == 0) = NaN;
  SeasFlags(iF,:,:,:,:) = repmat(ThisFlightSeasFlags,[sz(2),1,sz(4),sz(5)]);

  clear ThisFlightSeasFlags
end; 
Cost   = Cost   .* SeasFlags;
CostSE = CostSE .* SeasFlags;
clear iF sz SeasFlags


%finally, scale each value to the actual value of the index at the time of that flight
%remember, the indices ARE LAGGED DIFFERENTLY FOR EACH DIRECTION, so this has to be done by direction
for iDir=1:1:numel(Settings.Choices.Directions)



  for iIndex=1:1:numel(Settings.Indices.List)
    %flights, dir, season, index, cost
    c  = squeeze(Cost(  :,iDir,:,iIndex,:));
    cs = squeeze(CostSE(:,iDir,:,iIndex,:));
    I = FlightIndices.(Settings.Choices.Directions{iDir});
    I = I.(Settings.Indices.List{iIndex});
    I = repmat(I,1,numel(Settings.Seasons.List),numel(CostTypes));
    Cost(  :,iDir,:,iIndex,:) = c.*I;
    CostSE(:,iDir,:,iIndex,:) = cs.*I;

    clear I c cs

  end; clear iIndex
end; clear iDir


%remove 'time' from indices
idx = find(Settings.Indices.List ~= "Time");
Indices = Settings.Indices.List(idx);
Cost   = Cost(  :,:,:,idx,:);
CostSE = CostSE(:,:,:,idx,:);
clear idx

%add an index which is the sum of the daily values of all indices

Indices{end+1} = 'Sum';
Cost(  :,:,:,end+1,:) = nansum(Cost(  :,:,:,:,:),4);  Cost(    Cost == 0) = NaN;
CostSE(:,:,:,end+1,:) = nansum(CostSE(:,:,:,:,:),4);  CostSE(CostSE == 0) = NaN; %I *think* this is the right treatment?

%scale to millions of [USD or kg] per month
Cost   = Cost   ./ 1e6  .* (365/12);
CostSE = CostSE ./ 1e6  .* (365/12);



%remember: the order of these arrays is:
%[flights, directions, seasons, indices, cost types]
% CostTypes = {'Cost_Co2','Cost_Fuel','Cost_Co2_Fixed','Cost_Fuel_Fixed'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time plots - postprocessing and plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%for each month, find the sum of the daily effects of each index
%we're only going to use the all-seasons estimators again
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a monthly timescale and a daily timescale
[yy,~,~] = datevec(min(Settings.Choices.TimeRange));
MTime = datenum(yy,(1:1:100*12),1);
MTime(MTime > max(Settings.Choices.TimeRange)) = [];
DTime = Settings.Choices.TimeRange(1):1:Settings.Choices.TimeRange(2);
clear yy


%bin the flightwise cost estimates onto the daily timescale. Remember, these are scaled to represent MONTHLY values already
DailyCosts = NaN(numel(DTime),numel(Indices),numel(Settings.Choices.Directions),numel(CostTypes));
for iI=1:1:numel(Indices);
  for iDir=1:1:numel(Settings.Choices.Directions)
    for iType=1:1:numel(CostTypes)
      DailyCosts(:,iI,iDir,iType) = bin2matN(1,FlightData.Date,squeeze(Cost(:,iDir,find(contains(Settings.Seasons.List, 'All')),iI,iType)),DTime,'@nanmean');
      DailyCosts(:,iI,iDir,iType) = fillmissing(DailyCosts(:,iI,iDir,iType),'nearest');
    end; clear iType
  end; clear iDir
end; clear I

MonthlyCosts = NaN(numel(MTime),numel(Indices),numel(Settings.Choices.Directions),numel(CostTypes));
for iMonth=1:1:numel(MTime)-1
  idx = inrange(DTime,MTime([iMonth,iMonth+1]));
  MonthlyCosts(iMonth,:,:,:) = nanmean(DailyCosts(idx,:,:,:),1);
end; clear iMonth idx
clear DailyCosts

%prepare the data for a stacked histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first we need to drop the sum row
idx = find(Indices ~= "Sum");
MonthlyCosts = MonthlyCosts(:,idx,:,:);
clear idx

%ok. now for each month, generate a POSITIVE and NEGATIVE stack
PosCosts = MonthlyCosts; PosCosts(PosCosts <  0) = 0;
NegCosts = MonthlyCosts; NegCosts(NegCosts >= 0) = 0;


%overinterpolate in time to make it look (correctly) blocky
PosCosts = interp_1d_ndims(MTime,PosCosts,DTime,1,'nearest');
NegCosts = interp_1d_ndims(MTime,NegCosts,DTime,1,'nearest');


%% plot the data!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
subplot = @(m,n,p) subtightplot (m, n, p, 0.05, 0.08, 0.15);
hold on
k = 0;
Letters = 'abcdef';


for iDir=1:1:3
  % CostTypes = {'Cost_Co2','Cost_Fuel','Cost_Co2_Fixed','Cost_Fuel_Fixed'};
  for iCost = 2;%[1,2]
    k = k+1;
    subplot(3,1,k)


    %produce a background year pattern
    for iYear=1994:2:2024
      patch(datenum([0,1,1,0]+iYear,1,1),[-1,-1,1,1].*999,[1,1,1].*0.9,'edgecolor','none')
    end



    TheMax = 0;
    IndexOrder = [4,3,1,2]; % this forces the order TSI - QBO - ENSO - NAO, i.e. slowest-vaerying to fastest
    for iInd = numel(IndexOrder):-1:1
      for iPosNeg=1:2;
        switch iPosNeg
          case 1; Data = PosCosts;
          case 2; Data = NegCosts;
        end
        %colour for this index
        colour = Settings.Indices.Colours.(Indices{IndexOrder(iInd)});

        %sumcosts
        ThisSet = sum(Data(:,IndexOrder(1:iInd),iDir,iCost),2);
        if nanmax(abs(ThisSet)) > TheMax; TheMax = nanmax(abs(ThisSet)); end

        %generate and plot patch
        x = [min(DTime),DTime,max(DTime)];
        y = [0;ThisSet;0]';
        Good = find(~isnan(x+y));

        % patch(x(Good),y(Good),'w',   'edgecolor','none','facealpha',1)  %this is to stop the colours muddying together  
        patch(x(Good),y(Good),colour,'edgecolor',colour,'facealpha',1,'linewi',0.1)
        hold on


      end; clear iPosNeg

    end; clear iInd
    ylim([-1,1].*1.05.*TheMax)

    %vertical
    Ratio = FixedFuelPrice .* mean(MasterCost.GallonsPerMinute)./ mean(MasterCost.Co2PerMinute);
    switch iDir
      case 1; 
        set(gca,'ytick',[-75:25:75]);
        ylim([-1,1].*75)
        for iY=-75:25:75; text(datenum(2024,4,30),iY,num2str(round(iY.*Ratio.*10)./10),'HorizontalAlignment','left','VerticalAlignment','middle'); end
      case 2; 
        set(gca,'ytick',[-15: 5:15]);
        ylim([-1,1].*15)
        for iY=-15:5:15; text(datenum(2024,4,30),iY,num2str(round(iY.*Ratio.*10)./10),'HorizontalAlignment','left','VerticalAlignment','middle'); end
      case 3; 
        set(gca,'ytick',[-75:25:75]);
        ylim([-1,1].*75)
        for iY=-75:25:75; text(datenum(2024,4,30),iY,num2str(round(iY.*Ratio.*10)./10),'HorizontalAlignment','left','VerticalAlignment','middle'); end

    end



    %the first few months only have a small number of flights and hence a lot of filled missing values - truncate the graph
    xlim([datenum(1994,8,1),datenum(2024,4,1)])    
    % datetick('x','keeplimits')
    [yy,~,~] = datevec(datenum(1994:2:2024,6,20));
    set(gca,'xtick',datenum(1994:2:2024,6,20),'XTickLabel',yy)


    %overlay sum as a curve
    Sigma = PosCosts+NegCosts;
    Sigma = sum(Sigma(:,:,iDir,iCost),2);
    plot(DTime,smoothn(inpaint_nans(Sigma),31),'k-','LineWidth',0.5)

    % text(datenum(2023,6,1),0.9.*max(get(gca,'ylim')),['Total: ',num2str(round(nansum(Sigma./365.*12))),' ',CostUnits{iCost},' or ',num2str(round(nansum(Sigma./365.*12.*Ratio))),' M USD'],'HorizontalAlignment','right','fontsize',12)


    text(datenum(1994,10,1),0.9.*max(get(gca,'ylim')),['(',Letters(k),') ',Settings.Choices.Directions{iDir}])

    switch iCost
      case 1; ylabel('CO2 [kT]')
      case 2; ylabel('Fuel cost [million USD]')
    end

        %tidy up
        box on; grid on;
        set(gca,'layer','top')
        if     iDir == 1; set(gca,'xaxislocation','top'); 
        elseif iDir == 2; set(gca,'xticklabel',{});
        else              set(gca,'xaxislocation','bottom'); 
        end


    drawnow


  end
end





stop