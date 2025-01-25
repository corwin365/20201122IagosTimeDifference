function analysis_cost_histo(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cost analysis of delays implied by regression analysis
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/11/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('++++++++++++++++++++++++++++++++++++++++++++')
disp('Regression imputed cost analysis - histogram')
disp('++++++++++++++++++++++++++++++++++++++++++++')


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
%% violin plots - postprocessing and plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %choose a direction
% Direction = find(contains(Settings.Choices.Directions, 'R'));

%choose a season
Season = find(contains(Settings.Seasons.List, 'All'));

%choose a number of kdf fitting centres
%this is an initial value which wil be slightly adjusted to produce bins of round-number widths
nbins = 100;

%prepare figure
clf
subplot = @(m,n,p) subtightplot (m, n, p, [0.07,0.01], 0.06, 0.15);
k = 0;

%loop over types of cost
ctlist = [4] %[3,4];
for Direction=1:1:3;
for ict = 1:1:numel(ctlist)

  %get costs
  c = squeeze(Cost(:,Direction,Season,:,ctlist(ict)));

  % %rebase to zero median
  % for iInd=1:1:numel(Indices)
  %   c(:,iInd) = c(:,iInd) - nanmedian(c(:,iInd));
  % end


  %loop over indicex
  for iInd=[numel(Indices),1:1:numel(Indices)-1];

    %define x-axis, using a fixed number of bins across the actual data
    x = linspace(min(c(:,iInd)),max(c(:,iInd)),nbins);

    % %prepare KDF
    % kd = pdf(fitdist(c(:,iInd),'kernel','Kernel','normal'),x);
    kd = hist(c(:,iInd),x);
    kd = kd./max(kd);

    %prepare panel
    k = k+1;
    subplot(6,numel(Indices),k)
    grid on; hold on; box on;
    axis([[-1,1].*max(abs(x)) [-1,1].*1.25.*ceil(max(kd))])


    %get colour
    if strcmp(Indices{iInd},'Sum'); colour = [102,51,0]./255;%'k';
    else                            colour = Settings.Indices.Colours.(Indices{iInd});
    end

    if ctlist(ict) == 1 | ctlist(ict) == 3; set(gca,'XAxisLocation','top'); end


    %plot the violin
    patch([x,x(end:-1:1)],[kd,-kd(end:-1:1)],colour,'linewi',1,'facealpha',0.4,'edgecolor',colour);
    xlabel([CostUnits{ctlist(ict)},' per month'])


    %overlay a box plot
    prc = prctile(c(:,iInd),[2.5,18,50,82,97.5]);
    for iPrc=[2,4]; plot([1,1].*prc(iPrc),[-1,1].*0.25,'k-');end
    for iPrc=[1,5]; plot([1,1].*prc(iPrc),[-1,1].*0.15,'k-');end
    plot(prc([2,4]),[1,1].*0.25,'k-');  plot(prc([2,4]),-[1,1].*0.25,'k-')
    plot(prc([3,3]),[-1,1].*0.25,'k-','linewi',3);
    plot(prc([1,5]),[0,0],'k-')

    %print for numbers in text
    disp('================================')
    disp([Settings.Choices.Directions{Direction},'  ',Indices{iInd}])
    disp(['2stdev: ',num2str(round(prc([1,5]),2,'decimals')),'   ',num2str(round(range(prc([1,5])),2,'decimals'))])
    disp(['range : ',num2str(round(minmax(x),2, 'decimals')),'   ',num2str(round(range(minmax(x)),2,'decimals'))])



    %clean up axes
    set(gca,'yticklabel',{})
    % set(gca,'xticklabel',{})

    %name
    text(-max(abs(x)),-1.2,[' ',Indices{iInd}],'HorizontalAlignment','left','VerticalAlignment','bottom')


    % %label binsize
    % text(max(x),0.9.*ceil(max(kd)),[num2str(binsize),CostUnits{ctlist(ict)},' bins'],'HorizontalAlignment','right','FontSize',12)

    drawnow


  end; clear iIndex
end
end


stop