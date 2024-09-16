function process_optimal_lags(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute the best-fit time shift for each index
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/09/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++++')
disp('Computing optimal index lags')
disp('+++++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat'])

%load indices - arbitrarily load the first direction set only, as this must exist and all should be identical at zero lag
A = load([Settings.Paths.DataDir,'/',Settings.ID,'_indices.mat']);
RangeStore    = A.RangeStore.(   Settings.Choices.Directions{1});
DateIndices   = A.DateIndices.(  Settings.Choices.Directions{1});
FlightIndices = A.FlightIndices.(Settings.Choices.Directions{1});
clear A

%join indices onto flight data, and drop
Data = innerjoin(FlightData,FlightIndices);

clear FlightData FlightIndices RangeStore RouteData OverallMedianTimes 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% average flight delays by DAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TimeScale = DateIndices.Date;
TimeTakenDays = NaN(numel(TimeScale),numel(Settings.Choices.Directions));

for iDir=1:1:numel(Settings.Choices.Directions)

  %get flights and dates in this direction
  TimeTaken = Data.Delay.All(Data.Direction == Settings.Choices.Directions{iDir});
  Dates     = Data.Date(     Data.Direction == Settings.Choices.Directions{iDir});


  %average by day
  for iDay=1:1:numel(TimeScale)
    OnThisDay = find(floor(Dates) == TimeScale(iDay));
    TimeTakenDays(iDay,iDir) = nanmean(TimeTaken(OnThisDay));
  end

end
clear TimeTaken iDay TimeScale OnThisDay Data




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cross-correlate to find optimal lags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
subplot = @(m,n,p) subtightplot (m, n, p, 0.01, 0.2, 0.2);
Letters = 'abcdefgh';

Indices = Settings.Indices.List;
BestLag = NaN(numel(Indices),numel(Settings.Choices.Directions));
k = 0;
for iDir=1:1:numel(Settings.Choices.Directions)
  for iIndex = 1:1:numel(Indices)

    if strcmp(Indices{iIndex},'Time'); BestLag(iIndex,iDir) = 0; continue; end %time shouldn't show a lag!

    %load index
    I = DateIndices.(Indices{iIndex});

    Good = find(~isnan(I+TimeTakenDays(:,iDir)));
    I = I(Good); TimeTaken = TimeTakenDays(Good,iDir); clear Good

    %cross-correlate with flight data
    [c,lags] = xcorr(I,TimeTaken,Settings.Choices.MaxLag,'coeff');
    Bad = find(lags > 0);
    c(Bad) = []; lags(Bad) = [];
    clear Bad
    % c = abs(c);

    %find maximum absolute cross-correlation
    [~,idx] = max(abs(c));
    BestLag(iIndex,iDir) = lags(idx);


    % figure
    subplot(1,numel(Settings.Choices.Directions),iDir); hold on
    plot(lags,c,'color',Settings.Indices.Colours.(Settings.Indices.List{iIndex}),'linewi',1)
    plot(lags(idx),c(idx),'o','color',          Settings.Indices.Colours.(Settings.Indices.List{iIndex}), ...
                         'MarkerFaceColor',Settings.Indices.Colours.(Settings.Indices.List{iIndex}), ...
                         'MarkerSize',8)
    xlim([min(lags),0])
    ylim([-1,1].*0.38)
    box on; grid on

    switch Settings.Choices.Directions{iDir}
      case 'W'; Dir = 'Westwards';
      case 'R'; Dir = 'Round-trip';
      case 'E'; Dir = 'Eastwards';
      otherwise; Dir = '';
    end

    title(['(',Letters(iDir),') ',Dir])
    if iDir == 1; ylabel('Cross-correlation'); end
    if iDir == numel(Settings.Choices.Directions); set(gca,'yaxislocation','right'); end
    if iDir ~= 1 & iDir ~= numel(Settings.Choices.Directions); set(gca,'YTickLabel',{}); end
    xlabel('Lead time [days]')

    if iDir == numel(Settings.Choices.Directions)
      k = k+1;
      ypos = max(gca().YLim)-0.2.*range(gca().YLim).*k;
      text(50,ypos,Settings.Indices.List{iIndex},'color',Settings.Indices.Colours.(Settings.Indices.List{iIndex}),'fontsize',20)

    end

  end
end


%store
Lags = struct();
for iDir=1:1:numel(Settings.Choices.Directions)
  for iIndex=1:1:numel(Indices)
    Lags.(Settings.Choices.Directions{iDir}).(Indices{iIndex}) = BestLag(iIndex,iDir);
  end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save, then call the climate index generator with lag mode on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([Settings.Paths.DataDir,'/',Settings.ID,'_optimallags.mat'],'Lags')

disp('--------------------------')
disp('Lag data generated')
disp('--------------------------')

process_climate_indices_v3(Settings,1)
