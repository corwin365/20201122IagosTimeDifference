function analysis_frac_strat_trop(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%does the fraction of the flight spent in the stratosphere affect
%our results?
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/12/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Check stratosphere effects')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load flight metadata
Meta = load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat']);

%load track data
Z = load([Settings.Paths.DataDir,'/',Settings.ID,'_flighttracks.mat']);
L = Z.Store.Lat;
Z = Z.Store.Z;

%get tropopause height data
TP = load([Settings.Paths.DataDir,'/',Settings.ID,'_flightwinds.mat'],'TPStore');
TP = TP.TPStore;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the fraction of each flight in the stratosphere
%and also the maximum latitude reached
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FracS = NaN(numel(TP),1);
MaxL  = FracS;

textprogressbar('Finding fraction of each flight above stratopause and max lat ')
for iFlight=1:1:numel(TP)

  %get tropopause and flight altitude
  z = Z{iFlight}./1000;
  t = p2h(TP{iFlight});
  if numel(t) ~= numel(z); continue; end

  %hence, find fraction where flight is above the tropopause
  delta = z-t; ss = find(delta > 0);
  FracS(iFlight) = numel(ss)./numel(t);

  %also find the maximum latitude
  MaxL(iFlight) = nanmax(L{iFlight});

  if mod(iFlight,1000) == 1; textprogressbar(iFlight./numel(TP).*100); end
end; clear z t delta ss
textprogressbar('!')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now, for each direction and season, correlate them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Corrs = NaN(2,1,numel(Settings.Choices.Directions));


for iSeason=1 %we're interested in all data
  for iDir=1:1:numel(Settings.Choices.Directions)
    if Settings.Choices.Directions{iDir} == 'R'; continue; end

    %get flight times in this direction/season
    idx = find(Meta.FlightData.Direction == Settings.Choices.Directions{iDir} ...
             & Meta.FlightData.InSeasons.(Settings.Seasons.List{iSeason}) == 1);
    
    %hence, get delays for these flights
    delays = Meta.FlightData.Delay.(Settings.Seasons.List{iSeason});
    delays = delays(idx);

    %and correlate the delay with the fraction
    Good = find(~isnan(FracS(idx)+ delays));
    r = corrcoef(FracS(idx(Good)),delays(Good));
    Corrs(1,1,iDir) = r(2);

    %and the delay with the maxlat
    Good = find(~isnan(MaxL(idx)+ delays));
    r = corrcoef(MaxL(idx(Good)),delays(Good));
    Corrs(2,1,iDir) = r(2);




  end; clear iDir
end; clear iSeason
