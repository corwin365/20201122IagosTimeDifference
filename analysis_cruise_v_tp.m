function analysis_cruise_v_tp(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%is there a correlation between climate indices and cruise heights?
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/10/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Check for r(TP,z)')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load flight metadata
Meta = load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat']);

%load track data
Z = load([Settings.Paths.DataDir,'/',Settings.ID,'_flighttracks.mat']);
Z = Z.Store.Z;

%get tropopause height data
TP = load([Settings.Paths.DataDir,'/',Settings.ID,'_flightwinds.mat'],'TPStore');
TP = TP.TPStore;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for each flight, take the middle third, and correlate flight altitude
%with tropopause height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = NaN(numel(TP),1);

for iFlight=1:1:numel(r);

  %get TP and Z series
  T = TP{iFlight};
  z = Z{iFlight};

  %drop any bad data
  z(z < 8000 | z > 13000) = NaN; %out of usual height range for flights
  % T(T == 700) = NaN; %filler value for bad tropopause height
  if nansum(z) == 0; continue; end
  if nansum(T) == 0; continue; end
  


  %take middle third
  n = numel(z)./3;
  idx = floor(n):1:ceil(2.*n);
  T = p2h(T(idx));
  z = z(idx);

  %find good data
  Good = find(~isnan(T+z));
  if numel(Good) < 2; continue; end

  %correlate
  coeff = corrcoef(T(Good),z(Good));
  r(iFlight) = coeff(2);
 

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Good = find(~isnan(r));


x = -1:0.1:1;
y = hist(r(Good),x);

y = y./sum(y(:));

patch(x,y,[51,102,0]./255,'edgecolor','none')
xlabel('Correlation coefficient')
ylabel('Proportion of data in bin')
title('(d) r(plane height, tropopause height)')
box on
grid off

stop