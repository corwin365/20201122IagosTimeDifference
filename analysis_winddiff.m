function analysis_winddiff(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compare ERA5 and IAGOOS wind speed
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/09/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Wind speed comparison')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load flight metadata
Meta = load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat']);

%load indices
if Settings.Choices.ApplyLags == 0;
      Indices = load([Settings.Paths.DataDir,'/',Settings.ID,'_indices.mat'])
else; Indices = load([Settings.Paths.DataDir,'/',Settings.ID,'_laggedindices.mat']);
end

%load flight track and sampled wind data, and merge

%merge together
A = load([Settings.Paths.DataDir,'/',Settings.ID,'_flighttracks.mat']);
Tracks = A.Store; Tracks = rmfield(Tracks,{'U','V'});
Tracks.U.Planes = A.Store.U;
Tracks.V.Planes = A.Store.V;
clear A
B = load([Settings.Paths.DataDir,'/',Settings.ID,'_flightwinds.mat']);
Tracks.U.Model  = B.UStore;
Tracks.V.Model  = B.VStore;
clear B


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute a correlation coefficient and RMSD for every flight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Stats = NaN(numel(Tracks.Lat),2);
tic
for iFlight=1:1:numel(Tracks.Lat);

 %get track data
 a = Tracks.U.Planes{iFlight};
 b = Tracks.U.Model{ iFlight};
 if numel(a) ~= numel(b); continue; end %this shouldn't happen
 Good = find(~isnan(a+b));
 if numel(Good) < 0.95.*numel(a); continue; end
 
 %correlation
 r = corrcoef(a(Good),b(Good));
 Stats(iFlight,1) = r(2); 
 
 % if Stats(iFlight,1) < 0.2; stop; end
 %rmsd
 Stats(iFlight,2) = sqrt(sum((a-b).^2));

end
toc

stop
