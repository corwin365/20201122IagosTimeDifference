function analysis_index_stat_tests(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%assess independence of cliamte indices used
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/09/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Index statistical tests')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load indices. Arbitrarily take the first direction, as they are unlagged
A = load([Settings.Paths.DataDir,'/',Settings.ID,'_indices.mat']);
RangeStore    = A.RangeStore.(   Settings.Choices.Directions{1});
DateIndices   = A.DateIndices.(  Settings.Choices.Directions{1});
FlightIndices = A.FlightIndices.(Settings.Choices.Directions{1});
clear A

%load data
load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat'])

%convert time-based indices to an array
%first column is ignored because this is the date
%inpaint_nans is used because our date-based NAO series contains two NaN 
%values, which have no flight data nearby so don't affect any other analysed
Indices = inpaint_nans(table2array(DateIndices(:,2:end)));


%also get flight-based indices
FIndices= table2array(FlightIndices(:,3:end));
Times = FlightData.t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% corrplot on predictors
%requires matlab econometrics toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% corrplot(Indices)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% multicolinearity - VIF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%create a storage array for the VIFs
NIndices = size(Indices,2);
Store = NaN(NIndices,1);

%work out VIFs
for iIndex=1:1:NIndices

  idx = 1:1:NIndices; idx(idx == iIndex) = [];
  I = Indices(:,idx);
  J = Indices(:,iIndex);
  mdl = fitlm(I,J);
  
  VIF(iIndex) = 1./(1-mdl.Rsquared.Adjusted.^2);
end

disp(['VIFs lie in the range ',num2str(min(VIF)),' -- ',num2str(max(VIF))])


clear NIndices iIndex idx I J mdl VIF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% autocorrelation - Ljung-Box Q-test
%requires matlab econometrics toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Belsley correlation tests
%requires matlab econometrics toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,~,VarDecompTbl] = collintest(Indices);

disp(['Belsley tests lie in the range ',num2str(min(VarDecompTbl(:))),' -- ',num2str(max(VarDecompTbl(:)))])

