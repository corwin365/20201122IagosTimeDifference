function FlightData = import_ed_flights_v2(filename, dataLines)
%IMPORTFILE Import data from a text file
%  TFMS2007005 = IMPORTFILE(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  TFMS2007005 = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  tfms2007005 = importfile("C:\Data\corwin\ed_aircraft\tfms_2007005.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 23-Dec-2022 23:14:41

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
  dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 10);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "AIRCRAFT_ID", "UNIQUE_FLIGHT_ID", "DEPT_APRT", "DEPT_DATE_UTC", "DEPT_TIME_UTC", "ARR_APRT", "ARR_DATE_UTC", "ARR_TIME_UTC", "FLIGHT_TIME"];
opts.VariableTypes = ["double", "string", "double", "string", "double", "double", "string", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "AIRCRAFT_ID", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["AIRCRAFT_ID", "DEPT_APRT", "ARR_APRT"], "EmptyFieldRule", "auto");

% Import the data
FlightData = readtable(filename, opts);

end
