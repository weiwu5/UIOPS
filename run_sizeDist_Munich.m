%% Preamble
% This helper script populates variables to be input into the sizeDist.m
% script. For use with the Munich Data Processing Workshop, 7-9 July 2017.
%
% Copyright Joseph Finlon, Univ. Illinois 2017.
% 
%% Populate Input Variables
% Using 05 December 2015 case from OLYMPEX (partial flight)

clearvars;
fileDirectory = [pwd, '/files/'];
projectName = 'GPM'; % will use predefined probe settings already in sizeDist.m
date = '20151205';
probeName = {'2DS', 'HVPS'}; % use both '2DS' and 'HVPS' for this example
sizeMethod = 6; % will use D of minimum enclosing circle for this example
saMethod = 2; % will use option #2 (Heymsfield & Parrish correction) for this example

% Get flight variables
flightFile = [pwd, '/flightData_', date, '.cdf'];

flightTime = ncread(flightFile, 'Time'); % time (HHMMSS) since 0000 UTC
tas = ncread(flightFile, 'TAS'); % true air speed (m/s)
pressure = ncread(flightFile, 'Pres'); % static pressure (hPa)
temperature = ncread(flightFile, 'Temp'); % air temperature (deg C) corrected for dynamic heating

% Get inter-arrival threshold data filename
iaThreshFilename = [fileDirectory, 'intArrThreshold_', date, '.cdf'];

%% Run sizeDist.m Script

for iter=1:length(probeName)
    % The '.V' indicates using the vertical channel (orientation) of the 2DS
    % (HVPS)
    inFilename = [fileDirectory, 'proc', probeName{iter}, '.', date, '_subset.V.cdf'];
    outFilename = [fileDirectory, 'sd', probeName{iter}, '.', date, '_subset.V.cdf'];
    
    sizeDist(inFilename, outFilename, tas, flightTime, probeName{iter}, sizeMethod,...
        saMethod, pressure, temperature, projectName, date, iaThreshFilename)
end