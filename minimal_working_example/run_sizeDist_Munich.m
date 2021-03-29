%% Preamble
% This helper script populates variables to be input into the sizeDist.m
% script. For use with the Munich Data Processing Workshop, 7-9 July 2017.
%
% Copyright Joseph Finlon, Univ. Illinois 2017.
% 
%% Populate Input Variables
% Using 05 December 2015 case from OLYMPEX (partial flight)

clearvars;
addpath(fullfile(pwd, '..')) % directory for main scripts
fileDirectory = [pwd, '/files/'];
projectName = 'GPM'; % will use predefined probe settings already in sizeDist.m
date = '20151205';
probeName = {'2DS', 'HVPS'}; % use both '2DS' and 'HVPS' for this example
sizeMethod = 6; % will use D of minimum enclosing circle for this example
saMethod = 2; % will use option #2 (Heymsfield & Parrish correction) for this example
createAspectRatio = 0; % skipping processing of aspect ratio SDs
saveBadParticles = 0; % not saving data on rejected particles
saveIAandSV = 1; % not saving inter-arrival or sample volume info

% Get flight variables
flightFile = [pwd, '/flightData_', date, '.cdf'];

flightTime = ncread(flightFile, 'Time'); % time (HHMMSS) since 0000 UTC
tas = ncread(flightFile, 'TAS'); % true air speed (m/s)
pressure = ncread(flightFile, 'Pres'); % static pressure (hPa)
temperature = ncread(flightFile, 'Temp'); % air temperature (deg C) corrected for dynamic heating

% Apply other processing and output defaults


%% Run sizeDist.m Script

% 2DS Data (the '.V' indicates using the vertical channel of the 2DS)
inFilename = [fileDirectory, 'proc2DS.', date, '_subset.V.cdf'];
outFilename = [fileDirectory, 'sd2DS.', date, '_subset.V.cdf'];
iaThreshType = 1; % will apply time-dependent threshold
iaThreshFilename = [fileDirectory, '2DS_intArrThreshold_', date, '.cdf'];
sizeDist(inFilename, outFilename, tas, flightTime, probeName{1}, sizeMethod,...
    saMethod, pressure, temperature, iaThreshType, createAspectRatio, saveBadParticles,...
    saveIAandSV, projectName, date, iaThreshFilename)

% HVPS Data (the '.V' indicates using the vertical orientation of the HVPS)
inFilename = [fileDirectory, 'procHVPS.', date, '_subset.V.cdf'];
outFilename = [fileDirectory, 'sdHVPS.', date, '_subset.V.cdf'];
iaThreshType = 0; % will not apply threshold (per campaign/probe preset setting)
sizeDist(inFilename, outFilename, tas, flightTime, probeName{2}, sizeMethod,...
    saMethod, pressure, temperature, iaThreshType, createAspectRatio, saveBadParticles,...
    saveIAandSV, projectName, date)