%% Preamble
% This helper script populates variables to be input into the
% IntArrAnalysis.m script. For use with the Munich Data Processing
% Workshop, 7-9 July 2017.
%
% Copyright Joseph Finlon, Univ. Illinois 2017.
% 
%% Populate Input Variables
% Using 05 December 2015 case from OLYMPEX (partial flight)

clearvars;
addpath(fullfile(pwd, '..')) % directory for main scripts
fileDirectory = [pwd, '/files/'];
date = '20151205';
probeName = '2DS'; % will use only the 2DS for this example
runAnalysis = 1; % 0 for only plotting distributions; 1 for running fitting technique
numParticles = 3000; % # particles belonging to each bimodal fit
numCPU = 1; % will need only 1 CPU to process inter-arrival times

inFilename = [fileDirectory, 'proc', probeName, '.', date, '_subset.V.cdf'];

%% Run IntArrAnalysis.m Script

disp(['Processing ', probeName, ' inter-arrival times for ', date])

IntArrAnalysis(inFilename, fileDirectory, probeName, runAnalysis,...
    numParticles, date, numCPU)