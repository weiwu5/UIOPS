%% Preamble
% This helper script populates variables to be input into the imgProc_sm.m
% script. For use with the Munich Data Processing Workshop, 7-9 July 2017.
%
% Copyright Joseph Finlon, Univ. Illinois 2017.
% 
%% Stitch Together 2DS File (Original File too Large for GitHub)
fprintf('Initializing data files.\n')
fileDirectory = [pwd, '/'];
addpath(fullfile(pwd, '..')) % directory for main scripts
date = '20151205';

inFilename1 = [fileDirectory, '2DS.', date, '_subset1.V.cdf'];
inFilename2 = [fileDirectory, '2DS.', date, '_subset2.V.cdf'];

yr1 = ncread(inFilename1,'year'); yr2 = ncread(inFilename2,'year');
mo1 = ncread(inFilename1,'month'); mo2 = ncread(inFilename2,'month');
dy1 = ncread(inFilename1,'day'); dy2 = ncread(inFilename2,'day');
hr1 = ncread(inFilename1,'hour'); hr2 = ncread(inFilename2,'hour');
mn1 = ncread(inFilename1,'minute'); mn2 = ncread(inFilename2,'minute');
sec1 = ncread(inFilename1,'second'); sec2 = ncread(inFilename2,'second');
ms1 = ncread(inFilename1,'millisec'); ms2 = ncread(inFilename2,'millisec');
wk1 = ncread(inFilename1,'wkday'); wk2 = ncread(inFilename2,'wkday');
data1 = ncread(inFilename1,'data'); data2 = ncread(inFilename2,'data');

yr = [yr1; yr2]; mo = [mo1; mo2]; dy = [dy1; dy2]; hr = [hr1; hr2];
mn = [mn1; mn2]; sec = [sec1; sec2]; ms = [ms1; ms2]; wk = [wk1; wk2];
data = cat(3, data1, data2);

versionID = fopen('version.txt', 'r');
software_string = fscanf(versionID, '%s');
fclose(versionID);

f = netcdf.create([fileDirectory, '2DS.', date,...
    '_subset.V.cdf'], 'clobber');

NC_GLOBAL = netcdf.getConstant('NC_GLOBAL'); % Added file attributes ~ Joe Finlon 02/07/20
netcdf.putAtt(f, NC_GLOBAL, 'Software', software_string);
netcdf.putAtt(f, NC_GLOBAL, 'Creation Time', datestr(now, 'yyyy/mm/dd HH:MM:SS'));
netcdf.putAtt(f, NC_GLOBAL, 'Probe Channel_Orientation', 'Vertical');

dimid0 = netcdf.defDim(f,'time',size(data,3));
dimid1 = netcdf.defDim(f,'ImgRowlen',size(data,1));
dimid2 = netcdf.defDim(f,'ImgBlocklen',size(data,2));

varid0 = netcdf.defVar(f,'year','short',dimid0);
varid1 = netcdf.defVar(f,'month','byte',dimid0);
varid2 = netcdf.defVar(f,'day','byte',dimid0);
varid3 = netcdf.defVar(f,'hour','byte',dimid0);
varid4 = netcdf.defVar(f,'minute','byte',dimid0);
varid5 = netcdf.defVar(f,'second','byte',dimid0);
varid6 = netcdf.defVar(f,'millisec','short',dimid0);
varid7 = netcdf.defVar(f,'wkday','byte',dimid0);
varid8 = netcdf.defVar(f,'data','int',[dimid1 dimid2 dimid0]);
netcdf.endDef(f)

netcdf.putVar ( f, varid0, yr );
netcdf.putVar ( f, varid1, mo );
netcdf.putVar ( f, varid2, dy );
netcdf.putVar ( f, varid3, hr );
netcdf.putVar ( f, varid4, mn );
netcdf.putVar ( f, varid5, sec );
netcdf.putVar ( f, varid6, ms );
netcdf.putVar ( f, varid7, wk );
netcdf.putVar ( f, varid8, data );
netcdf.close(f)

%% Initialize Files Directory
if exist('files','dir')~=7
    mkdir files;
end

%% Populate Input Variables
% Using 05 December 2015 case from OLYMPEX (partial flight)

projectName = 'GPM'; % will use predefined probe settings already in sizeDist.m
probeName = {'2DS', 'HVPS'}; % use both '2DS' and 'HVPS' for this example
numCPU = 1; % will need only 1 CPU to process particle properties
framesPerCPU = 100000; % max # image records to process per CPU
createAspectRatio = 0; % do not output particle length/width using rectangle and ellipse fits
calcAllDiodeStats = 0; % do not output shadowed diode stats on a per-particle basis

flightFilename = [fileDirectory, 'flightData_', date, '.cdf'];
fltTime_HHMMSS = ncread(flightFilename, 'Time');
fltTAS = ncread(flightFilename, 'TAS');

%% Run imgProc_sm.m Script

for iter=1:length(probeName) % loop over desired probes
    disp(['Processing particles from the ', probeName{iter}, '.'])
    inFilename = [fileDirectory, probeName{iter}, '.', date, '_subset.V.cdf'];
    outFilename = [fileDirectory, 'files/proc', probeName{iter}, '.', date, '_subset.V.cdf'];

    imgProc_sm(inFilename, outFilename, probeName{iter}, numCPU, framesPerCPU,...
        projectName, createAspectRatio, calcAllDiodeStats, fltTime_HHMMSS, fltTAS)
end