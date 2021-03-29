%% Preamble
% This script combines SDs from two probes for plotting particle spectra &
% bulk properties. For use with the Munich Data Processing Workshop, 7-9
% July 2017.
%
% Copyright Joseph Finlon, Univ. Illinois 2017.
% 
%% Populate input parameters before reading the data
% Using 05 December 2015 case from OLYMPEX (partial flight)

clearvars;
addpath(fullfile(pwd, '..')) % directory for main scripts
fileDirectory = [pwd, '/files/'];
date = '20151205';
probeName = {'2DS', 'HVPS'}; % use both '2DS' and 'HVPS' for this example

startTime = '164130'; % start of time series plots
endTime = '164159'; % end of time series plots

inFilename1 = [fileDirectory, 'sd', probeName{1}, '.', date, '_subset.V.cdf'];
inFilename2 = [fileDirectory, 'sd', probeName{2}, '.', date, '_subset.V.cdf'];

%% Read and Trim the SD Data

% Read SD data from the higher resolution probe
time = ncread(inFilename1, 'time'); % flight time [UTC in HHMMSS]
binMin1 = ncread(inFilename1, 'bin_min'); % bin left endpoint [mm]
binMid1 = ncread(inFilename1, 'bin_mid'); % bin midpoint [mm]
binMax1 = ncread(inFilename1, 'bin_max'); % bin right endpoint [mm]
dD1 = ncread(inFilename1, 'bin_dD'); % bin width [mm]
ndf1 = ncread(inFilename1, 'conc_minR'); % N(D) using Dmax [cm^-4]
mdfHabit1 = ncread(inFilename1, 'mass'); % mass following Jackson et al. (2012) habit m-D relations [g cm^-4]
mdfBL1 = ncread(inFilename1, 'massBL'); % mass following Baker & Lawson (2006) m-A relationship [g cm^-4]
vt1 = ncread(inFilename1, 'vt'); % mass-weighted terminal velocity [g cm^-4]
pr1 = ncread(inFilename1, 'Prec_rate'); % precipitation rate [mm hr^-1]
area1 = ncread(inFilename1, 'Calcd_area'); % area following A-D relations [mm^2 cm^-4]
ar1 = ncread(inFilename1, 'mean_area_ratio'); % mean area ratio for each bin
perim1 = ncread(inFilename1, 'mean_perimeter'); % mean perimeter for each bin [um]
sv1 = ncread(inFilename1, 'sample_vol'); % sample volume for each bin [cm^3]

% Read SD data from the lower resolution probe
binMin2 = ncread(inFilename2, 'bin_min'); % bin left endpoint [mm]
binMid2 = ncread(inFilename2, 'bin_mid'); % bin midpoint [mm]
binMax2 = ncread(inFilename2, 'bin_max'); % bin right endpoint [mm]
dD2 = ncread(inFilename2, 'bin_dD'); % bin width [mm]
ndf2 = ncread(inFilename2, 'conc_minR'); % N(D) using Dmax [cm^-4]
mdfHabit2 = ncread(inFilename2, 'mass'); % mass following Jackson et al. (2012) habit m-D relations [g cm^-4]
mdfBL2 = ncread(inFilename2, 'massBL'); % mass following Baker & Lawson (2006) m-A relationship [g cm^-4]
vt2 = ncread(inFilename2, 'vt'); % mass-weighted terminal velocity [g cm^-4]
pr2 = ncread(inFilename2, 'Prec_rate'); % precipitation rate [mm hr^-1]
area2 = ncread(inFilename2, 'Calcd_area'); % area following A-D relations [mm^2 cm^-4]
ar2 = ncread(inFilename2, 'mean_area_ratio'); % mean area ratio for each bin
perim2 = ncread(inFilename2, 'mean_perimeter'); % mean perimeter for each bin [um]
sv2 = ncread(inFilename2, 'sample_vol'); % sample volume for each bin [cm^3]

% Find time indices to trim
timeInd1 = find(time==str2num(startTime));
timeInd2 = find(time==str2num(endTime));

% Find size bin indices to trim
startBinInd1 = find(binMin1==0.1); % ignore particles < 100 um
endBinInd1 = find(binMax1==1.0); % set probe cutoff at 1 mm
startBinInd2 = find(binMin2==1.0); % set probe cutoff at 1 mm
endBinInd2 = length(dD2);

% Trim data using time and particle size constraints
time = time(timeInd1:timeInd2);
binMin1 = binMin1(startBinInd1:endBinInd1);
binMid1 = binMid1(startBinInd1:endBinInd1);
binMax1 = binMax1(startBinInd1:endBinInd1);
dD1 = dD1(startBinInd1:endBinInd1);
ndf1 = ndf1(startBinInd1:endBinInd1, timeInd1:timeInd2);
mdfHabit1 = mdfHabit1(startBinInd1:endBinInd1, timeInd1:timeInd2);
mdfBL1 = mdfBL1(startBinInd1:endBinInd1, timeInd1:timeInd2);
vt1 = vt1(startBinInd1:endBinInd1, timeInd1:timeInd2);
pr1 = pr1(startBinInd1:endBinInd1, timeInd1:timeInd2);
area1 = area1(startBinInd1:endBinInd1, timeInd1:timeInd2);
ar1 = ar1(startBinInd1:endBinInd1, timeInd1:timeInd2);
perim1 = perim1(startBinInd1:endBinInd1, timeInd1:timeInd2);
sv1 = sv1(startBinInd1:endBinInd1, timeInd1:timeInd2);

binMin2 = binMin2(startBinInd2:endBinInd2);
binMid2 = binMid2(startBinInd2:endBinInd2);
binMax2 = binMax2(startBinInd2:endBinInd2);
dD2 = dD2(startBinInd2:endBinInd2);
ndf2 = ndf2(startBinInd2:endBinInd2, timeInd1:timeInd2);
mdfHabit2 = mdfHabit2(startBinInd2:endBinInd2, timeInd1:timeInd2);
mdfBL2 = mdfBL2(startBinInd2:endBinInd2, timeInd1:timeInd2);
vt2 = vt2(startBinInd2:endBinInd2, timeInd1:timeInd2);
pr2 = pr2(startBinInd2:endBinInd2, timeInd1:timeInd2);
area2 = area2(startBinInd2:endBinInd2, timeInd1:timeInd2);
ar2 = ar2(startBinInd2:endBinInd2, timeInd1:timeInd2);
perim2 = perim2(startBinInd2:endBinInd2, timeInd1:timeInd2);
sv2 = sv2(startBinInd2:endBinInd2, timeInd1:timeInd2);

binMin = [binMin1;binMin2]; binMid = [binMid1;binMid2]; binMax = [binMax1;binMax2];
dD = [dD1;dD2]; sv = [sv1;sv2];
ndf = [ndf1;ndf2]; mdfHabit = [mdfHabit1;mdfHabit2]; mdfBL = [mdfBL1;mdfBL2];
area = [area1;area2];

%% Generate 1-sec PSDs

timeArray = NaN(length(time),1);
for i=1:length(time)
    timeArray(i) = datenum(strcat(date, sprintf('%06d',time(i))),...
        'yyyymmddHHMMSS'); % convert time string to MATLAB date format  
end
timeArray_pcolor = [timeArray; timeArray(end)+1/86400]; % time array for pcolor plotting
ndf_pcolor = padarray(ndf, [1,1], 'post'); % N(D) padded for pcolor plotting
ndf_pcolor(ndf_pcolor==0) = NaN; % change 0 values to NaN

%% Generate 5-sec Average PSDs

timeArray_5sec = NaN(length(time)/5,1);
ndf_5sec = NaN(length(dD1)+length(dD2), floor(length(time)/5));
for i=1:floor(length(time)/5)
    timeArray_5sec(i) = datenum(strcat(date, sprintf('%06d',time(5*i-4))),...
        'yyyymmddHHMMSS'); % convert time string to MATLAB date format
    ndf_5sec(:,i) = nansum(ndf(:,5*i-4:5*i),2)/5;
end
timeArray_5sec_pcolor = [timeArray_5sec; timeArray_5sec(end)+5/86400]; % time array for pcolor plotting
ndf_5sec_pcolor = padarray(ndf_5sec, [1,1], 'post'); % N(D) padded for pcolor plotting
ndf_5sec_pcolor(ndf_5sec_pcolor==0) = NaN; % change 0 values to NaN

%% Generate Average SDs for Entire Period

ndf_avg = nansum(ndf,2) / length(time);
mdfHabit_avg = nansum(mdfHabit,2) / length(time);
mdfBL_avg = nansum(mdfBL,2) / length(time);
area_avg = nansum(area,2) / length(time);


%% Calculate Bulk Variables

% N for various particle size ranges
N = nansum(ndf.*repmat(dD/10,[1 length(time)]),1);

tempInds = find((binMin>=0.1) & (binMax<=0.4)); % indices for 100 < D < 400 um
N100 = nansum(ndf(tempInds,:).*repmat(dD(tempInds)/10,[1 length(time)]),1);

tempInds = find((binMin>=0.4) & (binMax<=1.0)); % indices for 400 < D < 1000 um
N400 = nansum(ndf(tempInds,:).*repmat(dD(tempInds)/10,[1 length(time)]),1);

tempInds = find(binMin>=1.0); % indices for D > 1000 um
N1000 = nansum(ndf(tempInds,:).*repmat(dD(tempInds)/10,[1 length(time)]),1);

% IWC following Jackson et al. (2012) habit m-D relations [g m^-3]
iwcHabit = (10^6)*nansum(mdfHabit.*repmat(dD,[1 length(time)])/10, 1);

% IWC following Baker & Lawson (2006) m-A relation [g m^-3]
iwcBL = (10^6)*nansum(mdfBL.*repmat(dD,[1 length(time)])/10,1);

% Dmm following Baker & Lawson (2006) m-A relation [mm]
massCumul = cumsum(mdfBL.*sv.*repmat(dD,[1 length(time)])/10, 1); % cumulative mass [g]
Dmm = NaN(1, length(time));
for i=1:length(time)
    if massCumul(end,i)>0
        if massCumul(1,i)>=0.5*massCumul(end,i)
            Dmm(i) = binMid(1);
        else
            Dmm(i) = binMid(find(massCumul(:,i)>0.5*massCumul(end,i),1,'first')-1);
        end
    end
end

%% Plot the data

% Time series N(D)
fig1 = figure('Position', [100, 100, 900, 600]); set(gcf, 'Color', 'w');
p(1) = subplot(2,2,[1 3]);
h = pcolor(timeArray_5sec_pcolor, [binMin; binMax2(end)], log10(ndf_5sec_pcolor)); hold on;
ln1 = plot(timeArray, Dmm, 'k', 'LineWidth', 2);
set(h,'edgecolor','none'); datetick('x', 'keeplimits'); ylim([10^-1 10]);
set(gca,'TickLength', 2*get(gca,'TickLength'));
set(gca, 'YScale', 'log', 'XMinorTick', 'on', 'YMinorTick', 'on', 'Layer', 'top');
colormap(jet); hcb = colorbar; caxis([-6 -1]);
xlabel('Time (UTC)'); ylabel('D_{max} (mm)'); ylabel(hcb,'log_{10}[N(D_{max})] (cm^{-4})')
leg1 = legend(ln1, {'D_{mm}'}, 'Location', 'northwest'); set(leg1, 'box', 'off');

% Time series N
p(2) = subplot(2,2,2);
ln2 = plot(timeArray, 1000*N, 'k', 'LineWidth', 2); hold on;
ln3 = plot(timeArray, 1000*N100, 'r', 'LineWidth', 2); hold on;
ln4 = plot(timeArray, 1000*N400, 'b', 'LineWidth', 2); hold on;
ln5 = plot(timeArray, 1000*N1000, 'g', 'LineWidth', 2);
datetick('x', 'keeplimits'); ylim([0 10]);
xlabel('Time (UTC)'); ylabel('N (L^{-1})');
leg2 = legend({'N_{>100}', 'N_{100-400}', 'N_{400-1000}', 'N_{>1000}'},...
    'Location', 'northeast', 'fontsize', 8); set(leg2, 'box', 'off');

% Averaged N(D)
p(3) = subplot(2,2,4);
st1 = stairs(binMin, ndf_avg, 'k', 'LineWidth', 2);
xlim([0.1 10]);
set(gca, 'XScale', 'log', 'YScale', 'log', 'XMinorTick', 'on', 'YMinorTick', 'on');
xlabel('D_{max} (mm)'); ylabel('N(D_{max}) (cm^{-4})');

print([fileDirectory, 'ConcentrationPlots_' , date, '.jpg'],'-djpeg','-r300')

% % Time series IWC
% subplot(2,2,4);
% plot(timeArray, iwcHabit, 'r'); hold on; plot(timeArray, iwcBL, 'b');
% xlabel('Time (UTC)'); ylabel('IWC (g m^{-3})');
% datetick('x', 'keeplimits'); ylim([0.2 0.8]);
% leg3 = legend({'Habit-dependent m-D', 'm-A following BL06'},...
%     'Location', 'southwest', 'fontsize', 10); set(leg3, 'box', 'off');

% Averaged M(D)
fig2 = figure('Position', [100, 100, 900, 600]); set(gcf, 'Color', 'w');
p(4) = subplot(1,2,1);
st2 = stairs(binMin, mdfHabit_avg, 'r', 'LineWidth', 2); hold on;
st3 = stairs(binMin, mdfBL_avg, 'b', 'LineWidth', 2);
xlim([0.1 10]);
set(gca, 'XScale', 'log', 'YScale', 'log', 'XMinorTick', 'on', 'YMinorTick', 'on');
xlabel('D_{max} (mm)'); ylabel('M(D_{max}) (g cm^{-4})');
leg3 = legend({'Habit-dependent m-D', 'm-A following BL06'},...
    'Location', 'northeast', 'fontsize', 10); set(leg3, 'box', 'off');

% Averaged A(D)
p(5) = subplot(1,2,2);
st4 = stairs(binMin, area_avg, 'k', 'LineWidth', 2); hold on;
xlim([0.1 10]);
set(gca, 'XScale', 'log', 'YScale', 'log', 'XMinorTick', 'on', 'YMinorTick', 'on');
xlabel('D_{max} (mm)'); ylabel('A(D_{max}) (mm^2 cm^{-4})');

print([fileDirectory, 'MassAreaPlots_' , date, '.jpg'],'-djpeg','-r300')
