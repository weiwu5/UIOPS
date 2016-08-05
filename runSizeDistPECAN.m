function runSizeDistPECAN(ddate)

rCIP = 0;

%% Run CIP
if rCIP == 1 
clearvars -except ddate
inst=num2str(ddate);
%inst = '20150617';
tasfilename =  '01CIP20150702015858.csv'; %ls('--color=none',['../PECAN/' num2str(ddate) '/01CIP*.csv']);
dirpath = pwd;
tasfilename = strtrim([dirpath '/' tasfilename]);
loadTASinfo

tas = True_Air_Speed;
timehhmmss=insec2hhmmss(Time);
pres=Static_Press+Diff_Press;
temp1=Ambient_Temp;

probe='CIP';
% Specify the desination file name infilename and diodesize ds in millimeters
outFile = ['sdistCI.' inst '.' probe '.cdf'];                  	% Input uncompressed netCDF file name
% Output file name
inFile =  ['proc2.' inst '.' probe '.cdf'];
inst
sizeDist(inFile,outFile, tas, floor(timehhmmss),probe, 6, 0, pres, temp1);
end
%% Run PIP
clearvars -except ddate
inst=num2str(ddate);
%inst = '20150617';
tasfilename = '00PIP20150620005009.csv'; %ls('--color=none',['../PECAN/' num2str(ddate) '/00PIP*.csv']); %'00PIP20150617005743.csv';
dirpath = pwd;
tasfilename = strtrim([dirpath '/' tasfilename]);
loadTASinfo

tas = True_Air_Speed;
timehhmmss=insec2hhmmss(Time);
pres=Static_Press+Diff_Press;
temp1=Ambient_Temp;

probe='PIP';
% Specify the desination file name infilename and diodesize ds in millimeters
outFile = ['sdistCI.' inst '.' probe '.cdf'];                  	% Input uncompressed netCDF file name
% Output file name
inFile =  ['proc2.TEST.' probe '.cdf'];

sizeDist(inFile,outFile, tas, floor(timehhmmss),probe, 6, 0, pres, temp1,'PECAN1',ddate);


end

