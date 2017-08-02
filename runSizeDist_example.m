function runSizeDistPECAN(ddate,runCIP,runPIP)
	
	% Added this line as MATLAB would not ls properly in -nodisplay mode without it...
	cd /data/pecan/a/stechma2/pecan/mp-data/IProcessingRelease
	
	%runCIP = 1;
	%runPIP = 1;
	
	projectname = 'PECAN';
		
	datadate=num2str(ddate);
	
	prmsF = ['/data/pecan/a/stechma2/pecan/' datadate '_PECANparams.nc'];
	
	%% Run CIP
	if runCIP
		clearvars -except ddate runCIP runPIP projectname datadate prmsF

		tasfilename =  ls('--color=none',['../' num2str(ddate) '/01CIP*.csv']);
		dirpath = pwd;
		cipTASf = strtrim([dirpath '/' tasfilename]);
		loadCIPcsv
		
		tas = CIP_True_Air_Speed;
		timehhmmss=insec2hhmmss(CIP_Time);
		pres=CIP_Static_Press+CIP_Diff_Press;
		temp1=CIP_Ambient_Temp;
		
		probe='CIP';
		
		% Specify the input filename, output filename, true air speed, time (hhmmss), probe name,
		% Dmax definition, surface area method, pressure and temperature
		outFile = ['sdistCI.' datadate '.' probe '.cdf']; % Uncompressed netCDF file name
		inFile =  ['proc2.' datadate '.' probe '.cdf'];
		
		sizeDist(inFile,outFile, tas, floor(timehhmmss),probe, 6, 0, pres, temp1, projectname, datadate, 'NONE', prmsF);
	end
	
	%% Run PIP
	if runPIP
		clearvars -except ddate runCIP runPIP projectname datadate prmsF

		inst=num2str(ddate);
		tasfilename = ls('--color=none',['../' num2str(ddate) '/00PIP*.csv']); %'00PIP20150617005743.csv';
		dirpath = pwd;
		pipTASf = strtrim([dirpath '/' tasfilename]);
		loadPIPcsv
		
		tas = PIP_True_Air_Speed;
		timehhmmss=insec2hhmmss(PIP_Time);
		pres=PIP_Static_Press+PIP_Diff_Press;
		temp1=PIP_Ambient_Temp;
		
		probe='PIP';
		
		% Specify the input filename, output filename, true air speed, time (hhmmss), probe name,
		% Dmax definition, surface area method, pressure and temperature
		outFile = ['sdistCI.' inst '.' probe '.cdf'];                  	% Input uncompressed netCDF file name
		inFile =  ['proc2.' inst '.' probe '.cdf'];
		
		sizeDist(inFile,outFile, tas, floor(timehhmmss),probe, 6, 0, pres, temp1, projectname, datadate, 'NONE', prmsF);
	end
end