%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Derive the area and size distribution for entire-in particles 
%  Include the IWC calculation
%  Include the effective radius 
%               Created by Will Wu, 09/18/2013
%
%	**************************
%	*** Modification Notes ***
%	**************************
%	* Modified to use the new maximum size and derive both maximum size distribution
%    and area-equivalent size distribution.  
%			Will Wu, 10/26/2013  
%	* Modified to calculate terminal velocity using Heymsfield and Westbrook (2010) method
%    and precipitation rate.     
%			Will Wu, 01/15/2014
%	* Modified to include mass size distribution with habit info. 
%			Will Wu, 02/09/2014
%	* Modified to include particle area using A-D relations. 
%			Will Wu, 02/14/2014 
%	* Special Edition for Boston Cloud workshop. 
%			Wei Wu, 04/01/2014
%	* Gneralized as a new sorting function for all probes. 
%			Wei Wu, 07/25/2014
%	* Modified to allow the option to ingest/use interarrival time dynamic threshold
%			Dan Stechman, 05/06/2016
%	* Added project and date specific capabilities (including spiral-dependent interarrival
%    thresholding). Also cleaned up code and improved efficiency in places.
%			Dan Stechman, 06/03/2016
%	* Added shatter removal using array of interarrival time thresholds (either constant or varrying [e.g., different threshold for
%	 each spiral in PECAN project]). Also added experimental shatter reacceptance option to allow for potential diffraction fringes
%    originally flagged as shattered to be reaccepted.
%			Dan Stechman, 06/09/2016
%
%  Usage: 
%    infile:   Input filename, string
%    outfile:       Output filename, string
%    tas:           True air speed, double array
%    timehhmmss:    Time in hhmmss format, double array
%    probename:     Should be one of 'HVPS', 'CIP', 'PIP', '2DC', '2DP', 'F2DC' 
%    d_choice:      the definition of Dmax, should use 6 usually. [1-6] 
%    SAmethod:      0: Center in; 1: Entire in; 2: With Correction
%    Pres:          1 second pressure data
%    Temp:          1 second temperature data
%    projectname:   Project name, string
%    ddate:         Date to be analyzed, string (YYYYMMDD)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sizeDist(infile, outfile, tas, timehhmmss, probename, d_choice, SAmethod, Pres, Temp, projectname, ddate)
iCreateBad = 0; % Default not to output bad particles PSDs and other info
iCreateAspectRatio = 0; % Default not to process aspect ratio info
%% Interarrival threshold file specification
% Can be implemented if a time-dependent threshold is required - add 'varargin' to arguments in function header above
%{
% if length(varargin) == 1
% 	iaThreshFile = varargin{1};
% else
% 	display(['You screwed up'])
% 	iaThreshFile = 'NONE';
% end
%}

%% Define input and output files and initialize time variable
f = netcdf.open(infile,'nowrite');
mainf = netcdf.create(outfile, 'clobber');

% tas_char = num2str(timehhmmss); %Unused
tas_time = floor(timehhmmss/10000)*3600+floor(mod(timehhmmss,10000)/100)*60+floor(mod(timehhmmss,100));
% averaging_time = 1;

%% Project-, probe-, and date-specific information
switch projectname
    case 'PECAN'
		switch probename
            case 'CIP'
                num_diodes =64;
                diodesize = 0.025; % units of mm
                armdst=100.;
%                 num_bins = 64;
%                 kk=diodesize/2:diodesize:(num_bins+0.5)*diodesize;
                num_bins=19;
                kk=[50.0   100.0   150.0   200.0   250.0   300.0   350.0   400.0   475.0   550.0   625.0 ...
                    700.0   800.0   900.0  1000.0  1200.0  1400.0  1600.0  1800.0  2000.0]/1000; %Array in microns - converted to mm
                probetype=1;
                tasMax=200; % Max airspeed that can be sampled without under-sampling (images would appear skewed)
                
				applyIntArrThresh = 1;
					defaultIntArrThresh = 1e-5;
				reaccptShatrs = 1;
					reaccptD = 0.5; % Diammeter (in mm) to reaccept if initially flagged as shattered
					reaccptMaxIA = 2.5e-7;	% Max interarrival time in seconds a particle can have to be reaccepted if 
											% size criteria are met. Possible definition of this is the time of one slice, so in
											% this case, with an airspeed of ~100 m/s and a slice of 25 um, this would be 2.5e-7.
                

				% Get start and end times (in seconds) of spirals; interarrival time thresholds for each spiral
				[startT, endT, ~, ~, intar_threshold_spirals] = getPECANparams(ddate, probename);
				
				intar_threshold = ones(size(tas_time))*defaultIntArrThresh;
				for ix = 1:length(tas_time)
					for iz = 1:length(startT)
						if (tas_time(ix) >= startT(iz) && tas_time(ix) < endT(iz))
							intar_threshold(ix) = intar_threshold_spirals(iz);
						end
					end
				end
                
            case 'PIP'
                num_diodes =64;
                diodesize = 0.1; %units of mm
                armdst=260.;
                num_bins = 64;
%                 kk=diodesize/2:diodesize:(num_bins+0.5)*diodesize;
                kk=diodesize/2:diodesize:(num_bins+0.6)*diodesize;
%                 num_bins=19;
%                 kk=[50.0   100.0   150.0   200.0   250.0   300.0   350.0   400.0   475.0   550.0   625.0 ...
%                    700.0   800.0   900.0  1000.0  1200.0  1400.0  1600.0  1800.0  2000.0]*4/1000;
                probetype=1;
                tasMax=200; 
                
				applyIntArrThresh = 1;
					defaultIntArrThresh = 1e-5;
				reaccptShatrs = 1;
					reaccptD = 0.5; % Diammeter (in mm) to reaccept if initially flagged as shattered
					reaccptMaxIA = 1e-6; % (Slice size [m])/(avg. airspeed [m/s])
                
				% Get start and end times (in seconds) of spirals; interarrival time thresholds for each spiral
				[startT, endT, ~, ~, intar_threshold_spirals] = getPECANparams(ddate, probename);
				
				intar_threshold = ones(size(tas_time))*defaultIntArrThresh;
				for ix = 1:length(tas_time)
					for iz = 1:length(startT)
						if (tas_time(ix) >= startT(iz) && tas_time(ix) < endT(iz))
							intar_threshold(ix) = intar_threshold_spirals(iz);
						end
					end
				end
		end

    otherwise
        switch probename
            case 'HVPS'
                % For the HVPS
                num_diodes =128;
                diodesize = .150;
                armdst=161.;
                num_bins = 28;
                kk=[200.0   400.0   600.0   800.0  1000.0  1200.0  1400.0  1600.0  1800.0  2200.0  2600.0 ...
                     3000.0  3400.0  3800.0  4200.0  4600.0  5000.0  6000.0  7000.0  8000.0  9000.0 10000.0 ...
                     12000.0 14000.0 16000.0 18000.0 20000.0 25000.0 30000.0]/1000;
                %num_bins =128;
                %kk=diodesize/2:diodesize:(num_bins+0.5)*diodesize;
                probetype=2;
                tasMax=170; 
				
				% Interarrival threshold and reaccept max interarrival time are often flight-/instrument-specific
				% **Values here may not be correct** 
				% The interarrival threshold can be modifided to change second-by-second if desired
                applyIntArrThresh = 0;
					defaultIntArrThresh = 4e-6;
				reaccptShatrs = 0;
					reaccptD = 0.5; 
					reaccptMaxIA = 1e-6; % (Slice size [m])/(avg. airspeed [m/s])
                
				intar_threshold = ones(size(tas_time))*defaultIntArrThresh;

            case '2DS'
                % For the HVPS
                num_diodes =128;
                diodesize = .010;
                armdst=63.;
                %num_bins = 28;
                %kk=[200.0   400.0   600.0   800.0  1000.0  1200.0  1400.0  1600.0  1800.0  2200.0  2600.0 ...
                %     3000.0  3400.0  3800.0  4200.0  4600.0  5000.0  6000.0  7000.0  8000.0  9000.0 10000.0 ...
                %     12000.0 14000.0 16000.0 18000.0 20000.0 25000.0 30000.0]/1000/15;
                num_bins =128;
                %kk=diodesize/2:diodesize:(num_bins+0.5)*diodesize;
                kk=diodesize/2:diodesize:(num_bins+0.6)*diodesize;
                probetype=2;
                tasMax=170;  
				
				% Interarrival threshold and reaccept max interarrival time are often flight-/instrument-specific
				% **Values here may not be correct** 
				% The interarrival threshold can be modifided to change second-by-second if desired
                applyIntArrThresh = 0;
					defaultIntArrThresh = 4e-6;
				reaccptShatrs = 0;
					reaccptD = 0.5; 
					reaccptMaxIA = 1e-6; % (Slice size [m])/(avg. airspeed [m/s])
                
				intar_threshold = ones(size(tas_time))*defaultIntArrThresh;

            case 'CIP'
                % For the CIP 
                num_diodes =64;
                diodesize = .025; %units of mm
                armdst=100.;
                %num_bins = 64;
                %kk=diodesize/2:diodesize:(num_bins+0.5)*diodesize;
                num_bins=19;
                kk=[50.0   100.0   150.0   200.0   250.0   300.0   350.0   400.0   475.0   550.0   625.0 ...
                    700.0   800.0   900.0  1000.0  1200.0  1400.0  1600.0  1800.0  2000.0]/1000; %Array in microns - converted to mm
                probetype=1;
                tasMax=200; % Max airspeed that can be sampled without under-sampling (images would appear skewed)
				
				% Interarrival threshold and reaccept max interarrival time are often flight-/instrument-specific
				% **Values here may not be correct** 
				% The interarrival threshold can be modifided to change second-by-second if desired
                applyIntArrThresh = 0;
					defaultIntArrThresh = 1e-5;
				reaccptShatrs = 0;
					reaccptD = 0.5; 
					reaccptMaxIA = 1e-6; % (Slice size [m])/(avg. airspeed [m/s])
                
				intar_threshold = ones(size(tas_time))*defaultIntArrThresh;

            case 'PIP'
                num_diodes =64;
                diodesize = .1; %units of mm
                armdst=260.;
                num_bins = 64;
%                 kk=diodesize/2:diodesize:(num_bins+0.5)*diodesize;
                kk=diodesize/2:diodesize:(num_bins+0.6)*diodesize;
%                 num_bins=19;
%                 kk=[50.0   100.0   150.0   200.0   250.0   300.0   350.0   400.0   475.0   550.0   625.0 ...
%                    700.0   800.0   900.0  1000.0  1200.0  1400.0  1600.0  1800.0  2000.0]*4/1000;
                probetype=1;
                tasMax=200; 
                
				% Interarrival threshold and reaccept max interarrival time are often flight-/instrument-specific
				% **Values here may not be correct** 
				% The interarrival threshold can be modifided to change second-by-second if desired
                applyIntArrThresh = 0;
					defaultIntArrThresh = 1e-5;
				reaccptShatrs = 0;
					reaccptD = 0.5; 
					reaccptMaxIA = 1e-6; % (Slice size [m])/(avg. airspeed [m/s])
                
				intar_threshold = ones(size(tas_time))*defaultIntArrThresh;


            case '2DC'
                % For the 2DC
                num_diodes =32;
                diodesize = .03; %.025;
                armdst=61.;
                %num_bins = 32;
                %kk=diodesize/2:diodesize:(num_bins+0.5)*diodesize;
                num_bins=19;
                kk=[50.0   100.0   150.0   200.0   250.0   300.0   350.0   400.0   475.0   550.0   625.0 ...
                    700.0   800.0   900.0  1000.0  1200.0  1400.0  1600.0  1800.0  2000.0]/1000;
                probetype=0;
                tasMax=125;  
				
				% Interarrival threshold and reaccept max interarrival time are often flight-/instrument-specific
				% **Values here may not be correct** 
				% The interarrival threshold can be modifided to change second-by-second if desired
                applyIntArrThresh = 0;
					defaultIntArrThresh = 4e-6;
				reaccptShatrs = 0;
					reaccptD = 0.5; 
					reaccptMaxIA = 1e-6; % (Slice size [m])/(avg. airspeed [m/s])
                
				intar_threshold = ones(size(tas_time))*defaultIntArrThresh;

            case '2DP'
                % For the 2DP
                num_diodes =32;
                diodesize = .200; %.025;
                armdst=260.; %75.77; %61.;
                %num_bins = 32;
                %kk=diodesize/2:diodesize:(num_bins+0.5)*diodesize;
                num_bins=19;
                kk=[50.0   100.0   150.0   200.0   250.0   300.0   350.0   400.0   475.0   550.0   625.0 ...
                    700.0   800.0   900.0  1000.0  1200.0  1400.0  1600.0  1800.0  2000.0]*8/1000;
                probetype=0;
				
				% Interarrival threshold and reaccept max interarrival time are often flight-/instrument-specific
				% **Values here may not be correct** 
				% The interarrival threshold can be modifided to change second-by-second if desired
                applyIntArrThresh = 0;
					defaultIntArrThresh = 4e-6;
				reaccptShatrs = 0;
					reaccptD = 0.5; 
					reaccptMaxIA = 1e-6; % (Slice size [m])/(avg. airspeed [m/s])
                
				intar_threshold = ones(size(tas_time))*defaultIntArrThresh;

            case 'F2DC'
                % For the 2DC
                num_diodes =64;
                diodesize = .025; %.025;
                armdst=61.; %60; %
                %num_bins = 32;
                %kk=diodesize/2:diodesize:(num_bins+0.5)*diodesize;
                num_bins=19;
                kk=[50.0   100.0   150.0   200.0   250.0   300.0   350.0   400.0   475.0   550.0   625.0 ...
                    700.0   800.0   900.0  1000.0  1200.0  1400.0  1600.0  1800.0  2000.0]/1000;
                probetype=0;
				
				% Interarrival threshold and reaccept max interarrival time are often flight-/instrument-specific
				% **Values here may not be correct** 
				% The interarrival threshold can be modifided to change second-by-second if desired
                applyIntArrThresh = 0;
					defaultIntArrThresh = 4e-6;
				reaccptShatrs = 0;
					reaccptD = 0.5; 
					reaccptMaxIA = 1e-6; % (Slice size [m])/(avg. airspeed [m/s])
                
				intar_threshold = ones(size(tas_time))*defaultIntArrThresh;
        end
end

if applyIntArrThresh && ~reaccptShatrs
	fprintf('Beginning sizeDist_Paris.m for %s %s - %s probe\n\t**Optional parameters active:\n\t- Shatter removal\n\n',projectname,ddate,probename);
elseif applyIntArrThresh && reaccptShatrs
	fprintf('Beginning sizeDist_Paris.m for %s %s - %s probe\n\t**Optional parameters active:\n\t- Shatter removal\n\t- Shatter reacceptance\n\n',...
		projectname,ddate,probename);
else
	fprintf('Beginning sizeDist_Paris.m for %s %s - %s probe\n\n',projectname,ddate,probename);
end

res=diodesize*1000;
binwidth=diff(kk);
% SAmethod = 2;
% for i=1:num_bins+1
%     kk(i)=  (diodesize*i)^2*3.1415926/4; 
% end


%% Define Variables

% Good particles (not rejected)
particle_dist_minR  = zeros(length(tas),num_bins)*NaN;
particle_dist_AreaR = zeros(length(tas),num_bins)*NaN;
particle_aspectRatio = zeros(length(tas),num_bins)*NaN;
particle_aspectRatio1 = zeros(length(tas),num_bins)*NaN;
particle_areaRatio1 = zeros(length(tas),num_bins)*NaN;
particle_area = zeros(length(tas),num_bins)*NaN;
cip2_meanp = zeros(length(tas),num_bins)*NaN;
cip2_iwc = zeros(length(tas),num_bins)*NaN;
cip2_iwcbl = zeros(length(tas),num_bins)*NaN;
cip2_vt = zeros(length(tas),num_bins)*NaN;
cip2_pr = zeros(length(tas),num_bins)*NaN;
cip2_partarea = zeros(length(tas),num_bins)*NaN;

cip2_re = zeros(1,length(tas))*NaN;
good_partpercent=zeros(length(tas),1);
goodintpercent=zeros(length(tas),1);
numGoodparticles=zeros(length(tas),1);
one_sec_ar=zeros(length(tas),1);


cip2_habitsd = zeros(length(tas),num_bins,10)*NaN;
cip2_habitmsd = zeros(length(tas),num_bins,10)*NaN;
area_dist2 = zeros(length(tas),num_bins,10)*NaN;

rejectpercentbycriterion=zeros(length(tas),14);


% Bad particles (rejected)
bad_particle_dist_minR  = zeros(length(tas),num_bins)*NaN;
bad_particle_dist_AreaR = zeros(length(tas),num_bins)*NaN;
bad_particle_aspectRatio = zeros(length(tas),num_bins)*NaN;
bad_particle_aspectRatio1 = zeros(length(tas),num_bins)*NaN;
bad_particle_areaRatio1 = zeros(length(tas),num_bins)*NaN;
bad_particle_area = zeros(length(tas),num_bins)*NaN;
bad_cip2_meanp = zeros(length(tas),num_bins)*NaN;
bad_cip2_iwc = zeros(length(tas),num_bins)*NaN;
bad_cip2_iwcbl = zeros(length(tas),num_bins)*NaN;
bad_cip2_vt = zeros(length(tas),num_bins)*NaN;
bad_cip2_pr = zeros(length(tas),num_bins)*NaN;
bad_cip2_partarea = zeros(length(tas),num_bins)*NaN;

bad_cip2_re = zeros(1,length(tas))*NaN;
badintpercent=zeros(length(tas),1);
numBadparticles=zeros(length(tas),1);
bad_one_sec_ar=zeros(length(tas),1);

bad_cip2_habitsd = zeros(length(tas),num_bins,10)*NaN;
bad_cip2_habitmsd = zeros(length(tas),num_bins,10)*NaN;
bad_area_dist2 = zeros(length(tas),num_bins,10)*NaN;


% particle_dist2 = zeros(length(tas),num_bins)*NaN; %Unused
% time_interval1 = zeros(length(tas), 1); %Unused
% cip2_ar = zeros(1,length(tas))*NaN; %Unused
% throwoutpercent=zeros(length(tas),1); %Used in legacy interarrival time analysis
% totalint=zeros(length(tas),1); %Used in legacy interarrival time analysis
% intsum=zeros(length(tas),1); %Used in legacy interarrival time analysis

area_bins = 0:.1:1.01;
one_sec_times = tas_time;
one_sec_dur = length(one_sec_times);
total_one_sec_locs(1:one_sec_dur) = 0;
start_time = floor(tas_time(1));
end_time = ceil(tas_time(end));
one_sec_tas(1:one_sec_dur) = 0;
one_sec_tas_entire(1:one_sec_dur) = 0;
deadtime(1:one_sec_dur) = 0;

warning off all

one_sec_times=[one_sec_times;one_sec_times(one_sec_dur)+1];
time_interval2 = zeros(one_sec_dur,1);

TotalPC1 = zeros(one_sec_dur,1)';
TotalPC2 = zeros(one_sec_dur,1)';

% Used for debugging of interarrival time analysis
shatrReject_times = [];
shatrReject_intArr = [];
shatrReject_diam = [];

rccptReject_times = [];
rccptReject_intArr = [];
rccptReject_diam = [];

loopedTimes = [];
loopedIntArr = [];
loopedDiam = [];
loopedAutoRej = [];

%% Load particles for each second, and then process them
% Only for large files cannot be processed at once
[~, NumofPart] = netcdf.inqDim(f,0); % Check the number of particles

if 1==probetype
   image_time_hhmmssall = netcdf.getVar(f,netcdf.inqVarID(f,'particle_time'));
elseif 2==probetype
   image_time_hhmmssallbuffer = netcdf.getVar(f,netcdf.inqVarID(f,'Time'));
%    image_time_hhmmssallbuffer(image_time_hhmmssallbuffer<10000 & image_time_hhmmssallbuffer>=0)=image_time_hhmmssallbuffer(image_time_hhmmssallbuffer<10000 & image_time_hhmmssallbuffer>=0)+240000;
   alltimeinseconds = netcdf.getVar(f,netcdf.inqVarID(f,'Time_in_seconds'));
   time_msec_all = netcdf.getVar(f,netcdf.inqVarID(f,'msec'),0,1);
  
   indexRollback=find(diff(alltimeinseconds)<-250)+1;
   
   for i=1:length(indexRollback)
       alltimeinseconds(indexRollback(i):end)=alltimeinseconds(indexRollback(i):end)+(2^32-1)*(res/10^6/tasMax);
   end
   
%    alltimeinsecondsstart=alltimeinseconds(indexBuffert);
%    increaseAllinseconds= alltimeinseconds-alltimeinseconds(1);
%    increaseAllinseconds(increaseAllinseconds<0)=increaseAllinseconds(increaseAllinseconds<0)+(2^32-1)*(res/10^6/170);
%    image_time_hhmmssall = insec2hhmmss(floor(47069+time_msec_all(1)/1000.0+increaseAllinseconds*170/110));
   image_time_hhmmssall = image_time_hhmmssallbuffer;
else
   image_time_hhmmssall = netcdf.getVar(f,netcdf.inqVarID(f,'Time'));
end

% image_time_hhmmssall = netcdf.getVar(f,netcdf.inqVarID(f,'Time'));
% image_time_hhmmssall(image_time_hhmmssall<50000 & image_time_hhmmssall>=0)=image_time_hhmmssall(image_time_hhmmssall<50000 & image_time_hhmmssall>=0)+120000;



% Find all indices (true/1) with a unique time in hhmmss - in other words, we're getting the particle index where each new
% one-second period starts
startindex=[true;(diff(hhmmss2insec(image_time_hhmmssall))>0)]; % & diff(hhmmss2insec(image_time_hhmmssall))<5)]; % Simplified (tested/changed by DS)

% startindex=int8(image_time_hhmmssall*0);
% for i=1:length(timehhmmss)
%    indexofFirstTime =  find(image_time_hhmmssall==timehhmmss(i),1);
%    if ( ~isempty(indexofFirstTime) )
%        startindex(indexofFirstTime)=1;
%    end
%    disp([i,length(timehhmmss)]);
% end


% Get the start time for each new one-second period
starttime=image_time_hhmmssall(startindex); % Simplified (tested/changed by DS)


% Find all instances where startindex is true (where image_time_hhmmssall changes by more than 0) and shift indices back by one to
% facilitate proper particle counts for each one-second period
start_all=find(startindex)-1; % Simplified (tested/changed by DS)

% Sort the particle one-second time array in the event it is out of order and redefine the start_all variable as needed
[starttime,newindexofsort]=sort(starttime);
start_all=start_all(newindexofsort);


%% Add 240000 to particle when their time is past midnight
%starttime(starttime<10000 & starttime>=0)=starttime(starttime<10000 & starttime>=0)+240000;

%% Remove times when there is no tas data available
% nNoTAS=0;
% for i=1:length(starttime)
%     if isempty(timehhmmss(timehhmmss == starttime(i)))
%         starttime(i)=500000;
%         nNoTAS=nNoTAS+1;
%     end
%     
%     if i>5 & i<length(starttime)-5 & hhmmss2insec(starttime(i))>mean(hhmmss2insec(starttime(i-5:i+5)))+5
%         starttime(i)=500000;
%     end
% end
% nNoTAS
% start_all = start_all(starttime<500000);
% count_all = count_all(starttime<500000);
% starttime = starttime(starttime<500000);

%% Remove any duplicate times and determine how many particles are present in each one-second period
fprintf('Number of duplicate times = %d\n\n',(length(starttime)-length(unique(starttime))));

[starttime, ia, ~] = unique(starttime,'first');
start_all = start_all(ia);
count_all= [diff(start_all); NumofPart-start_all(end)];
count_all(count_all<0)=1;

%% Remove times when there are less than 10 particles in one second
% starttime = starttime(count_all>10);
% start_all = start_all(count_all>10);
% count_all = count_all(count_all>10);

%if (int32(timehhmmss(1))>int32(starttime(2)))
%    error('Watch Out for less TAS time from begining!')
%end

%% Main loop over the length of the true air speed variable (1-sec resolution)
jjj=1;

sumIntArrGT1 = 0;
intArrGT1 = [];

% nThrow11=0; % Used in legacy interarrival time analysis
% maxRecNum=1; % Used in legacy interarrival time analysis

fprintf('Beginning size distribution calculations and sorting %s\n\n',datestr(now));

for i=1:length(tas) 
    
    if (int32(timehhmmss(i))>=int32(starttime(jjj)))
        
        % Attempt to sync TAS file time (timehhmmss) with particle time
        if (int32(timehhmmss(i))>int32(starttime(jjj)))
            jjj=jjj+1;
            
            if (jjj>length(start_all))
                break;
            end
        end
        
        start=start_all(jjj);
        count=count_all(jjj);
        jjj=min(jjj+1,length(start_all));
        
        % Load autoanalysis parameters. Start at beginning (start) of some one-second period and load the values for every
        % particle in that period (count)
        msec = netcdf.getVar(f,netcdf.inqVarID(f,'particle_millisec'),start,count);
        microsec = netcdf.getVar(f,netcdf.inqVarID(f,'particle_microsec'),start,count);
        auto_reject = netcdf.getVar(f,netcdf.inqVarID(f,'image_auto_reject'),start,count);
        im_width = netcdf.getVar(f,netcdf.inqVarID(f,'image_width'),start,count);
        im_length = netcdf.getVar(f,netcdf.inqVarID(f,'image_length'),start,count);
        area = netcdf.getVar(f,netcdf.inqVarID(f,'image_area'),start,count);
        perimeter = netcdf.getVar(f,netcdf.inqVarID(f,'image_perimeter'),start,count);
%         rec_nums = netcdf.getVar(f,netcdf.inqVarID(f,'parent_rec_num'),start,count); %Used in legacy interarrival time analysis
%         top_edges = netcdf.getVar(f,netcdf.inqVarID(f,'image_max_top_edge_touching'),start,count); %Unused
%         bot_edges = netcdf.getVar(f,netcdf.inqVarID(f,'image_max_bottom_edge_touching'),start,count); %Unused
%         longest_y = netcdf.getVar(f,netcdf.inqVarID(f,'image_longest_y'),start,count); %Unused
        size_factor = netcdf.getVar(f,netcdf.inqVarID(f,'size_factor'),start,count);
        habit1 = netcdf.getVar(f,netcdf.inqVarID(f,'holroyd_habit'),start,count);
        centerin = netcdf.getVar(f,netcdf.inqVarID(f,'image_center_in'),start,count);
        entirein = netcdf.getVar(f,netcdf.inqVarID(f,'image_touching_edge'),start,count);
        
        particle_diameter_AreaR = netcdf.getVar(f,netcdf.inqVarID(f,'image_diam_AreaR'),start,count);
        particle_diameter_AreaR = particle_diameter_AreaR * diodesize;

        Time_in_seconds = netcdf.getVar(f,netcdf.inqVarID(f,'Time_in_seconds'),start,count);
%         SliceCount = netcdf.getVar(f,netcdf.inqVarID(f,'SliceCount'),start,count); %Unused
        DMT_DOF_SPEC_OVERLOAD = netcdf.getVar(f,netcdf.inqVarID(f,'DMT_DOF_SPEC_OVERLOAD'),start,count);
        Particle_count = netcdf.getVar(f,netcdf.inqVarID(f,'Particle_number_all'),start,count);
        
        if 1==probetype
            auto_reject(DMT_DOF_SPEC_OVERLOAD~=0)='D';
        end
        
        if iCreateAspectRatio == 1
        aspectRatio = netcdf.getVar(f,netcdf.inqVarID(f,'image_RectangleW'),start,count)./netcdf.getVar(f,netcdf.inqVarID(f,'image_RectangleL'),start,count);
        aspectRatio1 = netcdf.getVar(f,netcdf.inqVarID(f,'image_EllipseW'),start,count)./netcdf.getVar(f,netcdf.inqVarID(f,'image_EllipseL'),start,count);
        end
        
        TotalPC1(i)=length(Particle_count);        
        TotalPC2(i)=Particle_count(end)-Particle_count(1);
			
        
        if 0==probetype
            int_arr=Time_in_seconds;
		else
			if start-1 <= 0
				int_arr = [0;diff(Time_in_seconds)];
				int_arr2 = []; %Won't bother reaccepting particles at the beginning or end of dataset				
			else
				Time_in_seconds2 = netcdf.getVar(f,netcdf.inqVarID(f,'Time_in_seconds'),start-1,count+1);
				int_arr = diff(Time_in_seconds2);
				
				if start ~= start_all(end)
					Time_in_seconds3 = netcdf.getVar(f,netcdf.inqVarID(f,'Time_in_seconds'),(start+count)-1,2);						
					int_arr2 = diff(Time_in_seconds3); %Single value describing interarrival time of first particle of next 1-sec period
				else
					int_arr2 = [];
				end
			end
			
			int_arr2(int_arr2<0)=0;
			
			if reaccptShatrs
				if start ~= start_all(end)
					Time_in_seconds4 = netcdf.getVar(f,netcdf.inqVarID(f,'Time_in_seconds'),start,count+1);
					int_arr3 = diff(Time_in_seconds4);  
				else
					Time_in_seconds4 = Time_in_seconds;
					int_arr3 = diff(Time_in_seconds4);  
					int_arr3 = [int_arr3;int_arr3(end)];
				end
				int_arr3(int_arr3<0)=0;
			end
			
        end
%         if 2==probetype 
%           int_arr=int_arr*(res/10^6/170); 
%         end
        
        if 2==probetype
            int_arr(int_arr<-10)=int_arr(int_arr<-10)+(2^32-1)*(res/10^6/tasMax);
        elseif 0==probetype
            int_arr(int_arr<0)=int_arr(int_arr<0)+(2^24-1)*(res/10^6/tasMax);            
        end
            
        if sum(int_arr<0)>0
			fprintf(2,'\nAt index %d number of int_arr < 0: %d\n',i,sum(int_arr<0));
            disp([int_arr(int_arr<0),int_arr(int_arr<0)+(2^32-1)*(res/10^6/tasMax)]);
        elseif sum(int_arr>1)>0
			sumIntArrGT1 = sumIntArrGT1 + sum(int_arr > 1);
			tempLocs = find(int_arr > 1); 
			intArrGT1 = vertcat(intArrGT1,int_arr(tempLocs));
			
%             fprintf(2,'\nAt index %d number of int_arr > 1: %d\n',i,sum(int_arr>1));
%             disp([int_arr(int_arr>1)-(2^32-1)*(res/10^6/tasMax), int_arr(int_arr>1), Time_in_seconds(int_arr>1)/(0.15/(10^3)/170), Time_in_seconds((int_arr>1))/(0.15/(10^3)/170)]);
        end
        
%         auto_reject(int_arr<0 | int_arr>1)='I';
		auto_reject(int_arr<0)='I';
        int_arr(int_arr<0)=0;
%         int_arr(int_arr>1)=0;


%         max_dimension = im_width;
%         max_dimension(im_length>im_width)=im_length(im_length>im_width);

        
        % Size definition chosen based on the d_choice given in the function call
        if 1==d_choice
            particle_diameter_minR = im_length * diodesize; %(im_length+
        elseif 2==d_choice
            particle_diameter_minR = im_width * diodesize; %(im_length+
        elseif 3==d_choice
            particle_diameter_minR = (im_length + im_width)/2 * diodesize; %(im_length+
        elseif 4==d_choice
            particle_diameter_minR = sqrt(im_width.^2+im_length.^2) * diodesize; %(im_length+
        elseif 5==d_choice
            particle_diameter_minR = max(im_width, im_length) * diodesize; %(im_length+
        elseif 6==d_choice
            particle_diameter_minR = netcdf.getVar(f,netcdf.inqVarID(f,'image_diam_minR'),start,count); % * diodesize
        end
        
%         if 1==strcmp('2DC',probename)  % Adjust resolution from 25 to 30
%            particle_diameter_minR=particle_diameter_minR*1.2;
%            area = area*1.44;
%         end



        % Legacy: Added for Paris meeting, 08/25/2014
        % Used in intercomparison with Environment Canada and University of Blaise Pascal
        %{
        diffPartCount=[1;diff(Particle_count)];
        time_interval22(i) = (Time_in_seconds(end)-Time_in_seconds(1));
        time_interval32(i) = sum(int_arr(diffPartCount==1));
        time_interval42(i) = sum(int_arr);
        time_interval52(i) = sum(int_arr(diffPartCount~=1));
        time_interval62(i) = sum(int_arr(DMT_DOF_SPEC_OVERLOAD==0));
        lengthForTemp = im_length * diodesize; 
        particle_diameter_minR(entirein~=0)=lengthForTemp(entirein~=0);
        
        if time_interval22(i)<0
            time_interval22(i)=time_interval22(i)+(2^32-1)*(res/10^6/tasMax); %#ok<*AGROW>
        end
        
        if RejectCriterier==1
           particle_diameter_minR = particle_diameter_minR .* size_factor;
        end
        

        if 1==probetype
            image_time_hhmmss    = netcdf.getVar(f,netcdf.inqVarID(f,'particle_time'),start,count);
            image_time_hhmmssnew = netcdf.getVar(f,netcdf.inqVarID(f,'particle_time'),start,count); 

        elseif 2==probetype
            alltimeinseconds = netcdf.getVar(f,netcdf.inqVarID(f,'Time_in_seconds'),start,count);
            increaseAllinseconds= alltimeinseconds-alltimeinseconds(1);
            increaseAllinseconds(increaseAllinseconds<0)=increaseAllinseconds(increaseAllinseconds<0)+(2^32-1)*(res/10^6/170);
            image_time_hhmmss = floor(hhmmss2insec(netcdf.getVar(f,netcdf.inqVarID(f,'Time'),start,count))+netcdf.getVar(f,netcdf.inqVarID(f,'msec'),start,count)/1000+increaseAllinseconds); % 'Time'?
            image_time_hhmmss = insec2hhmmss(image_time_hhmmss);
            image_time_hhmmssnew = image_time_hhmmss;
        end
        %}
        time_interval72(i) = sum(int_arr(DMT_DOF_SPEC_OVERLOAD~=0));
        
        % Simplified by DS - Removed image_time_hhmmssnew as it was defined by and never changed from image_time_hhmmss
        image_time_hhmmss = image_time_hhmmssall(start+1:start+count);
        
        % If image time crosses midnight, add 240000 to all times past midnight
%         image_time_hhmmss(image_time_hhmmss<10000)=image_time_hhmmss(image_time_hhmmss<10000)+240000;

        
        
        % Save an intermediate output file every 8000 steps through the loop
        if i==8000
            save([outfile(1:end-3) 'tempComp.mat']);
        end
        
        %% Calculate area of particle according to image reconstruction and airspeed (if tasMax exceeded)

        % Correct for airspeeds exceeding the max airspeed for the probe
        if(tas(i) > tasMax) % Set to threshold as necessary - stretch area of particle
			fprintf(2,'TAS at tas index %d exceeds tasMax (%.1f) of probe. Reconstructing area...\n\n',...
				i,tasMax);
            area = area*tas(i)/tasMax;
        end

        particle_mass = area*0;
        calcd_area = area*0;
        for iiii=1:length(area)
           particle_mass(iiii)=single_mass(particle_diameter_minR(iiii)/10, habit1(iiii));  % in unit of gram
           calcd_area(iiii)=single_area(particle_diameter_minR(iiii)/10, habit1(iiii));  % in unit of mm^2         
        end
        particle_massbl=0.115/1000*area.^(1.218); % in unit of gram


        %% Added by Robert Jackson -- old version did not have area ratio code
        area_ratio = area./(pi/4.*particle_diameter_minR.^2);
        auto_reject(area_ratio < .2) = 'z';
        
        %% Added by Will to calculate terminal velocity and precipitation rate
        particle_vt = area*0;
        for iiii=1:length(area)
           particle_vt(iiii)=single_vt(particle_diameter_minR(iiii)/1000, area_ratio(iiii), particle_mass(iiii)/1000,Pres(i),Temp(i));  % in unit of gram
        end
        particle_pr=particle_mass.*particle_vt;

        
        %% Time-dependent threshold for interarrival time - Added by Dan Stechman - 5/10/16
		% Enable this section to use a time-dependent threshold for interarrival time. Also need to enable section at top of
		% script allowing for threshold file to be pulled in
		
        % Ingest previously calculated interarrival time threshold and flag in auto_reject appropriately to remove particle
		% flagged with short inter arrv time, and the one immediately before it
        %{
        if strcmp(iaThreshFile,'NONE') == 0
			auto_reject_preIAT = auto_reject;
			iaThresh_ncid=netcdf.open(iaThreshFile,'nowrite');
			iaThresh=netcdf.getVar(iaThresh_ncid,netcdf.inqVarID(iaThresh_ncid,'threshold'));
			netcdf.close(iaThresh_ncid);
		
			if ((length(int_arr) == 1) && (int_arr(1) <= iaThresh(1)))
				auto_reject(1) = 'S';
			else
				if int_arr(1) <= iaThresh(1)
					auto_reject(1) = 'S';
				end
		
				for ix = 2:length(int_arr)
					if int_arr(ix) <= iaThresh(ix)
						auto_reject((ix-1):ix) = 'S';
					end
				end
			end
		end
        %}
   
        %% Legacy interarrival time integrity analysis 
        %{
        % Time and interarrival calculation. Modified by Will Wu 11/12/2013
        % Simplified (tested/changed by DS)
        if strcmp(probename,'2DC')==1 || strcmp(probename,'2DP')==1 || strcmp(probename,'F2DC')==1
            fracseccc= netcdf.getVar(f,netcdf.inqVarID(f,'msec'),start,count);
            image_timeia = hhmmss2insec(image_time_hhmmss)+fracseccc*1e-2; % for 2DC
        elseif strcmp(probename,'CIP')==1 || strcmp(probename,'PIP')==1
            image_timeia = hhmmss2insec(image_time_hhmmss)+msec*1e-3+microsec; % for CIP 
        else
            image_timeia = hhmmss2insec(image_time_hhmmss)+msec*1e-3+microsec/10^6; % for HVPS
        end
        
        disp('Checking Interarrival Times')

        nThrow=0;
        for(itemp=min(rec_nums):max(rec_nums))
           rec_particles = find(rec_nums == itemp);
           rej = auto_reject(rec_particles);
           arr = int_arr(rec_particles);
           sum_arr = sum(arr(2:end));
           if(~isempty(rec_particles) && length(rec_particles) > 1)
               int_arr(rec_particles(1)) = int_arr(rec_particles(2));
           elseif(length(rec_particles) == 1)
               int_arr(rec_particles(1)) = 0;
           end

           if (strcmp(probename,'CIP')==1 || strcmp(probename,'PIP')==1 || strcmp(probename,'HVPS')==1 ) % 2DC use the interarrival time for every particles, not absolute time 

               if(isempty(rec_particles))
                 sum_int_arr_good = 0;
               else
                 sum_int_arr_good = image_timeia(rec_particles(end))-image_timeia(rec_particles(1));
               end
               if ~(sum_int_arr_good >= .6*sum_arr && sum_int_arr_good <= 1.4*sum_arr)
                 auto_reject(rec_particles) = 'I';
                 %disp(['Record ' num2str(itemp) ' thrown out: Accepted time = ' num2str(sum_int_arr_good) ' total time = ' num2str(sum_arr)]);
                 nThrow=nThrow+1;  
                 nThrow11=nThrow11+1;
               end

           end
        end
        disp([num2str(100*nThrow/(max(rec_nums)-min(rec_nums)+1)),'% is thrown out']);
        throwoutpercent(i)=100*nThrow/(max(rec_nums)-min(rec_nums)+1);
        maxRecNum=max(max(rec_nums),maxRecNum);
        totalint(i)=sum_int_arr_good;
        intsum(i)=sum_arr;
        save('intarrhvps.mat','int_arr')
        %}
        
        %% Shatter identification and removal - Added by Dan Stechman on 5/31/2016
        % Currently this is spiral-dependent and uses a threshold defined in the header of this script
        % Flag particles as shattered if their interarrival time is less than or equal to the threshold. Also flag the particle
        % immediately before the target particle.
        if applyIntArrThresh
			% If the first particle in the next 1-sec period has a small interarrival time, we flag the last particle of
			% the current period as shattered as well
			if ~isempty(int_arr2)
				if int_arr2 <= intar_threshold(i)
					auto_reject(end) = 'S';
				end
			end

			if (length(int_arr) == 1 && int_arr(1) <= intar_threshold(i)) 
				auto_reject(1) = 'S';
			else
				if int_arr(1) <= intar_threshold(i)
					auto_reject(1) = 'S';
				end

				for ix = 2:length(int_arr)
					if int_arr(ix) <= intar_threshold(i)
						auto_reject(ix-1:ix) = 'S';
					end
				end
			end

			
			% Experimental option to reaccept particles flagged as shattered which may in fact be the result of diffraction
			% fringes
			% Added by Dan Stechman - 6/8/2015 - with base code by Wei Wu
			if reaccptShatrs
				% Start by defining the indices for the beginning and end of individual shattering events
				rBegin = ((int_arr > intar_threshold(i) & int_arr3 < intar_threshold(i)));
				rEnd = ((int_arr < intar_threshold(i) & int_arr3 > intar_threshold(i)));
				
				maxParticle = reaccptD;
				eIndex = [];
				
				% We search through each individual set of shattering events and check to see if any of the particles are both
				% larger than the reacceptance diameter and have an interarrival time less than the reacceptance threshold as we'd
				% expect diffraction fringes to be larger than shattered particles and to have a particularly small interarrival time
				for iEvent = find(rBegin):find(rEnd)
					if ((particle_diameter_minR(iEvent) > maxParticle) && (int_arr(iEvent) < reaccptMaxIA))
						maxParticle = particle_diameter_minR(iEvent);
						eIndex = iEvent;
					end
				end

				auto_reject(eIndex) = 'R';
			end
			
        
			% Following vars used for verifying shatter removal and reacceptance in external script - can be commented out if desired
			shatterLocs = find(auto_reject == 'S');
			shatterIA = int_arr(shatterLocs);
			shatterTimes = Time_in_seconds(shatterLocs);
            shatterDiam = particle_diameter_minR(shatterLocs);
            
			shatrReject_times = vertcat(shatrReject_times, shatterTimes);
			shatrReject_intArr = vertcat(shatrReject_intArr, shatterIA);
            shatrReject_diam = vertcat(shatrReject_diam, shatterDiam);
            
            rccptLocs = find(auto_reject == 'R');
			rccptIA = int_arr(rccptLocs);
			rccptTimes = Time_in_seconds(rccptLocs);
            rccptDiam = particle_diameter_minR(rccptLocs);
            
			rccptReject_times = vertcat(rccptReject_times, rccptTimes);
			rccptReject_intArr = vertcat(rccptReject_intArr, rccptIA);
            rccptReject_diam = vertcat(rccptReject_diam, rccptDiam);
            
                        
			loopedTimes = vertcat(loopedTimes, Time_in_seconds);
			loopedIntArr = vertcat(loopedIntArr, int_arr);
            loopedDiam = vertcat(loopedDiam, particle_diameter_minR);
			loopedAutoRej = vertcat(loopedAutoRej, auto_reject);
        end


        %% Apply rejection criteria and identify good and bad particles 
        % Modify the next line to include/exclude any particles you see fit. 
        
        good_particles = (auto_reject == '0' | auto_reject == 'H' | auto_reject == 'h' | auto_reject == 'u' | auto_reject == 'R');
        bad_particles = ~(auto_reject == '0' | auto_reject == 'H' | auto_reject == 'h' | auto_reject == 'u' | auto_reject == 'R');
% 		bad_particles = (auto_reject == 'S');
        
        % Legacy: Rejection criteria used in the past
        %{
        %if RejectCriterier==0
        %    good_particles = (auto_reject ~= 'c'); %  &  centerin==1; % & int_arr > 1e-5 int_arr > 1e-5 &
        %else
        %    good_particles = (auto_reject == '0' | auto_reject == 'H' | auto_reject == 'h' | auto_reject == 'u' & int_arr > intar_threshold) % | auto_reject == 'u'); % & centerin==1; % & int_arr > 1e-5; 
        %end
        %}
        
        
        if SAmethod==0
            good_particles = good_particles & centerin==1;
            bad_particles = bad_particles & centerin==1;
        elseif SAmethod==1
            good_particles = good_particles & entirein==0;
            bad_particles = bad_particles & centerin==0;
        end
        
        good_partpercent(i)=sum(good_particles)/length(good_particles);
        
        rejectpercentbycriterion(i,1)=sum(centerin==1)/length(good_particles);
        rejectpercentbycriterion(i,2)=sum(auto_reject == '0')/length(good_particles);
        rejectpercentbycriterion(i,3)=sum(auto_reject == 'H')/length(good_particles);
        rejectpercentbycriterion(i,4)=sum(auto_reject == 'h')/length(good_particles);
        rejectpercentbycriterion(i,5)=sum(auto_reject == 'u')/length(good_particles);
        rejectpercentbycriterion(i,6)=sum(auto_reject == 'a')/length(good_particles);
        rejectpercentbycriterion(i,7)=sum(auto_reject == 't')/length(good_particles);
        rejectpercentbycriterion(i,8)=sum(auto_reject == 'p')/length(good_particles);
        rejectpercentbycriterion(i,9)=sum(auto_reject == 's')/length(good_particles);
        rejectpercentbycriterion(i,10)=sum(auto_reject == 'z')/length(good_particles);
        rejectpercentbycriterion(i,11)=sum(auto_reject == 'i')/length(good_particles);
        rejectpercentbycriterion(i,12)=sum(auto_reject == 'A')/length(good_particles);
        rejectpercentbycriterion(i,13)=sum(auto_reject == 'S')/length(good_particles); %Shattered - Added DS
		rejectpercentbycriterion(i,14)=sum(auto_reject == 'R')/length(good_particles); %Reaccepted - Added DS
        
        numGoodparticles(i)=length(good_particles);
        numBadparticles(i)=length(bad_particles);
%         disp([int32(timehhmmss(i)), sum(good_particles),length(good_particles),length(good_particles)-sum(good_particles)]);

        image_time = hhmmss2insec(image_time_hhmmss);

        % Good (accepted) particles
        good_image_times = image_time(good_particles);
        good_particle_diameter_minR = particle_diameter_minR(good_particles);
        good_particle_diameter_AreaR = particle_diameter_AreaR(good_particles);
        good_int_arr=int_arr(good_particles);
        good_ar = area_ratio(good_particles);
        good_area = area(good_particles);
        good_perimeter = perimeter(good_particles);
        if iCreateAspectRatio == 1
        good_AspectRatio = aspectRatio(good_particles & entirein==0);        
        good_AspectRatio1 = aspectRatio1(good_particles & entirein==0);
        end
        good_ar1 = area_ratio(good_particles & entirein==0);
        good_image_times1 = image_time(good_particles  & entirein==0);
        good_iwc=particle_mass(good_particles);
        good_partarea=calcd_area(good_particles);
        good_iwcbl=particle_massbl(good_particles);
        good_vt=particle_vt(good_particles);
        good_pr=particle_pr(good_particles);        
        good_habit=habit1(good_particles);
        
        good_particle_diameter=good_particle_diameter_minR;
        good_particle_diameter1 = particle_diameter_minR(good_particles  & entirein==0);
        
        if iCreateBad == 1
        % Bad (rejected) particles
        bad_image_times = image_time(bad_particles);
        bad_particle_diameter_minR = particle_diameter_minR(bad_particles);
        bad_particle_diameter_AreaR = particle_diameter_AreaR(bad_particles);
        bad_int_arr=int_arr(bad_particles);
        bad_ar = area_ratio(bad_particles);
        bad_area = area(bad_particles);
        bad_perimeter = perimeter(bad_particles);
        bad_AspectRatio = aspectRatio(bad_particles & entirein==0);        
        bad_AspectRatio1 = aspectRatio1(bad_particles & entirein==0);
        bad_ar1 = area_ratio(bad_particles & entirein==0);
        bad_image_times1 = image_time(bad_particles  & entirein==0);
        bad_iwc=particle_mass(bad_particles);
        bad_partarea=calcd_area(bad_particles);
        bad_iwcbl=particle_massbl(bad_particles);
        bad_vt=particle_vt(bad_particles);
        bad_pr=particle_pr(bad_particles);        
        bad_habit=habit1(bad_particles);
        
        bad_particle_diameter=bad_particle_diameter_minR;
        bad_particle_diameter1 = particle_diameter_minR(bad_particles  & entirein==0);
        end
        %% Perform various status and error checks
        if mod(i,1000) == 0
			fprintf('%d/%d | %s\n',i,one_sec_dur,datestr(now));
        end
        
        total_one_sec_locs(i) = length(find(image_time >= one_sec_times(i) & image_time < one_sec_times(i+1)));
        time_interval2(i) = sum(int_arr(image_time >= one_sec_times(i) & image_time < one_sec_times(i+1)));

        if sum(image_time >= one_sec_times(i) & image_time < one_sec_times(i+1)) ~= length(image_time)
			fprintf(2,'%d / %d\tError on sizing at index %d\n',sum(image_time >= one_sec_times(i) & image_time < one_sec_times(i+1)),length(image_time),i);
        end

        if(total_one_sec_locs(i) == 0)
         time_interval2(i) = 1;
        end

        %% Sort good (accepted) particles into size distributions
        good_one_sec_locs = find(good_image_times >= one_sec_times(i) & good_image_times < one_sec_times(i+1));
        good_one_sec_locs1 = find(good_image_times1 >= one_sec_times(i) & good_image_times1 < one_sec_times(i+1));
        
        goodintpercent(i) = sum(good_int_arr(good_image_times >= one_sec_times(i) & good_image_times < one_sec_times(i+1)))/time_interval2(i);

        one_sec_ar(i) = mean(good_ar1(good_one_sec_locs1));
        
        if ~isempty(good_one_sec_locs)

           for j = 1:num_bins
               particle_dist_minR(i,j)  = length(find(good_particle_diameter_minR(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter_minR(good_one_sec_locs) < kk(j+1)));
               particle_dist_AreaR(i,j) = length(find(good_particle_diameter_AreaR(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter_AreaR(good_one_sec_locs) < kk(j+1)));

               % Create Habit Number Size Distribution 
               cip2_habitsd(i,j,1) = length(find(good_habit(good_one_sec_locs)=='s' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitsd(i,j,2) = length(find(good_habit(good_one_sec_locs)=='l' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitsd(i,j,3) = length(find(good_habit(good_one_sec_locs)=='o' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitsd(i,j,4) = length(find(good_habit(good_one_sec_locs)=='t' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitsd(i,j,5) = length(find(good_habit(good_one_sec_locs)=='h' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitsd(i,j,6) = length(find(good_habit(good_one_sec_locs)=='i' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitsd(i,j,7) = length(find(good_habit(good_one_sec_locs)=='g' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitsd(i,j,8) = length(find(good_habit(good_one_sec_locs)=='d' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitsd(i,j,9) = length(find(good_habit(good_one_sec_locs)=='a' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitsd(i,j,10) = length(find(good_habit(good_one_sec_locs)=='I' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));

               % Create Habit Mass Size Distribution 
               cip2_habitmsd(i,j,1) = sum(good_iwc(good_habit(good_one_sec_locs)=='s' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitmsd(i,j,2) = sum(good_iwc(good_habit(good_one_sec_locs)=='l' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitmsd(i,j,3) = sum(good_iwc(good_habit(good_one_sec_locs)=='o' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitmsd(i,j,4) = sum(good_iwc(good_habit(good_one_sec_locs)=='t' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitmsd(i,j,5) = sum(good_iwc(good_habit(good_one_sec_locs)=='h' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitmsd(i,j,6) = sum(good_iwc(good_habit(good_one_sec_locs)=='i' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitmsd(i,j,7) = sum(good_iwc(good_habit(good_one_sec_locs)=='g' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitmsd(i,j,8) = sum(good_iwc(good_habit(good_one_sec_locs)=='d' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitmsd(i,j,9) = sum(good_iwc(good_habit(good_one_sec_locs)=='a' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               cip2_habitmsd(i,j,10) = sum(good_iwc(good_habit(good_one_sec_locs)=='I' & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1)));


               particle_area(i,j) = nansum(good_area(good_one_sec_locs(good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1))));

               cip2_meanp(i,j) = nanmean(good_perimeter(good_one_sec_locs(good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1))));

               if iCreateAspectRatio == 1

               particle_aspectRatio(i,j) = nanmean(good_AspectRatio(good_one_sec_locs1(good_particle_diameter1(good_one_sec_locs1) >= kk(j) &...
                   good_particle_diameter1(good_one_sec_locs1) < kk(j+1))));

               particle_aspectRatio1(i,j) = nanmean(good_AspectRatio1(good_one_sec_locs1(good_particle_diameter1(good_one_sec_locs1) >= kk(j) &...
                   good_particle_diameter1(good_one_sec_locs1) < kk(j+1))));
               end
               particle_areaRatio1(i,j) = nanmean(good_ar1(good_one_sec_locs1(good_particle_diameter1(good_one_sec_locs1) >= kk(j) &...
                   good_particle_diameter1(good_one_sec_locs1) < kk(j+1))));


               cip2_iwc(i,j) = nansum(good_iwc(good_one_sec_locs(good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1))));

               cip2_partarea(i,j) = nansum(good_partarea(good_one_sec_locs(good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1))));

               cip2_iwcbl(i,j) = nansum(good_iwcbl(good_one_sec_locs(good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1))));

               cip2_vt(i,j) = nansum(good_vt(good_one_sec_locs(good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1))));

               cip2_pr(i,j) = nansum(good_pr(good_one_sec_locs(good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                   good_particle_diameter(good_one_sec_locs) < kk(j+1))));


               for k = 1:length(area_bins)-1
                   area_dist2(i,j,k) = length(find(good_ar(good_one_sec_locs) >= area_bins(k) & ...
                       good_ar(good_one_sec_locs) < area_bins(k+1) & good_particle_diameter(good_one_sec_locs) >= kk(j) &...
                       good_particle_diameter(good_one_sec_locs) < kk(j+1)));
               end
		   end
		   
           % Normalize by binwidth and convert from mm to cm
           particle_dist_minR(i,:)=particle_dist_minR(i,:)./binwidth*10;
           particle_dist_AreaR(i,:)=particle_dist_AreaR(i,:)./binwidth*10;
           cip2_iwc(i,:)=cip2_iwc(i,:)./binwidth*10; %g/cm
           cip2_iwcbl(i,:)=cip2_iwcbl(i,:)./binwidth*10;
           cip2_vt(i,:)=cip2_vt(i,:)./binwidth*10;
           cip2_pr(i,:)=cip2_pr(i,:)./binwidth*10;
           cip2_partarea(i,:)=cip2_partarea(i,:)./binwidth*10;
           particle_area(i,:)=particle_area(i,:)./binwidth*10;               
           
           for mmmmmm=1:10
               cip2_habitsd(i,:,mmmmmm)=cip2_habitsd(i,:,mmmmmm)./binwidth*10;
               cip2_habitmsd(i,:,mmmmmm)=cip2_habitmsd(i,:,mmmmmm)./binwidth*10;
           end
           
           for mmmmmm = 1:length(area_bins)-1
               area_dist2(i,:,mmmmmm) =area_dist2(i,:,mmmmmm)./binwidth*10 ;
           end
           
           % Generalized effective radius calculation from Fu (1996)
           cip2_re(i) = (sqrt(3)/(3*0.91))*1000*(sum(cip2_iwc(i,:)./binwidth,2)/sum(particle_area(i,:)./binwidth,2))*1000; % in unit of um
        
		else

           particle_dist_minR(i,1:num_bins) = 0;
           particle_dist_AreaR(i,1:num_bins) = 0;
           area_dist2(i,1:num_bins,1:length(area_bins)-1) = 0;
           cip2_partarea(i,:) = 0;
           cip2_iwc(i,:) = 0;
           cip2_iwcbl(i,:) = 0;
           cip2_vt(i,:) = 0;
           cip2_pr(i,:) = 0;
           cip2_re(i) = 0;
           cip2_habitsd(i,:,:) = 0;
           cip2_habitmsd(i,:,:) = 0;
           time_interval2(i) = 1;
                     
           % Legacy: used in Paris intercomparison
           %{
           time_interval22(i) = 1;
           time_interval32(i) = 1;
           time_interval42(i) = 1;
           time_interval52(i) = 0;
           time_interval62(i) = 1;
           %}
           time_interval72(i) = 0;

           TotalPC1(i)=1;        
           TotalPC2(i)=1;

        end
        
        if iCreateBad == 1
        %% Sort bad (rejected) particles into size distributions
        bad_one_sec_locs = find(bad_image_times >= one_sec_times(i) & bad_image_times < one_sec_times(i+1));
        bad_one_sec_locs1 = find(bad_image_times1 >= one_sec_times(i) & bad_image_times1 < one_sec_times(i+1));
        
        badintpercent(i) = sum(bad_int_arr(bad_image_times >= one_sec_times(i) & bad_image_times < one_sec_times(i+1)))/time_interval2(i);
        
        bad_one_sec_ar(i) = mean(bad_ar1(bad_one_sec_locs1));
        
        if ~isempty(bad_one_sec_locs)

           for j = 1:num_bins
               bad_particle_dist_minR(i,j)  = length(find(bad_particle_diameter_minR(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter_minR(bad_one_sec_locs) < kk(j+1)));
               bad_particle_dist_AreaR(i,j) = length(find(bad_particle_diameter_AreaR(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter_AreaR(bad_one_sec_locs) < kk(j+1)));

               % Create Habit Number Size Distribution 
               bad_cip2_habitsd(i,j,1) = length(find(bad_habit(bad_one_sec_locs)=='s' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitsd(i,j,2) = length(find(bad_habit(bad_one_sec_locs)=='l' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitsd(i,j,3) = length(find(bad_habit(bad_one_sec_locs)=='o' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitsd(i,j,4) = length(find(bad_habit(bad_one_sec_locs)=='t' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitsd(i,j,5) = length(find(bad_habit(bad_one_sec_locs)=='h' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitsd(i,j,6) = length(find(bad_habit(bad_one_sec_locs)=='i' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitsd(i,j,7) = length(find(bad_habit(bad_one_sec_locs)=='g' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitsd(i,j,8) = length(find(bad_habit(bad_one_sec_locs)=='d' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitsd(i,j,9) = length(find(bad_habit(bad_one_sec_locs)=='a' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitsd(i,j,10) = length(find(bad_habit(bad_one_sec_locs)=='I' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));

               % Create Habit Mass Size Distribution 
               bad_cip2_habitmsd(i,j,1) = sum(bad_iwc(bad_habit(bad_one_sec_locs)=='s' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitmsd(i,j,2) = sum(bad_iwc(bad_habit(bad_one_sec_locs)=='l' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitmsd(i,j,3) = sum(bad_iwc(bad_habit(bad_one_sec_locs)=='o' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitmsd(i,j,4) = sum(bad_iwc(bad_habit(bad_one_sec_locs)=='t' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitmsd(i,j,5) = sum(bad_iwc(bad_habit(bad_one_sec_locs)=='h' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitmsd(i,j,6) = sum(bad_iwc(bad_habit(bad_one_sec_locs)=='i' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitmsd(i,j,7) = sum(bad_iwc(bad_habit(bad_one_sec_locs)=='g' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitmsd(i,j,8) = sum(bad_iwc(bad_habit(bad_one_sec_locs)=='d' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitmsd(i,j,9) = sum(bad_iwc(bad_habit(bad_one_sec_locs)=='a' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               bad_cip2_habitmsd(i,j,10) = sum(bad_iwc(bad_habit(bad_one_sec_locs)=='I' & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));


               bad_particle_area(i,j) = nansum(bad_area(bad_one_sec_locs(bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1))));

               bad_cip2_meanp(i,j) = nanmean(bad_perimeter(bad_one_sec_locs(bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1))));


               bad_particle_aspectRatio(i,j) = nanmean(bad_AspectRatio(bad_one_sec_locs1(bad_particle_diameter1(bad_one_sec_locs1) >= kk(j) &...
                   bad_particle_diameter1(bad_one_sec_locs1) < kk(j+1))));

               bad_particle_aspectRatio1(i,j) = nanmean(bad_AspectRatio1(bad_one_sec_locs1(bad_particle_diameter1(bad_one_sec_locs1) >= kk(j) &...
                   bad_particle_diameter1(bad_one_sec_locs1) < kk(j+1))));

               bad_particle_areaRatio1(i,j) = nanmean(bad_ar1(bad_one_sec_locs1(bad_particle_diameter1(bad_one_sec_locs1) >= kk(j) &...
                   bad_particle_diameter1(bad_one_sec_locs1) < kk(j+1))));


               bad_cip2_iwc(i,j) = nansum(bad_iwc(bad_one_sec_locs(bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1))));

               bad_cip2_partarea(i,j) = nansum(bad_partarea(bad_one_sec_locs(bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1))));

               bad_cip2_iwcbl(i,j) = nansum(bad_iwcbl(bad_one_sec_locs(bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1))));

               bad_cip2_vt(i,j) = nansum(bad_vt(bad_one_sec_locs(bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1))));

               bad_cip2_pr(i,j) = nansum(bad_pr(bad_one_sec_locs(bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                   bad_particle_diameter(bad_one_sec_locs) < kk(j+1))));


               for k = 1:length(area_bins)-1
                   bad_area_dist2(i,j,k) = length(find(bad_ar(bad_one_sec_locs) >= area_bins(k) & ...
                       bad_ar(bad_one_sec_locs) < area_bins(k+1) & bad_particle_diameter(bad_one_sec_locs) >= kk(j) &...
                       bad_particle_diameter(bad_one_sec_locs) < kk(j+1)));
               end
           end
           
           % Normalize by binwidth and convert from mm to cm
           bad_particle_dist_minR(i,:)=bad_particle_dist_minR(i,:)./binwidth*10;
           bad_particle_dist_AreaR(i,:)=bad_particle_dist_AreaR(i,:)./binwidth*10;
           bad_cip2_iwc(i,:)=bad_cip2_iwc(i,:)./binwidth*10; %g/cm
           bad_cip2_iwcbl(i,:)=bad_cip2_iwcbl(i,:)./binwidth*10;
           bad_cip2_vt(i,:)=bad_cip2_vt(i,:)./binwidth*10;
           bad_cip2_pr(i,:)=bad_cip2_pr(i,:)./binwidth*10;
           bad_cip2_partarea(i,:)=bad_cip2_partarea(i,:)./binwidth*10;
           bad_particle_area(i,:)=bad_particle_area(i,:)./binwidth*10;               
           
           for mmmmmm=1:10
               bad_cip2_habitsd(i,:,mmmmmm)=bad_cip2_habitsd(i,:,mmmmmm)./binwidth*10;
               bad_cip2_habitmsd(i,:,mmmmmm)=bad_cip2_habitmsd(i,:,mmmmmm)./binwidth*10;
           end
           
           for mmmmmm = 1:length(area_bins)-1
               bad_area_dist2(i,:,mmmmmm)=bad_area_dist2(i,:,mmmmmm)./binwidth*10 ;
           end
           
           % Generalized effective radius calculation from Fu (1996)
           bad_cip2_re(i) = (sqrt(3)/(3*0.91))*1000*(sum(bad_cip2_iwc(i,:)./binwidth,2)/sum(bad_particle_area(i,:)./binwidth,2))*1000; % in unit of um
        
        else
           bad_particle_dist_minR(i,1:num_bins) = 0;
           bad_particle_dist_AreaR(i,1:num_bins) = 0;
           bad_area_dist2(i,1:num_bins,1:length(area_bins)-1) = 0;
           bad_cip2_partarea(i,:) = 0;
           bad_cip2_iwc(i,:) = 0;
           bad_cip2_iwcbl(i,:) = 0;
           bad_cip2_vt(i,:) = 0;
           bad_cip2_pr(i,:) = 0;
           bad_cip2_re(i) = 0;
           bad_cip2_habitsd(i,:,:) = 0;
           bad_cip2_habitmsd(i,:,:) = 0;

        end
        end
        warning on all
    elseif (int32(timehhmmss(i))<int32(starttime(jjj)))

       particle_dist_minR(i,1:num_bins) = NaN;
       particle_dist_AreaR(i,1:num_bins) = NaN;
       area_dist2(i,1:num_bins,1:length(area_bins)-1) = NaN;
       cip2_partarea(i,:) = NaN;
       cip2_iwc(i,:) = NaN;
       cip2_iwcbl(i,:) = NaN;
       cip2_vt(i,:) = NaN;
       cip2_pr(i,:) = NaN;
       cip2_re(i) = NaN;
       cip2_habitsd(i,:,:) = NaN;
       cip2_habitmsd(i,:,:) = NaN;
       one_sec_ar(i) = NaN;
       good_partpercent(i)=1;
       rejectpercentbycriterion(i,:)=NaN;
       numGoodparticles(i)=NaN;
       time_interval2(i) = 1;
       
       % Legacy: used in Paris intercomparison
       %{
       time_interval22(i) = 1;
       time_interval32(i) = 1;
       time_interval42(i) = 1;    
       time_interval52(i) = 0;    
       time_interval62(i) = 1;       
       %}
       time_interval72(i) = 0;
       
       TotalPC1(i)=1;        
       TotalPC2(i)=1;
       
    end

end

% Finished Sorting and close input file.
netcdf.close(f);

fprintf('int_arr > 1 mean: %.4f, max: %.4f\nNumber of particles with int_arr > 1: %d\n\n',...
	mean(intArrGT1),max(intArrGT1),sumIntArrGT1);

fprintf('Size distribution calculations and sorting completed %s\n\n', datestr(now));


%% Check TAS length, should be the same
% if (jjj~=length(start_all))
%     disp([jjj, length(start_all)])
%     %error('Watch Out for less TAS time at the end!')
% end

%disp([num2str(100*nThrow11/maxRecNum),'% is thrown out IN TOTAL']); 

%% Combine - calculate sample volumes, and divide by sample volumes 
% Modified by Will, Nov 27th, 2013. For flexible bins
cip2_binmin = kk(1:end-1); 
cip2_binmax = kk(2:end); 
cip2_binmid = (cip2_binmin+cip2_binmax)/2; 
cip2_bindD = diff(kk);

% Legacy bin and surface area calculations
%{
% cip2_binmin = diodesize/2:diodesize:(num_bins-0.5)*diodesize; %(12.5:25:(num_bins-0.5)*25);
% cip2_binmax = 3*diodesize/2:diodesize:(num_bins+0.5)*diodesize; %(37.5:25:(num_bins+0.5)*25);
% cip2_binmid = diodesize:diodesize:num_bins*diodesize; %(25:25:num_bins*25);
% cip2_bindD = diodesize*ones(1,num_bins);

% sa2 = calc_sa(num_bins,res,armdst,num_bins);  %mm2
% switch probename
%     case 'PIP'
%         sa2 = calc_sa_randombins_PIP(cip2_binmid,res,armdst,num_diodes, SAmethod); %(bins_mid,res,armdst,num_diodes)
%     case '2DS'
%         sa2 = calc_sa_randombins(cip2_binmid,res,armdst,num_diodes, SAmethod); %(bins_mid,res,armdst,num_diodes)
% end
%}

sa2 = calc_sa_randombins(cip2_binmid,res,armdst,num_diodes,SAmethod, probetype); %(bins_mid,res,armdst,num_diodes)

% Clocking problem correction
vol_scale_factor = tas/tasMax;  
vol_scale_factor(vol_scale_factor < 1) = 1;

TotalPC2_pre = TotalPC2;

if probetype==2
    time_interval200=1-time_interval72';

elseif probetype==1
	% Correct offset in probe particle count (TotalPC2) when we have negative values
    TotalPC2(TotalPC2<0)=TotalPC2(TotalPC2<0)+2^16;
	
	% Derive a linear scale factor based on the difference between number of images (TotalPC1)
	% and number of particles counted by the probe (TotalPC2).
    time_interval199=(TotalPC1./TotalPC2)';

elseif 0==probetype
    time_interval200=1-time_interval72';
end

% Experimental - Use with care!
% It was discovered that for data collected during the PECAN project, there were quite
% a few periods of time when the number of images we had for a 1-sec period of time was
% up to twice that of the number of particles the probe counted.
% This next if-statement contains code to find and change these instances to 1, resolving 
% the far exaggerated concentrations that resulted otherwise.
if probetype==1
	TotalPCerrIx = find(time_interval199 > 1);
	time_interval200 = time_interval199;
	time_interval200(TotalPCerrIx) = 1;
end

fprintf(['Total image count exceeded probe particle count %d times\ntime_interval200',...
	' was set to 1 in these cases. See TotalPCerrIx variable for indices of occurence.\n\n'],...
	length(TotalPCerrIx));


for j=1:num_bins
    % Sample volume is in m-3
%     svol_old(j,:)=dof/100.*sa/100.*tas;
    svol2(j,:) = sa2(j)*(1e-3)^2*time_interval200.*tas; %m3 .*vol_scale_factor
end

svol2 = svol2*100^3; %cm3

for j = 1:10
    svol2a(:,:,j) = svol2';
end

% Good (accepted) particles
cip2_conc_minR  = particle_dist_minR./svol2';
cip2_conc_AreaR = particle_dist_AreaR./svol2';
cip2_area = particle_area./svol2';
cip2_partarea = cip2_partarea./svol2';
cip2_iwc = cip2_iwc./svol2';
cip2_iwcbl = cip2_iwcbl./svol2';
cip2_vt = cip2_vt./svol2';
cip2_pr = cip2_pr./svol2';

cip2_countP_no  = particle_dist_minR;
cip2_conc_areaDist = permute(double(area_dist2)./svol2a, [3 2 1]);
cip2_n = nansum(cip2_conc_minR,2);
cip2_lwc = lwc_calc(cip2_conc_minR,cip2_binmid);

% Bad (rejected) particles
bad_cip2_conc_minR  = bad_particle_dist_minR./svol2';
bad_cip2_conc_AreaR = bad_particle_dist_AreaR./svol2';
bad_cip2_area = bad_particle_area./svol2';
bad_cip2_partarea = bad_cip2_partarea./svol2';
bad_cip2_iwc = bad_cip2_iwc./svol2';
bad_cip2_iwcbl = bad_cip2_iwcbl./svol2';
bad_cip2_vt = bad_cip2_vt./svol2';
bad_cip2_pr = bad_cip2_pr./svol2';

bad_cip2_countP_no  = bad_particle_dist_minR;
bad_cip2_conc_areaDist = permute(double(bad_area_dist2)./svol2a, [3 2 1]);
bad_cip2_n = nansum(bad_cip2_conc_minR,2);
bad_cip2_lwc = lwc_calc(bad_cip2_conc_minR,cip2_binmid);



%% Output results into NETCDF files (mainf)

fprintf('Now writing output files %s\n\n',datestr(now));

if applyIntArrThresh
	save([outfile(1:end-3) 'noShatters.mat']);
else
	save([outfile(1:end-3) 'withShatters.mat']);
end


% Define Dimensions
dimid0 = netcdf.defDim(mainf,'CIPcorrlen',num_bins);
dimid1 = netcdf.defDim(mainf,'CIParealen',10);
dimid2 = netcdf.defDim(mainf,'Time',length(timehhmmss));
dimid3 = netcdf.defDim(mainf,'Habit',10);

% Define Variables
varid0 = netcdf.defVar(mainf,'time','double',dimid2); 
netcdf.putAtt(mainf, varid0,'units','HHMMSS');
netcdf.putAtt(mainf, varid0,'name','Time');

varid1 = netcdf.defVar(mainf,'bin_min','double',dimid0); 
netcdf.putAtt(mainf, varid1,'units','millimeter');
netcdf.putAtt(mainf, varid1,'long_name','bin minimum size');
netcdf.putAtt(mainf, varid1,'short_name','bin min');

varid2 = netcdf.defVar(mainf,'bin_max','double',dimid0); 
netcdf.putAtt(mainf, varid2,'units','millimeter');
netcdf.putAtt(mainf, varid2,'long_name','bin maximum size');
netcdf.putAtt(mainf, varid2,'short_name','bin max');

varid3 = netcdf.defVar(mainf,'bin_mid','double',dimid0); 
netcdf.putAtt(mainf, varid3,'units','millimeter');
netcdf.putAtt(mainf, varid3,'long_name','bin midpoint size');
netcdf.putAtt(mainf, varid3,'short_name','bin mid');

varid4 = netcdf.defVar(mainf,'bin_dD','double',dimid0); 
netcdf.putAtt(mainf, varid4,'units','millimeter');
netcdf.putAtt(mainf, varid4,'long_name','bin size');
netcdf.putAtt(mainf, varid4,'short_name','bin size');

% Good (accepted) particles
varid5 = netcdf.defVar(mainf,'conc_minR','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid5,'units','cm-4');
netcdf.putAtt(mainf, varid5,'long_name','Size distribution using Dmax');
netcdf.putAtt(mainf, varid5,'short_name','N(Dmax)');

varid6 = netcdf.defVar(mainf,'area','double',[dimid1 dimid0 dimid2]); 
netcdf.putAtt(mainf, varid6,'units','cm-4');
netcdf.putAtt(mainf, varid6,'long_name','binned area ratio');
netcdf.putAtt(mainf, varid6,'short_name','binned area ratio');

varid7 = netcdf.defVar(mainf,'conc_AreaR','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid7,'units','cm-4');
netcdf.putAtt(mainf, varid7,'long_name','Size distribution using area-equivalent Diameter');
netcdf.putAtt(mainf, varid7,'short_name','N(Darea)');

varid8 = netcdf.defVar(mainf,'n','double',dimid2); 
netcdf.putAtt(mainf, varid8,'units','cm-3');
netcdf.putAtt(mainf, varid8,'long_name','number concentration');
netcdf.putAtt(mainf, varid8,'short_name','N');

varid9 = netcdf.defVar(mainf,'total_area','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid9,'units','mm2/cm4');
netcdf.putAtt(mainf, varid9,'long_name','projected area (extinction)');
netcdf.putAtt(mainf, varid9,'short_name','Ac');

varid10 = netcdf.defVar(mainf,'mass','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid10,'units','g/cm4');
netcdf.putAtt(mainf, varid10,'long_name','mass using m-D relations');
netcdf.putAtt(mainf, varid10,'short_name','mass');

varid11 = netcdf.defVar(mainf,'habitsd','double',[dimid3 dimid0 dimid2]); 
netcdf.putAtt(mainf, varid11,'units','cm-4');
netcdf.putAtt(mainf, varid11,'long_name','Size Distribution with Habit');
netcdf.putAtt(mainf, varid11,'short_name','habit SD');

varid12 = netcdf.defVar(mainf,'re','double',dimid2); 
netcdf.putAtt(mainf, varid12,'units','mm');
netcdf.putAtt(mainf, varid12,'long_name','effective radius');
netcdf.putAtt(mainf, varid12,'short_name','Re');

varid13 = netcdf.defVar(mainf,'ar','double',dimid2); 
netcdf.putAtt(mainf, varid13,'units','100/100');
netcdf.putAtt(mainf, varid13,'long_name','Area Ratio');
netcdf.putAtt(mainf, varid13,'short_name','AR');

varid14 = netcdf.defVar(mainf,'massBL','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid14,'units','g/cm4');
netcdf.putAtt(mainf, varid14,'long_name','mass using Baker and Lawson method');
netcdf.putAtt(mainf, varid14,'short_name','mass_BL');

varid15 = netcdf.defVar(mainf,'Reject_ratio','double',dimid2); 
netcdf.putAtt(mainf, varid15,'units','100/100');
netcdf.putAtt(mainf, varid15,'long_name','Reject Ratio');

varid16 = netcdf.defVar(mainf,'vt','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid16,'units','g/cm4');
netcdf.putAtt(mainf, varid16,'long_name','Mass-weighted terminal velocity');

varid17 = netcdf.defVar(mainf,'Prec_rate','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid17,'units','mm/hr');
netcdf.putAtt(mainf, varid17,'long_name','Precipitation Rate');

varid18 = netcdf.defVar(mainf,'habitmsd','double',[dimid3 dimid0 dimid2]); 
netcdf.putAtt(mainf, varid18,'units','g/cm-4');
netcdf.putAtt(mainf, varid18,'long_name','Mass Size Distribution with Habit');
netcdf.putAtt(mainf, varid18,'short_name','Habit Mass SD');

varid19 = netcdf.defVar(mainf,'Calcd_area','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid19,'units','mm^2/cm4');
netcdf.putAtt(mainf, varid19,'long_name','Particle Area Calculated using A-D realtions');
netcdf.putAtt(mainf, varid19,'short_name','Ac_calc');

varid20 = netcdf.defVar(mainf,'count','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid20,'units','1');
netcdf.putAtt(mainf, varid20,'long_name','number count for partial images without any correction');

if iCreateAspectRatio == 1
varid21 = netcdf.defVar(mainf,'mean_aspect_ratio_rectangle','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid21,'units','1');
netcdf.putAtt(mainf, varid21,'long_name','Aspect Ratio by Rectangle fit');

varid22 = netcdf.defVar(mainf,'mean_aspect_ratio_ellipse','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid22,'units','1');
netcdf.putAtt(mainf, varid22,'long_name','Aspect Ratio by Ellipse fit');
end
varid23 = netcdf.defVar(mainf,'mean_area_ratio','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid23,'units','1');
netcdf.putAtt(mainf, varid23,'long_name','Area Ratio');

varid24 = netcdf.defVar(mainf,'mean_perimeter','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid24,'units','um');
netcdf.putAtt(mainf, varid24,'long_name','mean perimeter');

if iCreateBad == 1

% Bad (rejected) particles
varid25 = netcdf.defVar(mainf,'REJ_conc_minR','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid25,'units','cm-4');
netcdf.putAtt(mainf, varid25,'long_name','Size distribution of rejected particles using Dmax');
netcdf.putAtt(mainf, varid25,'short_name','N(Dmax) rejected');

varid26 = netcdf.defVar(mainf,'REJ_area','double',[dimid1 dimid0 dimid2]); 
netcdf.putAtt(mainf, varid26,'units','cm-4');
netcdf.putAtt(mainf, varid26,'long_name','binned area ratio of rejected particles');
netcdf.putAtt(mainf, varid26,'short_name','binned area ratio of rejected particles');

varid27 = netcdf.defVar(mainf,'REJ_conc_AreaR','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid27,'units','cm-4');
netcdf.putAtt(mainf, varid27,'long_name','Size distribution of rejected particles using area-equivalent Diameter');
netcdf.putAtt(mainf, varid27,'short_name','N(Darea) rejected');

varid28 = netcdf.defVar(mainf,'REJ_n','double',dimid2); 
netcdf.putAtt(mainf, varid28,'units','cm-3');
netcdf.putAtt(mainf, varid28,'long_name','number concentration of rejected particles');
netcdf.putAtt(mainf, varid28,'short_name','N_rejected');

varid29 = netcdf.defVar(mainf,'REJ_total_area','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid29,'units','mm2/cm4');
netcdf.putAtt(mainf, varid29,'long_name','projected area (extinction) of rejected particles');
netcdf.putAtt(mainf, varid29,'short_name','Ac_rejected');

varid30 = netcdf.defVar(mainf,'REJ_mass','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid30,'units','g/cm4');
netcdf.putAtt(mainf, varid30,'long_name','mass of rejected particles using m-D relations');
netcdf.putAtt(mainf, varid30,'short_name','mass_rejected');

varid31 = netcdf.defVar(mainf,'REJ_habitsd','double',[dimid3 dimid0 dimid2]); 
netcdf.putAtt(mainf, varid31,'units','cm-4');
netcdf.putAtt(mainf, varid31,'long_name','Size Distribution with Habit of rejected particles');
netcdf.putAtt(mainf, varid31,'short_name','habit SD rejected');

varid32 = netcdf.defVar(mainf,'REJ_re','double',dimid2); 
netcdf.putAtt(mainf, varid32,'units','mm');
netcdf.putAtt(mainf, varid32,'long_name','effective radius of rejected particles');
netcdf.putAtt(mainf, varid32,'short_name','Re_rejected');

varid33 = netcdf.defVar(mainf,'REJ_ar','double',dimid2); 
netcdf.putAtt(mainf, varid33,'units','100/100');
netcdf.putAtt(mainf, varid33,'long_name','Area Ratio of rejected particles');
netcdf.putAtt(mainf, varid33,'short_name','AR_rejected');

varid34 = netcdf.defVar(mainf,'REJ_massBL','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid34,'units','g/cm4');
netcdf.putAtt(mainf, varid34,'long_name','mass of rejected particles using Baker and Lawson method');
netcdf.putAtt(mainf, varid34,'short_name','mass_BL_rejected');

varid35 = netcdf.defVar(mainf,'REJ_vt','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid35,'units','g/cm4');
netcdf.putAtt(mainf, varid35,'long_name','Mass-weighted terminal velocity of rejected particles');

varid36 = netcdf.defVar(mainf,'REJ_Prec_rate','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid36,'units','mm/hr');
netcdf.putAtt(mainf, varid36,'long_name','Precipitation Rate of rejected particles');

varid37 = netcdf.defVar(mainf,'REJ_habitmsd','double',[dimid3 dimid0 dimid2]); 
netcdf.putAtt(mainf, varid37,'units','g/cm-4');
netcdf.putAtt(mainf, varid37,'long_name','Mass Size Distribution with Habit of rejected particles');
netcdf.putAtt(mainf, varid37,'short_name','Habit Mass SD rejected');

varid38 = netcdf.defVar(mainf,'REJ_Calcd_area','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid38,'units','mm^2/cm4');
netcdf.putAtt(mainf, varid38,'long_name','Particle Area of rejected particles Calculated using A-D realtions');
netcdf.putAtt(mainf, varid38,'short_name','Ac_calc_rejected');

varid39 = netcdf.defVar(mainf,'REJ_count','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid39,'units','1');
netcdf.putAtt(mainf, varid39,'long_name','number count of rejected particles for partial images without any correction');

varid40 = netcdf.defVar(mainf,'REJ_mean_aspect_ratio_rectangle','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid40,'units','1');
netcdf.putAtt(mainf, varid40,'long_name','Aspect Ratio of rejected particles by Rectangle fit');

varid41 = netcdf.defVar(mainf,'REJ_mean_aspect_ratio_ellipse','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid41,'units','1');
netcdf.putAtt(mainf, varid41,'long_name','Aspect Ratio of rejected particles by Ellipse fit');

varid42 = netcdf.defVar(mainf,'REJ_mean_area_ratio','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid42,'units','1');
netcdf.putAtt(mainf, varid42,'long_name','Area Ratio of rejected particles');

varid43 = netcdf.defVar(mainf,'REJ_mean_perimeter','double',[dimid0 dimid2]); 
netcdf.putAtt(mainf, varid43,'units','um');
netcdf.putAtt(mainf, varid43,'long_name','mean perimeter of rejected particles');
end
netcdf.endDef(mainf)

% Output Variables
netcdf.putVar ( mainf, varid0, timehhmmss );
netcdf.putVar ( mainf, varid1, cip2_binmin );
netcdf.putVar ( mainf, varid2, cip2_binmax );
netcdf.putVar ( mainf, varid3, cip2_binmid );
netcdf.putVar ( mainf, varid4, cip2_bindD );

% Good (accepted) particles
netcdf.putVar ( mainf, varid5, cip2_conc_minR' );
netcdf.putVar ( mainf, varid6, cip2_conc_areaDist);
netcdf.putVar ( mainf, varid7, cip2_conc_AreaR' );
netcdf.putVar ( mainf, varid8, cip2_n);
netcdf.putVar ( mainf, varid9, cip2_area');
netcdf.putVar ( mainf, varid10, cip2_iwc');
netcdf.putVar ( mainf, varid11, permute(double(cip2_habitsd)./svol2a, [3 2 1]) );
netcdf.putVar ( mainf, varid12, cip2_re );
netcdf.putVar ( mainf, varid13, one_sec_ar );
netcdf.putVar ( mainf, varid14, cip2_iwcbl' );
netcdf.putVar ( mainf, varid15, 1-good_partpercent );
netcdf.putVar ( mainf, varid16, cip2_vt' );
netcdf.putVar ( mainf, varid17, cip2_pr' );
netcdf.putVar ( mainf, varid18, permute(double(cip2_habitmsd)./svol2a, [3 2 1]) );
netcdf.putVar ( mainf, varid19, cip2_partarea');
netcdf.putVar ( mainf, varid20, cip2_countP_no');
if iCreateAspectRatio == 1
netcdf.putVar ( mainf, varid21, particle_aspectRatio);
netcdf.putVar ( mainf, varid22, particle_aspectRatio1);
end
netcdf.putVar ( mainf, varid23, particle_areaRatio1);
netcdf.putVar ( mainf, varid24, cip2_meanp');

if iCreateBad == 1

% Bad (rejected) particles
netcdf.putVar ( mainf, varid25, bad_cip2_conc_minR' );
netcdf.putVar ( mainf, varid26, bad_cip2_conc_areaDist);
netcdf.putVar ( mainf, varid27, bad_cip2_conc_AreaR' );
netcdf.putVar ( mainf, varid28, bad_cip2_n);
netcdf.putVar ( mainf, varid29, bad_cip2_area');
netcdf.putVar ( mainf, varid30, bad_cip2_iwc');
netcdf.putVar ( mainf, varid31, permute(double(bad_cip2_habitsd)./svol2a, [3 2 1]) );
netcdf.putVar ( mainf, varid32, bad_cip2_re );
netcdf.putVar ( mainf, varid33, bad_one_sec_ar );
netcdf.putVar ( mainf, varid34, bad_cip2_iwcbl' );
netcdf.putVar ( mainf, varid35, bad_cip2_vt' );
netcdf.putVar ( mainf, varid36, bad_cip2_pr' );
netcdf.putVar ( mainf, varid37, permute(double(bad_cip2_habitmsd)./svol2a, [3 2 1]) );
netcdf.putVar ( mainf, varid38, bad_cip2_partarea');
netcdf.putVar ( mainf, varid39, bad_cip2_countP_no');
netcdf.putVar ( mainf, varid40, bad_particle_aspectRatio);
netcdf.putVar ( mainf, varid41, bad_particle_aspectRatio1);
netcdf.putVar ( mainf, varid42, bad_particle_areaRatio1);
netcdf.putVar ( mainf, varid43, bad_cip2_meanp');
end

netcdf.close(mainf) % Close output NETCDF file 

fprintf('sizeDist_Paris.m script completed %s\n',datestr(now));

end
