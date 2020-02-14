function read_binary_SPEC(infilename, outfilename, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Read the raw base*.2DS file, and then write into NETCDF file 
%% Follow the SPEC manual 
%%  by Will Wu, 08/01/2014
%%  **************************
%   *** Modification Notes ***
%   **************************
%   * Modified to prevent termination during rare instance when attempting
%   to decode beyond the image buffer
%               Joe Finlon, 03/08/2017
%   * Added netCDF4 support and data compression
%               Joe Finlon, 02/13/2019
%   * Improvements in handling instances when prob TAS is incorrect
%               Joe Finlon & Wei Wu, 11/3/2019
%   * Support for applying time offset correction (e.g., SOCRATES)
%               Joe Finlon, 11/3/2019
%   * Changed how the optional flight time and TAS are read in
%               Joe Finlon, 02/07/2020
%   * Ensures the hour timestamp is less than 24 when saving to file
%               Joe Finlon, 02/10/2020
%	* Small fix to global attributes
%				Joe Finlon, 02/12/2020
%   Interface:
%       infilename: The input file name outfilename: The output file name
%       (e.g., '/path/imgData.date.2DS') varargin: Optional arguments,
%       separated by commas, containing an array of 1-HZ flight times and
%       true airspeed values
%               flt_time: Array of times in HHMMSS (00 UTC will have value
%                       of 0)
%               flt_tas: Array of true airspeeds [m/s]
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Added support for probe TAS ratio calculation ~ Joe Finlon 11/3/19,
% 02/07/20
if ~isempty(varargin)
    fltTime = varargin{1}; % aircraft time in HHMMSS ~ Joe Finlon 02/07/20
    fltTAS = varargin{2}; % aircraft TAS in m/s
    fltTAS_fixed = interp_tas(fltTAS); % fix bad TAS values through linear interpolation
    %tasFile = varargin{1}; % path to aircraft TAS information
    %fltTAS = ncread(tasFile, 'TAS'); % aircraft TAS in m/s
    %fltTime1 = ncread(tasFile, 'Time'); % aircraft time in HHMMSS
end

% Added support for probe time offset (e.g., SOCRATES) ~ Joe Finlon 11/3/19
if length(varargin) == 3
   timeOffset = varargin{3}; % user-specified probe time offset [seconds]
end

starpos = find(infilename == '*',1,'last');

if ~isempty(starpos)
    files = dir(infilename);
    filenums = length(files);
    filedir = infilename(1:starpos-1);
else
    filenums = 1;
end

for i = 1:filenums
    if filenums > 1
        infilename = [filedir,files(i).name];
    end
    disp(['Raw File: ', infilename])
    
    if outfilename == '1'
        slashpos = find(infilename == '.',1,'last');
        outfilename = ['DIMG.',infilename(1:slashpos-1),'.cdf'];
    end
    
    outfilename1=[outfilename, '.H.nc'];
    outfilename2=[outfilename, '.V.nc'];
    
    fid=fopen(infilename,'r','l');
    
    % Added dynamic software version to netCDF metadata - Joe Finlon
    % 02/07/20
    versionID = fopen('version.txt', 'r');
    software_string = fscanf(versionID, '%s');
    fclose(versionID);

    f = netcdf.create(outfilename1, 'NETCDF4'); % netCDF-4/HDF5 compression support - Added by Joe Finlon 02/13/19
    
    dimid0 = netcdf.defDim(f,'time',netcdf.getConstant('NC_UNLIMITED'));
    dimid1 = netcdf.defDim(f,'ImgRowlen',8);
    dimid2 = netcdf.defDim(f,'ImgBlocklen',1700);
    
    NC_GLOBAL = netcdf.getConstant('NC_GLOBAL'); % Added file attributes ~ Joe Finlon 02/07/20
	netcdf.putAtt(f, NC_GLOBAL, 'Software', software_string);
	netcdf.putAtt(f, NC_GLOBAL, 'Creation Time', datestr(now, 'yyyy/mm/dd HH:MM:SS'));
	%netcdf.putAtt(f, NC_GLOBAL, 'Project', projectname);
	netcdf.putAtt(f, NC_GLOBAL, 'Probe Channel_Orientation', 'Horizontal');
    
    % Added data compression using 'defVarDeflate' argument ~ Joe Finlon 02/13/19
    varid0 = netcdf.defVar(f,'year','short',dimid0); netcdf.defVarDeflate(f, varid0, true, true, 9);
    varid1 = netcdf.defVar(f,'month','byte',dimid0); netcdf.defVarDeflate(f, varid1, true, true, 9);
    varid2 = netcdf.defVar(f,'day','byte',dimid0); netcdf.defVarDeflate(f, varid2, true, true, 9);
    varid3 = netcdf.defVar(f,'hour','byte',dimid0); netcdf.defVarDeflate(f, varid3, true, true, 9);
    varid4 = netcdf.defVar(f,'minute','byte',dimid0); netcdf.defVarDeflate(f, varid4, true, true, 9);
    varid5 = netcdf.defVar(f,'second','byte',dimid0); netcdf.defVarDeflate(f, varid5, true, true, 9);
    varid6 = netcdf.defVar(f,'millisec','short',dimid0); netcdf.defVarDeflate(f, varid6, true, true, 9);
    varid7 = netcdf.defVar(f,'wkday','byte',dimid0); netcdf.defVarDeflate(f, varid7, true, true, 9);
    varid8 = netcdf.defVar(f,'data','int',[dimid1 dimid2 dimid0]); netcdf.defVarDeflate(f, varid8, true, true, 9);
    varid9 = netcdf.defVar(f,'tas','double',dimid0); netcdf.defVarDeflate(f, varid9, true, true, 9); netcdf.defVarFill(f , varid9, false, -999.);
    varid10 = netcdf.defVar(f,'tasRatio','double',dimid0); netcdf.defVarDeflate(f, varid10, true, true, 9);
    netcdf.endDef(f)
        
    f1 = netcdf.create(outfilename2, 'NETCDF4');
    
    dimid01 = netcdf.defDim(f1,'time',netcdf.getConstant('NC_UNLIMITED'));
    dimid11 = netcdf.defDim(f1,'ImgRowlen',8);
    dimid21 = netcdf.defDim(f1,'ImgBlocklen',1700);
    
    NC_GLOBAL = netcdf.getConstant('NC_GLOBAL'); % Added file attributes ~ Joe Finlon 02/07/20
	netcdf.putAtt(f1, NC_GLOBAL, 'Software', software_string);
	netcdf.putAtt(f1, NC_GLOBAL, 'Creation Time', datestr(now, 'yyyy/mm/dd HH:MM:SS'));
	%netcdf.putAtt(f1, NC_GLOBAL, 'Project', projectname);
	netcdf.putAtt(f1, NC_GLOBAL, 'Probe Channel_Orientation', 'Vertical');
    
    % Added data compression using 'defVarDeflate' argument ~ Joe Finlon 02/13/19
    varid01 = netcdf.defVar(f1,'year','short',dimid01); netcdf.defVarDeflate(f1, varid01, true, true, 9);
    varid11 = netcdf.defVar(f1,'month','byte',dimid01); netcdf.defVarDeflate(f1, varid11, true, true, 9);
    varid21 = netcdf.defVar(f1,'day','byte',dimid01); netcdf.defVarDeflate(f1, varid21, true, true, 9);
    varid31 = netcdf.defVar(f1,'hour','byte',dimid01); netcdf.defVarDeflate(f1, varid31, true, true, 9);
    varid41 = netcdf.defVar(f1,'minute','byte',dimid01); netcdf.defVarDeflate(f1, varid41, true, true, 9);
    varid51 = netcdf.defVar(f1,'second','byte',dimid01); netcdf.defVarDeflate(f1, varid51, true, true, 9);
    varid61 = netcdf.defVar(f1,'millisec','short',dimid01); netcdf.defVarDeflate(f1, varid61, true, true, 9);
    varid71 = netcdf.defVar(f1,'wkday','byte',dimid01); netcdf.defVarDeflate(f1, varid71, true, true, 9);
    varid81 = netcdf.defVar(f1,'data','int',[dimid11 dimid21 dimid01]); netcdf.defVarDeflate(f1, varid81, true, true, 9);
    varid91 = netcdf.defVar(f1,'tas','double',dimid0); netcdf.defVarDeflate(f1, varid91, true, true, 9); netcdf.defVarFill(f1 , varid91, false, -999.);
    varid101 = netcdf.defVar(f1,'tasRatio','double',dimid0); netcdf.defVarDeflate(f1, varid101, true, true, 9);
    netcdf.endDef(f1)
    
    kk1=1;
    kk2=1;
    previousTAS = -999; % Joe Finlon 11/3/19
    endfile = 0;
    nNext=1;
    daystart = 99999999;
    dataprev=zeros(2048,1);
    recordTime = []; recordTime1 = []; % for use in TAS ratio calculation ~ Added by Joe Finlon 11/3/19
    
    %fseek(fid,4114*233191,'bof');
    while feof(fid)==0 && endfile == 0 
        %tic
        [year,month, wkday,day, hour, minute, second, millisec, data, discard]=readRecord(fid);
        % Apply record time correction if user-specified offset exists ~ Added by Joe Finlon 11/3/19
        if exist('timeOffset', 'var')
            [year, month, day] = ymd(datetime(year, month, day, hour, minute, second)+ seconds(timeOffset));
            [hour, minute, second] = hms(datetime(year, month, day, hour, minute, second)+ seconds(timeOffset));
        end
        % correct record time by the time offset (sec; if user-specified)
        
        %timebuffer = [year,month,day, hour, minute, second, millisec]
        %fprintf(formatSpec,A1,A2)
        [year1,month1, wkday1,day1, hour1, minute1, second1, millisec1, data1, discard1]=readRecord(fid);
        %fseek(fid,-4114,'cof');
        datan=[data' data1'];
        datan=datan';
        
        
        [imgH, imgV, nNext, imgTAS, previousTAS]=get_img(datan, hour*10000+minute*100+second+millisec/1000, outfilename, previousTAS);
        sizeimg= size(imgH);
        if sizeimg(2)>1700
            imgH=imgH(:,1:1700);
        end
        
        sizeimg= size(imgV);
        if sizeimg(2)>1700
            imgV=imgV(:,1:1700);
        end
        
        % Compute ratio between probe TAS and aircraft TAS ~ Joe Finlon
        % 11/3/19
        if ~isempty(varargin)
            if (hour<15)
                hour = hour + 24;
            end

            fltTimeInd = find(fltTime==hour*10000+minute*100+second);
            if (~isempty(fltTimeInd)) && (~isnan(fltTAS_fixed(fltTimeInd)))
                if imgTAS==-999. % uncorrectable probe TAS for current record
                    tasRatio = 1.0;
                else
                    tasRatio = imgTAS ./ fltTAS_fixed(fltTimeInd); % TAS ratio between SPEC probe and aircraft info
                end
            else
                tasRatio = 1.0; % assume SPEC probe TAS is correct as aircraft TAS does not exist for current time or is NaN
            end
        end
        
        if sum(sum(imgH))~=0
            for  mmm=1:8
                img1(mmm,1:1700)=sixteen2int(imgH((mmm-1)*16+1:mmm*16,1:1700));
            end
            % Ensures the hour timestamp is less than 24 when saving to file ~ Added by Joe Finlon 02/10/20
            if hour>=24
                hour = hour - 24;
            end
            netcdf.putVar ( f, varid0, kk1-1, 1, year );
            netcdf.putVar ( f, varid1, kk1-1, 1, month );
            netcdf.putVar ( f, varid2, kk1-1, 1, day );
            netcdf.putVar ( f, varid3, kk1-1, 1, hour );
            netcdf.putVar ( f, varid4, kk1-1, 1, minute );
            netcdf.putVar ( f, varid5, kk1-1, 1, second );
            netcdf.putVar ( f, varid6, kk1-1, 1, millisec );
            netcdf.putVar ( f, varid7, kk1-1, 1, wkday );
            netcdf.putVar ( f, varid8, [0, 0, kk1-1], [8,1700,1], img1 );
            netcdf.putVar ( f, varid9, kk1-1, 1, imgTAS); % Added variable ~ Joe Finlon 11/3/19
            if ~isempty(varargin)
                netcdf.putVar ( f, varid10, kk1-1, 1, tasRatio); % Added variable ~ Joe Finlon 11/3/19
            end
            
            kk1=kk1+1;
            if mod(kk1,1000) == 0
                 ['kk1=' num2str(kk1) ', ' datestr(now)]
            end
        end
        
        if sum(sum(imgV))~=0
            for  mmm=1:8
                img2(mmm,1:1700)=sixteen2int(imgV((mmm-1)*16+1:mmm*16,1:1700));
            end
            % Ensures the hour timestamp is less than 24 when saving to file ~ Added by Joe Finlon 02/10/20
            if hour>=24
                hour = hour - 24;
            end
            netcdf.putVar ( f1, varid01, kk2-1, 1, year );
            netcdf.putVar ( f1, varid11, kk2-1, 1, month );
            netcdf.putVar ( f1, varid21, kk2-1, 1, day );
            netcdf.putVar ( f1, varid31, kk2-1, 1, hour );
            netcdf.putVar ( f1, varid41, kk2-1, 1, minute );
            netcdf.putVar ( f1, varid51,  kk2-1, 1, second );
            netcdf.putVar ( f1, varid61, kk2-1, 1, millisec );
            netcdf.putVar ( f1, varid71, kk2-1, 1, wkday );
            netcdf.putVar ( f1, varid81, [0, 0, kk2-1], [8,1700,1], img2 );
            netcdf.putVar ( f1, varid91, kk2-1, 1, imgTAS); % Added variable ~ Joe Finlon 11/3/19
            if ~isempty(varargin)
                netcdf.putVar ( f1, varid101, kk2-1, 1, tasRatio); % Added variable ~ Joe Finlon 11/3/19
            end

            kk2=kk2+1;
            if mod(kk2,1000) == 0
                 ['kk2=' num2str(kk2) ', ' datestr(now)]
            end
        end
        
        %for j=1:4115
        bb=fread(fid,1,'int8');
        if feof(fid) == 1
            endfile=1;
            break
        end
        %end
        fseek(fid,-4115,'cof');
        %kk
        %ftell(fid)
        %toc
    end
    
    fclose(fid);
    
    netcdf.close(f);  
    netcdf.close(f1);      
end


end

function [year,month, wkday,day, hour, minute, second, millisec, data, discard, daystart]=readRecord(fid, daystart)

        year=fread(fid,1,'uint16');
        month=fread(fid,1,'uint16');
        wkday=fread(fid,1,'uint16');
        day=fread(fid,1,'uint16');
        hour=fread(fid,1,'uint16');
        minute=fread(fid,1,'uint16');
        second=fread(fid,1,'uint16');
        millisec=fread(fid,1,'uint16');
        data = fread(fid,2048,'uint16');
        discard=fread(fid,1,'uint16');
end

function [imgH, imgV, nNext, tas, tasPrev]=get_img(buf, timehhmmss, outfilename, tasPrev)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Decompress the image 
%% Follow the SPEC manual 
%%  by Will Wu, 06/20/2013
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imgH=zeros(128,1700);
    imgV=zeros(128,1700);
    tas = tasPrev; % TAS will be last known value if not overwritten in this subroutine ~ Added by Joe Finlon 11/3/19
    iSlice=0;
    iii=1;
    while iii<=2048 
        if 12883==buf(iii) %'0011001001010011'
              nH=bitand(buf(iii+1), 4095); %bin2dec('0000111111111111'));
              bHTiming=bitand(buf(iii+1), 4096); %bin2dec('0001000000000000'))/2^13;
              nV=bitand(buf(iii+2), 4095); %bin2dec('0000111111111111'));
              bVTiming=bitand(buf(iii+2), 4096); %bin2dec('0001000000000000'))/2^13;
              PC = buf(iii+3);
              nS = buf(iii+4);
              NHWord=buf(iii+1);
              NVWord=buf(iii+2);
              
            
              myformatNV = '%f, %d\n';
              fid = fopen([outfilename, '.NV.csv'],'a');
              fprintf(fid, myformatNV, [timehhmmss buf(iii+2)]);
              fclose(fid);

              if nH~=0 && nV~=0
                  system(['echo ' num2str(nH)  ' ' num2str(nV) ' >> output.txt']);
              end
              
              iii=iii+5;
              if bHTiming~=0 || bVTiming~=0     
                  iii=iii+nH+nV;
%               elseif nH~=0
              elseif nH~=0 && iii+nH-1<=length(buf) % Watch for index outside of image buffer - Added by Joe Finlon - 03/08/17
                  jjj=1;
                  kkk=0;
%                   while jjj<=nS && kkk<nH-2 % Last two slice is time
                  while iii+kkk<=length(buf) && jjj<=nS && kkk<nH-2 % Added by Joe Finlon - Last 2 slices is the time - 03/08/17
                      aa=bitand(buf(iii+kkk),16256)/2^7;  %bin2dec('0011111110000000')
                      bb=bitand(buf(iii+kkk),127); %bin2dec('0000000001111111')
                      imgH(min(128,bb+1):min(aa+bb,128),iSlice+jjj)=1;
                      bBase=min(aa+bb,128);
                      kkk=kkk+1;
%                       while( bitand(buf(iii+kkk),16384)==0  && kkk<nH-2) % bin2dec('1000000000000000')
                      while( iii+kkk<=length(buf) && bitand(buf(iii+kkk),16384)==0  && kkk<nH-2) % Added by Joe Finlon - 03/08/17
                          aa=bitand(buf(iii+kkk),16256)/2^7;
                          bb=bitand(buf(iii+kkk),127);
                          imgH(min(128,bBase+bb+1):min(bBase+aa+bb,128),iSlice+jjj)=1;
                          bBase=min(bBase+aa+bb,128);
                          kkk=kkk+1;
                      end
                      jjj=jjj+1;
                  end
                  iSlice=iSlice+nS+2;

                  imgH(:,iSlice-1)='10101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010';
                  imgH(:,iSlice-1)=imgH(:,iSlice-1)-48;
                  imgH(:,iSlice)=1;
                  tParticle=buf(iii+nH-2)*2^16+buf(iii+nH-1);
                  imgH(:,iSlice)=dec2bin(tParticle,128)-48;
                  imgH(97:128,iSlice)=dec2bin(tParticle,32)-48;
                  imgH(49:64,iSlice)=dec2bin(NHWord,16)-48;
                  imgH(65:80,iSlice)=dec2bin(PC,16)-48;
                  imgH(81:96,iSlice)=dec2bin(nS,16)-48;
                  iii=iii+nH;

%               elseif nV~=0
              elseif nV~=0 && iii+nV-1<=length(buf) % Watch for index outside of image buffer - Added by Joe Finlon - 03/08/17
                  jjj=1;
                  kkk=0;
%                   while jjj<=nS && kkk<nV-2 % Last two slice is time
                  while iii+kkk<=length(buf) && jjj<=nS && kkk<nV-2 % Added by Joe Finlon - Last 2 slices is the time - 03/08/17
                      aa=bitand(buf(iii+kkk),16256)/2^7;  %bin2dec('0011111110000000')
                      bb=bitand(buf(iii+kkk),127); %bin2dec('0000000001111111')
                      imgV(min(128,bb+1):min(aa+bb,128),iSlice+jjj)=1;
                      bBase=min(aa+bb,128);
                      kkk=kkk+1;
%                       while( bitand(buf(iii+kkk),16384)==0  && kkk<nV-2) % bin2dec('1000000000000000')
                      while( iii+kkk<=length(buf) && bitand(buf(iii+kkk),16384)==0  && kkk<nV-2) % Added by Joe Finlon - 03/08/17
                          aa=bitand(buf(iii+kkk),16256)/2^7;
                          bb=bitand(buf(iii+kkk),127);
                          imgV(min(128,bBase+bb+1):min(bBase+aa+bb,128),iSlice+jjj)=1;
                          bBase=min(bBase+aa+bb,128);
                          kkk=kkk+1;
                      end
                      jjj=jjj+1;
                  end
                  iSlice=iSlice+nS+2;

                  imgV(:,iSlice-1)='10101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010';
                  imgV(:,iSlice-1)=imgV(:,iSlice-1)-48;
                  imgV(:,iSlice)=1;
                  tParticle=buf(iii+nV-2)*2^16+buf(iii+nV-1);
                  imgV(:,iSlice)=dec2bin(tParticle,128)-48;
                  imgV(97:128,iSlice)=dec2bin(tParticle,32)-48;
                  imgV(49:64,iSlice)=dec2bin(NVWord,16)-48;
                  imgV(65:80,iSlice)=dec2bin(PC,16)-48;
                  imgV(81:96,iSlice)=dec2bin(nS,16)-48;
                  iii=iii+nV;
              end
              
        elseif 18507==buf(iii)
            tasTemp = typecast( uint32(bin2dec([dec2bin(buf(iii-1+50),16) dec2bin(buf(iii-1+51),16)])) ,'single');
            time = buf(iii-1+52)*2^16+buf(iii-1+53);
            nStereo = buf(iii-1+41);
            nTWM = buf(iii-1+42);
            nSCM = buf(iii-1+43);
            nHOL = buf(iii-1+44);
            nVOL = buf(iii-1+45);
            nVPD = buf(iii-1+34);
            nHPD = buf(iii-1+35);
            % Update TAS value for current record ~ Joe Finlon 11/3/19
            if tasTemp~=tasPrev
                if (tasTemp==170) || (nHOL+nVOL>0)
                    tas = tasPrev;
                else
                    tas = tasTemp;
                end
                tasPrev = tas; % update TAS for use in next iteration
            end
            % End of added code block ~ Joe Finlon 11/3/19
            myformat1 = '%f, %f, %f, %d, %d, %d, %d, %d, %d\n';
            fid = fopen([outfilename, '.tas.csv'],'a');
            fprintf(fid, myformat1, [timehhmmss tas time nTWM nSCM nHOL nVOL nVPD nHPD]);
            fclose(fid);
            
            iii = iii + 53;
        
        elseif 19787==buf(iii)
            timeWord = buf(iii-1+2)*2^16+buf(iii-1+3);
            MaskBits = [buf(iii-1+4:iii-1+19)]';
            timeStart = buf(iii-1+20)*2^16+buf(iii-1+21);
            timeEnd = buf(iii-1+22)*2^16+buf(iii-1+23);
            
            myformat2 = '%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d \n';
            fid = fopen([outfilename, '.MK.csv'],'a');
            fprintf(fid, myformat2, [timeWord MaskBits timeStart timeEnd]);
            fclose(fid);
            
            iii = iii + 23;
            
        else
            iii=iii+1;
            %disp(['Nowhere, move foreward...', num2str(iii)]);
        end
    end
    nNext=iii-2048;
end

function intres=sixteen2int(original)

intres=zeros(1,1700);
for i=1:16
    temp=original(i,:)*2^(16-i);
    intres=intres+temp;
end
end

function tasNew = interp_tas(tasValues) % interpolate TAS if bad values exist ~ Added by Joe Finlon 11/3/19

tasNew = tasValues;
good = find((tasValues>0) & ~isnan(tasValues)); bad = find((tasValues==0) | isnan(tasValues));
if ~isempty(bad)
    tasNew(bad) = interp1(good, tasValues(good), bad);
end

end
