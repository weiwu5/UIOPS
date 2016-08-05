function read_binary_SEA(infilename,outfilename)
%% Function to decompress SEA raw files 
% Need to double check the file format and code for each probes
% This only works for MC3E filed campaign
% * July 11, 2016, Created this new interface function, Wei Wu


starpos = find(infilename == '*',1,'last');
nWierdTotal = 0;

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
    
    if outfilename == '1'
        slashpos = find(infilename == '/',1,'last');
        outfilename = ['DIMG.',infilename(slashpos+1:end),'.cdf'];
    end
    
    fid=fopen(infilename,'r','l');
    infilename
    
    %%% Updated for new MATLAB NETCDF interface
    f = netcdf.create([outfilename, '.CIP.cdf'], 'clobber');
    
    dimid0 = netcdf.defDim(f,'time',netcdf.getConstant('NC_UNLIMITED'));
    dimid1 = netcdf.defDim(f,'ImgRowlen',8);
    dimid2 = netcdf.defDim(f,'ImgBlocklen',1700);
    
    varid0 = netcdf.defVar(f,'year','double',dimid0);
    varid1 = netcdf.defVar(f,'month','double',dimid0);
    varid2 = netcdf.defVar(f,'day','double',dimid0);
    varid3 = netcdf.defVar(f,'hour','double',dimid0);
    varid4 = netcdf.defVar(f,'minute','double',dimid0);
    varid5 = netcdf.defVar(f,'second','double',dimid0);
    varid6 = netcdf.defVar(f,'millisec','double',dimid0);
    varid7 = netcdf.defVar(f,'wkday','double',dimid0);
    varid8 = netcdf.defVar(f,'data','double',[dimid1 dimid2 dimid0]);
    netcdf.endDef(f)
    
    f1 = netcdf.create([outfilename, '.2DC.cdf'], 'clobber');
    
    dimid01 = netcdf.defDim(f1,'time',netcdf.getConstant('NC_UNLIMITED'));
    dimid11 = netcdf.defDim(f1,'ImgRowlen',8);
    dimid21 = netcdf.defDim(f1,'ImgBlocklen',1700);
    
    varid01 = netcdf.defVar(f1,'year','double',dimid01);
    varid11 = netcdf.defVar(f1,'month','double',dimid01);
    varid21 = netcdf.defVar(f1,'day','double',dimid01);
    varid31 = netcdf.defVar(f1,'hour','double',dimid01);
    varid41 = netcdf.defVar(f1,'minute','double',dimid01);
    varid51 = netcdf.defVar(f1,'second','double',dimid01);
    varid61 = netcdf.defVar(f1,'millisec','double',dimid01);
    varid71 = netcdf.defVar(f1,'wkday','double',dimid01);
    varid81 = netcdf.defVar(f1,'data','double',[dimid11 dimid21 dimid01]);
    netcdf.endDef(f)
    
    kk=1;
    wkday = 1;
    datatemp = [0 0 0 0 0 0 0 0 0 0];
    numFilename = 0;
    nread =0;
    endfile = 0;
    
   
    % while feof(fid)==0 & kk <= 3000
    while feof(fid)==0 & endfile == 0

       [datalast,datatemp]=readDir(fid,datatemp);
       doffset1 = -datatemp(2);
       if datatemp(1)==999
           nTagNext=1; 
       else
           nTagNext=0;
       end
       
       while nTagNext==0 
           
           [datalast,datatemp]=readDir(fid,datatemp);
           
           if datatemp(1)==33000 & datatemp(3)==4098
                %datalast
                datatemp
                [datalast,datatemp]=readDir(fid,datatemp);
                %datatemp
                ttt = readTime(fid);
                
                year =ttt(1);
                month=ttt(2);
                day  =ttt(3);
                hour =ttt(4);
                minute=ttt(5);
                second=ttt(6);
                millisec=ttt(7);
                wkday=1;
                
                data1 = fread(fid,4098,'uchar');
            
                fseek(fid,-4150,0);
                [datalast,datatemp]=readDir(fid,datatemp);
                data = data1(1:4096);
                %     data=reshape(fread(fid,4096*8,'ubit1'),4096,8);
                %     b1 = [num2str(data(:,1)),num2str(data(:,2)),num2str(data(:,3)),num2str(data(:,4)),num2str(data(:,5)),...
                %         num2str(data(:,6)),num2str(data(:,7)),num2str(data(:,8))];

                bytes=dec2hex(data,2);
                kk;

                i=1;
                ii=1;
                b1full=dec2bin(hex2dec(bytes(:,:)),8);
                b2 = bin2dec(b1full(:,4:8));


                while i<4096
                    b1 = b1full(i,:);
                    curi = i;
                    i=i+1;
                    if b1(3) == '1'
                        %             i=i+1;
                    elseif b1(1) == '0' & b1(2) == '0'
                        %            b2=bin2dec(b1(4:8));
                        for k=1:b2(curi)+1;
                            if i < length(bytes)
                                decomp(ii,:)=bytes(i,:);
                            else break
                            end
                            ii=ii+1;
                            i=i+1;
                        end
                    elseif b1(1) == '1' & b1(2) == '0'
                        %            b2=bin2dec(b1(4:8));
                        for k=1:b2(curi)+1;
                            decomp(ii,:)='00';
                            ii=ii+1;
                        end
                    elseif b1(2) == '1' & b1(1) == '0'
                        %            b2=bin2dec(b1(4:8));
                        for k=1:b2(curi)+1;
                            decomp(ii,:)='FF';
                            ii=ii+1;
                        end
                    else
                        kk;
                    end
                end

                found = 0;
                i=1;
                count=0;
                while found == 0
                    if decomp(i)=='AA'
                        count=count+1;
                    else
                        count=0;
                    end

                    if count == 8 
                        found=1;
                        dd=i+1:8:length(decomp)-7;
                        nWierd=0;
                    end
                    
                    if i==length(decomp) % Add to avoid no 'AA' even though wierd to have no 'AA'...
                        found =1;
                        nWierd=1;
                        nWierdTotal =nWierdTotal +1; 
                    end
                    i=i+1;

                end


                if nWierd ==0
                %
                %     decomp_convert=[hex2dec(decomp(dd,:)),hex2dec(decomp(dd+1,:)),hex2dec(decomp(dd+2,:)),hex2dec(decomp(dd+3,:)),...
                %         hex2dec(decomp(dd+4,:)),hex2dec(decomp(dd+5,:)),hex2dec(decomp(dd+6,:)),hex2dec(decomp(dd+7,:))];
                decomp_convert=[hex2dec(decomp(dd+7,:)),hex2dec(decomp(dd+6,:)),hex2dec(decomp(dd+5,:)),hex2dec(decomp(dd+4,:)),...
                    hex2dec(decomp(dd+3,:)),hex2dec(decomp(dd+2,:)),hex2dec(decomp(dd+1,:)),hex2dec(decomp(dd,:))];
                k2=[decomp(dd,:),decomp(dd+1,:),decomp(dd+2,:),decomp(dd+3,:),decomp(dd+4,:),decomp(dd+5,:),decomp(dd+6,:),decomp(dd+7,:)];

                %             length_diff=length(decomp_convert) - length(handles.matrix(kk-1,:,:));
                %             matrix_size(kk)=length(decomp_convert);
                %             if length_diff > 0
                %                 handles.matrix(1:kk-1,length(handles.matrix(kk-1,:,:)):length(decomp_convert),:)=-1;
                %             elseif length_diff < 0
                %                 decomp_convert(length(decomp_convert):length(handles.matrix(kk-1,:,:)),:)=-1;
                %             end
                if length(decomp_convert) < 1700
                    decomp_convert(length(decomp_convert):1700,:)=-1;
                end


                netcdf.putVar ( f, varid0, kk-1, 1, year )
                netcdf.putVar ( f, varid1, kk-1, 1, month );
                netcdf.putVar ( f, varid2, kk-1, 1, day );
                netcdf.putVar ( f, varid3, kk-1, 1, hour );
                netcdf.putVar ( f, varid4, kk-1, 1, minute );
                netcdf.putVar ( f, varid5, kk-1, 1, second );
                netcdf.putVar ( f, varid6, kk-1, 1, millisec );
                netcdf.putVar ( f, varid7, kk-1, 1, wkday );
                netcdf.putVar ( f, varid8, [0, 0, kk-1], [8,1700,1], decomp_convert' );

                kk=kk+1;
                if mod(kk,100) == 0
                    kk
                    datestr(now)
                end
                
                end
                
           elseif datatemp(1)==5000 && datatemp(3)==4096
                [datalast,datatemp]=readDir(fid,datatemp);
                [datalast,datatemp]=readDir(fid,datatemp); 
                [datalast,datatemp]=readDir(fid,datatemp);
                [datalast,datatemp]=readDir(fid,datatemp);
                [datalast,datatemp]=readDir(fid,datatemp);
                
                ttt = readTime(fid);
                year =ttt(1);
                month=ttt(2);
                day  =ttt(3);
                hour =ttt(4);
                minute=ttt(5);
                second=ttt(6);
                millisec=ttt(7);
                wkday=1;
                %dataother = fread(fid,14,'uchar');

                temp0011 = fread(fid,1,'int16');
                temp0012 = fread(fid,1,'int16');
                fread(fid,10,'char');
%                 temp003 = fread(fid,1,'char');
%                 temp004 = fread(fid,1,'short');
%                 temp003/temp0011*temp0012*2*0.001
%                 temp002*25*0.000001
                data1 = fread(fid,4096,'uchar');
            
                fseek(fid,-4162,0);
                [datalast,datatemp]=readDir(fid,datatemp);

                tas=temp0011/temp0012*50*1000*25*0.000001;

                
                datafinal = reshape(data1,4,1024);
                temp1234 = datafinal(1,:);
                datafinal(1,:)=datafinal(4,:);
                datafinal(4,:)=temp1234;
                temp1234 = datafinal(2,:);
                datafinal(2,:)=datafinal(3,:);
                datafinal(3,:)=temp1234;                    
                
                netcdf.putVar ( f1, varid01, kk-1, 1, year )
                netcdf.putVar ( f1, varid11, kk-1, 1, month );
                netcdf.putVar ( f1, varid21, kk-1, 1, day );
                netcdf.putVar ( f1, varid31, kk-1, 1, hour );
                netcdf.putVar ( f1, varid41, kk-1, 1, minute );
                netcdf.putVar ( f1, varid51, kk-1, 1, second );
                netcdf.putVar ( f1, varid61, kk-1, 1, millisec );
                netcdf.putVar ( f1, varid71, kk-1, 1, wkday );
                netcdf.putVar ( f1, varid91, kk-1, 1, tas );
                netcdf.putVar ( f1, varid81, [0, 0, kk-1], [4,1024,1], datafinal );

                kk=kk+1;
                if mod(kk,100) == 0
                    kk
                    datestr(now)
                end
           end
           
           clear decomp dd k2 b1 b2

           if datatemp(1)==999 
               doffset2 = datatemp(2);
               nTagNext = 1;
               fseek(fid,doffset2+doffset1,0);
           end
       end
       
       
       for j=1:16
           bb=fread(fid,1,'int8');
           if feof(fid) == 1
               endfile=1;
               break
           end
       end
       fseek(fid,-16,'cof');
       
    end
    
    
end

fclose(fid);
% close(f);
netcdf.close(f);  % New interface by Will
netcdf.close(f1);  % New interface by Will

nWierdTotal
end

function [ldata,tdata]=readDir(fid,tdata)

tagNumber=fread(fid,1,'uint16');
dataOffset=fread(fid,1,'uint16');
numberBytes=fread(fid,1,'uint16');
samples=fread(fid,1,'uint16');
bytesPerSample=fread(fid,1,'uint16');
type=fread(fid,1,'uint8');
param1=fread(fid,1,'uint8');
param2=fread(fid,1,'uint8');
param3=fread(fid,1,'uint8');
address=fread(fid,1,'uint16');
ldata = tdata;
tdata = [tagNumber dataOffset numberBytes samples bytesPerSample type  param1 param2 param3 address];
    
end

function time=readTime(fid)

for i=1:2
    year=fread(fid,1,'uint16');
    month=fread(fid,1,'uint16');
    day=fread(fid,1,'uint16');
    hour=fread(fid,1,'uint16');
    minute=fread(fid,1,'uint16');
    second=fread(fid,1,'uint16');
    fracsec=fread(fid,1,'uint16');
    maxfreq=fread(fid,1,'uint16');
    bls=fread(fid,1,'uint16');
    time=[year month day hour minute second fracsec maxfreq bls];
end
end

function readTable(fid)

filename = fread(fid,8,'uint8');
filename = char(filename);
filename=filename';

tfiles = fread(fid,datalast(3),'uint8');
abc = char(tfiles);
abc'
end
