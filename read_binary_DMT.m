function read_binary_DMT(infilename,outfilename)

decomp_convert1 = zeros(1700,8);

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
    
    if outfilename == '1'
        slashpos = find(infilename == '/',1,'last');
        outfilename = ['caps.',infilename(slashpos+1:end),'.cdf'];
    end
    
    fid=fopen(infilename,'r','l');
    
    %%% Updated for new MATLAB NETCDF interface
    f = netcdf.create(outfilename, 'clobber');
    
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
    
%     f = netcdf(outfilename,'clobber');
%     %f = netcdf('clobber');
%     
%     f('time') = 0;
%     f('ImgRowlen') = 8;
%     f('ImgBlocklen')= 1700;
%     
%     f{'year'} = 'time';
%     f{'month'} = 'time';
%     f{'day'} = 'time';
%     f{'hour'} = 'time';
%     f{'minute'} = 'time';
%     f{'second'} = 'time';
%     f{'millisec'} = 'time';
%     f{'wkday'} = 'time';
%     
%     f{'data'} = {'time','ImgBlocklen','ImgRowlen'};
    
    % recsize = f('time');
    
%     autoscale(f{'data'},1);
%     autonan(f{'data'},1);
    
    
    %
    % bytes=['00';'C0';'43';'00';'01';'81';'00';'F0';'41';'0F';'AA';'AA';'AA';'AA';'AA';'AA';'AA';'AA';'89';'9E';'91';'AA';'3C';'66';'6C';'67';'41';'03';'7F'...
    %     ;'00';'FC';'FF';'43';'03'];
    
    kk=1;
    
    endfile = 0;
    
    %if project == 1
    %    discard = fread(fid,1,'uint16'); %
    %end 
    
    % while feof(fid)==0 & kk <= 3000
    while feof(fid)==0 & endfile == 0

        year=fread(fid,1,'uint16');
        month=fread(fid,1,'uint16');
        day=fread(fid,1,'uint16');
        hour=fread(fid,1,'uint16');
        minute=fread(fid,1,'uint16');
        second=fread(fid,1,'uint16');
        millisec=fread(fid,1,'uint16');
        wkday=fread(fid,1,'uint16');
        
        data = fread(fid,4096,'uchar');
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
                kk
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
            end
            i=i+1;
            
        end        
        
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
            decomp_convert1(1:length(decomp_convert),:) = decomp_convert;
            decomp_convert1(length(decomp_convert)+1,:)=[170   170   170   170   170   170   170   170];
            decomp_convert1(length(decomp_convert)+2:1700,:)=-1;
        end
        %     recsize(:) = kk;
%         f{'year'}(kk) = year;
%         f{'month'}(kk) = month;
%         f{'day'}(kk) = day;
%         f{'hour'}(kk) = hour;
%         f{'minute'}(kk) = minute;
%         f{'second'}(kk) = second;
%         f{'millisec'}(kk) = millisec;
%         f{'wkday'}(kk) = wkday;
%         
%         
%         f{'data'}(kk,:,:)=decomp_convert(:,:);
        
        netcdf.putVar ( f, varid0, kk-1, 1, year );
        netcdf.putVar ( f, varid1, kk-1, 1, month );
        netcdf.putVar ( f, varid2, kk-1, 1, day );
        netcdf.putVar ( f, varid3, kk-1, 1, hour );
        netcdf.putVar ( f, varid4, kk-1, 1, minute );
        netcdf.putVar ( f, varid5, kk-1, 1, second );
        netcdf.putVar ( f, varid6, kk-1, 1, millisec );
        netcdf.putVar ( f, varid7, kk-1, 1, wkday );
        netcdf.putVar ( f, varid8, [0, 0, kk-1], [8,1700,1], decomp_convert1' );
    
        kk=kk+1;
        if mod(kk,100) == 0
            kk
            datestr(now)
        end
        clear decomp dd k2 b1 b2
        
        %% Read more bytes to test if the end of file is reached
        for j=1:20  % 4132
            bb=fread(fid,1,'int8');
            if feof(fid) == 1
                endfile=1;
                break
            end
        end
        fseek(fid,-20,'cof');
        
    end
    
    fclose(fid);
    % close(f);
    netcdf.close(f);  % New interface by Will
end
