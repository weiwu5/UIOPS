function read_binary_PMS(infilename,outfilename)

%%%
% NCAR OAP format; Read XML first and then the data
% May need to double check the XML output to make sure the probetye is
% correct

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
        slashpos = find(infilename == '.',1,'last');
        outfilename  = ['DIMG.',infilename(1:slashpos-1),'.2dc.cdf']; %(slashpos+1:end)
        outfilename1 = ['DIMG.',infilename(1:slashpos-1),'.f2dc.cdf'];
        outfilename2 = ['DIMG.',infilename(1:slashpos-1),'.2dp.cdf'];
    else % append probe type and extension to output file path and user-specified file name ~ Joe Finlon 02/09/18
        ouputDirectory = outfilename;
        outfilename = [ouputDirectory, '.2dc.cdf'];
        outfilename1 = [ouputDirectory, '.f2dc.cdf'];
        outfilename2 = [ouputDirectory, '.2dp.cdf'];
    end
    
    fid=fopen(infilename,'r','b');
    infilename
    
    n=4;
    n1=8;
    %%% Updated for new MATLAB NETCDF interface
    f = netcdf.create(outfilename, 'clobber');
    f1 = netcdf.create(outfilename1, 'clobber');
    f2 = netcdf.create(outfilename2, 'clobber');
    
    dimid0 = netcdf.defDim(f,'time',netcdf.getConstant('NC_UNLIMITED'));
    dimid1 = netcdf.defDim(f,'ImgRowlen',n);
    dimid2 = netcdf.defDim(f,'ImgBlocklen',4096/n);
    
    % changed data types to reduce file footprint ~ Joe Finlon 02/09/18
    varid0 = netcdf.defVar(f,'year','short',dimid0);
    varid1 = netcdf.defVar(f,'month','byte',dimid0);
    varid2 = netcdf.defVar(f,'day','byte',dimid0);
    varid3 = netcdf.defVar(f,'hour','byte',dimid0);
    varid4 = netcdf.defVar(f,'minute','byte',dimid0);
    varid5 = netcdf.defVar(f,'second','byte',dimid0);
    varid6 = netcdf.defVar(f,'millisec','short',dimid0);
    varid7 = netcdf.defVar(f,'wkday','short',dimid0);
    varid9 = netcdf.defVar(f,'tas','double',dimid0);
    varid8 = netcdf.defVar(f,'data','int',[dimid1 dimid2 dimid0]);

    dimid02 = netcdf.defDim(f1,'time',netcdf.getConstant('NC_UNLIMITED'));
    dimid12 = netcdf.defDim(f1,'ImgRowlen',n1);
    dimid22 = netcdf.defDim(f1,'ImgBlocklen',4096/n1);
    
    varid02 = netcdf.defVar(f1,'year','short',dimid02);
    varid12 = netcdf.defVar(f1,'month','byte',dimid02);
    varid22 = netcdf.defVar(f1,'day','byte',dimid02);
    varid32 = netcdf.defVar(f1,'hour','byte',dimid02);
    varid42 = netcdf.defVar(f1,'minute','byte',dimid02);
    varid52 = netcdf.defVar(f1,'second','byte',dimid02);
    varid62 = netcdf.defVar(f1,'millisec','short',dimid02);
    varid72 = netcdf.defVar(f1,'wkday','short',dimid02);
    varid92 = netcdf.defVar(f1,'tas','double',dimid02);
    varid82 = netcdf.defVar(f1,'data','int',[dimid12 dimid22 dimid02]);

    
    dimid01 = netcdf.defDim(f2,'time',netcdf.getConstant('NC_UNLIMITED'));
    dimid11 = netcdf.defDim(f2,'ImgRowlen',n);
    dimid21 = netcdf.defDim(f2,'ImgBlocklen',4096/n);
    
    varid01 = netcdf.defVar(f2,'year','short',dimid01);
    varid11 = netcdf.defVar(f2,'month','byte',dimid01);
    varid21 = netcdf.defVar(f2,'day','byte',dimid01);
    varid31 = netcdf.defVar(f2,'hour','byte',dimid01);
    varid41 = netcdf.defVar(f2,'minute','byte',dimid01);
    varid51 = netcdf.defVar(f2,'second','byte',dimid01);
    varid61 = netcdf.defVar(f2,'millisec','short',dimid01);
    varid71 = netcdf.defVar(f2,'wkday','short',dimid01);
    varid91 = netcdf.defVar(f2,'tas','double',dimid01);
    varid81 = netcdf.defVar(f2,'data','int',[dimid11 dimid21 dimid01]);

    netcdf.endDef(f)
    netcdf.endDef(f1)
    netcdf.endDef(f2)
    
    kk=1;
    kk1=1;
    kk2=1;
    kk3=1;
    
    endfile = 0;
    
    xmldoc = '<>';
    while ~isequal('</OAP>',strtrim(xmldoc))
       xmldoc=fgetl(fid);
    end
    
    % while feof(fid)==0 & kk <= 3000
    while feof(fid)==0 && endfile == 0

             
        probetype=fread(fid,2,'uchar');
        hour=fread(fid,1,'uint16');
        minute=fread(fid,1,'uint16');
        second=fread(fid,1,'uint16');
        year=fread(fid,1,'uint16');
        month=fread(fid,1,'uint16');
        day=fread(fid,1,'uint16');
        tas=fread(fid,1,'uint16');
        millisec=fread(fid,1,'uint16');
        wkday=fread(fid,1,'uint16');
        %disp([probetype(1) probetype(2) year month day hour minute second millisec tas]);
        data = fread(fid,4096,'uchar');

        %disp(probetype')

        if ( (probetype(1)==67) && (probetype(2)==49) ) % 67 % 80 % 2DC
            
            netcdf.putVar ( f, varid0, kk-1, 1, year );
            netcdf.putVar ( f, varid1, kk-1, 1, month );
            netcdf.putVar ( f, varid2, kk-1, 1, day );
            netcdf.putVar ( f, varid3, kk-1, 1, hour );
            netcdf.putVar ( f, varid4, kk-1, 1, minute );
            netcdf.putVar ( f, varid5, kk-1, 1, second );
            netcdf.putVar ( f, varid6, kk-1, 1, millisec );
            netcdf.putVar ( f, varid7, kk-1, 1, wkday );
            netcdf.putVar ( f, varid9, kk-1, 1, tas );

            netcdf.putVar ( f, varid8, [0, 0, kk-1], [n,4096/n,1], reshape(data,n,4096/n) );

            kk=kk+1;
            if mod(kk,100) == 0
                kk
                datestr(now)
            end
            clear decomp dd k2 b1 b2
        elseif ( (probetype(1)==67) && (probetype(2)==54) ) % 67 % 80 % F2DC
            
            netcdf.putVar ( f, varid02, kk2-1, 1, year );
            netcdf.putVar ( f, varid12, kk2-1, 1, month );
            netcdf.putVar ( f, varid22, kk2-1, 1, day );
            netcdf.putVar ( f, varid32, kk2-1, 1, hour );
            netcdf.putVar ( f, varid42, kk2-1, 1, minute );
            netcdf.putVar ( f, varid52, kk2-1, 1, second );
            netcdf.putVar ( f, varid62, kk2-1, 1, millisec );
            netcdf.putVar ( f, varid72, kk2-1, 1, wkday );
            netcdf.putVar ( f, varid92, kk2-1, 1, tas );

            netcdf.putVar ( f1, varid82, [0, 0, kk2-1], [n1,4096/n1,1], reshape(data,n1,4096/n1) );

            kk2=kk2+1;
            if mod(kk2,100) == 0
                kk2
                datestr(now)
            end
            clear decomp dd k2 b1 b2
            
        elseif ( (probetype(1)==67) && (probetype(2)==52) ) % 67 % 80 % F2DC
            % changed start index [kk2] to kk3 to match this condition ~
            % Joe Finlon 02/09/18
            netcdf.putVar ( f1, varid02, kk3-1, 1, year );
            netcdf.putVar ( f1, varid12, kk3-1, 1, month );
            netcdf.putVar ( f1, varid22, kk3-1, 1, day );
            netcdf.putVar ( f1, varid32, kk3-1, 1, hour );
            netcdf.putVar ( f1, varid42, kk3-1, 1, minute );
            netcdf.putVar ( f1, varid52, kk3-1, 1, second );
            netcdf.putVar ( f1, varid62, kk3-1, 1, millisec );
            netcdf.putVar ( f1, varid72, kk3-1, 1, wkday );
            netcdf.putVar ( f1, varid92, kk3-1, 1, tas );

            netcdf.putVar ( f1, varid82, [0, 0, kk3-1], [n1,4096/n1,1], reshape(data,n1,4096/n1) );

            kk3=kk3+1;
            if mod(kk3,100) == 0
                kk3
                datestr(now)
            end
            clear decomp dd k3 b1 b2     
            
            
        elseif probetype(1)==80 % 2DP
            netcdf.putVar ( f2, varid01, kk1-1, 1, year );
            netcdf.putVar ( f2, varid11, kk1-1, 1, month );
            netcdf.putVar ( f2, varid21, kk1-1, 1, day );
            netcdf.putVar ( f2, varid31, kk1-1, 1, hour );
            netcdf.putVar ( f2, varid41, kk1-1, 1, minute );
            netcdf.putVar ( f2, varid51, kk1-1, 1, second );
            netcdf.putVar ( f2, varid61, kk1-1, 1, millisec );
            netcdf.putVar ( f2, varid71, kk1-1, 1, wkday );
            netcdf.putVar ( f2, varid91, kk1-1, 1, tas );

            netcdf.putVar ( f2, varid81, [0, 0, kk1-1], [n,4096/n,1], reshape(data,n,4096/n) );

            kk1=kk1+1;
            if mod(kk1,100) == 0
                kk1
                datestr(now)
            end
            clear decomp dd k2 b1 b2
        else
            disp(probetype)
               
        end
        
        for j=1:4116
            bb=fread(fid,1,'int8');
            if feof(fid) == 1
                endfile=1;
                break
            end
        end
        fseek(fid,-4116,'cof');
        
    end
    
    fclose(fid);
    % close(f);
    netcdf.close(f);  % New interface by Will
    netcdf.close(f1);  % New interface by Will
    netcdf.close(f2);  % New interface by Will
end
