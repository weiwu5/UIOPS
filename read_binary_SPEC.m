function read_binary_SPEC(infilename,outfilename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Read the raw base*.2DS file, and then write into NETCDF file 
%% Follow the SPEC manual 
%%  by Will Wu, 08/01/2014
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        outfilename = ['DIMG.',infilename(1:slashpos-1),'.cdf'];
    end
    
    outfilename1=[outfilename, '.H.cdf'];
    outfilename2=[outfilename, '.V.cdf'];
    
    fid=fopen(infilename,'r','l');

    f = netcdf.create(outfilename1, 'clobber');
    
    dimid0 = netcdf.defDim(f,'time',netcdf.getConstant('NC_UNLIMITED'));
    dimid1 = netcdf.defDim(f,'ImgRowlen',8);
    dimid2 = netcdf.defDim(f,'ImgBlocklen',1700);
    
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
        
    f1 = netcdf.create(outfilename2, 'clobber');
    
    dimid01 = netcdf.defDim(f1,'time',netcdf.getConstant('NC_UNLIMITED'));
    dimid11 = netcdf.defDim(f1,'ImgRowlen',8);
    dimid21 = netcdf.defDim(f1,'ImgBlocklen',1700);
    
    varid01 = netcdf.defVar(f1,'year','short',dimid01);
    varid11 = netcdf.defVar(f1,'month','byte',dimid01);
    varid21 = netcdf.defVar(f1,'day','byte',dimid01);
    varid31 = netcdf.defVar(f1,'hour','byte',dimid01);
    varid41 = netcdf.defVar(f1,'minute','byte',dimid01);
    varid51 = netcdf.defVar(f1,'second','byte',dimid01);
    varid61 = netcdf.defVar(f1,'millisec','short',dimid01);
    varid71 = netcdf.defVar(f1,'wkday','byte',dimid01);
    varid81 = netcdf.defVar(f1,'data','int',[dimid11 dimid21 dimid01]);
    netcdf.endDef(f1)
    
    kk1=1;
    kk2=1;
    endfile = 0;
    nNext=1;
    dataprev=zeros(2048,1);   
    
    %fseek(fid,4114*233191,'bof');
    while feof(fid)==0 && endfile == 0 
        %tic
        [year,month, wkday,day, hour, minute, second, millisec, data, discard]=readRecord(fid);         
        timebuffer = [year,month,day, hour, minute, second, millisec]
        %fprintf(formatSpec,A1,A2)
        [year1,month1, wkday1,day1, hour1, minute1, second1, millisec1, data1, discard1]=readRecord(fid);
        %fseek(fid,-4114,'cof');
        datan=[data' data1'];
        datan=datan';
        
        
        [imgH, imgV, nNext]=get_img(datan, hour*10000+minute*100+second+millisec/1000,outfilename);
        sizeimg= size(imgH);
        if sizeimg(2)>1700
            imgH=imgH(:,1:1700);
            sizeimg(2)
        end
        
        sizeimg= size(imgV);
        if sizeimg(2)>1700
            imgV=imgV(:,1:1700);
            sizeimg(2)
        end
        
        if sum(sum(imgH))~=0
            for  mmm=1:8
                img1(mmm,1:1700)=sixteen2int(imgH((mmm-1)*16+1:mmm*16,1:1700));
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
            
            kk1=kk1+1;
            if mod(kk1,1000) == 0
                 ['kk1=' num2str(kk1) ', ' datestr(now)]
            end
        end
        
        if sum(sum(imgV))~=0
            for  mmm=1:8
                img2(mmm,1:1700)=sixteen2int(imgV((mmm-1)*16+1:mmm*16,1:1700));
            end

            netcdf.putVar ( f1, varid01, kk2-1, 1, year );
            netcdf.putVar ( f1, varid11, kk2-1, 1, month );
            netcdf.putVar ( f1, varid21, kk2-1, 1, day );
            netcdf.putVar ( f1, varid31, kk2-1, 1, hour );
            netcdf.putVar ( f1, varid41, kk2-1, 1, minute );
            netcdf.putVar ( f1, varid51, kk2-1, 1, second );
            netcdf.putVar ( f1, varid61, kk2-1, 1, millisec );
            netcdf.putVar ( f1, varid71, kk2-1, 1, wkday );
            netcdf.putVar ( f1, varid81, [0, 0, kk2-1], [8,1700,1], img2 );

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

function [year,month, wkday,day, hour, minute, second, millisec, data, discard]=readRecord(fid)

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

function [imgH, imgV, nNext]=get_img(buf, timehhmmss,outfilename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Decompress the image 
%% Follow the SPEC manual 
%%  by Will Wu, 06/20/2013
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imgH=zeros(128,1700);
    imgV=zeros(128,1700);
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
              elseif nH~=0
                  jjj=1;
                  kkk=0;
                  while jjj<=nS && kkk<nH-2 % Last two slice is time
                      aa=bitand(buf(iii+kkk),16256)/2^7;  %bin2dec('0011111110000000')
                      bb=bitand(buf(iii+kkk),127); %bin2dec('0000000001111111')
                      imgH(min(128,bb+1):min(aa+bb,128),iSlice+jjj)=1;
                      bBase=min(aa+bb,128);
                      kkk=kkk+1;
                      while( bitand(buf(iii+kkk),16384)==0  && kkk<nH-2) % bin2dec('1000000000000000')
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

              elseif nV~=0
                  jjj=1;
                  kkk=0;
                  while jjj<=nS && kkk<nV-2 % Last two slice is time
                      aa=bitand(buf(iii+kkk),16256)/2^7;  %bin2dec('0011111110000000')
                      bb=bitand(buf(iii+kkk),127); %bin2dec('0000000001111111')
                      imgV(min(128,bb+1):min(aa+bb,128),iSlice+jjj)=1;
                      bBase=min(aa+bb,128);
                      kkk=kkk+1;
                      while( bitand(buf(iii+kkk),16384)==0  && kkk<nV-2) % bin2dec('1000000000000000')
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
            tas = typecast( uint32(bin2dec([dec2bin(buf(iii-1+50),16) dec2bin(buf(iii-1+51),16)])) ,'single');
            time = buf(iii-1+52)*2^16+buf(iii-1+53);
            nStereo = buf(iii-1+41);
            nTWM = buf(iii-1+42);
            nSCM = buf(iii-1+43);
            nHOL = buf(iii-1+44);
            nVOL = buf(iii-1+45);
            nVPD = buf(iii-1+34);
            nHPD = buf(iii-1+35);
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
