function imgProc_sm(infile, outfile, probename, n, nEvery, projectname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%  This function is the image processing part of OAP processing using 
%  distributed memory parallisation. The function use one simple interface
%  for all probes. 
%
%  Interface:
%    infile   :   The input file name
%    outfile  :   The output file name
%    probetype:   One of the following: '2DC','2DP','CIP','PIP','HVPS' and '2DS'
%    n        :   The nth chuck to be processed.  
%    nEvery   :   The individual chuck size. nChuck*nEvery shoudl equal the
%                 total frame number 
%    projectname: The name of project so that you can write the specific
%                 code for you data
%
%  Note other important variables used in the program
%    handles:  a structure to store information. It is convinient to use a
%          struture to store the global information rather than using
%          various varibles
%
%  Update Dates:
%    * Initially Written by Will Wu, 06/24/2013 
%          imgprocdm(File,probetype,n)
%    * Updated by Will Wu, 10/11/2013   
%          New function interface 
%          imgprocdm(infile,outfile,probetype,n, nEvery) and updated documentation. 
%          This version is a major update to include all probes and simplify
%          the function interface significantly
%    * Updated by Will Wu, 07/10/2013
%          New function interface imgProc_dm(infile,outfile,probetype,n, nEvery)
%          Output perimeter, rectangle length/width/angle and eclispe
%          length/width/angle
%    * Added by Wei Wu, May 11, 2016
%          Add the project specific code with projectname in the following format:
%            if strcmp(projectname, 'PECAN')  % For example for PECAN dataset 
%               ...
%            end 
%    * Updated by Will Wu, 07/11/2016
%          New function name with the option to turn CGAL on and off for
%          speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%% Setting probe information according to probe type
%    use ProbeType to indicate three type of probes:
%       0: 2DC/2DP, 32 doides, boundary 85, 
%       1: CIP/PIP, 64 doides, boundary 170
%       2: HVPS/2DS, 128 doides, boundary 170
 

iRectEllipse = 0;  % Set defualt to no Rectangle fit and Ellipse fit
switch probename
    case '2DC'
        boundary=[255 255 255 255];
        boundarytime=85;

        ds = 0.025;			     % Size of diode in millimeters
        handles.diodesize = ds;  
        handles.diodenum  = 32;  % Diode number
        handles.current_image = 1;
        probetype=0;

    case '2DP'
        boundary=[255 255 255 255];
        boundarytime=85;

        ds = 0.200;			     % Size of diode in millimeters
        handles.diodesize = ds;  
        handles.diodenum  = 32;  % Diode number
        handles.current_image = 1;
        probetype=0;

    case 'CIP'
        boundary=[170, 170, 170, 170, 170, 170, 170, 170];
        boundarytime=NaN;

        ds = 0.025;			     % Size of diode in millimeters
        handles.diodesize = ds;
        handles.diodenum  = 64;  % Diode number
        handles.current_image = 1;
        probetype=1;

    case 'PIP'
        boundary=[170, 170, 170, 170, 170, 170, 170, 170];
        boundarytime=NaN;

        ds = 0.100;			     % Size of diode in millimeters
        handles.diodesize = ds;
        handles.diodenum  = 64;  % Diode number
        handles.current_image = 1;
        probetype=1;
        
    case 'HVPS'
        boundary=[43690, 43690, 43690, 43690, 43690, 43690, 43690, 43690];
        boundarytime=0;

        ds = 0.150;			     % Size of diode in millimeters
        handles.diodesize = ds;
        handles.diodenum  = 128; % Diode number
        handles.current_image = 1;
        probetype=2;

    case '2DS'
        boundary=[43690, 43690, 43690, 43690, 43690, 43690, 43690, 43690];
        boundarytime=0;

        ds = 0.010;			     % Size of diode in millimeters
        handles.diodesize = ds;
        handles.diodenum  = 128; % Diode number
        handles.current_image = 1;
        probetype=2;
end

diodenum = handles.diodenum;
byteperslice = diodenum/8;  
handles.disagree = 0;

%% Read the particle image files
handles.f = netcdf.open(infile,'nowrite');
[~, dimlen] = netcdf.inqDim(handles.f,2);
[~, handles.img_count] = netcdf.inqDim(handles.f,0);
size_mat = dimlen; 
warning off all
diode_stats = zeros(1,diodenum);

if strcmp(projectname, 'PECAN')  % For example for PECAN dataset 
    disp('Testing...')  %% Add project specific code if you like
end
%% Create output NETCDF file and variables
f = netcdf.create(outfile, 'clobber');
dimid0 = netcdf.defDim(f,'time',netcdf.getConstant('NC_UNLIMITED'));
dimid1 = netcdf.defDim(f,'pos_count',2);
dimid2 = netcdf.defDim(f,'bin_count',diodenum);

varid1 = netcdf.defVar(f,'Date','double',dimid0);
varid0  = netcdf.defVar(f,'Time','double',dimid0);
varid2  = netcdf.defVar(f,'msec','double',dimid0);
varid101  = netcdf.defVar(f,'Time_in_seconds','double',dimid0);

varid102  = netcdf.defVar(f,'SliceCount','double',dimid0);
varid103  = netcdf.defVar(f,'DMT_DOF_SPEC_OVERLOAD','double',dimid0);
varid104  = netcdf.defVar(f,'Particle_number_all','double',dimid0);

%varid3 = netcdf.defVar(f,'wkday','double',dimid0);
varid4  = netcdf.defVar(f,'position','double',[dimid1 dimid0]);
varid5  = netcdf.defVar(f,'particle_time','double',dimid0);
varid6  = netcdf.defVar(f,'particle_millisec','double',dimid0);
varid7  = netcdf.defVar(f,'particle_microsec','double',dimid0);
varid8  = netcdf.defVar(f,'parent_rec_num','double',dimid0);
varid9  = netcdf.defVar(f,'particle_num','double',dimid0);
varid10 = netcdf.defVar(f,'image_length','double',dimid0);                                
varid11 = netcdf.defVar(f,'image_width','double',dimid0);                                 
varid12 = netcdf.defVar(f,'image_area','double',dimid0);                                  
varid13 = netcdf.defVar(f,'image_longest_y','double',dimid0);                             
varid14 = netcdf.defVar(f,'image_max_top_edge_touching','double',dimid0);                 
varid15 = netcdf.defVar(f,'image_max_bottom_edge_touching','double',dimid0);              
varid16 = netcdf.defVar(f,'image_touching_edge','double',dimid0);                         
varid17 = netcdf.defVar(f,'image_auto_reject','double',dimid0);                           
varid18 = netcdf.defVar(f,'image_hollow','double',dimid0);                                
varid19 = netcdf.defVar(f,'image_center_in','double',dimid0);                             
varid20 = netcdf.defVar(f,'image_axis_ratio','double',dimid0);                            
varid21 = netcdf.defVar(f,'image_diam_circle_fit','double',dimid0);                       
varid22 = netcdf.defVar(f,'image_diam_horiz_chord','double',dimid0);                      
varid23 = netcdf.defVar(f,'image_diam_horiz_chord_corr','double',dimid0);                 
varid24 = netcdf.defVar(f,'image_diam_following_bamex_code','double',dimid0);             
varid25 = netcdf.defVar(f,'image_diam_vert_chord','double',dimid0);                       
varid26 = netcdf.defVar(f,'image_diam_minR','double',dimid0);                       
varid27 = netcdf.defVar(f,'image_diam_AreaR','double',dimid0);     
varid45 = netcdf.defVar(f,'image_perimeter','double',dimid0);
if 1==iRectEllipse 
    varid46 = netcdf.defVar(f,'image_RectangleL','double',dimid0);                       
    varid47 = netcdf.defVar(f,'image_RectangleW','double',dimid0);                         
    varid67 = netcdf.defVar(f,'image_RectangleAngle','double',dimid0);                         
    varid48 = netcdf.defVar(f,'image_EllipseL','double',dimid0);                         
    varid49 = netcdf.defVar(f,'image_EllipseW','double',dimid0);                            
    varid69 = netcdf.defVar(f,'image_EllipseAngle','double',dimid0);   
end
varid28 = netcdf.defVar(f,'percent_shadow_area','double',dimid0);                         
varid29 = netcdf.defVar(f,'edge_at_max_hole','double',dimid0);                            
varid30 = netcdf.defVar(f,'max_hole_diameter','double',dimid0);                           
varid31 = netcdf.defVar(f,'part_z','double',dimid0);                                      
varid32 = netcdf.defVar(f,'size_factor','double',dimid0);                                 
varid33 = netcdf.defVar(f,'holroyd_habit','double',dimid0);                               
varid34 = netcdf.defVar(f,'area_hole_ratio','double',dimid0);                             
varid35 = netcdf.defVar(f,'inter_arrival','double',dimid0);                               
varid36 = netcdf.defVar(f,'bin_stats','double',dimid2);                                   
netcdf.endDef(f)

%% Variabels initialization 
kk=1;
w=-1;
wstart = 0;

time_offset_hr = 0;
time_offset_mn = 0;
time_offset_sec = 0;
time_offset_ms = 0;
timeset_flag = 0;



%% Processing nth chuck. Every chuck is nEvery frame
%% Analyze each individual particle images and Outpu the particle by particle information
for i=((n-1)*nEvery+1):min(n*nEvery,handles.img_count)

    handles.year     = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'year'    ),i-1,1);
    handles.month    = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'month'  ),i-1,1);
    handles.day      = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'day'  ),i-1,1);
    handles.hour     = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'hour'    ),i-1,1);
    handles.minute   = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'minute'  ),i-1,1);
    handles.second   = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'second'  ),i-1,1);
    handles.millisec = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'millisec'),i-1,1);
  
    if mod(i,100) == 0
        [num2str(i),'/',num2str(handles.img_count), ', ',datestr(now)]
    end
    varid = netcdf.inqVarID(handles.f,'data');
    
    if probetype==0
        temp = netcdf.getVar(handles.f,varid,[0, 0, i-1], [4,1024,1]);
    else
        temp = netcdf.getVar(handles.f,varid,[0, 0, i-1], [8,1700,1]);
    end
    data(:,:) = temp';  
    
    j=1;
    start=0;
    firstpart = 1;
    
    %c=[dec2bin(data(:,1),8),dec2bin(data(:,2),8),dec2bin(data(:,3),8),dec2bin(data(:,4),8)];
    while data(j,1) ~= -1 && j < size(data,1)
        % Calculate every particles
        if (isequal(data(j,:), boundary) && ( (isequal(data(j+1,1), boundarytime) || probetype==1) ) )
           if start ==0
               if 1 == probetype 
                   start = 2;
               elseif 0 == probetype
                   start = 2;
               else
                   start = 1;
               end
           end
            
               if probetype==0
                   if start+1 > (j-1)  % Remove Corrupted Data
                    break;
                   end
               else
                   if start > (j-1)  % Remove Corrupted Data
                    break;
                   end
               end 
                
                header_loc = j+1;
                w=w+1;
                %% Create binary image according to probe type
                   
                if probetype==0    
                    ind_matrix(1:j-start-1,:) = data(start+1:j-1,:);  % 2DC has 3 slices between particles (sync word timing word and end of particle words)
                    c=[dec2bin(ind_matrix(:,1),8),dec2bin(ind_matrix(:,2),8),dec2bin(ind_matrix(:,3),8),dec2bin(ind_matrix(:,4),8)];
                elseif probetype==1
                    ind_matrix(1:j-start,:) = data(start:j-1,:);
                    c=[dec2bin(ind_matrix(:,1),8), dec2bin(ind_matrix(:,2),8),dec2bin(ind_matrix(:,3),8),dec2bin(ind_matrix(:,4),8), ...
                    dec2bin(ind_matrix(:,5),8), dec2bin(ind_matrix(:,6),8),dec2bin(ind_matrix(:,7),8),dec2bin(ind_matrix(:,8),8)];
                elseif probetype==2
                    ind_matrix(1:j-start,:) = 65535 - data(start:j-1,:); % I used 1 to indicate the illuminated doides for HVPS
                    c=[dec2bin(ind_matrix(:,1),16), dec2bin(ind_matrix(:,2),16),dec2bin(ind_matrix(:,3),16),dec2bin(ind_matrix(:,4),16), ...
                    dec2bin(ind_matrix(:,5),16), dec2bin(ind_matrix(:,6),16),dec2bin(ind_matrix(:,7),16),dec2bin(ind_matrix(:,8),16)];
                end
                
                % Just to test if there is bad images, usually 0 area images
                figsize = size(c);
                if figsize(2)~=diodenum
                    disp('Not equal to doide number');
                    return
                end
                
                
                images.position(kk,:) = [start, j-1];
                parent_rec_num(kk)=i;
                particle_num(kk) = mod(kk,66536); %hex2dec([dec2hex(data(start-1,7)),dec2hex(data(start-1,8))]);
                
                %  Get the particle time 
                if probetype==0 
                    bin_convert = [dec2bin(data(header_loc,2),8),dec2bin(data(header_loc,3),8),dec2bin(data(header_loc,4),8)];
                    part_time = bin2dec(bin_convert);       % Interarrival time in tas clock cycles
                    tas2d = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'tas'),i-1, 1);
                    part_time = part_time/tas2d*handles.diodesize/(10^3);                    
                    time_in_seconds(kk) = part_time;

                    images.int_arrival(kk) = part_time;
                    
                    if(firstpart == 1)
                        firstpart = 0;
                        start_hour = handles.hour;
                        start_minute = handles.minute;
                        start_second = handles.second;
                        start_msec = handles.millisec*10;
                        % First, we get the hours....
                        start_msec = start_msec;
                        start_microsec = 0;
                        time_offset_hr = 0;
                        time_offset_mn = 0;
                        time_offset_sec = 0;
                        time_offset_ms = 0;

                        part_hour(kk) = start_hour;
                        part_min(kk) = start_minute;
                        part_sec(kk) = start_second;
                        part_mil(kk) = start_msec;
                        part_micro(kk) = 0;
                    else
                        frac_time = part_time - floor(part_time);
                        frac_time = frac_time * 1000;
                        part_micro(kk) = part_micro(kk-1) + (frac_time - floor(frac_time))*1000;
                        part_mil(kk) = part_mil(kk-1) + floor(frac_time);
                        part_sec(kk) = part_sec(kk-1) + floor(part_time);
                        part_min(kk) = part_min(kk-1);
                        part_hour(kk) = part_hour(kk-1);
                    end
                    
                    part_sec(part_mil >= 1000) = part_sec(part_mil >= 1000) + 1;
                    part_mil(part_mil >= 1000) = part_mil(part_mil >= 1000) - 1000;

                    part_min(part_sec >= 60) = part_min(part_sec >= 60) + 1;
                    part_sec(part_sec >= 60) = part_sec(part_sec >= 60) - 60;

                    part_hour(part_min >= 60) = part_hour(part_min >= 60) + 1;
                    part_min(part_min >= 60) = part_min(part_min >= 60) - 60;
                    part_hour(part_hour >= 24) = part_hour(part_hour >= 24) - 24;
                elseif probetype==1
                    bin_convert = [dec2bin(data(start-1,2),8),dec2bin(data(start-1,3),8),dec2bin(data(start-1,4),8), ...
                        dec2bin(data(start-1,5),8), dec2bin(data(start-1,6),8)];

                    part_hour(kk) = bin2dec(bin_convert(1:5));
                    part_min(kk) = bin2dec(bin_convert(6:11));
                    part_sec(kk) = bin2dec(bin_convert(12:17));
                    part_mil(kk) = bin2dec(bin_convert(18:27));
                    part_micro(kk) = bin2dec(bin_convert(28:40))*125e-9;
                
                    particle_sliceCount(kk)=bitand(data(start-1,1),127);
                    particle_DOF(kk)=bitand(data(start-1,1),128);
                    particle_partNum(kk)=bin2dec([dec2bin(data(start-1,7),8),dec2bin(data(start-1,8),8)]);

                    time_in_seconds(kk) = part_hour(kk) * 3600 + part_min(kk) * 60 + part_sec(kk) + part_mil(kk)/1000 + part_micro(kk);
                    if kk > 1
                        images.int_arrival(kk) = time_in_seconds(kk) - time_in_seconds(kk-1);
                    else
                        images.int_arrival(kk) = time_in_seconds(kk);
                    end


                elseif probetype==2

                    particle_DOF(kk)=bitand(data(header_loc,4), 32768);
                    particle_partNum(kk)=double(data(header_loc,5));
                    particle_sliceCount(kk)=double(data(header_loc,6));

                    part_time = double(data(header_loc,7))*2^16+double(data(header_loc,8));       % Interarrival time in tas clock cycles
                    part_micro(kk) = part_time;
                    part_mil(kk)   = 0;
                    part_sec(kk)   = 0;
                    part_min(kk)   = 0;
                    part_hour(kk)  = 0;
                    time_in_seconds(kk) = part_time*(handles.diodesize/(10^3)/170);
                    if(kk>1)
                        images.int_arrival(kk) = part_time-part_micro(kk-1); 
                    else
                        images.int_arrival(kk) = 0;
                    end
                end
                
                temptimeinhhmmss = part_hour(kk) * 10000 + part_min(kk) * 100 + part_sec(kk);
                %if (temptimeinhhmmss<200000 || temptimeinhhmmss>240000)
                %    temptimeinhhmmss
                %end
                
                slices_ver = length(start:j-1);
                rec_time(kk)=double(handles.hour)*10000+double(handles.minute)*100+double(handles.second);
                rec_date(kk)=double(handles.year)*10000+double(handles.month)*100+double(handles.day);
                rec_millisec(kk)=handles.millisec;
                %                 rec_wkday(kk)=handles.wkday(i);
   
                %% Determine the Particle Habit
                %  We use the Holroyd algorithm here
                handles.bits_per_slice = diodenum;
                diode_stats = diode_stats + sum(c=='1',1);
                csum = sum(c=='1',1);

                images.holroyd_habit(kk) = holroyd(handles,c);
                
                %% Determine if the particle is rejected or not
                %  Calculate the Particle Length, Width, Area, Auto Reject 
                %  Status And more... See calculate_reject_unified()
                %  funtion for more information
                
                [images.image_length(kk),images.image_width(kk),images.image_area(kk), ...
                    images.longest_y_within_a_slice(kk),images.max_top_edge_touching(kk),images.max_bottom_edge_touching(kk),...
                    images.image_touching_edge(kk), images.auto_reject(kk),images.is_hollow(kk),images.percent_shadow(kk),images.part_z(kk),...
                    images.sf(kk),images.area_hole_ratio(kk),handles]=calculate_reject_unified(c,handles,images.holroyd_habit(kk));

                images.max_hole_diameter(kk) = handles.max_hole_diameter;
                images.edge_at_max_hole(kk) = handles.edge_at_max_hole;

                max_horizontal_length = images.image_length(kk);
                max_vertical_length = images.longest_y_within_a_slice(kk);
                image_area = images.image_area(kk);

                diode_size= handles.diodesize;
                corrected_horizontal_diode_size = handles.diodesize;
                largest_edge_touching  = max(images.max_top_edge_touching(kk), images.max_bottom_edge_touching(kk));
                smallest_edge_touching = min(images.max_top_edge_touching(kk), images.max_bottom_edge_touching(kk));

                %% Calculate more size deciptor using more advanced techniques
                %  See dropsize for more information
                [images.center_in(kk),images.axis_ratio(kk),images.diam_circle_fit(kk),images.diam_horiz_chord(kk),images.diam_vert_chord(kk),...
                    images.diam_horiz_mean(kk), images.diam_spheroid(kk)]=dropsize(max_horizontal_length,max_vertical_length,image_area...
                    ,largest_edge_touching,smallest_edge_touching,diode_size,corrected_horizontal_diode_size, diodenum);
                
                %% Calculate size deciptor using bamex code
                %  See dropsize_new for more information
                % images.diam_bamex(kk) = dropsize_new(c, largest_edge_touching, smallest_edge_touching, diodenum, corrected_horizontal_diode_size, handles.diodesize, max_vertical_length);
                
                %% Using OpenCV C program to calculate length, width and radius. This                 
                %% Get diameter of the smallest-enclosing circle, rectangle and ellipse
                %images.minR(kk)=particlesize_cgal(c);
                images.minR(kk)=CGAL_minR(c);
                images.AreaR(kk)=2*sqrt(images.image_area(kk)/3.1415926);  % Calculate the Darea (area-equivalent diameter)
                images.Perimeter(kk)=ParticlePerimeter(c);
                
                if 1==iRectEllipse 
                    [images.RectangleL(kk), images.RectangleW(kk), images.RectangleAngle(kk)] = CGAL_RectSize(c);
                    [images.EllipseL(kk), images.EllipseW(kk), images.EllipseAngle(kk)]       = CGAL_EllipseSize(c);
                end
                %% Get the area ratio using the DL=max(DT,DP), only observed area are used
                if images.image_length(kk) > images.image_width(kk)
                    images.percent_shadow(kk) = images.image_area(kk) / (pi * images.image_length(kk).^ 2 / 4);
                elseif images.image_width(kk) ~= 0
                    images.percent_shadow(kk) = images.image_area(kk) / (pi * images.image_width(kk).^ 2 / 4);
                else
                    images.percent_shadow(kk) = 0;
                end

                start = j + 2;
                kk = kk + 1;
                clear c ind_matrix
           %end
        end

        j = j + 1;
    end

    %% Write out the processed information on NETCDF
    if kk > 1
        
       
        netcdf.putVar ( f, varid0, wstart, w-wstart+1, rec_time(:) );
        netcdf.putVar ( f, varid1, wstart, w-wstart+1, rec_date(:) );
        
        netcdf.putVar ( f, varid101, wstart, w-wstart+1, time_in_seconds(:) );
        netcdf.putVar ( f, varid102, wstart, w-wstart+1, particle_sliceCount );
        netcdf.putVar ( f, varid103, wstart, w-wstart+1, particle_DOF );
        netcdf.putVar ( f, varid104, wstart, w-wstart+1, particle_partNum );

        
        netcdf.putVar ( f, varid2, wstart, w-wstart+1, rec_millisec(:) );
        %netcdf.putVar ( f, varid3, wstart, w-wstart+1, rec_wkday(:) );
        netcdf.putVar ( f, varid4, [0 wstart], [2 w-wstart+1], images.position' );
        netcdf.putVar ( f, varid5, wstart, w-wstart+1, part_hour(:)*10000+part_min(:)*100+part_sec(:) );
        netcdf.putVar ( f, varid6, wstart, w-wstart+1, part_mil(:) );
        netcdf.putVar ( f, varid7, wstart, w-wstart+1, part_micro(:) );
        netcdf.putVar ( f, varid8, wstart, w-wstart+1, parent_rec_num );    
        netcdf.putVar ( f, varid9, wstart, w-wstart+1, particle_num(:) );
        netcdf.putVar ( f, varid10, wstart, w-wstart+1, images.image_length);                         
        netcdf.putVar ( f, varid11, wstart, w-wstart+1, images.image_width);                          
        netcdf.putVar ( f, varid12, wstart, w-wstart+1, images.image_area*diode_size*diode_size);                           
        netcdf.putVar ( f, varid13, wstart, w-wstart+1, images.longest_y_within_a_slice);             
        netcdf.putVar ( f, varid14, wstart, w-wstart+1, images.max_top_edge_touching);                
        netcdf.putVar ( f, varid15, wstart, w-wstart+1, images.max_bottom_edge_touching); 
        netcdf.putVar ( f, varid16, wstart, w-wstart+1, images.image_touching_edge-'0');                  
        netcdf.putVar ( f, varid17, wstart, w-wstart+1, double(images.auto_reject));                  
        netcdf.putVar ( f, varid18, wstart, w-wstart+1, images.is_hollow);                            
        netcdf.putVar ( f, varid19, wstart, w-wstart+1, images.center_in);                            
        netcdf.putVar ( f, varid20, wstart, w-wstart+1, images.axis_ratio);                           
        netcdf.putVar ( f, varid21, wstart, w-wstart+1, images.diam_circle_fit);                      
        netcdf.putVar ( f, varid22, wstart, w-wstart+1, images.diam_horiz_chord);                     
        netcdf.putVar ( f, varid23, wstart, w-wstart+1, images.diam_horiz_chord ./ images.sf);        
        netcdf.putVar ( f, varid24, wstart, w-wstart+1, images.diam_horiz_mean);              
        netcdf.putVar ( f, varid25, wstart, w-wstart+1, images.diam_vert_chord);                           
        netcdf.putVar ( f, varid26, wstart, w-wstart+1, images.minR*diode_size);                      
        netcdf.putVar ( f, varid27, wstart, w-wstart+1, images.AreaR*diode_size);       
        netcdf.putVar ( f, varid45, wstart, w-wstart+1, images.Perimeter*diode_size); 
        if 1==iRectEllipse 
            netcdf.putVar ( f, varid46, wstart, w-wstart+1, images.RectangleL*diode_size);                      
            netcdf.putVar ( f, varid47, wstart, w-wstart+1, images.RectangleW*diode_size);         
            netcdf.putVar ( f, varid67, wstart, w-wstart+1, images.RectangleAngle);         
            netcdf.putVar ( f, varid48, wstart, w-wstart+1, images.EllipseL*diode_size);                      
            netcdf.putVar ( f, varid49, wstart, w-wstart+1, images.EllipseW*diode_size); 
            netcdf.putVar ( f, varid69, wstart, w-wstart+1, images.EllipseAngle); 
        end
        netcdf.putVar ( f, varid28, wstart, w-wstart+1, images.percent_shadow);                       
        netcdf.putVar ( f, varid29, wstart, w-wstart+1, images.max_hole_diameter);                    
        netcdf.putVar ( f, varid30, wstart, w-wstart+1, images.edge_at_max_hole);                     
        netcdf.putVar ( f, varid31, wstart, w-wstart+1, images.part_z);                               
        netcdf.putVar ( f, varid32, wstart, w-wstart+1, images.sf);                                   
        netcdf.putVar ( f, varid33, wstart, w-wstart+1, double(images.holroyd_habit));                
        netcdf.putVar ( f, varid34, wstart, w-wstart+1, images.area_hole_ratio);                      
        netcdf.putVar ( f, varid35, wstart, w-wstart+1, images.int_arrival);                          
        netcdf.putVar ( f, varid36, diode_stats );
        
        wstart = w+1;
        kk = 1;
        clear rec_time rec_date rec_millisec part_hour part_min part_sec part_mil part_micro parent_rec_num particle_num images time_in_seconds particle_sliceCount particle_DOF particle_partNum

    end
    clear images
end
warning on all

netcdf.close(f);
end
