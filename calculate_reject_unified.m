function [p_length,width,area,longest_y,max_top,max_bottom,touching_edge,reject_status,is_hollow,percent_shadow_area,part_z,size_factor,area_hole_ratio,handles]=calculate_reject_unified(image_buffer,handles,habit)

% /*  RETURN CODE                                               */
% /* 0 = not rejected                                           */
% /* 'a' = reject max. aspect ratio							  */
% /* 't' = reject max. aspect ratio touch edg					  */
% /* 'p' = reject percent shadowed area      					  */
% /*	'h' = reject Hollow image								  */
% /*	's' = reject split image								  */
% /*	'z' = reject 0 area image                                 */
% /*	'f' = reject fake 0 area image                            */


z_d = 0 : .05 : 8.15;


part_z = -1;
size_factor = 1;


edge_0 = [1.000, 1.054, 1.083, 1.101, 1.095, 1.110, 1.148, 1.162, 1.155, 1.123, ...
    1.182, 1.121, 1.162, 1.210, 1.242, 1.134, 1.166, 1.202, 1.238, 1.270, ...
    1.294, 1.278, 1.130, 1.148, 1.170, 1.194, 1.218, 1.242, 1.265, 1.288, ...
    1.310, 1.331, 1.351, 1.369, 1.386, 1.400, 1.411, 1.416, 1.407, 1.074, ...
    1.080, 1.087, 1.096, 1.106, 1.117, 1.127, 1.139, 1.150, 1.162, 1.173, ...
    1.185, 1.197, 1.208, 1.220, 1.232, 1.243, 1.255, 1.266, 1.277, 1.289, ...
    1.300, 1.311, 1.322, 1.333, 1.344, 1.355, 1.366, 1.376, 1.387, 1.397, ...
    1.407, 1.418, 1.428, 1.438, 1.448, 1.458, 1.467, 1.477, 1.486, 1.496, ...
    1.505, 1.515, 1.524, 1.533, 1.542, 1.551, 1.559, 1.568, 1.577, 1.585, ...
    1.594, 1.602, 1.610, 1.618, 1.626, 1.634, 1.642, 1.650, 1.657, 1.665, ...
    1.673, 1.680, 1.687, 1.694, 1.702, 1.709, 1.716, 1.722, 1.729, 1.736, ...
    1.742, 1.749, 1.755, 1.761, 1.768, 1.774, 1.780, 1.786, 1.791, 1.797, ...
    1.803, 1.808, 1.813, 1.819, 1.824, 1.829, 1.834, 1.839, 1.843, 1.848, ...
    1.852, 1.857, 1.861, 1.865, 1.869, 1.872, 1.876, 1.880, 1.883, 1.886, ...
    1.889, 1.892, 1.895, 1.897, 1.899, 1.901, 1.903, 1.905, 1.906, 1.907, ...
    1.908, 1.908, 1.908, 1.908, 1.907, 1.905, 1.903, 1.900, 1.897, 1.892, ...
    1.885, 1.877, 1.865, 1.845];


spot_edge = [0.003, 0.008, 0.017, 0.024, 0.033, 0.040, 0.047, 0.054, 0.062, 0.072, ...
    0.076, 0.088, 0.093, 0.096, 0.101, 0.119, 0.123, 0.127, 0.130, 0.134, ...
    0.139, 0.148, 0.175, 0.180, 0.184, 0.188, 0.192, 0.195, 0.199, 0.202, ...
    0.206, 0.209, 0.213, 0.217, 0.221, 0.225, 0.230, 0.235, 0.243, 0.327, ...
    0.334, 0.340, 0.345, 0.351, 0.355, 0.360, 0.365, 0.369, 0.373, 0.377, ...
    0.381, 0.385, 0.389, 0.393, 0.397, 0.400, 0.404, 0.408, 0.411, 0.415, ...
    0.419, 0.422, 0.426, 0.429, 0.433, 0.436, 0.439, 0.443, 0.446, 0.450, ...
    0.453, 0.457, 0.460, 0.463, 0.467, 0.470, 0.473, 0.477, 0.480, 0.484, ...
    0.487, 0.490, 0.494, 0.497, 0.501, 0.504, 0.507, 0.511, 0.514, 0.518, ...
    0.521, 0.525, 0.528, 0.532, 0.535, 0.539, 0.543, 0.547, 0.550, 0.554, ...
    0.558, 0.562, 0.566, 0.569, 0.572, 0.575, 0.578, 0.581, 0.584, 0.587, ...
    0.590, 0.593, 0.596, 0.598, 0.601, 0.605, 0.610, 0.614, 0.618, 0.623, ...
    0.627, 0.631, 0.635, 0.640, 0.644, 0.648, 0.653, 0.657, 0.662, 0.666, ...
    0.671, 0.676, 0.680, 0.685, 0.690, 0.695, 0.700, 0.705, 0.711, 0.716, ...
    0.721, 0.727, 0.733, 0.738, 0.744, 0.751, 0757, 0.763, 0.770, 0.777, ...
    0.784, 0.792, 0.800, 0.808, 0.817, 0.826, 0.836, 0.846, 0.858, 0.870, ...
    0.884, 0.901, 0.921, 0.950];




% temp = [dec2bin(image_buffer(:,1),16),dec2bin(image_buffer(:,2),16),dec2bin(image_buffer(:,3),16),...
%     dec2bin(image_buffer(:,4),16),dec2bin(image_buffer(:,5),16),dec2bin(image_buffer(:,6),16),...
%     dec2bin(image_buffer(:,7),16),dec2bin(image_buffer(:,8),16)];
% clear image_buffer
% image_buffer(:,:)=temp(:,:);

n_size=size(image_buffer);
n_slices=n_size(1);


handles.rej_zero_area = 1;
handles.rej_split = 1;
handles.rej_hollow = 1;
handles.bits_per_slice = n_size(2);
handles.shadowed_area = 25;
handles.max_edge_img_ar = 6;
handles.max_comp_img_ar = 5;

handles.max_hole_diameter = 0;
handles.edge_at_max_hole = 0;



min_length=-1;
max_length=-1;
max_width=1;
min_width=n_size(2);

total_area=0;

touch=0;
width = 0;
ndrops=0;
split=0;
hollow=0;
met_image=0;
is_hollow=0;

aspect_ratio=0;
percent_shadow_area=0;
area_hole_ratio = 0;


area=0;
max_top=0;
max_bottom=0;
longest_y=0;
touching_edge=0;

top_min_x=-1;
top_max_x=-1;
bottom_min_x=-1;
bottom_max_x=-1;

if n_slices==0
    p_length=0;
    area=0;



    if handles.rej_zero_area==1
        reject_status='z';
        return;
    else
        reject_status='0';
    end

else

    for i=1:n_slices

        [min_pos_lite,max_pos_lite,n_lite]=scan_slice(image_buffer(i,:),handles);

        if longest_y < n_lite
            longest_y=n_lite;
        end

        if i>1
            vertical_split=vertical_split & image_buffer(i,:);
        else
            vertical_split=image_buffer(i,:);
        end

        if n_lite>0
            if touch==0
                if max_pos_lite==handles.bits_per_slice || min_pos_lite==1
                    touch=1;
                end
            end

            if min_pos_lite==1
                if bottom_min_x==-1
                    bottom_min_x=i;
                end
                if bottom_max_x<i
                    bottom_max_x=i;
                end
            end

            if max_pos_lite==handles.bits_per_slice
                if top_min_x==-1
                    top_min_x=i;
                end
                if top_max_x<i
                    top_max_x=i;
                end

            end

            if max_pos_lite > max_width
                max_width=max_pos_lite;
            end
            if min_pos_lite < min_width
                min_width=min_pos_lite;
            end

            if min_length == -1
                min_length = i;
                max_length = i;
            else
                max_length = i;
            end


            total_area=n_lite+total_area;
        end

        if met_image == 0 & n_lite > 0
            met_image=1;
            ndrops=ndrops+1;
            if ndrops > 1 & handles.rej_split==1
                split=1;
            end
        elseif met_image==1 & n_lite == 0
            met_image=0;
        end

    end


    area=total_area;
    if top_min_x == -1
        max_top = 0;
    else
        max_top = (top_max_x - top_min_x) + 1;
    end
    if bottom_min_x == -1
        max_bottom = 0;
    else
        max_bottom = (bottom_max_x - bottom_min_x) + 1;
    end
    if touch == 1
        touching_edge = 't';
    else
        touching_edge = '0';
    end

    

    if total_area == 0;
        p_length = 0;
        width = 0;
    else
        p_length = max_length - min_length + 1;
        width = max_width - min_width + 1;
    end
    
    if total_area > .8 * handles.bits_per_slice * n_slices
        reject_status = 'A';
        return;
    end

    if split == 1
        reject_status = 's';
        return;
    end

    if exist('vertical_split') == 1
        [min_pos_lite,max_pos_lite,n_lite]=scan_slice(vertical_split,handles);
    else
        min_pos_lite=0;
        max_pos_lite=0;
        n_lite=0;
    end
    if n_lite == 1 & n_lite ~= max_pos_lite - min_pos_lite + 1 & handles.rej_split == 0
        reject_status = 's';
        return;
    end

    if total_area > 0
        if p_length > 0 & width > 0
            aspect_ratio =  p_length / width;
            percent_shadow_area = total_area / (p_length * width ) * 100;
        else
            aspect_ratio = 0;
            percent_shadow_area = 0;
        end
    else
        aspect_ratio = 0;
        percent_shadow_area = 0;
    end

    if total_area == 0 && handles.rej_zero_area == 1
        reject_status = 'z';
        return;
    elseif total_area == 0 && handles.rej_zero_area == 0
        reject_status = '0';
        return;
    elseif ( aspect_ratio > handles.max_comp_img_ar || aspect_ratio < 1/handles.max_comp_img_ar ) % Second critirion added on Dec 2nd, 2013 by Will for small aspect ratio
        reject_status = 'a';
        return;
    elseif touch == 1 && aspect_ratio > handles.max_edge_img_ar
        reject_status = 't';
        return;
    elseif percent_shadow_area < handles.shadowed_area
        reject_status = 'p';
        return;
    elseif handles.rej_hollow == 1
        [hollow_status,edge_at_max_hole,max_hole_diameter]=is_it_hollow(image_buffer(1:n_slices,:),n_slices,handles);
        [hollow_status2,edge_at_max_hole2,max_hole_diameter2]=is_it_hollow(image_buffer(n_slices:-1:1,:),n_slices,handles);

        [hollow_status_side1,edge_at_max_hole_side1,max_hole_diameter_side1]=is_it_hollow_sidescan(image_buffer(1:n_slices,:)',n_slices,handles);
        [hollow_status_side2,edge_at_max_hole_side2,max_hole_diameter_side2]=is_it_hollow_sidescan(image_buffer(1:n_slices,32:-1:1)',n_slices,handles);

        
        if hollow_status ~= hollow_status2
            hollow_status;
%             handles.disagree = handles.disagree + 1;
        end
% 
%         if hollow_status == 1 & hollow_status2 == 0
%             hollow_status = 0;
%         end



if hollow_status + hollow_status2 == 1 & (habit == 's' | habit == 'h' | habit == 'i' | habit == 't')
%if hollow_status + hollow_status2 == 1
    hollow_status = 1;
elseif habit == 'd' & percent_shadow_area < 35 & hollow_status + hollow_status2 == 1
    hollow_status = 1;
elseif hollow_status == 1 & hollow_status2 == 0
    hollow_status = 0;
    
end

if percent_shadow_area < 30
    hollow_status = 0;
end


        if hollow_status == 1
            %         ratio = max_hole_diameter./(max_width-min_width+1);
            if edge_at_max_hole <= 0
                ratio = 0;
            else
                ratio = max_hole_diameter./edge_at_max_hole;
            end
            if ratio == 0
                part_z = 0;
                size_factor = 1;
                reject_status = 'h';
                area_hole_ratio = 0;
            elseif max_hole_diameter <= 1
                part_z = 0;
                size_factor = 1;
                area_hole_ratio = 0;
                reject_status = '0';
            else
                part_z = z_d(find(spot_edge < ratio,1,'last'));
                size_factor = edge_0(find(z_d <= part_z,1,'last'));
                reject_status = 'H';
                area_hole_ratio = area/max_hole_diameter;
                
                if hollow_status_side1 + hollow_status_side2 < 1
                    part_z = 0;
                    size_factor = 1;
                    reject_status = 'i';
                end
                
                if area_hole_ratio > 20 & habit == 'i'
                    part_z = 0;
                    size_factor = 1;
                    reject_status = 'u';
                elseif area_hole_ratio > 35 & habit == 'h'
                    part_z = 0;
                    size_factor = 1;
                    reject_status = 'u';
                elseif area_hole_ratio > 40
                    part_z = 0;
                    size_factor = 1;
                    reject_status = 'u';
                end
                
            end
            handles.edge_at_max_hole = edge_at_max_hole;
            handles.max_hole_diameter = max_hole_diameter;
            return
        else
            ratio = -1;
            part_z = -1;
            size_factor = 1;
            area_hole_ratio = 0;

        end
%         if hollow_status ==1
%             reject_status='h';
%             return
%         end
    end
    reject_status = '0';
    return


end



function [min_pos_lite,max_pos_lite,n_lite]=scan_slice(image_buf,handles)

n_lite=0;
max_pos_lite=0;
min_pos_lite=0;


zeros = find(image_buf == '0');
n_lite = length(zeros);
if n_lite == 0
    return
else
    min_pos_lite = zeros(1);
    max_pos_lite = zeros(n_lite);
end

% for i=1:handles.bits_per_slice
%     if image_buf(i) == '0'
%         n_lite=n_lite+1;
%         if min_pos_lite==0
%             min_pos_lite=i;
%             max_pos_lite=i;
%         else
%             max_pos_lite=i;
%         end
%     end
% end
return

function [status,edge_at_max_hole,max_hole_diameter] = is_it_hollow(image_buf,slices,handles)

current = 0;
old = 0;
new = 0;

possibly_hollow = 0;

max_hole_diameter = 0;
edge_at_max_hole = 0;
status = 0;

start_img = 0;
end_img = 0;
i = 1;
while end_img == 0
    zero_amt = sum(image_buf(i,:) == '0');
    if zero_amt > 0 & start_img == 0
        start_img = i;
    end
    if zero_amt == 0 & start_img > 0
        end_img = i;
    end
    i = i + 1;
    if i > slices
        if start_img == 0
            start_img = 1;
        end
        if end_img == 0
            end_img = slices;
        end
    end

end

slices = end_img-start_img+1;

for i=start_img:end_img

    [min_pos_lite,max_pos_lite,n_lite]=scan_slice(image_buf(i,:),handles);

    num_empty = max_pos_lite - min_pos_lite + 1 - n_lite;
      
    if slices > 6
        slices_third = floor(slices/3);
    else
        slices_third = 1;
    end
    
    
    if num_empty > max_hole_diameter & i > slices_third & i < slices - slices_third
        max_hole_diameter = num_empty;
        edge_at_max_hole = max_pos_lite - min_pos_lite + 1;
    end

    if possibly_hollow == 1 & status == 0
        if n_lite > 0 & n_lite ~= max_pos_lite - min_pos_lite + 1
            new = bin2dec(image_buf(i,17:32)) + bitshift(bin2dec(image_buf(i,1:16)),16);
            olddec = bin2dec(old(17:32)) + bitshift(bin2dec(old(1:16)),16);

            newandold = bitand(new , olddec);
            if newandold == zeros

                status=1;
%                 return
            else
                old = mask_start_end(max_pos_lite, min_pos_lite, image_buf(i,:),handles.bits_per_slice);
            end

        elseif n_lite > 0
            bufdec = bin2dec(image_buf(i,17:32)) + bitshift(bin2dec(image_buf(i,1:16)),16);
            olddec = bin2dec(old(17:32)) + bitshift(bin2dec(old(1:16)),16);

            bufdec1 = bin2dec(image_buf(i,1:16));
            olddec1 = bin2dec(old(1:16));

            bufdec2 = bin2dec(image_buf(i,17:32));
            olddec2 = bin2dec(old(17:32));
            
            bufandold1 = bitand(bufdec1,olddec1);
            bufandold2 = bitand(bufdec2,olddec2);
            
            hole_size = length(find(old == '1'));
            cover_size = length(find(dec2bin(bufandold1) == '1')) + length(find(dec2bin(bufandold2) == '1'));
            
            
%             
%             bufandold = bitand(bufdec , olddec);
%             hole_size = length(find(old == '1'));
%             cover_size = length(find(dec2bin(bufandold) == '1'));
%             
%             
            

            if bufandold1 + bufandold2 > 0
%                 if cover_size <= 2 & hole_size ~=1
if cover_size <= .65*hole_size
                    status = 1;
%                     return;
                end
                possibly_hollow = 0;
                old = 0;
            elseif i > 1
                status = 1;
%                 return;
            else
                possibly_hollow = 0;
                old = 0;
            end

        else
            possibly_hollow = 0;

        end
    elseif status == 0
        if n_lite > 0 & n_lite ~= max_pos_lite - min_pos_lite + 1
            old = mask_start_end(max_pos_lite, min_pos_lite, image_buf(i,:),handles.bits_per_slice);

            possibly_hollow = 1;
        end
    end
end
return;

function [status,edge_at_max_hole,max_hole_diameter] = is_it_hollow_sidescan(image_buf,slices,handles)

current = 0;
old = 0;
new = 0;

possibly_hollow = 0;

max_hole_diameter = 0;
edge_at_max_hole = 0;
status = 0;

im_width = size(image_buf);

if im_width(2) > 32
    im_width(2) = 32;
end


start_img = 0;
end_img = 0;
i = 1;

    
while end_img == 0
    zero_amt = sum(image_buf(i,:) == '0');
    if zero_amt > 0 && start_img == 0
        start_img = i;
    end
    if zero_amt == 0 && start_img > 0
        end_img = i;
    end
    i = i + 1;
    if i > 32
        if start_img == 0
            start_img = 1;
        end
        if end_img == 0
            end_img = 32;
        end
    end

end

slices = end_img-start_img+1;

for i=start_img:end_img

    [min_pos_lite,max_pos_lite,n_lite]=scan_slice(image_buf(i,:),handles);

    num_empty = max_pos_lite - min_pos_lite + 1 - n_lite;
      
    if slices > 6
        slices_third = floor(slices/3);
    else
        slices_third = 1;
    end
    
    
    if num_empty > max_hole_diameter & i > slices_third & i < slices - slices_third
        max_hole_diameter = num_empty;
        edge_at_max_hole = max_pos_lite - min_pos_lite + 1;
    end

    if possibly_hollow == 1 & status == 0
        if n_lite > 0 & n_lite ~= max_pos_lite - min_pos_lite + 1
            if im_width(2) <= 16
                new = bin2dec(image_buf(i,1:im_width(2)));
                olddec = bin2dec(old(1:im_width(2)));
            
            else
                new = bin2dec(image_buf(i,17:im_width(2))) + bitshift(bin2dec(image_buf(i,1:16)),16);
                olddec = bin2dec(old(17:im_width(2))) + bitshift(bin2dec(old(1:16)),16);
            end
            
            newandold = bitand(new , olddec);
            if newandold == zeros

                status=1;
%                 return
            else
                old = mask_start_end(max_pos_lite, min_pos_lite, image_buf(i,:),im_width(2));
            end

        elseif n_lite > 0
            
            %             bufdec = bin2dec(image_buf(i,33:64)) + bitshift(bin2dec(image_buf(i,1:32)),32);
            %             olddec = bin2dec(old(33:64)) + bitshift(bin2dec(old(1:32)),32);


            if im_width(2) <= 16
                
                bufdec1 = bin2dec(image_buf(i,1:im_width(2)));
                olddec1 = bin2dec(old(1:im_width(2)));
                
                bufdec2 = 0;
                olddec2 = 0;

            else
                bufdec1 = bin2dec(image_buf(i,1:16));
                olddec1 = bin2dec(old(1:16));
                
                bufdec2 = bin2dec(image_buf(i,17:im_width(2)));
                olddec2 = bin2dec(old(17:im_width(2)));
            end
            
            bufandold1 = bitand(bufdec1,olddec1);
            bufandold2 = bitand(bufdec2,olddec2);
            
            hole_size = length(find(old == '1'));
            cover_size = length(find(dec2bin(bufandold1) == '1')) + length(find(dec2bin(bufandold2) == '1'));
            
            
%             
%             bufandold = bitand(bufdec , olddec);
%             hole_size = length(find(old == '1'));
%             cover_size = length(find(dec2bin(bufandold) == '1'));
%             
%             
            

            if bufandold1 + bufandold2 > 0
%                 if cover_size <= 2 & hole_size ~=1
if cover_size <= .65*hole_size
                    status = 1;
%                     return;
                end
                possibly_hollow = 0;
                old = 0;
            elseif i > 1
                status = 1;
%                 return;
            else
                possibly_hollow = 0;
                old = 0;
            end

        else
            possibly_hollow = 0;

        end
    elseif status == 0
        if n_lite > 0 & n_lite ~= max_pos_lite - min_pos_lite + 1
            old = mask_start_end(max_pos_lite, min_pos_lite, image_buf(i,:),im_width(2));

            possibly_hollow = 1;
        end
    end
end
return;

function old = mask_start_end(end_mask, start_mask, to_mask,bits_per_slice)

result=0;

if start_mask == 0 & end_mask == 0
    result = to_mask;

else
    %     to_mask_dec = bin2dec(to_mask(33:64)) + bitshift(bin2dec(to_mask(1:32)),32);
    %     result = bitshift(bitshift(to_mask_dec,start_mask),-start_mask);
    %     result = bitshift(bitshift(result,-(bits_per_slice - end_mask) + 1),(bits_per_slice - end_mask)+1);
    %     result = dec2bin(result,bits_per_slice);

    result(1:bits_per_slice) = '0';
    result(start_mask:end_mask) = to_mask(start_mask:end_mask);

end

old=char(result);
return;




