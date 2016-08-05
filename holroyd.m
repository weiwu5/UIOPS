% holroyd - identified particle habit according to Holroyd  (1987)
% inputs:
%    handles - handles structure outlined in run_img_processing.m
%    image_buffer - n x photodiodes/8 raw image buffer without timestamps
% outputs:
%    holroyd_habit - habit code as listed below
function [holroyd_habit] = holroyd(handles, image_buffer)

%/***************************************************************/%
%/*	Return code                                                */
%/*	                                                           */
%/*	                                                           */
%/*	reference: J. Atmos. and Oceanic Tech. Vol 4, Sept. '87    */
%/*	           pages 498- 511.                                 */
%/*	                                                           */
%/*	'M' = not calculated, zero image                           */
%/*	'C' = not calculated, center is out                        */
%/*	't' = tiny                                                 */
%/*	'o' = oriented                                             */
%/*	'l' = linear                                               */
%/*	'a' = aggregate                                            */
%/*	'g' = graupel                                              */
%/*	's' = spherel                                              */
%/*	'h' = hexagonal                                            */
%/*	'i' = irregular                                            */
%/*	'd' = dendrite                                             */
%/*	                                                           */
%/***************************************************************/


	image_size = size(image_buffer);
    probe_resolution = .025;
	n_slices  = image_size(1);
	

	if (n_slices == 0)
		holroyd_habit = 'M';
		return;
	else
		if (parabola_fit_center_is_in(image_buffer, n_slices) == 1) 
			[x_length, y_length, d_length, w_width,a_angle,area,r2_correlation, F_fine_detail, S_ratio] = calc_stat(handles,image_buffer, n_slices);

			

			if (area == 0 )
				holroyd_habit = 'M';
				return;
			elseif (area < 25)
				holroyd_habit = 't';
				return;
			elseif  (r2_correlation >= .4) | ( (d_length < 64) & ( (x_length >= 4*y_length) |  (y_length >= 4*x_length)))
				
				if ((a_angle> 30.0) & (a_angle < 60.0))
					holroyd_habit = 'o';
					return;
				else
					holroyd_habit = 'l';
					return;
				end
			elseif ( (d_length * probe_resolution > 6.4 ) |	(d_length > 160.0)) 
				holroyd_habit = 'a';
				return;
			elseif (S_ratio >= .7)
				holroyd_habit = 'g';
				return;
            elseif (d_length >= 64)
				if (F_fine_detail <= 13)
					holroyd_habit = 'g';
					return;
				else
					holroyd_habit = 'a';
					return;
				end
			elseif (F_fine_detail < 5.5) 
				holroyd_habit = 's';
				return;
			elseif (F_fine_detail < 10.0)
				if (d_length >= 32)
					holroyd_habit = 'g';
					return;
				else
					holroyd_habit = 'h';
					return;
				end
			elseif ((F_fine_detail < 16.0) | (x_length <= 7.0))
				holroyd_habit = 'i';
				return;
			else
				holroyd_habit = 'd';
				return;
			end
				
		else
			holroyd_habit = 'C';
			return;
        end
    end
end

%/*************************************************************************/
		


%/*************************************************************************/
function [x_length, y_length, d_length, w_width,a_angle,area,r2_correlation, F_fine_detail, S_ratio] = calc_stat(handles, image_buffer, n_slices)

	BITS_PER_SLICE = handles.bits_per_slice;
    MAX_TWOD_DATA_LENGTH = 6000;



	area = 0.0;
	n_count = 0;
	sum_x2= 0.0;
	sum_y2= 0.0;
	sum_x = 0.0;
	sum_y = 0.0;
	sum_xy= 0.0;
	cross_x2= 0.0;
	cross_y2= 0.0;
	cross_xy= 0.0;
	p_perimeter_change = 0;
	min_x = MAX_TWOD_DATA_LENGTH*3;
	min_y = BITS_PER_SLICE;
	max_x = 0;
	max_y = 0;
    
	spot_on_off = 0;
	fully_on_count = 0;
	partial_on_count =0;

	if (n_slices <= 0)
		return;
	end
	
	for i=1:n_slices
		fully_on_temp = 0;
		for j=1:BITS_PER_SLICE
			if ((image_buffer(i,j)) == '0')     

				tx = i;
				ty =  j;
				if (tx > max_x) 
					max_x = tx;
				end
				
				if (tx < min_x) 
					min_x = tx;
				end
				if (ty > max_y) 
					max_y = ty;
				end
				if (ty < min_y) 
					min_y = ty;
				end
				sum_x2 = sum_x2 + tx * tx;
				sum_y2 = sum_y2 + ty * ty;
				sum_x  = sum_x  + tx;
				sum_y  = sum_y  + ty;
				sum_xy = sum_xy + tx * ty;

                n_count = n_count + 1;
				p(n_count).x = tx;
				p(n_count).y = ty;



				fully_on_temp = fully_on_temp + 1;

				if (spot_on_off == 0)
					spot_on_off = 1;
					p_perimeter_change = p_perimeter_change + 1;
				end
			else
				if spot_on_off == 1 
					spot_on_off = 0;
					p_perimeter_change = p_perimeter_change + 1;
				end
			end
		end

		if (fully_on_temp == BITS_PER_SLICE)
			fully_on_count = fully_on_count + 1;
		end
		if (fully_on_temp ~= 0)
			partial_on_count = partial_on_count + 1;
		end
    end
	area = n_count;

%/*** scan the other way for perimeter change ****/

	spot_on_off = 0;
	for j=1:BITS_PER_SLICE
		for i=1:n_slices
			if ((image_buffer(i,j)) == '0')     
				if (spot_on_off == 0)
					spot_on_off = 1;
					p_perimeter_change = p_perimeter_change + 1;
				end
			else
				if (spot_on_off == 1)
					spot_on_off = 0;
					p_perimeter_change = p_perimeter_change + 1;
				end
			end
		end
	end

	if (max_x >= min_x) 
		x_length = max_x - min_x +1;
	else
		x_length = 0.0;
	end

	if (max_y >= min_y) 
		y_length = max_y - min_y +1;
	else
		y_length = 0.0;
	end
		
	cross_xy = sum_xy - (sum_x * sum_y / area);
	cross_x2 = sum_x2 - (sum_x * sum_x / area);
	cross_y2 = sum_y2 - (sum_y * sum_y / area);

	slope = cross_xy / cross_x2;
	intercept = (sum_y/(area)) - slope * (sum_x/(area));

	angle_radian = atan(slope);
	a_angle = atan(slope) * (180.0/pi);

	if (a_angle < 0) 
		a_angle = a_angle + 180.0;
		angle_radian = angle_radian + pi;
	end



	dmin_x = MAX_TWOD_DATA_LENGTH*3;
	dmin_y = BITS_PER_SLICE;
	dmax_x = 0;
	dmax_y = 0;

	if ( (angle_radian > (pi/2.0)) & (angle_radian <= (pi))) 
		angle_radian = (pi - angle_radian);
	elseif ( angle_radian > pi) 
		['HEY: something is wrong here  a_angle = ', num2str(a_angle)]; 
        return
	end
	for i=1:n_count
		new_x = (p(i).x * cos(angle_radian)) + (p(i).y * sin(angle_radian));
		new_y = (p(i).y * cos(angle_radian)) - (p(i).x * sin(angle_radian));
		if (new_x > dmax_x) 
			dmax_x = new_x;
		end
		if (new_y > dmax_y) 
			dmax_y = new_y;
		end
		if (new_x < dmin_x) 
			dmin_x = new_x;
		end
		if (new_y < dmin_y) 
			dmin_y = new_y;
		end
	end
	
	d_length = (dmax_x - dmin_x) +1;
	w_width = (dmax_y - dmin_y) +1;

	r2_correlation = (cross_xy) / (sqrt( cross_x2 * cross_y2));

	F_fine_detail = p_perimeter_change * (d_length)/ area;
	if (partial_on_count ~=0 )
		S_ratio = fully_on_count / partial_on_count;
	else
		S_ratio = 0.0;
	end
end
%/**************************************************************************/


%/**************************************************************************/
function result = parabola_fit_center_is_in(image_buffer, n_slices) 


    result = 1;

	return;

end






				



