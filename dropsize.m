function [center_in,axis_ratio,diameter_circle_fit,diameter_horiz_chord,diameter_vert_chord,diameter_horiz_mean, diameter_spheroid]=...
    dropsize(max_horizontal_length,max_vertical_length,image_area,largest_edge_touching,...
    smallest_edge_touching,diode_size,corrected_horizontal_diode_size,number_diodes_in_array)

global max_vertical_chords max_horizontal_chords vertical_chords vertical_chord_equivalent_spherical_diameters...
    horizontal_chords horizontal_chord_equivalent_spherical_diameters
% 
% 
% 	DROP SIZING FOR 2DP PROBES (equiv sph vol diam)
% 		one diode added to max_horizontal_length_in_pixels for missing first slice
% 		uncalculated size appears as zero; oversize as 9.1
% 
% 
% 	diameter_horz_chord
% 		d0 from horz chord (max_horizontal_length_in_pixels + 1)
% 		designed for sideways-looking probe                
% 		but can be used for any probe orientation 
% 		with center-in image of equil shape
% 
% 
% 	diameter_circle_fit
% 		Heimsfield-Parish CIRCLE FIT SIZES FOR 2-EDGE & 1-EDGE(CENTER OUT) IMAGES:
% 		d0 from horz chord (circle fit)
% 			designed for downward-looking probe
% 
% 
% 	diameter_vert_chord
% 		max_vertical_length_in_pixels & 2-CHORD SIZES FOR ENTIRE-IN IMAGES:
% 		d0 from vert chord (max_vertical_length_in_pixels)
% 			designed for sideways-looking probe
% 			optional size for entire-in images
% 
% 
% 	[If diam_vchord differs from diam_hchord, then drop is not equil shape]
% 	
% 		diameter_horz_mean
% 			d0 from mean horz chord [(xmax+1)*ymax]^0.5
% 			designed for downward-looking probe
% 			optional size for distorted entire-in images
% 
% 		diameter_spheroid
% 			d0 from spheroid assumption (hchord*hchord*vchord)^1/3
% 			designed for sideways-looking probe
% 			optional size for distorted entire-in images
% 
% 
% 	AXIS RATIO
% 		axis_ratio = max_vertical_length_in_pixels / max_horizontal_length_in_pixels
% 		for entire-in image (no smoothing)
%       center_in = (0 if particle is not center in, 1 if particle is center
%   in)

max_horizontal_chords = 117;
max_vertical_chords = 55;

diodes_added_to_length = 0;
diodes_added_to_height = 1.0;

horizontal_chords = [...
	 0.0000,  0.1000,  0.2000,  0.3000,  0.4000,  0.5000,  0.6000,  0.7000,  0.8000,  0.9000,...
	 1.0000,  1.1000,  1.2000,  1.3000,  1.4000,  1.5000,  1.6000,  1.7000,  1.8000,  1.9000,...
	 2.0000,  2.1000,  2.2000,  2.3000,  2.4000,  2.5000,  2.6000,  2.7000,  2.8000,  2.9000,...
	 3.0000,  3.1000,  3.2000,  3.3000,  3.4000,  3.5000,  3.6000,  3.7000,  3.8000,  3.9000,...
	 4.0000,  4.1000,  4.2000,  4.3000,  4.4000,  4.5000,  4.6000,  4.7000,  4.8000,  4.9000,...
	 5.0000,  5.1000,  5.2000,  5.3000,  5.4000,  5.5000,  5.6000,  5.7000,  5.8000,  5.9000,...
	 6.0000,  6.1000,  6.2000,  6.3000,  6.4000,  6.5000,  6.6000,  6.7000,  6.8000,  6.9000,...
	 7.0000,  7.1000,  7.2000,  7.3000,  7.4000,  7.5000,  7.6000,  7.7000,  7.8000,  7.9000,...
	 8.0000,  8.1000,  8.2000,  8.3000,  8.4000,  8.5000,  8.6000,  8.7000,  8.8000,  8.9000,...
	 9.0000,  9.1000,  9.2000,  9.3000,  9.4000,  9.5000,  9.6000,  9.7000,  9.8000,  9.9000,...
	10.0000, 10.1000, 10.2000, 10.3000, 10.4000, 10.5000, 10.6000, 10.7000, 10.8000, 10.9000,...
	11.0000, 11.1000, 11.2000, 11.3000, 11.4000, 11.5000, 11.6000];

horizontal_chord_equivalent_spherical_diameters = [...
	0.0000, 0.1000, 0.2000, 0.3000, 0.3998, 0.4996, 0.5992, 0.6986, 0.7976, 0.8964,...
	0.9947, 1.0927, 1.1903, 1.2875, 1.3842, 1.4804, 1.5761, 1.6711, 1.7657, 1.8597,...
	1.9531, 2.0460, 2.1385, 2.2304, 2.3217, 2.4124, 2.5026, 2.5921, 2.6812, 2.7696,...
	2.8574, 2.9448, 3.0319, 3.1185, 3.2046, 3.2903, 3.3755, 3.4603, 3.5446, 3.6286,...
	3.7120, 3.7950, 3.8776, 3.9598, 4.0416, 4.1229, 4.2036, 4.2841, 4.3642, 4.4440,...
	4.5233, 4.6024, 4.6811, 4.7597, 4.8378, 4.9152, 4.9924, 5.0690, 5.1452, 5.2210,...
	5.2961, 5.3711, 5.4457, 5.5201, 5.5942, 5.6681, 5.7416, 5.8148, 5.8877, 5.9602,...
	6.0323, 6.1041, 6.1756, 6.2467, 6.3175, 6.3879, 6.4580, 6.5278, 6.5973, 6.6664,...
	6.7353, 6.8040, 6.8720, 6.9399, 7.0075, 7.0744, 7.1411, 7.2074, 7.2733, 7.3388,...
	7.4040, 7.4688, 7.5332, 7.5973, 7.6610, 7.7241, 7.7864, 7.8490, 7.9117, 7.9739,...
	8.0359, 8.0976, 8.1591, 8.2203, 8.2813, 8.3421, 8.4027, 8.4631, 8.5233, 8.5834,...
	8.6433, 8.7030, 8.7625, 8.8220, 8.8814, 8.9407, 8.9998];

vertical_chords = [...
	0.0000,  0.1000,  0.2000,  0.3000,  0.4000,  0.5000,  0.6000,  0.7000,  0.8000,  0.9000,...
	1.0000,  1.1000,  1.2000,  1.3000,  1.4000,  1.5000,  1.6000,  1.7000,  1.8000,  1.9000,...
	2.0000,  2.1000,  2.2000,  2.3000,  2.4000,  2.5000,  2.6000,  2.7000,  2.8000,  2.9000,...
	3.0000,  3.1000,  3.2000,  3.3000,  3.4000,  3.5000,  3.6000,  3.7000,  3.8000,  3.9000,...
	4.0000,  4.1000,  4.2000,  4.3000,  4.4000,  4.5000,  4.6000,  4.7000,  4.8000,  4.9000,...
	5.0000,  5.1000,  5.2000,  5.3000,  5.4000];
vertical_chord_equivalent_spherical_diameters = [...
	0.0000,  0.1000,  0.2000,  0.3001,  0.4003,  0.5008,  0.6016,  0.7028,  0.8048,  0.9075,...
	1.0110,  1.1155,  1.2208,  1.3271,  1.4348,  1.5438,  1.6545,  1.7668,  1.8809,  1.9967,...
	2.1142,  2.2336,  2.3553,  2.4792,  2.6059,  2.7353,  2.8677,  3.0023,  3.1394,  3.2792,...
	3.4219,  3.5678,  3.7169,  3.8695,  4.0259,  4.1865,  4.3512,  4.5201,  4.6922,  4.8687,...
	5.0525,  5.2438,  5.4432,  5.6474,  5.8580,  6.0775,  6.3073,  6.5485,  6.8018,  7.0729,...
	7.3699,  7.7007,  8.0781,  8.4868,  8.9230];


diameter_circle_fit = 0.0;
diameter_horiz_chord = 0.0;
diameter_vert_chord = 0.0;
diameter_horiz_mean = 0.0;
diameter_spheroid = 0.0;
axis_ratio = 0.0;
center_in = 0;

scaling_factor_for_horizontal_lengths = corrected_horizontal_diode_size / diode_size;


corrected_diodes_added_to_length = 0;

if(image_area < 1)
    diameter_circle_fit = corrected_diodes_added_to_length;
    diameter_horiz_chord = corrected_diodes_added_to_length;
    diameter_vert_chord = corrected_diodes_added_to_length;
    diameter_horiz_mean = corrected_diodes_added_to_length;
    diameter_spheroid = corrected_diodes_added_to_length;
    center_in=0;
    axis_ratio = 1.0;
    return;
end
    

    
largest_edge_touching_length =  largest_edge_touching * scaling_factor_for_horizontal_lengths;
smallest_edge_touching_length = smallest_edge_touching * scaling_factor_for_horizontal_lengths;
max_horizontal_length = max_horizontal_length * scaling_factor_for_horizontal_lengths;
corrected_diodes_added_to_length = diodes_added_to_length * scaling_factor_for_horizontal_lengths;



% determine no. of edges
number_edges_touching = 0;
if largest_edge_touching > 0
    number_edges_touching = number_edges_touching + 1;
    if smallest_edge_touching > 0
        number_edges_touching = number_edges_touching + 1;
    end
elseif smallest_edge_touching > 0
    number_edges_touching = number_edges_touching + 1;
end


center_in = 1;

if max_horizontal_length <= largest_edge_touching
    center_in = 0;
end

if number_edges_touching == 2 & center_in == 1
    temp = number_diodes_in_array + (largest_edge_touching_length^2 - smallest_edge_touching_length^2 ) / ( 4 * number_diodes_in_array); % + is replaced by -, Will 10/17/2013
    horizontal_size = sqrt(temp^2 + smallest_edge_touching_length^2);
    horizontal_chord = horizontal_size * diode_size;
    diameter_circle_fit = horizontal_chord_to_spherical_dia(horizontal_chord);
    diameter_horiz_chord = horizontal_chord_to_spherical_dia(horizontal_chord);
elseif number_edges_touching == 2
    largest_edge_touching_length = largest_edge_touching_length + corrected_diodes_added_to_length;
    temp = number_diodes_in_array + (largest_edge_touching_length^2 - smallest_edge_touching_length^2 ) / ( 4 * number_diodes_in_array);  % + is replaced by -, Will 10/17/2013
    horizontal_size = sqrt(temp^2 + smallest_edge_touching_length^2);
    horizontal_chord = horizontal_size * diode_size;
    diameter_circle_fit = horizontal_chord_to_spherical_dia(horizontal_chord);
elseif number_edges_touching == 1 & center_in == 1
    horizontal_chord = (max_horizontal_length + diodes_added_to_length) * corrected_horizontal_diode_size;
    diameter_horiz_chord = horizontal_chord_to_spherical_dia(horizontal_chord);
    diameter_circle_fit = horizontal_chord_to_spherical_dia(horizontal_chord);
elseif number_edges_touching == 1
    largest_edge_touching_length = largest_edge_touching_length + corrected_diodes_added_to_length;
    max_vertical_length = max_vertical_length + diodes_added_to_height * 0.5;
    horizontal_size = (0.25 * largest_edge_touching_length^2 + max_vertical_length^2)/(max_vertical_length);
    horizontal_chord = horizontal_size * diode_size;
    diameter_circle_fit = horizontal_chord_to_spherical_dia(horizontal_chord);
else
    horizontal_chord = (max_horizontal_length + diodes_added_to_length) * corrected_horizontal_diode_size;
    diameter_horiz_chord = horizontal_chord_to_spherical_dia(horizontal_chord);
    diameter_circle_fit = horizontal_chord_to_spherical_dia(horizontal_chord);
    vertical_chord = (max_vertical_length + diodes_added_to_length) * diode_size;
    axis_ratio = vertical_chord / horizontal_chord;
    diameter_vert_chord = vertical_chord_to_spherical_dia(vertical_chord);
    horizontal_mean_chord = sqrt(horizontal_chord * vertical_chord);
    diameter_horiz_mean = horizontal_chord_to_spherical_dia(horizontal_mean_chord);
    diameter_spheroid = exp(log(horizontal_chord^2 * vertical_chord)/3);
end
end

    
    

function diameter=vertical_chord_to_spherical_dia(vertical_chord)
global max_vertical_chords max_horizontal_chords vertical_chords vertical_chord_equivalent_spherical_diameters...
    horizontal_chords horizontal_chord_equivalent_spherical_diameters

    delta_vertical_chord = .1;
    
    i = round(vertical_chord * 10);
    if i+1 < max_vertical_chords & i ~=0
        delta_diameter = vertical_chord_equivalent_spherical_diameters(i+1) - vertical_chord_equivalent_spherical_diameters(i);
        diameter = vertical_chord_equivalent_spherical_diameters(i) + (delta_diameter / delta_vertical_chord) * (vertical_chord - vertical_chords(i));
    elseif i == 0
        diameter = 0;
    else
        diameter = 9.1;
    end
    
   
end


function diameter = horizontal_chord_to_spherical_dia(horizontal_chord)
global max_vertical_chords max_horizontal_chords vertical_chords vertical_chord_equivalent_spherical_diameters...
    horizontal_chords horizontal_chord_equivalent_spherical_diameters


    delta_horizontal_chord = .1;

    i = round(horizontal_chord * 10);
    if i+1 < max_horizontal_chords & i ~= 0
        delta_diameter = horizontal_chord_equivalent_spherical_diameters(i+1) - horizontal_chord_equivalent_spherical_diameters(i);
        diameter = horizontal_chord_equivalent_spherical_diameters(i) + (delta_diameter / delta_horizontal_chord) * (horizontal_chord - horizontal_chords(i));
    elseif i == 0
        diameter = 0;
    else
        diameter = 9.1;
    end

end










        

    
