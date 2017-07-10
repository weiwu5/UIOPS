# UIOPS: University of Illinois Optical Array Probe (OAP) Processing Software

This is the software developed in University of Illinois at Urbana-Champaign to process the Optical Array Probe Datasets.
The documentation below contains detailed instructions to run the software, and also information on optional parameters to run.

## Dependencies
**Supported OSes:** Linux, Mac (Windows support coming in the future)  
**MATLAB Versions:** R2015a and higher (see https://www.mathworks.com/products/matlab.html for download and more info; student version from $99 USD; check w/ university for possible university-wide license)  
**MATLAB Toolboxes:** Parallel Computing, Statistics and Machine Learning, Curve Fitting, Image Processing  
**Supported Probes:** 2DS, HVPS, CIP, PIP, 2DC, 2DP  
**Supported Platforms:** SPEC, DMT, PMS, NCAR  

## Included Programs
- **calculate\_reject\_unified.m:** Determines the reason particles are rejected (or contains separate flag if not rejected)
- **calc\_sa\_randombins.m:** Calculates probe sample area of particle as a function of Dmax
- **CGAL\_RecSize.mexa64:** Calculates minimum rectangle around particle using the CGAL library (https://www.cgal.org/)
- **CGAL\_EllipseSize.mexa64:** Calculates minimum ellipse around particle using the CGAL library
- **dropsize.m:** Generates Dmax for particle using different definitions
- **holroyd.m:** Runs automated habit identification scheme of Holroyd (1987) on particle
- **lwc\_calc.m:** Calculates liquid water content from a size distribution
- **Imageview.m:** Views the image file (no longer supported - see http://www.github.com/joefinlon/HydroViewer for replacement)
- **imgProc\_sm.m:** Generates particle-by-particle (PBP) parameters from an image netCDF file using shared-memory parallelization
- **ImgView.m:** Code for the ImgView program (no longer supported)
- **ImgView.fig:** ImgView program (no longer supported)
- **intArrAnalysis.m:** Determines the inter-arrival time threshold for populations of predefined particle count 
- **read_binary.m:** Wrapper function for the following probes/platforms
   **read\_binary\_SPEC.m:** Converts SPEC raw data to netCDF format  
   **read\_binary\_PMS.m:** Converts PMS data to netCDF format  
   **read\_binary\_SEA.m:** Converts SEA raw data to netCDF format  
   **read\_binary\_DMT.m:** Convert DMT raw data to netCDF format  
- **runImgProc.m:** Function to set variables to run imgProc\_sm.m
- **run\_imgProc\_Munich.m:** Helper function to run imgProc\_sm on included sample dataset
- **run\_intArrAnalysis\_Munich.m:** Helper function to run intArrAnalysis on included sample dataset
- **runSizeDistPECAN.m:** Example function to do size distribution
- **run\_sdPlots\_Munich.m:** Plots example distribution results on included sample dataset
- **run\_sizeDist\_Munich.m:** Helper function to run sizeDist on included sample dataset
- **single_area.m:** Calculates the area of a paricle using the A-D relationship from Mitchell (1996) (used in calculation of particle fallspeed)
- **single_mass.m:** Calculates the mass of a single particle
- **single_vt.m:** Calculates the particle fallspeed of a single particle
- **sizeDist.m:** Generates size distributions


## Minimum Working Example
The UIOPS repository contains a sample dataset for you to test UIOPS on your system. Included are sample image files and helper scripts that allow the individual scripts to be executed. To save file space, the first step in processing (read\_binary\_\*) is skipped. Output files and plots are saved to a 'files' directory within the downloaded repository. To run this example execute the following scripts below.

1. run\_imgProc\_Munich
2. run\_intArrAnalysis\_Munich (only analyzes 2DS inter-arrival times)
3. run\_sizeDist\_Munich
4. run\_sdPlots\_Munich

## Reading Binary Files
UIOPS contains separate scripts to decompress the binary file for different probe types. Each variant contains the same inputs and outputs, and are included below.

**read\_binary\_SPEC(** inFilename [string], outFilename [string] **)**  
**read\_binary\_DMT(** inFilename [string], outFilename [string] **)**  
**read\_binary\_PMS(** inFilename [string], outFilename [string] **)**   
**read\_binary\_SEA(** inFilename [string], outFilename [string] **)** [NOTE: Ensure that the tags contained in your .sea file match the ones in this script]

## Processing Individual Particles
The imgProc\_sm script parses an image record to compute various size, morphological, and diagnistic properites for individual particles.

**imgProc_sm(** inFilename [string], outFilename [string], probeName [string], nthDataChunk [num], frameChunkSize [num], projectName [string], flightFilename [optional string] **)**  
   **inFilename:** netCDF file generated from _read\_binary\_\*_ script   
   **probeName:** "2DS", "HVPS", "CIP", "PIP", "2DC", "2DP"  
   **nthDataChunk:** If not running in parallel set this value to 1; if running in parallel set up parfor loop with this variable changing within the parallel loop  
   **frameChunkSize:** Maximum # of image records to process for each output file created  
   **projectName:** Insert an empty string if you haven't added modified probe presets specific to a field project   
   **flightFilename:** Only needed for SPEC probes; requires netCDF file generated from _formatFlightData_ script

### Notable Techniques Applied
- Dmax following a minimum enclosing circle, rectangle, and ellipse follows the CGAL library (https://www.cgal.org/) using C++ code compiled within MATLAB wrappers
- Korolev (2007) correction applied to hollow, spherical particles (not to be used with ice crystals)
### Output Variables
- **Date:** Date of image record in which the particle resides [type: double; units: YYYYMMDD]
- **Time:** Time (UTC) of image record in which the particle resides [type: double; units: HHMMSS]
- **msec:** Sub-second time of image record in which the particle resides [type: double; units: ms]
- **Time_in_seconds:** Time since probe was activated [type: double; units: s]
- **SliceCount:** Number of slices containing particle [SPEC/DMT only; type: double]
- **DMT_DOF_SPEC_OVERLOAD:** Flag denoting out-of-focus (DMT) or overloaded (SPEC) particles [SPEC/DMT only; type: double]
- **Particle_number_all:** Index of particle in 2-D buffer [SPEC/DMT only; type: double]
- **position:** Slice within image record where particle is first sampled [type: double]
- **particle_time:** Particle time (not available for 2DS/HVPS) [type: double; units: HHMMSS]
- **particle_millisec:** Sub-second particle time (not available for 2DS/HVPS) [type: double; units: ms]
- **particle_microsec:** Sub-second particle time [type: double; units: microsec (DMT/PMS), unitless clock count (SPEC)]
- **parent_rec_num:** Index of image record in which the particle resides [type: double]
- **particle_num:** Particle index within current image record [type: double]
- **image_length:** Particle length in time direction using # photodiodes [type: double]
- **image_width:** Particle length in photodiode direction using # photodiodes [type: double]
- **image_area:** Projected area using the # shadowed photodiodes [type: double]
- **image_longest_y:** Longest vertically-oriented chord through particle in time direction [type: double]
- **image_max_top_edge_touching:** Maximum # of times the top diode is shadowed in succession [type: double]
- **image_max_bottom_edge_touching:** Maximum # of times the bottom diode is shadowed in succession [type: double]
- **image_touching_edge:** 0 denotes image entirely in array; 1 denotes image touching edge [type: double]
- **image_auto_reject:** ASCII reject code ("0" or 48: accepted; "a" or 97: aspect ratio > 6; "t" or 116: aspect ratio > 5 + image touching edge; "p" or 112: < 25% shadowed diodes in rectangle encompassing particle; "h" or 104 or 72 or 117: hollow particle; "s" or 115: split image; "z" or 122: zero area image; "f" or 102: zero area image) [type: double]
- **image_hollow:** 0 denotes not hollow; 1 denotes a hollow image [type: double]
- **image_center_in:** 0 denotes center of particle outside array; 1 denotes center is inside [type: double]
- **image_axis_ratio:** Ratio between maximum vertical length and maximum horizontal length [type: double]
- **image_diam_circle_fit:** Dmax following Heymsfield & Parrish (1978) [type: double; units: mm]
- **image_diam_horiz_chord:** Dmax from # slices+1; best for sideways-looking probes [type: double; units: mm]
- **image_diam_horiz_chord_corr:** Dmax from # slices+1, with Korolev (2007) correction applied to hollow spherical particles [type: double; units: mm]
- **image_diam_following_bamex_code:** Dmax from maximum length chord through the particle [type: double; units: mm]
- **image_diam_vert_chord:** Dmax from maximum length in photodiode direction; best for sideways-looking probes [type: double; units: mm]
- **image_diam_minR:** Dmax of smallest circle enclosing the particle [type: double; units: mm]
- **image_diam_AreaR:** Area equivalent diameter [type: double; units: mm]
- **image_perimeter:** Perimeter following the particle boundary [type: double; units: mm]
- **percent_shadow_area:** Ratio between the projected area and L\*W [type: double; units: %]
- **edge_at_max_hole:** # diodes between edges of the particle for the slice containing the largest gap inside the particle [type: double]
- **max_hole_diameter:** Diameter of the largest hole inside the particle [type: double]
- **part_z:** Particle depth in object plane calculated via lookup table [type: double]
- **size_factor:** Dmax reduction factor folowing Korolev (2007) correction [type: double]
- **holroyd_habit:** ASCII habit code following the Holroyd (1987) algorithm ("M" or 77: zero image; "C" or 67: center-out image; "t" or 116: tiny; "o" or 111: oriented; "l" or 108: linear; "a" or 97: aggregate; "g" or 103: graupel; "s" or 115: sphere; "h" or 104: hexagonal; "i" or 105: irregular; "d" or 100: dendrite) [type: double]
- **area_hole_ratio:** Ratio between the projected area and the area of hole inside particle [type: double]
- **inter_arrival:** Inter-arrival time between particles (use diff(time_in_seconds) for 2DS/HVPS)[type: double; units: s]
- **bin_stats:** # times specified photodiode is shadowed for particles in this file [type: double]
