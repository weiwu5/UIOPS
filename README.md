# UIOPS: University of Illinois Optical Array Probe (OAP) Processing Software

This is the software developed in University of Illinois at Urbana-Champaign to process the Optical Array Probe Datasets.
The documentation below contains detailed instructions to run the software, and also information on optional parameters to run.

## Current Developers & Contributers
- **Joe Finlon** (finlon2@illinois.edu)
- **Dan Stechman** (stechma2@illinois.edu)
- **Adam Majewski** (amajewsk@uwyo.edu)

## What's New
**February 26, 2018**
- **sizeDist:** Added support for adjusting the TAS when NaN values are encountered  

**February 09, 2018**
- **read\_binary\_PMS:** Changed netCDF data types to save space and improve write times  
- **read\_binary\_PMS:** Fixed rare instance where data was not saved to file  
- **read\_binary\_DMT:** Changed netCDF data types to save space and improve write times  
- **imgProc\_sm:** Minor corrections to Fast2DC implementation from previous update  
- **imgProc\_sm:** Addresses handling of 2DS clock time when the TAS is a NaN value  
- **IntArrAnalysis:** Added support for the Fast2DC probe  
- **IntArrAnalysis:** Improved filename conventions for 2DS H and V channels  
- **IntArrAnalysis:** Improved handling of time for 2D histogram plots  
- **sizeDist:** Added support for the max TAS values w/ Fast2DC  
- **sizeDist:** Added support for SOCRATES defaults; improved handling of inter-arrival times for Fast2DC  
- **sizeDist:** Addresses rare issue w/ particle reacceptance at very end of flight  
- **sizeDist:** Addresses issue where a flight data file spans multiple days but sampled particles do not  

**January 05, 2018**
- **read\_binary\_PMS:** Improved image filename creation if a directory is specified in the ‘outfilename’ input argument  

**December 28, 2017**
- **imgProc\_sm:** Added support for Fast2DC probe (64 diode array)  
- **imgProc\_sm:** Moved iRectEllipse and iCalcAllDiodeStats to input parameters  
- **imgProc\_sm:** Added software preamble printed statements  

**November 17, 2017**
- **sizeDist:** Implemented user argument for handling how inter-arrival time thresholds are implemented  
- **sizeDist:** Moved toggles for processing aspect ratio info and saving of info for rejected particles, inter-arrival times, and sample volume up to the user input arguments  
- **hhmmss2insec** and **insec2hhmmss:** Changed ingested arguments to double precision to fix rare issue where "Error on sizing" message pops up when running _sizeDist_  

**November 9, 2017**
- **imgProc\_sm:** Fixed PECAN-specific corrupt record identification that prevented code to run for other projects  
- **imgProc\_sm:** Romoved unnecessary 'else' statement immediately following the data integrity statement block  
- **imgProc\_sm:** Changed netCDF type for overload variable to work with the 2DC/2DP  

**November 6, 2017**
- **imgProc\_sm:** Improved handling of image buffers when zero-image particles are detected  
- **imgProc\_sm:** Fixed millisec conversion and handling of millisec/microsec values > 1000  
- **imgProc\_sm:** Buffer overload times now saved for PMS platform data for subsequent sample volume correction in sizeDist  
- **sizeDist:** Improved sample volume treatment for 2DC/2DP using the buffer overload time  

**August 8, 2017**
- **imgProc\_sm:** Improved variable types when saving data to netCDF to allow up to 40% smaller file sizes!  
- **imgProc\_sm:** Improvements to handling corrupted records with CIP/PIP  
- **sizeDist:** Cleaned up some code and made minor code formatting changes  
- Added a few example helper scripts to demonstrate how base UIOPS scripts are run (will be moved to separate 'Examples' folder in future)  

## Dependencies
**Supported OSes:** Linux, Mac (MATLAB Windows support coming in the future)  
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
- **runSizeDistPECAN.m:** Example function to do size distribution
- **single_area.m:** Calculates the area of a paricle using the A-D relationship from Mitchell (1996) (used in calculation of particle fallspeed)
- **single_mass.m:** Calculates the mass of a single particle
- **single_vt.m:** Calculates the particle fallspeed of a single particle
- **sizeDist.m:** Generates size distributions
- **minimal\_working\_example/:** -- DIRECTORY CONTAINING HELPER SCRIPTS FOR A SAMPLE DATASET --
   **../run\_imgProc\_Munich.m:** Helper function to run imgProc\_sm on included sample dataset
   **../run\_intArrAnalysis\_Munich.m:** Helper function to run intArrAnalysis on included sample dataset
   **../run\_sizeDist\_Munich.m:** Helper function to run sizeDist on included sample dataset
   **../run\_sdPlots\_Munich.m:** Plots example distribution results on included sample dataset

## Minimum Working Example
The UIOPS repository contains a sample dataset for you to test UIOPS on your system. Included are sample image files and helper scripts that allow the individual scripts to be executed. To save file space, the first step in processing (_read\_binary\_\*_) is skipped. Output files and plots are saved to a 'files' directory within the downloaded repository. To run this example execute the following scripts below within the _minimum\_working\_example/_ directory.

1. run\_imgProc\_Munich
2. run\_intArrAnalysis\_Munich (only analyzes 2DS inter-arrival times)
3. run\_sizeDist\_Munich
4. run\_sdPlots\_Munich

## Reading Binary Files
UIOPS contains separate scripts to decompress the binary file for different probe types. Each variant contains the same inputs and outputs, and are included below. The image buffer, denoted by the _data_ variable in the output (image) file, is saved as a double representation as follows: # image blocks (4 for PMS probes; 8 for other probes) by the # of slices (512 for PMS probes; 1700 for other probes). Each image block can later be decoded using an 8-bit representation (16-bit for SPEC probes) to represent each diode in the slice as a binary string (e.g., 16 x 8 = 128 diodes for SPEC probes).

**read\_binary\_SPEC(** inFilename [string], outFilename [string] **)**  
**read\_binary\_DMT(** inFilename [string], outFilename [string] **)**  
**read\_binary\_PMS(** inFilename [string], outFilename [string] **)**   
**read\_binary\_SEA(** inFilename [string], outFilename [string] **)** [NOTE: Ensure that the tags contained in your .sea file match the ones in this script]

## Processing Individual Particles
The imgProc\_sm script parses an image record to compute various size, morphological, and diagnistic properites for individual particles.

**imgProc_sm(** inFilename [string], outFilename [string], probeName [string], nthDataChunk [num], frameChunkSize [num], projectName [string], iRectEllipse [num], iCalcAllDiodeStats [num], flightFilename [optional string] **)**  
   **inFilename:** netCDF file generated from _read\_binary\_\*_ script   
   **probeName:** "2DS", "HVPS", "CIP", "PIP", "2DC", "2DP", "Fast2DC"  
   **nthDataChunk:** If not running in parallel set this value to 1; if running in parallel set up parfor loop with this variable changing within the parallel loop  
   **frameChunkSize:** Maximum # of image records to process for each output file created  
   **projectName:** Insert an empty string if you haven't added modified probe presets specific to a field project  
   **iRectEllipse:** 0 - Do not process rectangle/ellipse fit dimensions; 1 - Process this info  
   **iCalcAllDiodeStats:** 0 - Do not save diode stats for every particle (large files); 1 - Save  
   **flightFilename:** Only needed for SPEC probes; requires netCDF file generated from _formatFlightData_ script
### Particle Techniques Applied
- Dmax following a minimum enclosing circle, rectangle, and ellipse follows the CGAL library (https://www.cgal.org/) using C++ code compiled within MATLAB wrappers
- Korolev (2007) correction applied to hollow, spherical particles (not to be used with ice crystals)
### Particle Output Variables
Variables outputted from sizeDist are saved with a combination of dimensions: time (variable based on # of particles sampled), position (length of 2), and diode count (# probe photodiodes). Global attributes of the software used, institution name, time of file creation, project dataset, image file source path, probe type, and active optional parameters are outlined in the file's metadata.
- **Date:** Date of image record in which the particle resides [dimensions: time; units: YYYYMMDD]
- **Time:** Time (UTC) of image record in which the particle resides [dimensions: time; units: HHMMSS]
- **msec:** Sub-second time of image record in which the particle resides [dimensions: time; units: ms]
- **Time_in_seconds:** Time since probe was activated [dimensions: time; units: s]
- **SliceCount:** Number of slices containing particle [SPEC/DMT only; dimensions: time]
- **DMT_DOF_SPEC_OVERLOAD:** Flag denoting out-of-focus (DMT) or overloaded (SPEC) particles. Overload time for each particle (PMS) in encompassing record. [dimensions: time; units (PMS): ms]
- **Particle_number_all:** Index of particle in 2-D buffer [SPEC/DMT only; dimensions: time]
- **position:** Slice within image record where particle is first/last sampled [dimensions: position x time]
- **particle_time:** Particle time (not available for 2DS/HVPS) [dimensions: time; units: HHMMSS]
- **particle_millisec:** Sub-second particle time (not available for 2DS/HVPS) [dimensions: time; units: ms]
- **particle_microsec:** Sub-second particle time [dimensions: time; units: microsec (DMT/PMS), unitless clock count (SPEC)]
- **parent_rec_num:** Index of image record in which the particle resides [dimensions: time]
- **particle_num:** Particle index within current image record [dimensions: time]
- **image_length:** Particle length in time direction using # photodiodes [dimensions: time]
- **image_width:** Particle length in photodiode direction using # photodiodes [dimensions: time]
- **image_area:** Projected area using the # shadowed photodiodes [dimensions: time]
- **image_longest_y:** Longest vertically-oriented chord through particle in time direction [dimensions: time]
- **image_max_top_edge_touching:** Maximum # of times the top diode is shadowed in succession [dimensions: time]
- **image_max_bottom_edge_touching:** Maximum # of times the bottom diode is shadowed in succession [dimensions: time]
- **image_touching_edge:** 0 denotes image entirely in array; 1 denotes image touching edge [dimensions: time]
- **image_auto_reject:** ASCII reject code ("0" or 48: accepted; "a" or 97: aspect ratio > 6; "t" or 116: aspect ratio > 5 + image touching edge; "p" or 112: < 25% shadowed diodes in rectangle encompassing particle; "h" or 104 or 72 or 117: hollow particle; "s" or 115: split image; "z" or 122: zero area image; "f" or 102: zero area image) [dimensions: time]
- **image_hollow:** 0 denotes not hollow; 1 denotes a hollow image [dimensions: time]
- **image_center_in:** 0 denotes center of particle outside array; 1 denotes center is inside [dimensions: time]
- **image_axis_ratio:** Ratio between maximum vertical length and maximum horizontal length [dimensions: time]
- **image_diam_circle_fit:** Dmax following Heymsfield & Parrish (1978) [dimensions: time; units: mm]
- **image_diam_horiz_chord:** Dmax from # slices+1; best for sideways-looking probes [dimensions: time; units: mm]
- **image_diam_horiz_chord_corr:** Dmax from # slices+1, with Korolev (2007) correction applied to hollow spherical particles [dimensions: time; units: mm]
- **image_diam_following_bamex_code:** Dmax from maximum length chord through the particle [dimensions: time; units: mm]
- **image_diam_vert_chord:** Dmax from maximum length in photodiode direction; best for sideways-looking probes [dimensions: time; units: mm]
- **image_diam_minR:** Dmax of smallest circle enclosing the particle [dimensions: time; units: mm]
- **image_diam_AreaR:** Area equivalent diameter [dimensions: time; units: mm]
- **image_perimeter:** Perimeter following the particle boundary [dimensions: time; units: mm]
- **percent_shadow_area:** Ratio between the projected area and L\*W [dimensions: time; units: %]
- **edge_at_max_hole:** # diodes between edges of the particle for the slice containing the largest gap inside the particle [dimensions: time]
- **max_hole_diameter:** Diameter of the largest hole inside the particle [dimensions: time]
- **part_z:** Particle depth in object plane calculated via lookup table [dimensions: time]
- **size_factor:** Dmax reduction factor folowing Korolev (2007) correction [dimensions: time]
- **holroyd_habit:** ASCII habit code following the Holroyd (1987) algorithm ("M" or 77: zero image; "C" or 67: center-out image; "t" or 116: tiny; "o" or 111: oriented; "l" or 108: linear; "a" or 97: aggregate; "g" or 103: graupel; "s" or 115: sphere; "h" or 104: hexagonal; "i" or 105: irregular; "d" or 100: dendrite) [dimensions: time]
- **area_hole_ratio:** Ratio between the projected area and the area of hole inside particle [dimensions: time]
- **inter_arrival:** Inter-arrival time between particles (use diff(time_in_seconds) for all except 2DC/2DP)[dimensions: time; units: s]
- **bin_stats:** # times specified photodiode is shadowed for particles in this file [dimensions: diode count]
- **image_bin_stats:** # times specified photodiode is shadowed for each particle [dimensions: diode count x time]

## Analyzing Inter-arrival Times
The _intArrAnalysis_ script fits the inter-arrival times of a population of particles to a bimodal distribution following Field et al. (2006). Thresholds are applied based on (1) the minimum frequency in the (1) fit and (2) histogram between the peaks of both modes.

**intArrAnalysis(** inFilename [string], fileDirectory [string], probeName [string], runAnalysis [num], date [string], numChunks [num], nthDataChunk [optional num] **)**  
   **inFilename:** File path of PBP data file  
   **fileDirectory:** Directory path to save data and plots  
   **probeName:** "2DS", "HVPS", "CIP", "PIP", "2DC", "2DP"  
   **runAnalysis:** 1 runs bimodal fit analysis; 0 only plots histograms of inter-arrival times  
   **date:** Start date of flight [YYYYMMDD]  
   **numChunks:** If not running in parallel set this value to 1; if running in parallel assign the number of parallel jobs to run (maximum of 7)  
   **nthDataChunk:** If not running in parallel then skip this input; if running in parallel set up a parfor loop with this variable changing within the parallel loop
### Inter-arrival Techniques Applied
- tau_1 and tau_2 are determined following Field et al. (2006)
- The bimodal fit relies on initial (beta) coefficients, with up to 1000 iterations allowed for the distribution parameters to converge on a value
- **NOTE:** If the distribution does not represent a bimodal distribution well, then the fit will be ill-constrained and the determined threshold will need to be manually adjusted.
### Threshold Values from Output .mat File
- **threshold:** 2 \* Inter-arrival value from peak in smaller inter-arrival mode
- **threshold_ak:** Inter-arrival value from minimum histogram frequency between both mode peaks
- **threshold_ww:** Inter-arrival value from minimum fit frequency between both mode peaks
### Threshold Values from Output .cdf File
**NOTE:** This file contains the time and inter-arrival threshold to be used when generating the size distribuions. If threshold values need to be occasionally modified, or if processing was done in chunks and needs to be combined into one file, make sure the final netCDF contains the same variables listed below.
- **Time:** Particle time (length must equal the # of particles) [type: double; units: HHMMSS]
- **threshold:** Inter-arrival time threshold (length must equal the # of particles) [type: double; units: s]

## Creating Size Distributions
The _sizeDist_ script generates size distributions of number _N_, mass _M_, area _A_, area ratio, aspect ratio, etc. as well as multiple bulk variables for a given probe. The size bins in which to partition particles can be set to the probe's original diode resolution or a defined array with varying bin widths at the beginning of the script. This and other presets can be fine-tuned in the beginning by adding a project name to the project case clause.

**sizeDist(** inFilename [string], outFilename [string], tas [array], flightTime [array], probeName [string], sizeMethod [num], saMethod [num], pressure [array], temperature [array], iaThreshType [num], createAspectRatio [num], saveBadParticles [num], saveIAandSV [num], projectName [string], date [string], iaThreshFilename [optional string] **)**
   **inFilename:** File path of PBP data file  
   **outFilename:** Filename to save size distribution data  
   **tas:** True air speed loaded as an array (e.g., within a helper script) [m/s]  
   **flightTime:** Aircraft time loaded as an array (e.g., within a helper script) [HHMMSS]  
   **probeName:** "2DS", "HVPS", "CIP", "PIP", "2DC", "2DP"  
   **sizeMethod:** 1 (Lx), 2 (Ly), 3 (mean(Lx,Ly)), 4 (hypotenuse(Lx,Ly)), 5 (max(Lx,Ly)), 6 (D of minimum enclosing circle)  
   **saMethod:** 0 (center-in), 1 (entire-in), 2 (Heymsfield & Parrish (1978) correction)  
   **pressure:** Aircraft pressure loaded as an array (e.g., within a helper script) [hPa]  
   **temperature:** Aircraft temperature loaded as an array (e.g., within a helper script) [deg C]  
   **iaThreshType:** 0: Campaign/probe default; 1: Time-dependent; 2: Spiral-dependent  
   **createAspectRatio:** 0: Do not process aspect ratio info; 1: Process this info  
   **saveBadParticles:** 0: Do not save info on rejected particles, inc. PSDs; 1: Save info  
   **saveIAandSV:** 0: Do not save info on inter-arrival time and sample volume, inc. PSDs; 1: Save info  
   **projectName:** If project-specific presets are not used then supply an empty string  
   **date:** Start date of flight [HHMMSS]  
   **iaThreshFilename:** Filename for inter-arrival time threshold data (length must equal the # of particles; only needed if you're using a time-varying threshold)
### Distribution Techniques Applied
- Option to reject particles using a constant inter-arrival time threshold, one that is constant for provided time ranges, or one that is assigned to subsets of particles
- Mass distribution functions M(D) are computed two ways: (1) using the m-A relationship following Baker & Lawson (2006), and (2) using several habit-specific m-D relationships following McFarquhar et al. (2002) and subsequently applied in Jackson et al. (2012)
- Area distribution function A(D) is computed using several habit-specific A-D relationships for particles entirely in the photodiode array
- Terminal velocity distributions are computed one of two ways: (1- default) using the modified Best Number following Heymsfield and Westbrook (2010), or (2) using Mitchell (1996), a special case from Heymsfield and Westbrook (2010) with different parameter coefficients
- Mean aspect ratio distributions (elliptical and rectangular fits) and mean area ratio distributions are computed from particles entirely in the photodiode array
- Effective radius computed following Fu (1996)
- Particle precipitation rate computed using particle mass from habit-specific m-D relations and from the particle's terminal velocity (see methodology above)
### Distribution Output Variables
Variables outputted from sizeDist are saved with a combination of dimensions: number of size bins (varies based on configuration), number of area ratio bins (10; ranging between 0.1 and 1 by increments of 0.1), time (variable based on flight duration), and habit (10; ordered as: spherical, linear, oriented, tiny, hexagonal, irregular, graupel, dendrite, aggregate). Global attributes of the software used, institution name, time of file creation, project dataset, PBP file source path, probe type, sample area and maximum dimension methods used, and active optional parameters are outlined in the file's metadata.
- **time:** Aircraft time in UTC [dimensions: time; units: HHMMSS]
- **bin_min:** Bin minimum size [dimensions: size; units: mm]
- **bin_max:** Bin maximum size [dimensions: size; units: mm]
- **bin_mid:** Bin midpoint size [dimensions: size; units: mm]
- **bin_dD:** Bin size [dimensions: size; units: mm]
- **conc_minR:** Size distribution using the maximum dimension [dimensions: size x time; units: cm^-4]
- **area:** Size distribution partitioned by area ratio (ar=0.1:0.1:1.0) [dimensions: area ratio x size x time; units: cm^-4]
- **conc_AreaR:** Size distribution using the area-equivalent diameter [dimensions: size x time; units: cm^-4]
- **n:** Number concentration [dimensions: time; units: cm^-3]
- **total_area:** Binned projected area (extinction) [dimensions: size x time; units: mm^2/cm^-4]
- **mass:** Binned mass using habit-specific m-D relations [dimensions: size x time; units: g/cm^-4]
- **habitsd:** Size distribution partitioned by habit [dimensions: habit x size x time; units: cm^-4]
- **re:** Effective radius [dimensions: time; units: mm]
- **ar:** Mean area ratio following a minimum enclosing circle [dimensions: time]
- **massBL:** Binned mass using an m-A relation [dimensions: size x time; units: g/cm^-4]
- **Reject_ratio:** Percent of particles rejected for each time interval [dimensions: time; units: %]
- **vt:** Binned mass-weighted terminal velocity [dimensions: size x time; units: g/cm^-4]
- **Prec_rate:** Binned precipitation rate [dimensions: size x time; units: mm/hr]
- **habitmsd:** Binned mass partitioned by habit [dimensions: habit x size x time; units: g/cm^-4]
- **Calcd_area:** Binned mean particle area using habit-specific A-D relations [dimensions: size x time; units: mm^2/cm^-4]
- **count:** Binned particle count for all accepted particles [dimensions: size x time]
- **mean_aspect_ratio_rectangle:** Binned mean aspect ratio following a rectangular fit [dimensions: size x time]
- **mean_aspect_ratio_ellipse:** Binned mean aspect ratio following an elliptical fit [dimensions: size x time]
- **mean_area_ratio:** Binned mean area ratio following a minimum enclosing circle [dimensions: size x time]
- **mean_perimeter:** Binned mean perimeter [dimensions: size x time]
- _Variables with "REJ" appended retain the above definitions for particles that have been rejected_
- **sum_IntArr:** Sum of inter-arrival times, excluding the overload time for particles affected by saving of image data [dimensions: time; units: s]
- **sample_vol:** Binned sample volume [dimensions: size x time; units: cm^3]
