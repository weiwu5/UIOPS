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
   **read\_binary\_SPEC.m:** Converts SPEC raw data to netCDF format . 
   **read\_binary\_PMS.m:** Converts PMS data to netCDF format . 
   **read\_binary\_SEA.m:** Converts SEA raw data to netCDF format . 
   **read\_binary\_DMT.m:** Convert DMT raw data to netCDF format . 
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

**imgProc_sm(** inFilename [string], outFilename [string], probeName [string], nthDataChunk [num],
             frameChunkSize [num], projectName [string], flightFilename [optional string] **)**
   **inFilename:** netCDF file generated from _read\_binary\_\*_ script . 
   **probeName:** "2DS", "HVPS", "CIP", "PIP", "2DC", "2DP" . 
   **nthDataChunk:** If not running in parallel set this value to 1; if running in parallel set up parfor loop with this variable changing within the parallel loop . 
   **frameChunkSize:** Maximum # of image records to process for each output file created . 
   **projectName:** Insert an empty string if you haven't added modified probe presets specific to a field project . 
   **flightFilename:** Only needed for SPEC probes; requires netCDF file generated from _formatFlightData_ script

### Notable Techniques Applied
- Dmax following a minimum enclosing circle, rectangle, and ellipse follows the CGAL library (https://www.cgal.org/) using C++ code compiled within MATLAB wrappers
- Korolev (2007) correction applied to hollow, spherical particles (not to be used with ice crystals)
###
