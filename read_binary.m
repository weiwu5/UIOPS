function read_binary(inputfile, outputfile, platform)
%% Wrapper function for decompressing the raw files
%%
%% Parameters:
%% inputfile: the file name of input file
%% outputfile: the file name of output file, or the prefix of output file since multiple files maybe created
%% platform choice: 
%%   'DMT': For CIP and PIP probes 
%%   'SPEC': For 2DS and HVPS
%%   'PMS': For 2DC and 2DP
%%   'SEA': For all the probes in SEA platform. Need to read the code for each probe, 
%%          and the current one is based on MC3E using UND Citation
%%
%% Example Usage: read_binary('Imagefile_20150617011652','DIMG.150617.CIP.cdf','DMT') 
%%
%% Created by Wei Wu, July 1st 2016

switch platform
    case 'DMT'
        read_binary_DMT(inputfile, outputfile);
    case 'PMS'
        read_binary_PMS(inputfile, outputfile);
    case 'SPEC'
        read_binary_SPEC(inputfile, outputfile);
    case 'SEA'
        read_binary_SEA(inputfile, outputfile);
    otherwise
        disp('Wrong Platform input!')
end

end