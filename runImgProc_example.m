function runImgProc_sm_PECAN(nSlice, ddate, probe)

inOutDir = '/data/pecan/a/stechma2/pecan/mp-data/IProcessingRelease/';


%  nChuck*nEvery shoudl equal the total frame number 
nChucks=12; % Number of chucks (sepeate files)
nEvery =floor(nSlice/12)+1; % Size of every chucks. 
numb=11:10+nChucks;  % Start from 11 to avoid sigle numbers in file name for later convinience

% Assign the number of CPUs for this program
gcp=parpool(12)        % Assign n CPUs to process

% Choose the start and end of chucks to be processed. Remember you can
% split the chucks into different programs to process, since matlabpool can
% only use 8 CPUs at once
parfor iii=1:nChucks  % iiith chuck will be processed 
	infile  = [inOutDir 'DIMG.' probe '.' num2str(ddate) '.cdf'];  % Input file
    outfile = [inOutDir 'proc2.' num2str(ddate) '.' num2str(numb(iii)) '.' probe '.cdf'];	% Output image autoanalysis file
    imgProc_sm(infile,outfile, probe, iii, nEvery, 'PECAN');  % See imgProc_sm documentation for more information
end

delete(gcp)
%matlabpool close
exit

end

