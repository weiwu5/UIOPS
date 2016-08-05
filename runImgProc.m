function runImgProc


%  nChuck*nEvery shoudl equal the total frame number 
nChucks=1; % 48;          % Number of chucks (sepeate files)
nEvery =300; %9800;       % Size of every chucks. 
numb=11:10+nChucks;  % Start from 11 to avoid sigle numbers in file name for later convinience

% Assign the number of CPUs for this program
%parpool(2)        % Assign n CPUs to process

% Choose the start and end of chucks to be processed. Remember you can
% split the chucks into different programs to process, since matlabpool can
% only use 8 CPUs at once
for iii=1:nChucks % 33:40  % iiith chuck will be processed 
    infile  = 'DIMG.example.cdf';  % Input file
    outfile = ['proc2.TEST.PIP.cdf'];			% Output image autoanalysis file
    %outfile = ['./HVPSFAST/proc2.TEST_debug.HVPS' num2str(numb(iii)) '.cdf'];			% Output image autoanalysis file
    imgProc_sm(infile,outfile, 'PIP', iii, nEvery, 'PECAN');  % See imgProc_sm documentation for more information
end

%delete(gcp)

end

