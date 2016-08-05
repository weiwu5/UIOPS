import os, subprocess

os.chdir('/data/mcfarq/a/weiwu3/IProcessingRelease')
os.system('source ~/.profile')
#os.system('ls -alt')
ddate = "20150620"

# Decompress all binary files
#os.system('/sw/matlab-r2015a/bin/matlab -nodisplay -r "runalldecompCIP(' + ddate  + ')"')
#os.system('/sw/matlab-r2015a/bin/matlab -nodisplay -r "runalldecompPIP(' + ddate  + ')"')

#os.system('ncrcat -O DIMG.CIP.' + ddate + '* ' +'DIMG.CIP.' + ddate + '.cdf')
#os.system('ncrcat -O DIMG.PIP.' + ddate + '* ' +'DIMG.PIP.' + ddate + '.cdf')

# Image Processing Part. Produce the particle by particle information from the uncompressed images
#timeline=subprocess.Popen('ncdump -h DIMG.CIP.' + ddate + '.cdf | grep UNLIMITED', shell=True, stdout=subprocess.PIPE)
#nframes = timeline.stdout.read().split('(')[1].split()[0]
#os.system('/sw/matlab-r2015a/bin/matlab -nodisplay -r "runImgProc_dmPECAN_CIP(' + nframes + ',' + ddate + ')"')
#os.system('ncrcat proc2.' + ddate + '.*.CIP.cdf ' +'proc2.' + ddate + '.CIP.cdf')

#timeline=subprocess.Popen('ncdump -h DIMG.PIP.' + ddate + '.cdf | grep UNLIMITED', shell=True, stdout=subprocess.PIPE)
#nframes = timeline.stdout.read().split('(')[1].split()[0]
#os.system('/sw/matlab-r2015a/bin/matlab -nodisplay -r "runImgProc_dmPECAN_PIP(' + nframes + ',' + ddate + ')"')
#os.system('ncrcat proc2.' + ddate + '.*.PIP.cdf ' +'proc2.' + ddate + '.PIP.cdf')

# Produce the particle size distribution from the particle by particle information
os.system('/sw/matlab-r2015a/bin/matlab -nodisplay -r "runSizeDistPECAN(' + ddate +')"')

