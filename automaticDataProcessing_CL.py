import os, subprocess
import sys
from datetime import datetime as dt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("ddate", help="date of the flight to be processed in YYYYMMDD")
parser.add_argument("-C", "--runCIP", help="process CIP data", action="store_true")
parser.add_argument("-P", "--runPIP", help="process PIP data", action="store_true")
parser.add_argument("-D", "--decompress", help="run raw image file decompression", action="store_true")
parser.add_argument("-I", "--imgProc", help="run particle-by-particle processing", action="store_true")
parser.add_argument("-S", "--sizeDist", help="run size distribution processing", action="store_true")

args = parser.parse_args()

if not args.runCIP and not args.runPIP:
    sys.exit('You must specify at least one probe to process (--runCIP and/or --runPIP)')
if not args.decompress and not args.imgProc and not args.sizeDist:
    sys.exit('You must specify at least one processing step (options include: --decompress, --imgProc, --sizeDist)')

ddate = args.ddate

runCIP      = args.runCIP
runPIP      = args.runPIP

decompress  = args.decompress
imgProc     = args.imgProc
sizeDist    = args.sizeDist


inOutDir = '/data/pecan/a/stechma2/pecan/mp-data/UIOPS/'

os.chdir('/data/pecan/a/stechma2/pecan/mp-data/UIOPS')

os.system('source ~/.bashrc')

print('Starting UIOPS processing for PECAN ' + ddate)

if decompress:
    if runCIP:
        print('Starting decompression of raw CIP image data... ' + dt.strftime(dt.now(),'%H:%M:%S %Y%m%d'))
        os.system('/sw/matlab-r2015a/bin/matlab -nodisplay -r "runalldecompCIP(' + ddate  + ')"')
        print('CIP Image decompression complete. Starting netCDF concatenation... ' + dt.strftime(dt.now(),'%H:%M:%S %Y%m%d'))
        os.system('ncrcat -O ' + inOutDir + 'DIMG.CIP.' + ddate + '* ' +'DIMG.CIP.' + ddate + '.cdf')
    
    if runPIP:
        print('Starting decompression of raw PIP image data... ' + dt.strftime(dt.now(),'%H:%M:%S %Y%m%d'))
        os.system('/sw/matlab-r2015a/bin/matlab -nodisplay -r "runalldecompPIP(' + ddate  + ')"')
        print('PIP Image decompression complete. Starting netCDF concatenation... ' + dt.strftime(dt.now(),'%H:%M:%S %Y%m%d'))
        os.system('ncrcat -O ' + inOutDir + 'DIMG.PIP.' + ddate + '* ' +'DIMG.PIP.' + ddate + '.cdf')

if imgProc:
    if runCIP:
        print('Starting CIP particle-by-particle image processing... ' + dt.strftime(dt.now(),'%H:%M:%S %Y%m%d'))
        timeline=subprocess.Popen('ncdump -h DIMG.CIP.' + ddate + '.cdf | grep UNLIMITED', shell=True, stdout=subprocess.PIPE)
        nframes = timeline.stdout.read().split('(')[1].split()[0]
        os.system('/sw/matlab-r2016b/bin/matlab -nodisplay -r "runImgProc_sm_PECAN(' + nframes + ',' + ddate + ',' + '\'CIP\'' + ')"')
        print('CIP particle-by-particle image processing complete. Starting netCDF concatenation... ' + dt.strftime(dt.now(),'%H:%M:%S %Y%m%d'))
        os.system('ncrcat ' + inOutDir + 'proc2.' + ddate + '.*.CIP.cdf ' +'proc2.' + ddate + '.CIP.cdf')

    if runPIP:
        print('Starting PIP particle-by-particle image processing... ' + dt.strftime(dt.now(),'%H:%M:%S %Y%m%d'))
        timeline=subprocess.Popen('ncdump -h DIMG.PIP.' + ddate + '.cdf | grep UNLIMITED', shell=True, stdout=subprocess.PIPE)
        nframes = timeline.stdout.read().split('(')[1].split()[0]
        os.system('/sw/matlab-r2016b/bin/matlab -nodisplay -r "runImgProc_sm_PECAN(' + nframes + ',' + ddate + ',' + '\'PIP\'' + ')"')
        print('PIP particle-by-particle image processing complete. Starting netCDF concatenation... ' + dt.strftime(dt.now(),'%H:%M:%S %Y%m%d'))
        os.system('ncrcat ' + inOutDir + 'proc2.' + ddate + '.*.PIP.cdf ' +'proc2.' + ddate + '.PIP.cdf')

if sizeDist:
    if runCIP and runPIP:
        print('Starting CIP and PIP size distribution processing... ' + dt.strftime(dt.now(),'%H:%M:%S %Y%m%d'))
        os.system('/sw/matlab-r2016b/bin/matlab -nodisplay -r "runSizeDistPECAN(' + ddate + ',' + repr(1) + ',' + repr(1) +')"')
        print('CIP and PIP size distribution processing complete... ' + dt.strftime(dt.now(),'%H:%M:%S %Y%m%d'))
    elif runCIP and not runPIP:
        print('Starting CIP size distribution processing... ' + dt.strftime(dt.now(),'%H:%M:%S %Y%m%d'))
        os.system('/sw/matlab-r2016b/bin/matlab -nodisplay -r "runSizeDistPECAN(' + ddate + ',' + repr(1) + ',' + repr(0) +')"')
        print('CIP size distribution processing complete... ' + dt.strftime(dt.now(),'%H:%M:%S %Y%m%d'))
    elif runPIP and not runCIP:
        print('Starting PIP size distribution processing... ' + dt.strftime(dt.now(),'%H:%M:%S %Y%m%d'))
        os.system('/sw/matlab-r2016b/bin/matlab -nodisplay -r "runSizeDistPECAN(' + ddate + ',' + repr(0) + ',' + repr(1) +')"')
        print('PIP size distribution processing complete... ' + dt.strftime(dt.now(),'%H:%M:%S %Y%m%d'))
    