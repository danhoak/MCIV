#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is the wrapper for the MCIV/dMCIV channel mapping routines.
# Calls deriv_map or glitch_map as functions, with appropriate options.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys
from optparse import OptionParser

usage = """usage: %prog [options]

This script will measure the glitch rate across the phase space of the IFO and the associated 
glitch distribution over a given epoch and SNR range.  Required options:

--ifo
--scan_type
--channelX
--channelY
--start_time
--end_time
--frame type (R, RDS, T, etc)
--sample rate
--filter_rate
--lower_snr
--upper_snr
--lower_freq
--upper_freq
--use_science

Example command line:

python MCIV_start.py -i H1 -t deriv -s 940032015 -e 940464015 -x SUS-ETMX_COIL_UL -f RDS -r 128 -l 6 -u 15 -g 60 -j 180 -a 1 -m 1

python MCIV_start.py -i H1 -t chans -s 940032015 -e 940464015 -x ASC-QPDX_Y -y ASC-QPDX_P -f RDS -r 128 -l 6 -u 1000 -g 40 -j 2048 -m 1


"""

parser = OptionParser(usage=usage)

parser.add_option("-i", "--ifo", action="store", type="str",\
                      help="IFO (H1, L1)")

parser.add_option("-t", "--scan_type", action="store", type="str", \
                      help="type of scan: chans (X vs Y), deriv (X vs dX/dt)")

parser.add_option("-x", "--horizontal_channel", action="store", type="str",\
                      help="horizontal channel")

parser.add_option("-y", "--vertical_channel", action="store", type="str", default=None, \
                      help="vertical channel for scan (not required if using deriv scan type)")

parser.add_option("-s", "--start_time", action="store", type="str",\
                      help="start gps time for data")

parser.add_option("-e", "--end_time", action="store", type="str",\
                      help="end gps time for data")

parser.add_option("-f", "--frame_type", action="store", type="str",\
                      help="desired frame type for both channels")

parser.add_option("-r", "--rate", action="store", type="str", default='128', \
                      help="desired sample rate for both channels (sets AA filter rolloff)")

parser.add_option("-l", "--lower_snr", action="store", type="int", default=6, \
                      help="lower snr threshold, above which to count triggers (inclusive)")

parser.add_option("-u", "--upper_snr", action="store", type="int", default=1000, \
                      help="upper snr threshold, below which to count triggers (inclusive)")

parser.add_option("-g", "--lower_freq", action="store", type="int", default=40, \
                      help="lower freq threshold, above which to count triggers (inclusive)")

parser.add_option("-j", "--upper_freq", action="store", type="int", default=2048, \
                      help="upper freq threshold, below which to count triggers (inclusive)")

parser.add_option("-a", "--filter", action="store", type="int", default=None, \
                      help="rate for cheby1 rolloff filter")

parser.add_option("-m", "--use_science_mode", action="store", type="int", default='1', \
                      help="If true, consider science mode segments only")

(options, args) = parser.parse_args()

ifo = options.ifo
scan_type = options.scan_type
channelX = options.horizontal_channel
channelY = options.vertical_channel
start_time = options.start_time
end_time = options.end_time
sample_rate = options.rate
frame_choice = options.frame_type
snr_low = options.lower_snr
snr_high = options.upper_snr
freq_low = options.lower_freq
freq_high = options.upper_freq
filter_rate = options.filter
science_flag = options.use_science_mode

code_version = 'v0.1'

if not(filter_rate):
    filter_rate = int(sample_rate)

if scan_type == 'deriv':

    from deriv_map import deriv_map

    output_path = deriv_map(ifo,channelX,start_time,end_time,sample_rate,frame_choice,snr_low,snr_high,freq_low,freq_high,filter_rate,science_flag,code_version)

elif scan_type == 'chans':

    if not(channelY):
        print "Must specify a vertical channel for a channel map."
        sys.exit()

    else:

        from glitch_map import glitch_map

        output_path = glitch_map(ifo,channelX,channelY,start_time,end_time,sample_rate,frame_choice,snr_low,snr_high,freq_low,freq_high,filter_rate,science_flag,code_version)



from mask_gen import *

mask_gen(output_path,0.1)
