#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script calls seg_study as a function to collect efficiency/deadtime statistics for a set of veto segments
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import commands, sys, time, os
from optparse import OptionParser
#import h5py
#from MVIC_tools import *
#from glue import segmentsUtils
#from glue.segments import *
from seg_study import *

usage = """usage: %prog [options]

This script will load a list of veto segments and study their efficiency for a given epoch.

Required options:

--output_path
--ifo
--segment_file
--start_time
--end_time


Example command line: 

python study_vetoes.py -p detchar/MCIV/v0.1/H1/01Oct2009_01Nov2009/deriv_SUS-ETMX_COIL_UL/ -s 940032015 -e 940464015 -i H1 -f same

python study_vetoes.py -p detchar/MCIV/v0.0/H1/S6B/cat3/ -s 937785615 -e 947203215 -i H1 -f /home/dhoak/detchar/S6_segments/H1_S6B_sciencecat1_cat3veto.txt

"""

parser = OptionParser(usage=usage)

parser.add_option("-p", "--output_path", action="store", type="str",\
                      help="output path for plots and log file, assumes /home/dhoak/public_html/")

parser.add_option("-s", "--start_time", action="store", type="str",\
                      help="start gps time for epoch")

parser.add_option("-e", "--end_time", action="store", type="str",\
                      help="end gps time for epoch")

parser.add_option("-i", "--ifo", action="store", type="str",\
                      help="IFO (H1, L1, etc)")

parser.add_option("-f", "--segment_file", action="store", type="str", default='same',\
                      help="absolute path to file with veto segments; if 'same' then uses 'veto_segments_coalesced.txt' in output_path")


(options, args) = parser.parse_args()

data_path = options.output_path
ifo = options.ifo
start_time = options.start_time
end_time = options.end_time
segment_file = options.segment_file

root_directory = '/home/dhoak/public_html/'
web_directory = root_directory + data_path

if segment_file == 'same':
    veto_segfile = web_directory + 'veto_segments_coalesced.txt'
else:
    veto_segfile = segment_file


efficiency_per_SNR, deadtime, num_inj, num_vetoed_inj, safety_significance = seg_study(ifo,data_path,start_time,end_time,veto_segfile)


# Run web page generation
"""
if segment_file == 'same':

    mciv_web_gen(data_path)

else:

    study_web_gen(data_path)
"""
