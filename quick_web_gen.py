#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script will generate an MCIV results webpage
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import commands, sys, time, os
from optparse import OptionParser
import h5py
from MCIV_tools import *
from glue import segmentsUtils
from glue.segments import *

usage = """usage: %prog [options]

This script will generate an html webpage using plots & data in the user-defined web directory

Required options:

--path

Example command line: 

python mciv_web_gen.py -p detchar/MCIV/beta/H1/20Feb2010_01Apr2010/ASC-WFS1_QP_OMC-QPD1_P_OUT_DAQ/

"""

parser = OptionParser(usage=usage)

parser.add_option("-p", "--path", action="store", type="str",\
                      help="path to directory with veto mask file (assumes /home/dhoak/public_html root directory)")

parser.add_option("-i", "--ifo", action="store", type="str",\
                      help="path to directory with veto mask file (assumes /home/dhoak/public_html root directory)")

parser.add_option("-e", "--epoch", action="store", type="str",\
                      help="")

parser.add_option("-c", "--category", action="store", type="str",\
                      help="")


(options, args) = parser.parse_args()

data_path = options.path
ifo = options.ifo
epoch = options.epoch
category = options.category

# Function to get total time in a list of segments stored in array or list
# Kludgy but necessary...hdf5 stores glue.segments objects as arrays
def segmentlist_time(seglist):
    total_time = 0.0
    j=0
    for segment in seglist:
        j += 1
        #seg_start,seg_stop = segment.split()
        start = int(segment[0])
        stop = int(segment[1])
        total_time += stop-start

    return j, total_time

root_directory = '/home/dhoak/public_html/'
web_directory = root_directory + data_path

www_root = 'https://ldas-jobs.ligo.caltech.edu/~dhoak/'
www_directory = www_root + data_path

#input_datafile = web_directory + 'veto_inputs.hdf5'
output_datafile = web_directory + 'veto_outputs.hdf5'

#veto_segfile = web_directory + 'veto_segments_coalesced.txt'

html_file = 'veto_study.html'
webfilename = web_directory + html_file


f = h5py.File(output_datafile,'r')

# gpstimes of study epoch
dataset = f['study_epoch']
study_epoch = dataset[...]

# source file for vetoes in study
dataset = f['study_veto_file']
study_veto_file = dataset[...]

# deadtime of vetoes in study epoch, in seconds
dataset = f['study_deadtime']
study_deadtime = dataset[...]

# science segments used to study vetoes (science-cat1 time in study_epoch)
dataset = f['study_segments']
study_segments = dataset[...]

# Number of veto segments in study epoch
dataset = f['num_vetoes_in_epoch']
num_vetoes_in_epoch = dataset[...]

# Efficiency for SNR thresholds
dataset = f['efficiency_per_SNR']
efficiency_per_SNR = dataset[...]

# Number of used veto segments in study epoch
dataset = f['used_vetoes_in_epoch']
used_vetoes_in_epoch = dataset[...]

# Number of hardware injections in study epoch
dataset = f['num_inj']
num_inj = dataset[...]

# Number of vetoed hardware injections in study epoch
dataset = f['num_vetoed_inj']
num_vetoed_inj = dataset[...]

# Safety significance
dataset = f['safety_significance']
safety_significance = dataset[...]

# Number of glitches in study epoch
dataset = f['num_study_glitches']
num_study_glitches = dataset[...]

# Number of vetoed glitches in study epoch
dataset = f['num_study_vetoed_glitches']
num_study_vetoed_glitches = dataset[...]

f.close()

cmd = 'lalapps_tconvert ' + str(study_epoch[0])
glitchEpoch_start = commands.getoutput(cmd)

cmd = 'lalapps_tconvert ' + str(study_epoch[1])
glitchEpoch_end = commands.getoutput(cmd)



# Get the science time used to generate the veto segments, calculate the deadtime fraction in the veto epoch
num_study_segments, study_epoch_time = segmentlist_time(study_segments)
study_deadtime_fraction = study_deadtime / study_epoch_time

study_glitch_rate = num_study_glitches/study_epoch_time
study_veto_glitch_rate = num_study_vetoed_glitches/study_deadtime

# Calculate efficiency/deadtime values for SNR threshold (in study_epoch)

snr_level = [5, 6, 8, 10, 20, 50, 100]

EFF_DED = efficiency_per_SNR / study_deadtime_fraction
used_percentage = used_vetoes_in_epoch/num_vetoes_in_epoch


# Calculate overlap with cat2, cat3 in study epoch

# Generate a segment covering the epoch for the glitch study
study_epoch_start = study_epoch[0]
study_epoch_stop = study_epoch[1]
study_epoch_seg = segmentlist(segmentsUtils.segmentlist_range(study_epoch_start, study_epoch_stop, study_epoch_stop-study_epoch_start))

S6_epoch = getS6epoch(study_epoch_start,study_epoch_stop)

cat2_segfile = '/home/dhoak/detchar/S6_segments/' + ifo + '_' + S6_epoch + '_sciencecat1_cat2veto.txt'
cat3_segfile = '/home/dhoak/detchar/S6_segments/' + ifo + '_' + S6_epoch + '_sciencecat1_cat3veto.txt'

cat2_seglist = segmentsUtils.fromsegwizard(open(cat2_segfile), coltype = float, strict=False).coalesce()
cat3_seglist = segmentsUtils.fromsegwizard(open(cat3_segfile), coltype = float, strict=False).coalesce()

veto_seglist = segmentsUtils.fromsegwizard(open(study_veto_file), coltype = float, strict=False).coalesce()

# We're only interested in the vetoes that intersect the study epoch in question
vetoes_in_study_epoch = veto_seglist & study_epoch_seg

veto_in_study_epoch_time = abs(vetoes_in_study_epoch)

cat2_overlap_seglist = vetoes_in_study_epoch & cat2_seglist
cat3_overlap_seglist = vetoes_in_study_epoch & cat3_seglist

cat2_overlap = abs(cat2_overlap_seglist) / veto_in_study_epoch_time
cat3_overlap = abs(cat3_overlap_seglist) / veto_in_study_epoch_time


webfile = open(webfilename, 'w')

title = 'Veto Study for ' + ifo + ' ' + category + ' in ' + epoch
header_text = 'Veto Study Results'

webfile.write('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n')
webfile.write('<html>\n')
webfile.write('<head>\n')
webfile.write('  <title> ' + title + ' </title>\n')
webfile.write('  <style type="text/css">\n')
webfile.write('body {\n')
#webfile.write('font-family: Helvetica,Geneva,sans-serif;\n')
# Georgia is a nice font, too.
webfile.write('font-family: Garamond,Times,serif;\n')
webfile.write('color: black;\n')
webfile.write('background-color: white;\n')
webfile.write('}\n')
webfile.write('h1 {\n')
webfile.write('color: #ffffff;\n')
#webfile.write('background-color: #000088;\n')
webfile.write('background-color: #881c1c;\n')
webfile.write('padding: 0.35em;\n')
webfile.write('border: 1px solid black;\n')
webfile.write('}\n')
webfile.write('h2 {\n')
webfile.write('color: #000000;\n')
#webfile.write('background-color: #cccccc;\n')
webfile.write('background-color: #b7b6b6;\n')
webfile.write('padding: 0.35em;\n')
webfile.write('border: 1px solid black;\n')
webfile.write('}\n')
webfile.write('  </style>\n')
webfile.write(' <script type="text/javascript">\n')
webfile.write(' function toggleVisible(division) {\n')
webfile.write(' if (document.getElementById("div_" + division).style.display == "none") {\n')
webfile.write(' document.getElementById("div_" + division).style.display = "block";\n')
webfile.write(' document.getElementById("input_" + division).checked = true;\n')
webfile.write(' } else {\n')
webfile.write(' document.getElementById("div_" + division).style.display = "none";\n')
webfile.write(' document.getElementById("input_" + division).checked = false;\n')
webfile.write(' }\n')
webfile.write(' }\n')
webfile.write(' </script>\n')
webfile.write('</head>\n')
webfile.write('<body>\n')

webfile.write('<h1>Veto Study for ' + ifo + ' ' + category + ' in ' + epoch + '</h1>\n')

#webfile.write('<h2><font face="Verdana" size="5">Human added comments.</font></h2>\n')


webfile.write('<h2><font face="Times" size="5">Veto Segment Properties</font></h2>\n')

webfile.write('<br><table border=1 cellspacing="1" cellpadding="5">\n')

webfile.write('<tr>\n')
webfile.write('<td> IFO </td>\n')
webfile.write('<td> ' + ifo + '\n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Epoch </td>\n')
webfile.write('<td> ' + epoch + '\n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Start time for veto study </td>\n')
webfile.write('<td> ' + glitchEpoch_start + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Stop time for veto study </td>\n')
webfile.write('<td> ' + glitchEpoch_end + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Science time in study epoch </td>\n')
webfile.write('<td> ' + str(study_epoch_time) + ' seconds \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Number of veto segments</td>\n')
webfile.write('<td> ' + str(int(num_vetoes_in_epoch)) + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Veto deadtime in study epoch</td>\n')
webfile.write('<td> ' + '%.2f' % study_deadtime + ' seconds \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Deadtime fraction in study epoch</td>\n')
webfile.write('<td> ' + '%.4f' % study_deadtime_fraction + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Hardware injections in study epoch</td>\n')
webfile.write('<td> ' + str(int(num_inj)) + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Vetoed hardware injections</td>\n')
webfile.write('<td> ' + str(int(num_vetoed_inj)) + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Safety significance</td>\n')
webfile.write('<td> ' + '%.2f' % safety_significance + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Overlap with category 2</td>\n')
webfile.write('<td> ' + '%.4f' % cat2_overlap + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Overlap with category 3</td>\n')
webfile.write('<td> ' + '%.4f' % cat3_overlap + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Link to glitch study log file </td>\n')
webfile.write('<td> <A HREF="' + www_directory + 'segstudy_logfile.txt">segstudy_logfile.txt</A> \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('</table><br>\n')

webfile.write('<br><TABLE BORDER=0><TR>')
webfile.write('<TD ALIGN="left"><A HREF="studyvetoDurationHist.png"><IMG SRC="studyvetoDurationHist.png" WIDTH=300 ALT="Veto Duration ALT" TITLE="Veto Segment Duration Histogram"></A></TD>')
webfile.write('</TR></TABLE>')


webfile.write('<h2><font face="Times" size="5">Vetoed Glitch Properties</font></h2>\n')

webfile.write('<br><table border=1 cellspacing="1" cellpadding="5">\n')

webfile.write('<tr>\n')
webfile.write('<td> Glitch rate in study epoch, SNR > 5</td>\n')
webfile.write('<td> ' + '%.2f' % study_glitch_rate + ' Hz \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Glitch rate in veto segments</td>\n')
webfile.write('<td> ' + '%.2f' % study_veto_glitch_rate + ' Hz \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> SNR > 5 veto efficiency </td>\n')
webfile.write('<td> ' + '%.4f' % efficiency_per_SNR[0] + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Use percentage of veto segments </td>\n')
webfile.write('<td> ' + '%.4f' % used_percentage + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('</table><br>\n')

webfile.write('<br><TABLE BORDER=0><TR>')
webfile.write('<TD ALIGN="left"><A HREF="strengthFreq.png"><IMG SRC="strengthFreq.png" WIDTH=300 ALT="SNR FREQ ALT" TITLE="Freq vs SNR"></A></TD>')
webfile.write('<TD ALIGN="left"><A HREF="snrHist.png"><IMG SRC="snrHist.png" WIDTH=300 ALT="SNR Hist ALT" TITLE="SNR Histogram"></A></TD>')
webfile.write('<TD ALIGN="left"><A HREF="frequencyHist.png"><IMG SRC="frequencyHist.png" WIDTH=300 ALT="Freq Hist ALT" TITLE="Frequency Histogram"></A></TD>')
webfile.write('<TD ALIGN="left"><A HREF="vetoPercentFreqBin.png"><IMG SRC="vetoPercentFreqBin.png" WIDTH=300 ALT="Perc Freq ALT" TITLE="Percent Frequency Vetoed"></A></TD>')

webfile.write('<br><TABLE BORDER=0><TR>')
webfile.write('<TD ALIGN="left"><A HREF="segment_eff_dead_stats.png"><IMG SRC="segment_eff_dead_stats.png" WIDTH=300 ALT="Veto Performance ALT" TITLE="Veto Segment Statistics"></A></TD>')
webfile.write('</TR></TABLE>')


webfile.write('</TR></TABLE>')


webfile.write('<h2><font face="Times" size="5">Veto Efficiency vs SNR</font></h2>\n')


webfile.write('<br><table border=1 cellspacing="1" cellpadding="5">\n')

webfile.write('<tr>\n')
webfile.write('<td> SNR </td>\n')
webfile.write('<td> <b> ' + str(snr_level[0]) + ' </b> \n')
webfile.write('<td> <b> ' + str(snr_level[1]) + ' </b> \n')
webfile.write('<td> <b> ' + str(snr_level[2]) + ' </b> \n')
webfile.write('<td> <b> ' + str(snr_level[3]) + ' </b> \n')
webfile.write('<td> <b> ' + str(snr_level[4]) + ' </b> \n')
webfile.write('<td> <b> ' + str(snr_level[5]) + ' </b> \n')
webfile.write('<td> <b> ' + str(snr_level[6]) + ' </b> \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Efficiency </td>\n')
webfile.write('<td> ' + '%.3f' % efficiency_per_SNR[0] + ' \n')
webfile.write('<td> ' + '%.3f' % efficiency_per_SNR[1] + ' \n')
webfile.write('<td> ' + '%.3f' % efficiency_per_SNR[2] + ' \n')
webfile.write('<td> ' + '%.3f' % efficiency_per_SNR[3] + ' \n')
webfile.write('<td> ' + '%.3f' % efficiency_per_SNR[4] + ' \n')
webfile.write('<td> ' + '%.3f' % efficiency_per_SNR[5] + ' \n')
webfile.write('<td> ' + '%.3f' % efficiency_per_SNR[6] + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Eff / Deadtime </td>\n')
webfile.write('<td> ' + '%.2f' % EFF_DED[0] + ' \n')
webfile.write('<td> ' + '%.2f' % EFF_DED[1] + ' \n')
webfile.write('<td> ' + '%.2f' % EFF_DED[2] + ' \n')
webfile.write('<td> ' + '%.2f' % EFF_DED[3] + ' \n')
webfile.write('<td> ' + '%.2f' % EFF_DED[4] + ' \n')
webfile.write('<td> ' + '%.2f' % EFF_DED[5] + ' \n')
webfile.write('<td> ' + '%.2f' % EFF_DED[6] + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('</table><br>\n')

webfile.write('\n')
webfile.write('\n')
webfile.write('\n')
webfile.write('\n')
webfile.write('\n')
webfile.write('\n')
webfile.write('\n')
webfile.write('\n')
webfile.write('\n')
webfile.write('\n')
webfile.write('\n')
webfile.write('\n')



webfile.write('\n')

webfile.close()
