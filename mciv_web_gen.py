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

This script will generate an mciv.html webpage using plots & data in the user-defined web directory

Required options:

--path

Example command line: 

python mciv_web_gen.py -p detchar/MCIV/beta/H1/20Feb2010_01Apr2010/ASC-WFS1_QP_OMC-QPD1_P_OUT_DAQ/

"""

parser = OptionParser(usage=usage)

parser.add_option("-p", "--path", action="store", type="str",\
                      help="path to directory with veto mask file (assumes /home/dhoak/public_html root directory)")

(options, args) = parser.parse_args()

data_path = options.path

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

input_datafile = web_directory + 'veto_inputs.hdf5'
output_datafile = web_directory + 'veto_outputs.hdf5'

veto_segfile = web_directory + 'veto_segments_coalesced.txt'

html_file = 'mciv.html'
webfilename = web_directory + html_file


# Load hdf5 input file

# most of this is used for filling in text boxes on the webpage

g = h5py.File(input_datafile,'r')

dataset = g['ifo']
ifo = dataset[...]

# gpstimes of map epoch
dataset = g['map_epoch']
map_epoch = dataset[...]

# science segments used to generate map
#dataset = g['segments']
dataset = g['map_segments']
map_segments = dataset[...]

# glitch count map - to calculate glitch rate
dataset = g['glitch_count']
glitch_count = dataset[...]

# User threshold for user-defined veto mask
dataset = g['user_threshold']
user_threshold = dataset[...]

# Bin edges for time map histogram
#dataset = g['bin_edges']
#bin_edges = dataset[...]

dataset = g['channels']
channels = dataset[...]
channelX = channels[0]
channelY = channels[1]

dataset = g['spans']
spans = dataset[...]
xspan = spans[0]
yspan = spans[1]

dataset = g['means']
means = dataset[...]
xmean = means[0]
ymean = means[1]

dataset = g['snr_limits']
snr_limits = dataset[...]
snr_low = snr_limits[0]
snr_high = snr_limits[1]

dataset = g['freq_limits']
freq_limits = dataset[...]
freq_low = freq_limits[0]
freq_high = freq_limits[1]

dataset = g['frame_type']
frame_type = dataset[...]

# Sample rate used to generate time map
dataset = g['rate']
rate = dataset[...]

# Low-pass filter rate used to generate time map
dataset = g['filter_rate']
filter_rate = dataset[...]

dataset = g['scan_type']
scan_type = dataset[...]

g.close()

f = h5py.File(output_datafile,'r')

# gpstimes of veto epoch - dates for which vetoes were produced
dataset = f['veto_epoch']
veto_epoch = dataset[...]

# deadtime of vetoes in veto epoch in seconds
dataset = f['veto_deadtime']
veto_deadtime = dataset[...]

# science segments used to generate vetoes (science-cat1 time in veto_epoch)
dataset = f['veto_segments']
veto_segments = dataset[...]

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

# Veto mask used to generate veto segments
dataset = f['mask_label']
mask_label = dataset[...]

# Per-segment use percentage
#dataset = g['segment_UPV']
#segment_UPV = dataset[...]

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

# Deadtime in seconds, for veto segs in epoch over which glitch study was performed
#dataset = g['veto_in_epoch_time']
#veto_in_epoch_time = dataset[...]

#dataset = g['']
# = dataset[...]




# Save the science segments used to generate the time map to a txt file to link in the web page
science_segment_file = web_directory + 'timemap_scisegs.txt'
savetxt(science_segment_file,map_segments,fmt='%.0f', delimiter=' ')


# Convert import gpstimes into date strings

cmd = 'lalapps_tconvert ' + str(map_epoch[0])
timeMap_start = commands.getoutput(cmd)

cmd = 'lalapps_tconvert ' + str(map_epoch[1])
timeMap_end = commands.getoutput(cmd)


cmd = 'lalapps_tconvert ' + str(veto_epoch[0])
vetoEpoch_start = commands.getoutput(cmd)

cmd = 'lalapps_tconvert ' + str(veto_epoch[1])
vetoEpoch_end = commands.getoutput(cmd)


cmd = 'lalapps_tconvert ' + str(study_epoch[0])
glitchEpoch_start = commands.getoutput(cmd)

cmd = 'lalapps_tconvert ' + str(study_epoch[1])
glitchEpoch_end = commands.getoutput(cmd)



# Get the science time used to generate the map, calculate the mean glitch rate in the SNR and freq limits
num_map_segments, map_epoch_time = segmentlist_time(map_segments)
total_glitch_count = int(sum(sum(glitch_count)))
mean_rate = total_glitch_count/map_epoch_time


# Get the science time used to generate the veto segments, calculate the deadtime fraction in the veto epoch
num_veto_segments, veto_epoch_time = segmentlist_time(veto_segments)
veto_deadtime_fraction = veto_deadtime / veto_epoch_time

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
hveto_segfile = '/home/dhoak/detchar/S6_segments/' + ifo + '_' + S6_epoch + '_sciencecat1_hveto.txt'

cat2_seglist = segmentsUtils.fromsegwizard(open(cat2_segfile), coltype = float, strict=False).coalesce()
cat3_seglist = segmentsUtils.fromsegwizard(open(cat3_segfile), coltype = float, strict=False).coalesce()
hveto_seglist = segmentsUtils.fromsegwizard(open(hveto_segfile), coltype = float, strict=False).coalesce()

veto_seglist = segmentsUtils.fromsegwizard(open(study_veto_file), coltype = float, strict=False).coalesce()

# We're only interested in the vetoes that intersect the study epoch in question
vetoes_in_study_epoch = veto_seglist & study_epoch_seg

veto_in_study_epoch_time = abs(vetoes_in_study_epoch)

cat2_overlap_seglist = vetoes_in_study_epoch & cat2_seglist
cat3_overlap_seglist = vetoes_in_study_epoch & cat3_seglist
hveto_overlap_seglist = vetoes_in_study_epoch & hveto_seglist

cat2_overlap = abs(cat2_overlap_seglist) / veto_in_study_epoch_time
cat3_overlap = abs(cat3_overlap_seglist) / veto_in_study_epoch_time
hveto_overlap = abs(hveto_overlap_seglist) / veto_in_study_epoch_time


webfile = open(webfilename, 'w')

if scan_type == 'deriv':
    title = 'X vs dX/dt MCIV Output ' + ifo + ' ' + channelX[3:]
    header_text = 'X vs dX/dt MCIV Output'
elif scan_type == 'chans':
    title = 'X vs Y MCIV Output ' + ifo + ' ' + channelX[3:] + ' ' + channelY[3:]
    header_text = 'X vs Y MCIV Output'
else:
    print 'Scan type unrecognized!'
    sys.exit()

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

if scan_type=='deriv':
    webfile.write('<h1>X vs dX/dt Instrumental Veto</h1>\n')
else:
    webfile.write('<h1>X vs Y Instrumental Veto</h1>\n')

#webfile.write('<h2><font face="Verdana" size="5">Human added comments.</font></h2>\n')

webfile.write('<h2><font face="Times" size="5">Time Map Parameters</font></h2>\n')

webfile.write('<br><table border=1 cellspacing="1" cellpadding="5">\n')
#webfile.write('<tr>\n')
#webfile.write('<td><b> Codes </b></td>\n')
#webfile.write('<td><b> svnversion </b></td>\n')
#webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> IFO </td>\n')
webfile.write('<td> ' + ifo + '\n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Channel </td>\n')
webfile.write('<td> ' + channelX + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Channel Y </td>\n')
webfile.write('<td> ' + channelY + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Sample rate </td>\n')
webfile.write('<td> ' + str(rate) + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Low-pass filter rate </td>\n')
webfile.write('<td> ' + str(filter_rate) + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Frame type </td>\n')
webfile.write('<td> ' + frame_type + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Start time </td>\n')
webfile.write('<td> ' + timeMap_start + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Stop time </td>\n')
webfile.write('<td> ' + timeMap_end + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Lower SNR Limit </td>\n')
webfile.write('<td> ' + str(snr_low) + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Upper SNR Limit </td>\n')
webfile.write('<td> ' + str(snr_high) + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Lower Freq Limit </td>\n')
webfile.write('<td> ' + str(freq_low) + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Upper Freq Limit </td>\n')
webfile.write('<td> ' + str(freq_high) + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Channel X Mean </td>\n')
webfile.write('<td> ' + str(xmean) + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Channel X Span </td>\n')
webfile.write('<td> ' + str(xspan) + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Channel Y Mean </td>\n')
webfile.write('<td> ' + str(ymean) + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Channel Y Span </td>\n')
webfile.write('<td> ' + str(yspan) + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('</table><br>\n')



webfile.write('<br><TABLE BORDER=0><TR>')
webfile.write('<TD ALIGN="left"><A HREF="time_map.png"><IMG SRC="time_map.png" WIDTH=300 ALT="Time Map ALT" TITLE="Time Map"></A></TD>')
webfile.write('<TD ALIGN="left"><A HREF="log_glitch_count.png"><IMG SRC="log_glitch_count.png" WIDTH=300 ALT="Glitch Count ALT" TITLE="Glitch Count"></A></TD>')
webfile.write('<TD ALIGN="left"><A HREF="log_glitch_rate.png"><IMG SRC="log_glitch_rate.png" WIDTH=300 ALT="Glitch Rate ALT" TITLE="Glitch Rate"></A></TD>')
webfile.write('<TD ALIGN="left"><A HREF="glitch_probability_map.png"><IMG SRC="glitch_probability_map.png" WIDTH=300 ALT="Glitch Probability ALT" TITLE="Glitch Probability"></A></TD>')
webfile.write('</TR></TABLE>')

webfile.write('<br><TABLE BORDER=0><TR>')
webfile.write('<TD ALIGN="left"><A HREF="time_map_surface.png"><IMG SRC="time_map_surface.png" WIDTH=300 ALT="Time Map 3D ALT" TITLE="Time Map 3D"></A></TD>')
webfile.write('<TD ALIGN="left"><A HREF="time_map_surface_linear.png"><IMG SRC="time_map_surface_linear.png" WIDTH=300 ALT="Time Map 3D Linear ALT" TITLE="Time Map 3D Linear"></A></TD>')
webfile.write('<TD ALIGN="left"><A HREF="glitch_count_surface_linear.png"><IMG SRC="glitch_count_surface_linear.png" WIDTH=300 ALT="Glitch Count 3D ALT" TITLE="Glitch Count 3D 3D"></A></TD>')
webfile.write('<TD ALIGN="left"><A HREF="glitch_rate_surface_linear.png"><IMG SRC="glitch_rate_surface_linear.png" WIDTH=300 ALT="Glitch Rate 3D ALT" TITLE="Glitch Rate 3D"></A></TD>')
webfile.write('</TR></TABLE>')

webfile.write('<br><table border=1 cellspacing="1" cellpadding="5">\n')

webfile.write('<tr>\n')
webfile.write('<td> Total Science Time In Map </td>\n')
webfile.write('<td> ' + str(map_epoch_time) + ' seconds \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Number of Omega Events Included In Map </td>\n')
webfile.write('<td> ' + str(total_glitch_count) + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Mean Glitch Rate (in SNR and freq limits) </td>\n')
webfile.write('<td> ' + '%.2f' % mean_rate + ' Hz \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Link to Science Segments </td>\n')
webfile.write('<td> <A HREF="' + www_directory + 'timemap_scisegs.txt">timemap_scisegs.txt</A> \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Link to Time Map Log File </td>\n')
webfile.write('<td> <A HREF="' + www_directory + 'timemap_logfile.txt">timemap_logfile.txt</A> \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('</table><br>\n')

webfile.write('<h2><font face="Times" size="5">Veto Mask Generation</font></h2>\n')

webfile.write('<br>')
webfile.write('Veto masks at various rate thresholds:')
webfile.write('<br>')

webfile.write('<br><TABLE BORDER=0><TR>')
webfile.write('<TD ALIGN="left"><A HREF="peak_veto_mask.png"><IMG SRC="peak_veto_mask.png" WIDTH=300 ALT="Peak Mask ALT" TITLE="Peak Eff/Deadtime Mask"></A></TD>')
webfile.write('<TD ALIGN="left"><A HREF="user_veto_mask.png"><IMG SRC="user_veto_mask.png" WIDTH=300 ALT="User Mask ALT" TITLE="User Threshold Mask"></A></TD>')
webfile.write('<TD ALIGN="left"><A HREF="rate_veto_mask.png"><IMG SRC="rate_veto_mask.png" WIDTH=300 ALT="Rate Mask ALT" TITLE="Rate Mask"></A></TD>')
webfile.write('<TD ALIGN="left"><A HREF="ersatz_veto_mask.png"><IMG SRC="ersatz_veto_mask.png" WIDTH=300 ALT="Ersatz Mask ALT" TITLE="Ersatz Mask"></A></TD>')
webfile.write('</TR></TABLE>')

webfile.write('<br><TABLE BORDER=0><TR>')
webfile.write('<TD ALIGN="left"><A HREF="peak_veto_mask_simple.png"><IMG SRC="peak_veto_mask_simple.png" WIDTH=300 ALT="Peak Mask ALT" TITLE="Peak Eff/Deadtime Mask"></A></TD>')
webfile.write('<TD ALIGN="left"><A HREF="user_veto_mask_simple.png"><IMG SRC="user_veto_mask_simple.png" WIDTH=300 ALT="User Mask ALT" TITLE="User Threshold Mask"></A></TD>')
webfile.write('<TD ALIGN="left"><A HREF="rate_veto_mask_simple.png"><IMG SRC="rate_veto_mask_simple.png" WIDTH=300 ALT="Rate Mask ALT" TITLE="Rate Mask"></A></TD>')
webfile.write('</TR></TABLE>')

webfile.write('<br>')
webfile.write('<br>')
webfile.write('Approximate efficiency and deadtime for various rate thresholds:')
webfile.write('<br>')

webfile.write('<br><TABLE BORDER=0><TR>')
webfile.write('<TD ALIGN="left"><A HREF="eff_dead.png"><IMG SRC="eff_dead.png" WIDTH=300 ALT="Eff/Dead ALT" TITLE="Eff & Deadtime Curves for ' + mask_label + ' Mask"></A></TD>')
webfile.write('</TR></TABLE>')

webfile.write('<h2><font face="Times" size="5">Veto Segment Generation</font></h2>\n')

webfile.write('<br><table border=1 cellspacing="1" cellpadding="5">\n')

webfile.write('<tr>\n')
webfile.write('<td> Start time </td>\n')
webfile.write('<td> ' + vetoEpoch_start + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Stop time </td>\n')
webfile.write('<td> ' + vetoEpoch_end + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Veto mask </td>\n')
webfile.write('<td> ' + mask_label + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Science time in veto epoch </td>\n')
webfile.write('<td> ' + str(veto_epoch_time) + ' seconds \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Total veto deadtime </td>\n')
webfile.write('<td> ' + '%.2f' % veto_deadtime + ' seconds \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Veto deadtime fraction </td>\n')
webfile.write('<td> ' + '%.4f' % veto_deadtime_fraction + ' \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Link to veto segments </td>\n')
webfile.write('<td> <A HREF="' + www_directory + 'veto_segments_coalesced.txt">veto_segments_coalesced.txt</A> \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('<tr>\n')
webfile.write('<td> Link to veto generation log file </td>\n')
webfile.write('<td> <A HREF="' + www_directory + 'vetogen_logfile.txt">vetogen_logfile.txt</A> \n')
webfile.write(' </td>\n')
webfile.write('</tr>\n')

webfile.write('</table><br>\n')

webfile.write('<br><TABLE BORDER=0><TR>')
webfile.write('<TD ALIGN="left"><A HREF="vetoDurationHist.png"><IMG SRC="vetoDurationHist.png" WIDTH=300 ALT="Veto Duration ALT" TITLE="Veto Segment Duration Histogram"></A></TD>')
webfile.write('</TR></TABLE>')

webfile.write('<h2><font face="Times" size="5">Veto Segment Properties</font></h2>\n')

webfile.write('<br><table border=1 cellspacing="1" cellpadding="5">\n')

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
webfile.write('<td> Overlap with hveto</td>\n')
webfile.write('<td> ' + '%.4f' % hveto_overlap + ' \n')
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
