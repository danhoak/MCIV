#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script generates veto segments for a given veto mask and GPS interval.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
matplotlib.use("Agg")
import pylab
import commands, sys, time, os
from optparse import OptionParser
import h5py
from scipy import signal
from MCIV_tools import *
from glue import segmentsUtils
from glue.segments import *

usage = """usage: %prog [options]

This script will load a two-channel veto mask and generate veto segments, based on that mask, for science segments over a given interval.  

Required options:

--path
--start_time
--end_time
--scan_type
--mask_label
--plot_flag

IFO, frame_type, sample rate, filter rate, and channel list are taken from snr6_data.hdf5 file in 'path', to ensure that vetoes are generated from the same data 
that was used to generate the masks.

Example command line: 

python veto_gen.py -p detchar/dMCIV/beta/H1/01Oct2009_01Dec2009/SUS-ETMX_COIL_UL/ -s 940032015 -e 940464015 -m user -t deriv

"""

parser = OptionParser(usage=usage)

parser.add_option("-p", "--path", action="store", type="str",\
                      help="path to directory with veto mask file (assumes /home/dhoak/public_html root directory)")

parser.add_option("-s", "--start_time", action="store", type="str",\
                      help="start gps time for data")

parser.add_option("-e", "--end_time", action="store", type="str",\
                      help="end gps time for data")

parser.add_option("-m", "--mask_label", action="store", type="str", default='peak', \
                      help="veto mask to use for veto segment generation")

parser.add_option("-t", "--scan_type", action="store", type="str", \
                      help="type of scan: chans (X vs Y), deriv (X vs dX/dt)")

parser.add_option("-f", "--plot_flag", action="store", type="int", default=0, \
                      help="boolean variable, generates X-Y channel map with vetoed regions when true")


(options, args) = parser.parse_args()

data_path = options.path
start_time = options.start_time
end_time = options.end_time
mask_label = options.mask_label
plot_flag = options.plot_flag

scan_type = options.scan_type

root_directory = '/home/dhoak/public_html/'
web_directory = root_directory + data_path
inputfile = web_directory + 'veto_inputs.hdf5'
datafile = web_directory + 'veto_outputs.hdf5'
logfilename = web_directory + 'vetogen_logfile.txt'

num_days = 1 + (int(end_time) - int(start_time))/86400

command_string = 'python veto_gen.py -p ' + data_path + ' -s ' + start_time + ' -e ' + end_time + \
    ' -m ' + mask_label + ' -f ' + str(plot_flag) + ' -t ' + str(scan_type)


# Load the necessary data fields from the data file in /path
# Grab the simplified masks for vetoing

g = h5py.File(inputfile,'r')

dataset = g['ifo']
ifo = dataset[...]

dataset = g['peak_veto_mask_simple']
peak_veto_mask = dataset[...]

dataset = g['rate_veto_mask_simple']
rate_veto_mask = dataset[...]

dataset = g['user_veto_mask_simple']
user_veto_mask = dataset[...]

dataset = g['ersatz_veto_mask']
ersatz_veto_mask = dataset[...]

dataset = g['user_threshold']
user_threshold = dataset[...]

dataset = g['bin_edges']
bin_edges = dataset[...]

bins = [bin_edges,bin_edges]

dataset = g['channels']
channels = dataset[...]
channelX = channels[0]
channelY = channels[1]

dataset = g['scan_type']
input_scan = dataset[...]

if input_scan != scan_type:
    print 'Scan type is different from veto inputs!'
    sys.exit()

if scan_type == 'deriv':
    channel_list = [channelX]
elif scan_type=='chans':
    channel_list = [channelX, channelY]
else:
    print 'Scan type unsupported!'
    sys.exit()    

dataset = g['spans']
spans = dataset[...]
xspan = spans[0]
yspan = spans[1]

dataset = g['means']
means = dataset[...]
xmean = means[0]
ymean = means[1]

dataset = g['frame_type']
frame_type = dataset[...]

dataset = g['rate']
rate = dataset[...]

sample_rate = str(rate)

dataset = g['filter_rate']
filter_rate = dataset[...]

g.close()

observatory = ifo[0]

# Get the current time, start time & end times of the analysis

cmd = 'lalapps_tconvert now'
current_gpstime = commands.getoutput(cmd)
current_gpstime_string = ''.join(current_gpstime.rsplit('\n'))
cmd = 'lalapps_tconvert ' + current_gpstime_string
current_date = commands.getoutput(cmd)
current_date_string = ''.join(current_date)

cmd = 'lalapps_tconvert ' + start_time
start_date = commands.getoutput(cmd)

cmd = 'lalapps_tconvert ' + end_time
end_date = commands.getoutput(cmd)

start_date_string, month, day, year = getDateString(start_time)
end_date_string, month, day, year = getDateString(end_time)

file_prefix = ifo + '_VetoGenSegs_' + start_date_string + '_' + end_date_string + '_'

# For now just use one veto mask...eventually generalize to all three
if mask_label == 'ersatz':
    veto_mask = ersatz_veto_mask
elif mask_label == 'peak':
    veto_mask = peak_veto_mask
elif mask_label == 'user':
    veto_mask = user_veto_mask
elif mask_label == 'rate':
    veto_mask = rate_veto_mask
else:
    print 'Veto mask unrecognized!'
    sys.exit()


# Start a log file to store output
logfile = open(logfilename, 'w',0)
logfile.write(current_date_string + '\n\n')
logfile.write(command_string + '\n\n')
logfile.write('Generating veto segments for ' + ifo + '.\n\n')
logfile.write('~'*100 + '\n\n')
logfile.write('Parameters for veto segments:\n\n')

if scan_type=='deriv':
    logfile.write('Scan type is X vs dX/dt.\n')
    logfile.write('Veto segments produced for channel ' + channelX + ' and derivative.\n')
else:
    logfile.write('Scan type is X vs Y.\n')
    logfile.write('Veto segments produced for channel ' + channelX + ' and ' + channelY + '.\n')

logfile.write('Number of days is ' + str(num_days) + '\n')
logfile.write('Frame type is ' + frame_type + '\n')
logfile.write('Sample rate is ' + sample_rate + '\n')
logfile.write('Low-pass filter rate is ' + str(filter_rate) + '\n')
logfile.write('Veto mask used is ' + mask_label + ' mask in veto_inputs.hdf5.\n\n')
logfile.write('~'*100 + '\n\n')


# Get science segments

S6_epoch = getS6epoch(start_time,end_time)

logfile.write('Collecting science data from ' + S6_epoch + ': ' + ''.join(start_date) + ' to ' + ''.join(end_date) + '\n\n')

science_cat1_filename = '/home/dhoak/detchar/S6_segments/' + ifo + '_' + S6_epoch + '_sciencecat1.txt'

logfile.write('Science-minus-CAT1 segment file is ' + science_cat1_filename + '\n')


#cmd1 = 'ligolw_segment_query -dqa "' + ifo + ':DMT-SCIENCE" -s ' + start_time + ' -e ' + end_time + ' -o ' + file_prefix + 'segments.xml'
#cmd2 = 'ligolw_print -t segment -c start_time -c end_time ' + file_prefix + 'segments.xml -d " "'

#logfile.write('Collecting science data from ' + ''.join(start_date) + ' to ' + ''.join(end_date) + '\n\n')
#logfile.write('Finding science segments...\n')

#commands.getoutput(cmd1)
#segment_list = commands.getoutput(cmd2)
#cmd3 = 'rm '+ file_prefix + 'segments.xml'
#commands.getoutput(cmd3)
    
#segments = segment_list.rsplit('\n')

epoch_start = int(start_time)
epoch_stop = int(end_time)
epoch_seg = segmentlist(segmentsUtils.segmentlist_range(epoch_start,epoch_stop,epoch_stop-epoch_start))

epoch_science = segmentsUtils.fromsegwizard(open(science_cat1_filename), coltype = float, strict=False).coalesce()

segments = epoch_seg & epoch_science

if not segments:
    print 'No science data in specified epoch.'
    sys.exit()

total_seg_time = abs(segments)
#for segment in segments:
#    seg_start,seg_stop = segment.split()        
#    start = int(seg_start)
#    stop = int(seg_stop)
#    total_seg_time = total_seg_time + stop - start

logfile.write('...done.\n\n')
logfile.write('Number of segments is ' + str(len(segments)) + '.  Total time in segments is ' + str(total_seg_time) + ' seconds.\n\n')

logfile.write('~'*100 + '\n\n')
logfile.write('Looping through science segments...\n\n')

# Initialize list to hold veto segs
veto_segs = []

clock_big = time.time()

#print len(bin_edges)
#print shape(veto_mask)

j=0
for segment in segments:
    j = j+1
    seg_start = str(int(segment[0]))
    seg_stop = str(int(segment[1]))
    #seg_start,seg_stop = segment.split()
    start = int(seg_start)
    stop = int(seg_stop)

    logfile.write('Segment ' + str(j) + ' ' + seg_start + ' ' + seg_stop + ' ' + str(int(seg_stop)-int(seg_start)) + '\n')
    clock1 = time.time()

    # Get the data for this segment and normalize in the same way as glitch_map
    vector = getChannelVector(ifo,seg_start,seg_stop,frame_type,channel_list,rate)

    q = int(rate)/filter_rate

    if q > 1:  # build the filter
        n = 3
        b, a = signal.cheby1(n, 0.05, 0.8 / q)
    elif q < 1:
        print 'Low pass filter rate is larger than sample rate!'
        sys.exit()

    if scan_type=='deriv':

        x = ( vector[:,0] - mean(vector[:,0]) ) / xspan
        if q > 1:
            z = signal.filtfilt(b, a, x)
        elif q==1:
            z = x

        y = gradient(z)*40  # fudge factor to get the scales right

    else:

        xx = ( vector[:,0] - xmean ) / xspan
        yy = ( vector[:,1] - ymean ) / yspan

        # Filter the data

        if q > 1:
            x = signal.filtfilt(b, a, xx)
            y = signal.filtfilt(b, a, yy)
        elif q==1:
            x = xx
            y = yy


    #data_points = transpose(vstack((y,x)))
    #rows, columns = shape(data_points)
    
    data_time_vector = arange(start,stop,1.0/rate)
    if rate > 32:
        veto_time_vector = arange(start,stop,1.0/32.0)
    else:
        veto_time_vector = arange(start,stop,1.0/rate)

    # Throw an error if the amount of returned data is not what we expect
    if len(data_time_vector) != len(x):
        print 'Wrong number of data points!'
        sys.exit()

    time_idx = 1

    # Loop through the time series
    for idx in range(len(x)):

        #data_point = vstack((y[idx],x[idx]))
        #hist_test, edges = histogramdd(data_point.T, bins)
        #bin_y, bin_x = unravel_index(hist_test.argmax(), hist_test.shape)

        # digitize is faster than histogramdd -- need to quantify this!
        # using digitize will sometimes place the data point in the wrong pixel -- this should only matter on the margins of the
        # veto region, need to quantify how much it changes the veto efficiency/deadtime

        bin_x = digitize([x[idx]],bin_edges)
        bin_y = digitize([y[idx]],bin_edges)

        # If the data point is outside of the pre-determined range from glitch_map, assign it to the largest bin
        # This should only be a problem if we're using digitize?
        if bin_x[0] >= len(bin_edges)-1:
            bin_x[0] = len(bin_edges)-2

        if bin_y[0] >= len(bin_edges)-1:
            bin_y[0] = len(bin_edges)-2

        if not(veto_mask[bin_y[0],bin_x[0]]):
        #if not(veto_mask[bin_y,bin_x]):

            # Build the veto segment.  Find the 1/32nd second clock count that is closest to the time of the sample.

            veto_time = data_time_vector[idx]

            # The digitize(x,bins) function returns the bin index i such that bins[i-1] <= x < bins[i].
            #vidx = digitize([veto_time],veto_time_vector)[0]
            vidx = digitize([veto_time],veto_time_vector[time_idx:])[0] + time_idx
            time_idx = vidx - 1

            # Instead of digitize use abs(t-time).argmin(), this allows us to save some time
            
            #vidx = (abs(veto_time_vector[time_idx:]-veto_time)).argmin() + time_idx
            #time_idx = vidx

            # We want our veto segments to start on the 1/32nd sec boundary before the 1/32nd sec that contains the veto time,
            # and stop on the 1/32nd sec boundary after the veto time, so that all veto segments are 3/32 sec long or longer,
            # and are always multiples of 1/32 sec.

            # If the preferred channel rate is less than 32Hz, the veto segments are defined using that rate, and we define the 
            # veto segment only using the bin that contains the veto time

            if vidx<2:
                veto_start = veto_time_vector[0]
            elif rate < 32.0:
                veto_start = veto_time_vector[vidx-1]
            else:
                veto_start = veto_time_vector[vidx-2]

            if vidx > len(veto_time_vector)-2:
                veto_stop = veto_time_vector[-1]
            elif rate < 32.0:
                veto_stop = veto_time_vector[vidx]
            else:
                veto_stop = veto_time_vector[vidx+1]

            veto_segment = [veto_start, veto_stop]

            veto_segs.append(veto_segment)  # this should use the python list append, not the numpy append (list append doesn't copy the whole array)
    

    clock2 = time.time()
    elapsed = clock2-clock1
    logfile.write('...done.  Elapsed time is ' + str("{0:.2f}".format(elapsed)) + ' seconds.\n\n')



veto_segmentlist = vstack((veto_segs))

clock3 = time.time()
elapsed = clock3-clock_big
logfile.write('Segment generation finished.\n')
logfile.write('Total elapsed time is ' + str(elapsed) + ' seconds.\n\n')

logfile.write('~'*100 + '\n\n')

logfile.write('Coalescing segments and saving to "veto_segments_coalesced.txt".\n')
clock1 = time.time()

# Coalesce & save the data

raw_vetofile = web_directory + 'veto_segments_raw.txt'
coalesced_vetofile = web_directory + 'veto_segments_coalesced.txt'

savetxt(raw_vetofile,veto_segmentlist,fmt='%.6f', delimiter=' ')

seglist = segmentsUtils.fromsegwizard(open(raw_vetofile), coltype = float).coalesce()

segfile = open(coalesced_vetofile, 'w',0)
segmentsUtils.tosegwizard(segfile, seglist, header=False, coltype = float)
segfile.close()

clock2 = time.time()
elapsed = clock2-clock1
logfile.write('...done.  Elapsed time is ' + str("{0:.2f}".format(elapsed)) + ' seconds.\n\n')

# Get durations of veto segments for histogram
data=genfromtxt(coalesced_vetofile)
durations=data[:,3]

total_veto_time = abs(seglist)

logfile.write('Total time in science segments is ' + str(total_seg_time) + ' seconds.\n\n')

logfile.write('Number of veto segments is ' + str(len(durations)) + '.\n')
logfile.write('Total time in veto segments is ' + str(total_veto_time) + ' seconds.\n')
logfile.write('Median veto segment duration is ' + str(median(durations)) + ' seconds.\n')
logfile.write('Deadtime for this veto over this interval is ' + str("{0:.2f}".format(total_veto_time/total_seg_time)) + '.\n\n')

# Save list of science segments, mask identifier, and epoch of science segments to datafile

veto_epoch = [int(start_time), int(end_time)]

# Start the veto output file
f = h5py.File(datafile,'w')

if 'veto_segments' in f:
    del f['veto_segments']

if 'mask_label' in f:
    del f['mask_label']

if 'veto_epoch' in f:
    del f['veto_epoch']

if 'veto_deadtime' in f:
    del f['veto_deadtime']

vetoseg_dset = f.create_dataset('veto_segments', data=segments)
mlabel_dset = f.create_dataset('mask_label', data=mask_label)
vetotime_dset = f.create_dataset('veto_epoch', data=veto_epoch)
veto_deadtime_dset = f.create_dataset('veto_deadtime', data=total_veto_time)

f.close()


logfile.write('Generating plots...')

matplotlib.rcParams.update({'savefig.dpi':250})

fignum=0
fignum=fignum+1
pylab.figure(fignum)

dmin = 0.08
dmax = 10.0

nbins = 50
logdmin=log10(dmin)
logdmax=log10(dmax)
dbinwidth = (logdmax-logdmin)/nbins
dbins=[0.0]*(nbins+1)
for i in range(len(dbins)):
   dbins[i] = pow(10,logdmin+i*dbinwidth)

n, b = histogram(durations, dbins)

plotbins = zeros((len(b)*2)) + 0.001
plotcounts = zeros((len(b)*2)) + 0.001

for i in range(len(b)):
    plotbins[2*i] = b[i] - b[i]*0.0003
    plotbins[2*i+1] = b[i] + b[i]*0.0003

for i in range(len(n)-1):
    plotcounts[2*i+1] = n[i] + 0.001
    plotcounts[2*i+2] = n[i] + 0.001

pylab.plot(plotbins,plotcounts,'k-',linewidth=1.0)

pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('Duration [sec]',fontsize=14)
pylab.grid(True, which='major', linestyle=':')
pylab.grid(True, which='minor', linestyle=':')
pylab.xlim(dmin,dmax)
pylab.ylim(0.7,1.2*max(plotcounts))
pylab.title('Veto Segment Duration')
pylab.savefig(web_directory + 'vetoDurationHist.png')


# If plot_flag is true, generate channel map with vetoed segments
# Probably should only do this if requested interval is short (~100sec)

if plot_flag:

    fignum=fignum+1
    pylab.figure(fignum)

    pylab.plot(x,y,'k-',linewidth=0.8)

    lower = 50
    upper = 200

    pylab.plot([bin_edges[lower],bin_edges[lower]],[bin_edges[0],bin_edges[249]],'r-.')
    pylab.plot([bin_edges[upper],bin_edges[upper]],[bin_edges[0],bin_edges[249]],'r-.')
    pylab.plot([bin_edges[0],bin_edges[249]],[bin_edges[lower],bin_edges[lower]],'r-.')
    pylab.plot([bin_edges[0],bin_edges[249]],[bin_edges[upper],bin_edges[upper]],'r-.')

    for veto in veto_segmentlist:
        veto_start = veto[0]
        veto_stop = veto[1]

        start_idx = (abs(data_time_vector-veto_start)).argmin()
        stop_idx = (abs(data_time_vector-veto_stop)).argmin()
        
        vx = x[start_idx:stop_idx]
        vy = y[start_idx:stop_idx]

        pylab.plot(vx,vy,'r-',linewidth=1.2)

    pylab.axis([-0.51, 0.51, -0.51, 0.51])
    pylab.grid(True)
    pylab.xticks(visible=False)
    pylab.yticks(visible=False)
    pylab.xlabel(channelX,fontsize=14)
    pylab.ylabel(channelY,fontsize=14)

    pylab.title('Channel Map With Vetoed Segments',fontsize=12)

    pylab.savefig(web_directory + 'plot_2chan_vetomap.png')

logfile.write('done.\n\n')
logfile.write('~'*100 + '\n\n')

cmd = 'lalapps_tconvert now'
current_gpstime = commands.getoutput(cmd)
current_gpstime_string = ''.join(current_gpstime.rsplit('\n'))
cmd = 'lalapps_tconvert ' + current_gpstime_string
current_date = commands.getoutput(cmd)
current_date_string = ''.join(current_date)

logfile.write(current_date_string + '\n')
logfile.close()
